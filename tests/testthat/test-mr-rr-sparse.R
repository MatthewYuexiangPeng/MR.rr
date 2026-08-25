sparse_test_data <- function() {
  set.seed(20260825)

  list(
    Y = matrix(stats::rnorm(40), nrow = 20, ncol = 2),
    X = matrix(stats::rnorm(60), nrow = 20, ncol = 3),
    Sigma_X = diag(rep(1e-4, 3)),
    W = matrix(c(1.5, 0.2, 0.2, 1.0), nrow = 2, ncol = 2)
  )
}


skip_if_sparse_solver_unavailable <- function() {
  testthat::skip_if_not_installed("CVXR")
  testthat::skip_if_not_installed("osqp")
  testthat::skip_if(
    !"OSQP" %in% CVXR::installed_solvers(),
    "The OSQP solver is not available to CVXR."
  )
}


reference_is_psd <- function(A, tol = 1e-8) {
  A <- as.matrix(A)
  A <- (A + t(A)) / 2
  min(eigen(A, symmetric = TRUE, only.values = TRUE)$values) >= -tol
}


reference_nearest_psd <- function(M, epsilon = 1e-6) {
  M <- (M + t(M)) / 2
  eig <- eigen(M, symmetric = TRUE)
  eig_values <- pmax(eig$values, epsilon)
  out <- eig$vectors %*% (eig_values * t(eig$vectors))
  (out + t(out)) / 2
}


reference_construct_gamma_tilde <- function(Y, X, Sigma_X) {
  n <- nrow(X)
  P_Y <- Y %*% solve(t(Y) %*% Y) %*% t(Y)
  P_Y_complement <- diag(1, n) - P_Y
  Sigma_X_hat <- t(X) %*% X / n
  projected_X <- P_Y %*% X
  matrix_part1 <-
    Sigma_X_hat - t(projected_X) %*% projected_X / n
  matrix_part2 <- matrix_part1 - Sigma_X

  if (!reference_is_psd(matrix_part2)) {
    matrix_part2 <- reference_nearest_psd(
      matrix_part2,
      epsilon = 1e-6
    )
  }

  R <- chol(matrix_part1)
  Q <- chol(matrix_part2)
  L <- solve(R) %*% Q

  P_Y %*% X + (P_Y_complement %*% X) %*% L
}


reference_mr_rr_sparse <- function(
    Y,
    X,
    W,
    Sigma_X,
    lambda,
    r,
    max_iter,
    tol,
    threshold) {
  n <- nrow(X)
  px <- ncol(X)

  W_sqrt <- .sqrt_matrix(W)
  gamma_tilde <- reference_construct_gamma_tilde(Y, X, Sigma_X)

  init_result <- mr_rr(
    Y = Y,
    X = X,
    r = r,
    Sigma_X = Sigma_X,
    W = W
  )
  A_hat <- init_result$A

  for (iter in seq_len(max_iter)) {
    B_var <- .cvxr_matrix_variable(r, px)

    loss <-
      CVXR::sum_squares(
        t(A_hat) %*% W %*% t(Y) - B_var %*% t(gamma_tilde)
      ) / n +
      Reduce(
        `+`,
        lapply(seq_len(px), function(k) {
          lambda[k] *
            CVXR::sum_entries(
              CVXR::norm1(B_var[, k, drop = FALSE])
            )
        })
      )

    problem <- CVXR::Problem(CVXR::Minimize(loss))
    B_hat <- .solve_cvxr_variable(
      problem = problem,
      variable = B_var,
      solver = "OSQP"
    )
    B_hat <- matrix(B_hat, nrow = r, ncol = px)

    svd_result <- svd(
      B_hat %*% t(gamma_tilde) %*% Y %*% W_sqrt
    )
    A_hat_new <-
      solve(W_sqrt) %*%
      svd_result$v %*%
      t(svd_result$u)

    dist <-
      norm(A_hat_new %*% B_hat - A_hat %*% B_hat, "F") /
      max(1e-8, norm(A_hat %*% B_hat, "F"))

    A_hat <- A_hat_new
    if (dist < tol) {
      break
    }
  }

  B_raw <- B_hat
  AB_raw <- A_hat %*% B_raw
  B_hat[abs(B_hat) < threshold] <- 0

  list(
    A = A_hat,
    B = B_hat,
    AB = A_hat %*% B_hat,
    B_raw = B_raw,
    AB_raw = AB_raw
  )
}


test_that("mr_rr_sparse matches the frozen paper algorithm", {
  skip_if_sparse_solver_unavailable()
  dat <- sparse_test_data()

  fit <- mr_rr_sparse(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    lambda = rep(1e-3, ncol(dat$X)),
    W = dat$W,
    max_iter = 3,
    tol = 1e-12,
    threshold = 1e-2,
    solver = "OSQP"
  )

  reference <- reference_mr_rr_sparse(
    Y = dat$Y,
    X = dat$X,
    W = dat$W,
    Sigma_X = dat$Sigma_X,
    lambda = rep(1e-3, ncol(dat$X)),
    r = 1,
    max_iter = 3,
    tol = 1e-12,
    threshold = 1e-2
  )

  expect_equal(fit$AB_raw, reference$AB_raw, tolerance = 1e-5)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-5)
  expect_equal(abs(fit$B_raw), abs(reference$B_raw), tolerance = 1e-5)
  expect_identical(fit$B == 0, reference$B == 0)

  expect_equal(dim(fit$A), c(2L, 1L))
  expect_equal(dim(fit$B), c(1L, 3L))
  expect_equal(dim(fit$AB), c(2L, 3L))
  expect_equal(fit$AB, fit$A %*% fit$B, tolerance = 1e-10)
  expect_equal(fit$AB_raw, fit$A %*% fit$B_raw, tolerance = 1e-10)

  below_threshold <- abs(fit$B_raw) < 1e-2
  expect_true(all(fit$B[below_threshold] == 0))
  expect_true(all(is.finite(fit$A)))
  expect_true(all(is.finite(fit$B_raw)))
  expect_true(fit$iter >= 1L && fit$iter <= 3L)
  expect_length(fit$dist, 1L)
  expect_type(fit$converged, "logical")
})


test_that("mr_rr_sparse validates its tuning controls", {
  dat <- sparse_test_data()

  common_args <- list(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    W = dat$W
  )

  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(lambda = c(1e-3, 1e-3)))),
    "`lambda`"
  )
  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(lambda = -1e-3))),
    "`lambda`"
  )
  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(max_iter = 0))),
    "`max_iter`"
  )
  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(tol = 0))),
    "`tol`"
  )
  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(threshold = -1))),
    "`threshold`"
  )
  expect_error(
    do.call(mr_rr_sparse, c(common_args, list(solver = ""))),
    "`solver`"
  )
})
