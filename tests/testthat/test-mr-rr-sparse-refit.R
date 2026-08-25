sparse_refit_test_data <- function() {
  set.seed(20260826)

  list(
    Y = matrix(stats::rnorm(48), nrow = 24, ncol = 2),
    X = matrix(stats::rnorm(72), nrow = 24, ncol = 3),
    Sigma_X = diag(rep(1e-4, 3)),
    W = matrix(c(1.5, 0.2, 0.2, 1.0), nrow = 2, ncol = 2),
    support = matrix(c(TRUE, FALSE, TRUE), nrow = 1, ncol = 3)
  )
}


reference_mr_rr_sparse_refit <- function(
    Y,
    X,
    W,
    Sigma_X,
    support,
    A_init = NULL,
    B_init = NULL,
    max_iter = 100,
    tol = 1e-2) {
  px <- ncol(X)
  r <- nrow(support)

  gamma_tilde <- .construct_gamma_tilde(Y, X, Sigma_X)
  W_sqrt <- .sqrt_matrix(W)

  if (is.null(A_init) || is.null(B_init)) {
    init_result <- mr_rr(
      Y = Y,
      X = X,
      r = r,
      Sigma_X = Sigma_X,
      W = W
    )
    A_hat <- init_result$A
    B_hat <- init_result$B
  } else {
    A_hat <- A_init
    B_hat <- B_init
  }

  B_hat[!support] <- 0

  for (iter in seq_len(max_iter)) {
    B_hat <- matrix(0, nrow = r, ncol = px)

    for (k in seq_len(r)) {
      selected_idx <- which(support[k, ])

      if (length(selected_idx) > 0L) {
        latent_response <- as.vector(Y %*% W %*% A_hat[, k])
        B_hat[k, selected_idx] <- qr.solve(
          gamma_tilde[, selected_idx, drop = FALSE],
          latent_response
        )
      }
    }

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

  list(
    A = A_hat,
    B = B_hat,
    AB = A_hat %*% B_hat,
    support = support,
    iter = iter,
    dist = dist
  )
}


test_that("mr_rr_sparse_refit matches the frozen paper algorithm", {
  dat <- sparse_refit_test_data()

  init <- mr_rr(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    W = dat$W
  )

  fit <- mr_rr_sparse_refit(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    support = dat$support,
    W = dat$W,
    A_init = init$A,
    B_init = init$B,
    max_iter = 5,
    tol = 1e-12
  )

  reference <- reference_mr_rr_sparse_refit(
    Y = dat$Y,
    X = dat$X,
    W = dat$W,
    Sigma_X = dat$Sigma_X,
    support = dat$support,
    A_init = init$A,
    B_init = init$B,
    max_iter = 5,
    tol = 1e-12
  )

  expect_equal(fit$A, reference$A, tolerance = 1e-10)
  expect_equal(fit$B, reference$B, tolerance = 1e-10)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-10)
  expect_equal(fit$AB, fit$A %*% fit$B, tolerance = 1e-12)
  expect_identical(fit$support, dat$support)
  expect_true(all(fit$B[!dat$support] == 0))
  expect_true(all(is.finite(fit$A)))
  expect_true(all(is.finite(fit$B)))
  expect_true(fit$iter >= 1L && fit$iter <= 5L)
  expect_length(fit$dist, 1L)
  expect_type(fit$converged, "logical")
})


test_that("default and explicit MR-rr initialization agree", {
  dat <- sparse_refit_test_data()

  init <- mr_rr(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    W = dat$W
  )

  default_fit <- mr_rr_sparse_refit(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    support = dat$support,
    W = dat$W,
    max_iter = 3,
    tol = 1e-12
  )

  explicit_fit <- mr_rr_sparse_refit(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    support = dat$support,
    W = dat$W,
    A_init = init$A,
    B_init = init$B,
    max_iter = 3,
    tol = 1e-12
  )

  expect_equal(default_fit$A, explicit_fit$A, tolerance = 1e-12)
  expect_equal(default_fit$B, explicit_fit$B, tolerance = 1e-12)
  expect_equal(default_fit$AB, explicit_fit$AB, tolerance = 1e-12)
})


test_that("mr_rr_sparse_refit validates support and initialization", {
  dat <- sparse_refit_test_data()

  common_args <- list(
    Y = dat$Y,
    X = dat$X,
    Sigma_X = dat$Sigma_X,
    W = dat$W
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(common_args, list(support = c(TRUE, FALSE, TRUE)))
    ),
    "`support`"
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(common_args, list(support = matrix(TRUE, nrow = 1, ncol = 2)))
    ),
    "`support`"
  )

  invalid_support <- dat$support
  invalid_support[1, 1] <- NA
  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(common_args, list(support = invalid_support))
    ),
    "`support`"
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(
        common_args,
        list(support = matrix(c(1, 2, 0), nrow = 1, ncol = 3))
      )
    ),
    "`support`"
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(
        common_args,
        list(support = dat$support, A_init = matrix(1, nrow = 2, ncol = 1))
      )
    ),
    "both"
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(common_args, list(support = dat$support, max_iter = 0))
    ),
    "`max_iter`"
  )

  expect_error(
    do.call(
      mr_rr_sparse_refit,
      c(common_args, list(support = dat$support, tol = 0))
    ),
    "`tol`"
  )
})
