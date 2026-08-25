.frozen_mr_rr <- function(
    Y,
    X,
    r,
    Sigma_X,
    W = NULL,
    W_inv = NULL
) {
  frozen_sqrt_matrix <- function(mat) {
    eig <- eigen(mat)

    eig$vectors %*%
      diag(sqrt(eig$values)) %*%
      t(eig$vectors)
  }

  n <- nrow(Y)

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)

  if (is.null(W)) {
    sqrt_Gamma <- diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv <- sqrt_Gamma
  } else {
    sqrt_Gamma <- frozen_sqrt_matrix(W)

    if (is.null(W_inv)) {
      sqrt_Gamma_inv <- solve(sqrt_Gamma)
    } else {
      sqrt_Gamma_inv <- frozen_sqrt_matrix(W_inv)
    }
  }

  V <- sqrt_Gamma %*%
    t(Sigma_xy) %*%
    debiased_Sigma_xx_inv %*%
    Sigma_xy %*%
    sqrt_Gamma

  V <- eigen(V)$vectors[, seq_len(r), drop = FALSE]

  A <- sqrt_Gamma_inv %*% V

  B <- crossprod(
    V,
    sqrt_Gamma %*%
      t(Sigma_xy) %*%
      debiased_Sigma_xx_inv
  )

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}


.make_mr_rr_test_data <- function() {
  X <- matrix(
    c(
      0.2, -0.1,  0.4,
      0.5,  0.3, -0.2,
      -0.4,  0.6,  0.1,
      0.7, -0.5,  0.3,
      -0.2,  0.4,  0.8,
      0.1,  0.2, -0.6
    ),
    nrow = 6,
    byrow = TRUE
  )

  Y <- matrix(
    c(
      0.1,  0.3,
      0.4, -0.2,
      -0.3,  0.5,
      0.6, -0.4,
      0.2,  0.7,
      -0.1, -0.3
    ),
    nrow = 6,
    byrow = TRUE
  )

  list(
    X = X,
    Y = Y,
    Sigma_X = diag(c(0.01, 0.015, 0.02)),
    W = matrix(
      c(
        2.0, 0.3,
        0.3, 1.5
      ),
      nrow = 2,
      byrow = TRUE
    )
  )
}


test_that("mr_rr reproduces the frozen unweighted estimator", {
  data <- .make_mr_rr_test_data()

  fit <- mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X
  )

  reference <- .frozen_mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X
  )

  expect_equal(fit$A, reference$A, tolerance = 1e-12)
  expect_equal(fit$B, reference$B, tolerance = 1e-12)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-12)
})


test_that("mr_rr reproduces the frozen weighted estimator", {
  data <- .make_mr_rr_test_data()
  W_inv <- solve(data$W)

  fit <- mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X,
    W = data$W,
    W_inv = W_inv
  )

  reference <- .frozen_mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X,
    W = data$W,
    W_inv = W_inv
  )

  expect_equal(fit$A, reference$A, tolerance = 1e-12)
  expect_equal(fit$B, reference$B, tolerance = 1e-12)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-12)

  expect_equal(
    crossprod(fit$A, data$W %*% fit$A),
    diag(1),
    tolerance = 1e-10
  )
})


test_that("mr_rr reduces to mr_rr_naive when Sigma_X is zero", {
  data <- .make_mr_rr_test_data()
  zero_Sigma_X <- matrix(
    0,
    nrow = ncol(data$X),
    ncol = ncol(data$X)
  )

  corrected_fit <- mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = zero_Sigma_X,
    W = data$W
  )

  naive_fit <- mr_rr_naive(
    Y = data$Y,
    X = data$X,
    r = 1,
    W = data$W
  )

  expect_equal(
    corrected_fit$AB,
    naive_fit$AB,
    tolerance = 1e-10
  )
})


test_that("mr_rr requires Sigma_X", {
  data <- .make_mr_rr_test_data()

  expect_error(
    mr_rr(
      Y = data$Y,
      X = data$X,
      r = 1,
      Sigma_X = NULL
    ),
    "must be provided"
  )
})
