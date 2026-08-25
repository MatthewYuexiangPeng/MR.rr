regularized_test_data <- function() {
  X <- rbind(
    diag(3),
    c(1.0, 1.0, 1.0),
    c(2.0, -1.0, 0.5),
    c(-0.5, 2.0, 1.0)
  )

  Y <- rbind(
    c(0.5, -0.2),
    c(-0.3, 0.7),
    c(0.4, 0.1),
    c(0.8, -0.5),
    c(1.2, 0.3),
    c(-0.7, 1.1)
  )

  list(
    X = X,
    Y = Y,
    Sigma_X = diag(c(0.01, 0.015, 0.02)),
    W = matrix(c(2.0, 0.3, 0.3, 1.5), nrow = 2)
  )
}

reference_sqrt_matrix <- function(x) {
  decomposition <- eigen(x)

  decomposition$vectors %*%
    diag(
      sqrt(decomposition$values),
      nrow = length(decomposition$values)
    ) %*%
    t(decomposition$vectors)
}

reference_mr_rr_regularized <- function(
    Y, X, r, Sigma_X, regularization_rate, W = NULL, W_inv = NULL) {
  n <- nrow(Y)

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)

  stable_Sigma_xx <-
    debiased_Sigma_xx +
    regularization_rate * debiased_Sigma_xx_inv

  stable_Sigma_xx_inv <- solve(stable_Sigma_xx)

  if (is.null(W)) {
    sqrt_W <- diag(ncol(Y))
    sqrt_W_inv <- sqrt_W
  } else {
    sqrt_W <- reference_sqrt_matrix(W)

    if (is.null(W_inv)) {
      sqrt_W_inv <- solve(sqrt_W)
    } else {
      sqrt_W_inv <- reference_sqrt_matrix(W_inv)
    }
  }

  target_matrix <-
    sqrt_W %*%
    t(Sigma_xy) %*%
    stable_Sigma_xx_inv %*%
    Sigma_xy %*%
    sqrt_W

  eigenvectors <- eigen(target_matrix)$vectors[
    , seq_len(r), drop = FALSE
  ]

  A <- sqrt_W_inv %*% eigenvectors

  B <- crossprod(
    eigenvectors,
    sqrt_W %*% t(Sigma_xy) %*% stable_Sigma_xx_inv
  )

  list(A = A, B = B, AB = A %*% B)
}

test_that("mr_rr_regularized matches the frozen unweighted estimator", {
  dat <- regularized_test_data()

  observed <- mr_rr_regularized(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    regularization_rate = 0.05
  )

  expected <- reference_mr_rr_regularized(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    regularization_rate = 0.05
  )

  expect_equal(observed$A, expected$A, tolerance = 1e-12)
  expect_equal(observed$B, expected$B, tolerance = 1e-12)
  expect_equal(observed$AB, expected$AB, tolerance = 1e-12)

  expect_equal(dim(observed$A), c(2L, 1L))
  expect_equal(dim(observed$B), c(1L, 3L))
  expect_equal(dim(observed$AB), c(2L, 3L))
})

test_that("mr_rr_regularized matches the frozen weighted estimator", {
  dat <- regularized_test_data()

  observed <- mr_rr_regularized(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    regularization_rate = 0.05,
    W = dat$W
  )

  expected <- reference_mr_rr_regularized(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    regularization_rate = 0.05,
    W = dat$W
  )

  expect_equal(observed$A, expected$A, tolerance = 1e-12)
  expect_equal(observed$B, expected$B, tolerance = 1e-12)
  expect_equal(observed$AB, expected$AB, tolerance = 1e-12)

  expect_equal(
    crossprod(observed$A, dat$W %*% observed$A),
    diag(1),
    tolerance = 1e-10
  )
})

test_that("zero regularization reproduces mr_rr", {
  dat <- regularized_test_data()

  regularized_fit <- mr_rr_regularized(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X,
    regularization_rate = 0
  )

  corrected_fit <- mr_rr(
    Y = dat$Y,
    X = dat$X,
    r = 1,
    Sigma_X = dat$Sigma_X
  )

  expect_equal(
    regularized_fit$AB,
    corrected_fit$AB,
    tolerance = 1e-10
  )
})

test_that("regularization rate is validated", {
  dat <- regularized_test_data()

  expect_error(
    mr_rr_regularized(
      dat$Y, dat$X, 1, dat$Sigma_X,
      regularization_rate = -0.1
    ),
    "regularization_rate"
  )

  expect_error(
    mr_rr_regularized(
      dat$Y, dat$X, 1, dat$Sigma_X,
      regularization_rate = Inf
    ),
    "regularization_rate"
  )

  expect_error(
    mr_rr_regularized(
      dat$Y, dat$X, 1, dat$Sigma_X,
      regularization_rate = c(0.1, 0.2)
    ),
    "regularization_rate"
  )
})
