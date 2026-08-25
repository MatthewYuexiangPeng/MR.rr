.frozen_mr_rr_naive <- function(Y, X, r, W = NULL, W_inv = NULL) {
  frozen_sqrt_matrix <- function(mat) {
    eig <- eigen(mat)
    eig$vectors %*%
      diag(sqrt(eig$values)) %*%
      t(eig$vectors)
  }

  XtY <- crossprod(X, Y)
  XtX_inv <- solve(crossprod(X))

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
    t(XtY) %*%
    XtX_inv %*%
    XtY %*%
    sqrt_Gamma / nrow(Y)

  V <- eigen(V)$vectors[, seq_len(r), drop = FALSE]

  A <- sqrt_Gamma_inv %*% V
  B <- crossprod(V, sqrt_Gamma %*% t(XtY) %*% XtX_inv)

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}


test_that("mr_rr_naive reproduces the frozen unweighted estimator", {
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

  fit <- mr_rr_naive(Y = Y, X = X, r = 1)
  reference <- .frozen_mr_rr_naive(Y = Y, X = X, r = 1)

  expect_equal(fit$A, reference$A, tolerance = 1e-12)
  expect_equal(fit$B, reference$B, tolerance = 1e-12)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-12)
})


test_that("mr_rr_naive reproduces the frozen weighted estimator", {
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

  W <- matrix(
    c(
      2.0, 0.3,
      0.3, 1.5
    ),
    nrow = 2,
    byrow = TRUE
  )

  W_inv <- solve(W)

  fit <- mr_rr_naive(
    Y = Y,
    X = X,
    r = 1,
    W = W,
    W_inv = W_inv
  )

  reference <- .frozen_mr_rr_naive(
    Y = Y,
    X = X,
    r = 1,
    W = W,
    W_inv = W_inv
  )

  expect_equal(fit$A, reference$A, tolerance = 1e-12)
  expect_equal(fit$B, reference$B, tolerance = 1e-12)
  expect_equal(fit$AB, reference$AB, tolerance = 1e-12)

  expect_equal(
    crossprod(fit$A, W %*% fit$A),
    diag(1),
    tolerance = 1e-10
  )
})


test_that("mr_rr_naive returns matrices with the expected dimensions", {
  X <- rbind(
    diag(3),
    c(1.0,  1.0, 1.0),
    c(2.0, -1.0, 0.5),
    c(-0.5, 2.0, 1.0)
  )

  Y <- matrix(
    c(
      1,  2,
      2, -1,
      3,  4,
      4,  0,
      5,  3,
      6,  5
    ),
    nrow = 6,
    byrow = TRUE
  )

  fit <- mr_rr_naive(Y = Y, X = X, r = 1)

  expect_identical(dim(fit$A), c(2L, 1L))
  expect_identical(dim(fit$B), c(1L, 3L))
  expect_identical(dim(fit$AB), c(2L, 3L))
})
