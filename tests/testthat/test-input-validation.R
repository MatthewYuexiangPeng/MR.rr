test_that("valid estimator inputs are returned in normalized form", {
  X <- data.frame(
    x1 = c(1, 2, 3, 4),
    x2 = c(2, 1, 4, 3)
  )

  Y <- data.frame(
    y1 = c(1, 3, 2, 5),
    y2 = c(2, 4, 1, 3)
  )

  result <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = 1,
    Sigma_X = diag(2) / 10,
    W = diag(2),
    W_inv = diag(2),
    require_sigma_x = TRUE
  )

  expect_true(is.matrix(result$X))
  expect_true(is.matrix(result$Y))
  expect_type(result$X, "double")
  expect_type(result$Y, "double")
  expect_identical(result$r, 1L)
})


test_that("X and Y must have compatible dimensions", {
  X <- matrix(1:8, nrow = 4)
  Y <- matrix(1:6, nrow = 3)

  expect_error(
    .validate_estimator_inputs(Y, X, r = 1),
    "same number of rows"
  )
})


test_that("rank must be an integer in the allowed range", {
  X <- matrix(1:12, nrow = 4, ncol = 3)
  Y <- matrix(1:8, nrow = 4, ncol = 2)

  expect_error(
    .validate_estimator_inputs(Y, X, r = 0),
    "between 1"
  )

  expect_error(
    .validate_estimator_inputs(Y, X, r = 3),
    "between 1"
  )

  expect_error(
    .validate_estimator_inputs(Y, X, r = 1.5),
    "single integer"
  )
})


test_that("Sigma_X must have the exposure dimension", {
  X <- matrix(1:12, nrow = 4, ncol = 3)
  Y <- matrix(1:8, nrow = 4, ncol = 2)

  expect_error(
    .validate_estimator_inputs(
      Y,
      X,
      r = 1,
      Sigma_X = diag(2),
      require_sigma_x = TRUE
    ),
    "3 by 3"
  )

  expect_error(
    .validate_estimator_inputs(
      Y,
      X,
      r = 1,
      require_sigma_x = TRUE
    ),
    "must be provided"
  )
})


test_that("weight matrices must be symmetric and positive definite", {
  X <- matrix(1:12, nrow = 4, ncol = 3)
  Y <- matrix(1:8, nrow = 4, ncol = 2)

  nonsymmetric_W <- matrix(
    c(
      1, 2,
      0, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  indefinite_W <- matrix(
    c(
      1, 2,
      2, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_error(
    .validate_estimator_inputs(
      Y,
      X,
      r = 1,
      W = nonsymmetric_W
    ),
    "symmetric"
  )

  expect_error(
    .validate_estimator_inputs(
      Y,
      X,
      r = 1,
      W = indefinite_W
    ),
    "positive definite"
  )
})


test_that("W_inv requires W", {
  X <- matrix(1:12, nrow = 4, ncol = 3)
  Y <- matrix(1:8, nrow = 4, ncol = 2)

  expect_error(
    .validate_estimator_inputs(
      Y,
      X,
      r = 1,
      W_inv = diag(2)
    ),
    "cannot be supplied"
  )
})


test_that("non-finite values are rejected", {
  X <- matrix(c(1, 2, 3, NA), nrow = 2)
  Y <- matrix(1:4, nrow = 2)

  expect_error(
    .validate_estimator_inputs(Y, X, r = 1),
    "finite values"
  )
})
