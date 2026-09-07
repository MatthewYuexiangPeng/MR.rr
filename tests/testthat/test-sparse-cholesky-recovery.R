test_that("the sparse fallback preserves successful Cholesky factors", {
  x <- matrix(c(2, 0.1, 0.1, 1), 2, 2)
  result <- .sparse_corrected_chol(x)
  expect_identical(result$factor, chol(x))
  expect_identical(result$covariance, x)
  expect_false(result$diagnostics$corrected_covariance_projected)
})

test_that("PSD tolerance does not hide a failed sparse Cholesky factor", {
  x <- diag(c(2e-4, 1e-4, -3.990090310468490e-9))
  expect_true(.is_psd(x))
  expect_error(chol(x))
  result <- .sparse_corrected_chol(x)
  expect_identical(result$diagnostics$projection_reason, "chol_failed")
  expect_identical(result$diagnostics$eigenvalue_floor, 1e-6)
  expect_equal(crossprod(result$factor), diag(c(2e-4, 1e-4, 1e-6)), tolerance = 1e-12)
})

test_that("the historical projection rule still handles clearly indefinite inputs", {
  x <- diag(c(2e-4, -2e-5))
  result <- .sparse_corrected_chol(x)
  expect_identical(result$diagnostics$projection_reason, "not_psd")
  expect_equal(crossprod(result$factor), .nearest_psd(x, 1e-6), tolerance = 1e-12)
})
