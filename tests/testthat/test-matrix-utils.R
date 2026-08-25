test_that(".second_moment computes the uncentered second moment", {
  x <- matrix(
    c(
      1, 2,
      3, 4,
      5, 6
    ),
    nrow = 3,
    byrow = TRUE
  )

  expected <- crossprod(x) / nrow(x)

  expect_equal(.second_moment(x), expected)
})


test_that(".sqrt_matrix reconstructs a positive-definite matrix", {
  x <- matrix(
    c(
      4, 1,
      1, 3
    ),
    nrow = 2,
    byrow = TRUE
  )

  x_sqrt <- .sqrt_matrix(x)
  x_inv_sqrt <- .sqrt_matrix(x, inv = TRUE)

  expect_equal(
    x_sqrt %*% x_sqrt,
    x,
    tolerance = 1e-10
  )

  expect_equal(
    x_inv_sqrt %*% x %*% x_inv_sqrt,
    diag(2),
    tolerance = 1e-10
  )
})


test_that(".sqrt_matrix handles a one-dimensional matrix", {
  x <- matrix(4, nrow = 1, ncol = 1)

  expect_equal(
    .sqrt_matrix(x),
    matrix(2, nrow = 1, ncol = 1)
  )
})


test_that(".is_psd identifies positive-semidefinite matrices", {
  positive_definite <- matrix(
    c(
      4, 1,
      1, 3
    ),
    nrow = 2,
    byrow = TRUE
  )

  indefinite <- matrix(
    c(
      1, 2,
      2, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  expect_true(.is_psd(positive_definite))
  expect_false(.is_psd(indefinite))
  expect_false(.is_psd(matrix(1:6, nrow = 2)))
})


test_that(".nearest_psd corrects an indefinite matrix", {
  x <- matrix(
    c(
      1, 2,
      2, 1
    ),
    nrow = 2,
    byrow = TRUE
  )

  corrected <- .nearest_psd(x, epsilon = 1e-6)

  expect_true(.is_psd(corrected))

  expect_equal(
    corrected,
    t(corrected),
    tolerance = 1e-12
  )

  expect_gte(
    min(eigen(corrected, symmetric = TRUE)$values),
    1e-6 - 1e-10
  )
})


test_that("matrix square-root utilities reject non-square matrices", {
  x <- matrix(1:6, nrow = 2)

  expect_error(
    .sqrt_matrix(x),
    "square matrix"
  )

  expect_error(
    .nearest_psd(x),
    "square matrix"
  )
})
