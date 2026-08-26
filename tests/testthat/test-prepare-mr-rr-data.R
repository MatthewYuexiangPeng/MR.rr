mr_rr_preparation_fixture <- function() {
  instrument_names <- paste0("rs", 1:4)
  exposure_names <- c("X1", "X2")
  outcome_names <- c("Y1", "Y2")

  beta_exposure <- matrix(
    c(
      0.10, 0.04,
      0.12, 0.03,
      0.08, 0.05,
      0.11, 0.02
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(instrument_names, exposure_names)
  )

  se_exposure <- matrix(
    c(
      0.010, 0.020,
      0.012, 0.018,
      0.011, 0.019,
      0.013, 0.017
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(instrument_names, exposure_names)
  )

  beta_outcome <- matrix(
    c(
      0.030, 0.010,
      0.025, 0.012,
      0.035, 0.008,
      0.028, 0.011
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(instrument_names, outcome_names)
  )

  se_outcome <- matrix(
    c(
      0.008, 0.009,
      0.007, 0.010,
      0.009, 0.008,
      0.008, 0.011
    ),
    nrow = 4,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(instrument_names, outcome_names)
  )

  list(
    beta_exposure = beta_exposure,
    se_exposure = se_exposure,
    beta_outcome = beta_outcome,
    se_outcome = se_outcome,
    cor_exposure = matrix(
      c(1, 0.25, 0.25, 1),
      nrow = 2,
      ncol = 2,
      dimnames = list(exposure_names, exposure_names)
    ),
    cor_outcome = matrix(
      c(1, -0.15, -0.15, 1),
      nrow = 2,
      ncol = 2,
      dimnames = list(outcome_names, outcome_names)
    ),
    variant_variance = c(0.32, 0.40, 0.28, 0.36)
  )
}


reference_average_covariance <- function(standard_errors, correlation) {
  dimension <- ncol(standard_errors)
  covariance_sum <- matrix(0, nrow = dimension, ncol = dimension)

  for (j in seq_len(nrow(standard_errors))) {
    D_j <- diag(
      standard_errors[j, ],
      nrow = dimension,
      ncol = dimension
    )
    covariance_sum <- covariance_sum + D_j %*% correlation %*% D_j
  }

  covariance_sum / nrow(standard_errors)
}


test_that("prepare_mr_rr_data reproduces the paper covariance construction", {
  fixture <- mr_rr_preparation_fixture()

  prepared <- prepare_mr_rr_data(
    beta_exposure = fixture$beta_exposure,
    se_exposure = fixture$se_exposure,
    beta_outcome = fixture$beta_outcome,
    se_outcome = fixture$se_outcome,
    cor_exposure = fixture$cor_exposure,
    cor_outcome = fixture$cor_outcome,
    variant_variance = fixture$variant_variance
  )

  row_scale <- sqrt(fixture$variant_variance)
  expected_X <- fixture$beta_exposure * row_scale
  expected_se_X <- fixture$se_exposure * row_scale
  expected_Y <- fixture$beta_outcome * row_scale
  expected_se_Y <- fixture$se_outcome * row_scale

  expected_Sigma_X <- reference_average_covariance(
    expected_se_X,
    fixture$cor_exposure
  )
  expected_Sigma_Y <- reference_average_covariance(
    expected_se_Y,
    fixture$cor_outcome
  )

  dimnames(expected_Sigma_X) <- list(c("X1", "X2"), c("X1", "X2"))
  dimnames(expected_Sigma_Y) <- list(c("Y1", "Y2"), c("Y1", "Y2"))

  expect_s3_class(prepared, "mr_rr_data")
  expect_equal(prepared$X, expected_X, tolerance = 0)
  expect_equal(prepared$se_X, expected_se_X, tolerance = 0)
  expect_equal(prepared$Y, expected_Y, tolerance = 0)
  expect_equal(prepared$se_Y, expected_se_Y, tolerance = 0)
  expect_equal(prepared$Sigma_X, expected_Sigma_X, tolerance = 1e-15)
  expect_equal(prepared$Sigma_Y, expected_Sigma_Y, tolerance = 1e-15)
  expect_equal(
    unname(prepared$W %*% prepared$Sigma_Y),
    diag(2),
    tolerance = 1e-12
  )
  expect_identical(prepared$n_instruments, 4L)
  expect_identical(prepared$n_exposures, 2L)
  expect_identical(prepared$n_outcomes, 2L)
})


test_that("identity trait correlations are used by default", {
  fixture <- mr_rr_preparation_fixture()

  prepared <- prepare_mr_rr_data(
    beta_exposure = fixture$beta_exposure,
    se_exposure = fixture$se_exposure,
    beta_outcome = fixture$beta_outcome,
    se_outcome = fixture$se_outcome
  )

  expected_Sigma_X <- diag(colMeans(fixture$se_exposure^2))
  expected_Sigma_Y <- diag(colMeans(fixture$se_outcome^2))
  dimnames(expected_Sigma_X) <- list(c("X1", "X2"), c("X1", "X2"))
  dimnames(expected_Sigma_Y) <- list(c("Y1", "Y2"), c("Y1", "Y2"))

  expect_equal(prepared$Sigma_X, expected_Sigma_X, tolerance = 1e-15)
  expect_equal(prepared$Sigma_Y, expected_Sigma_Y, tolerance = 1e-15)
})


test_that("a user-supplied outcome weight matrix is retained", {
  fixture <- mr_rr_preparation_fixture()
  custom_W <- matrix(c(2, 0.1, 0.1, 1.5), nrow = 2, ncol = 2)

  prepared <- prepare_mr_rr_data(
    beta_exposure = fixture$beta_exposure,
    se_exposure = fixture$se_exposure,
    beta_outcome = fixture$beta_outcome,
    se_outcome = fixture$se_outcome,
    W = custom_W
  )

  expect_equal(unname(prepared$W), custom_W, tolerance = 0)
})


test_that("prepare_mr_rr_data validates aligned inputs", {
  fixture <- mr_rr_preparation_fixture()

  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = fixture$se_exposure[-1, , drop = FALSE],
      beta_outcome = fixture$beta_outcome,
      se_outcome = fixture$se_outcome
    ),
    "`se_exposure`"
  )

  invalid_se <- fixture$se_exposure
  invalid_se[1, 1] <- 0
  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = invalid_se,
      beta_outcome = fixture$beta_outcome,
      se_outcome = fixture$se_outcome
    ),
    "`se_exposure`"
  )

  misordered_outcomes <- fixture$beta_outcome[c(2, 1, 3, 4), , drop = FALSE]
  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = fixture$se_exposure,
      beta_outcome = misordered_outcomes,
      se_outcome = fixture$se_outcome[c(2, 1, 3, 4), , drop = FALSE]
    ),
    "row names"
  )

  invalid_correlation <- matrix(c(1, 0.2, 0.4, 1), nrow = 2)
  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = fixture$se_exposure,
      beta_outcome = fixture$beta_outcome,
      se_outcome = fixture$se_outcome,
      cor_exposure = invalid_correlation
    ),
    "symmetric"
  )

  misordered_correlation <- fixture$cor_exposure[c(2, 1), c(2, 1)]
  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = fixture$se_exposure,
      beta_outcome = fixture$beta_outcome,
      se_outcome = fixture$se_outcome,
      cor_exposure = misordered_correlation
    ),
    "trait order"
  )

  expect_error(
    prepare_mr_rr_data(
      beta_exposure = fixture$beta_exposure,
      se_exposure = fixture$se_exposure,
      beta_outcome = fixture$beta_outcome,
      se_outcome = fixture$se_outcome,
      variant_variance = c(0.2, 0.3)
    ),
    "`variant_variance`"
  )
})


test_that("one-exposure and one-outcome inputs are supported", {
  instrument_names <- paste0("rs", 1:3)

  prepared <- prepare_mr_rr_data(
    beta_exposure = matrix(
      c(0.10, 0.12, 0.08),
      ncol = 1,
      dimnames = list(instrument_names, "X1")
    ),
    se_exposure = matrix(
      c(0.01, 0.02, 0.015),
      ncol = 1,
      dimnames = list(instrument_names, "X1")
    ),
    beta_outcome = matrix(
      c(0.03, 0.02, 0.04),
      ncol = 1,
      dimnames = list(instrument_names, "Y1")
    ),
    se_outcome = matrix(
      c(0.008, 0.009, 0.01),
      ncol = 1,
      dimnames = list(instrument_names, "Y1")
    )
  )

  expect_equal(unname(prepared$Sigma_X), matrix(mean(c(0.01, 0.02, 0.015)^2)))
  expect_equal(unname(prepared$Sigma_Y), matrix(mean(c(0.008, 0.009, 0.01)^2)))
  expect_identical(dim(prepared$W), c(1L, 1L))
})
