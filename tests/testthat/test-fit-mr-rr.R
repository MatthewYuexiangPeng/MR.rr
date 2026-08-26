fit_mr_rr_test_data <- function() {
  set.seed(20260828)

  n <- 60L
  instrument_names <- paste0("rs", seq_len(n))
  exposure_names <- c("X1", "X2", "X3")
  outcome_names <- c("Y1", "Y2")

  X <- matrix(
    stats::rnorm(n * length(exposure_names)),
    nrow = n,
    ncol = length(exposure_names),
    dimnames = list(instrument_names, exposure_names)
  )

  causal_effect <- matrix(
    c(
      0.40, -0.20, 0.10,
      0.15, 0.25, -0.10
    ),
    nrow = length(outcome_names),
    ncol = length(exposure_names),
    byrow = TRUE
  )

  Y <- X %*% t(causal_effect) +
    matrix(stats::rnorm(n * length(outcome_names), sd = 0.10), nrow = n)
  dimnames(Y) <- list(instrument_names, outcome_names)

  se_X <- matrix(
    0.02,
    nrow = n,
    ncol = length(exposure_names),
    dimnames = list(instrument_names, exposure_names)
  )
  se_Y <- matrix(
    0.03,
    nrow = n,
    ncol = length(outcome_names),
    dimnames = list(instrument_names, outcome_names)
  )

  prepare_mr_rr_data(
    beta_exposure = X,
    se_exposure = se_X,
    beta_outcome = Y,
    se_outcome = se_Y,
    W = matrix(c(1.4, 0.2, 0.2, 1.1), nrow = 2, ncol = 2)
  )
}


skip_if_fit_wrapper_solver_unavailable <- function() {
  testthat::skip_if_not_installed("CVXR")
  testthat::skip_if(
    !"OSQP" %in% CVXR::installed_solvers(),
    "The OSQP solver is not available to CVXR."
  )
}


test_that("fit_mr_rr dispatches to the corrected estimator", {
  data <- fit_mr_rr_test_data()

  wrapped <- fit_mr_rr(
    data = data,
    method = "corrected",
    rank = 1
  )
  direct <- mr_rr(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X,
    W = data$W
  )

  expect_s3_class(wrapped, "mr_rr_fit")
  expect_identical(wrapped$method, "corrected")
  expect_identical(wrapped$rank, 1L)
  expect_false(wrapped$rank_selected)
  expect_null(wrapped$rank_selection)
  expect_equal(unname(wrapped$AB), unname(direct$AB), tolerance = 1e-12)
  expect_identical(rownames(wrapped$AB), colnames(data$Y))
  expect_identical(colnames(wrapped$AB), colnames(data$X))
  expect_identical(wrapped$data, data)
})


test_that("fit_mr_rr dispatches to naive and regularized estimators", {
  data <- fit_mr_rr_test_data()

  wrapped_naive <- fit_mr_rr(data, method = "naive", rank = 1)
  direct_naive <- mr_rr_naive(
    Y = data$Y,
    X = data$X,
    r = 1,
    W = data$W
  )

  expect_equal(
    unname(wrapped_naive$AB),
    unname(direct_naive$AB),
    tolerance = 1e-12
  )

  wrapped_regularized <- fit_mr_rr(
    data,
    method = "regularized",
    rank = 1,
    regularization_rate = 1e-6,
    regularized_implementation = "spectral"
  )
  direct_regularized <- mr_rr_regularized(
    Y = data$Y,
    X = data$X,
    r = 1,
    Sigma_X = data$Sigma_X,
    regularization_rate = 1e-6,
    W = data$W,
    implementation = "spectral"
  )

  expect_equal(
    unname(wrapped_regularized$AB),
    unname(direct_regularized$AB),
    tolerance = 1e-12
  )
  expect_identical(
    wrapped_regularized$control$implementation,
    "spectral"
  )
})


test_that("fit_mr_rr selects rank when rank is NULL", {
  data <- fit_mr_rr_test_data()

  expected_selection <- select_mr_rr_rank(
    Y = data$Y,
    X = data$X,
    Sigma_X = data$Sigma_X,
    W = data$W,
    alpha = 0.05,
    min_rank = 1
  )

  wrapped <- fit_mr_rr(
    data,
    method = "corrected",
    rank = NULL,
    rank_alpha = 0.05,
    rank_min = 1
  )

  expect_true(wrapped$rank_selected)
  expect_s3_class(wrapped$rank_selection, "mr_rr_rank_selection")
  expect_identical(wrapped$rank, expected_selection$selected_rank)
  expect_equal(
    wrapped$rank_selection$p_value,
    expected_selection$p_value,
    tolerance = 1e-12
  )
})


test_that("fit_mr_rr performs sparse selection and refitting", {
  skip_if_fit_wrapper_solver_unavailable()
  data <- fit_mr_rr_test_data()

  wrapped <- fit_mr_rr(
    data,
    method = "sparse",
    rank = 1,
    sparse_lambda = 0,
    sparse_refit = TRUE,
    sparse_threshold = 0,
    sparse_max_iter = 2,
    sparse_tol = 1e-10,
    sparse_solver = "OSQP"
  )

  expect_true(wrapped$sparse_refitted)
  expect_true(is.list(wrapped$sparse_selection))
  expect_true(is.matrix(wrapped$details$support))
  expect_identical(
    wrapped$details$support,
    wrapped$sparse_selection$B != 0
  )
  expect_equal(wrapped$AB, wrapped$A %*% wrapped$B, tolerance = 1e-10)
  expect_equal(wrapped$AB, wrapped$details$AB, tolerance = 0)
})


test_that("fit_mr_rr validates its interface", {
  data <- fit_mr_rr_test_data()

  expect_error(
    fit_mr_rr(list(), method = "corrected", rank = 1),
    "mr_rr_data"
  )
  expect_error(
    fit_mr_rr(data, method = "unknown", rank = 1),
    "arg"
  )
  expect_error(
    fit_mr_rr(data, method = "corrected", rank = 0),
    "`rank`"
  )
  expect_error(
    fit_mr_rr(data, method = "corrected", rank = NULL, rank_min = 0),
    "`rank_min`"
  )
  expect_error(
    fit_mr_rr(data, method = "corrected", rank = 1, sparse_refit = NA),
    "`sparse_refit`"
  )
})
