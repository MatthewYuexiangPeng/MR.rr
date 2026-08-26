fit_all_mr_rr_test_data <- function() {
  set.seed(20260829)

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


skip_if_fit_all_solver_unavailable <- function() {
  testthat::skip_if_not_installed("CVXR")
  testthat::skip_if(
    !"OSQP" %in% CVXR::installed_solvers(),
    "The OSQP solver is not available to CVXR."
  )
}


test_that("fit_all_mr_rr returns the requested estimator collection", {
  data <- fit_all_mr_rr_test_data()
  requested_methods <- c("naive", "corrected", "regularized")

  fits <- fit_all_mr_rr(
    data = data,
    methods = requested_methods,
    rank = 1,
    regularization_rate = 1e-6
  )

  expect_s3_class(fits, "mr_rr_fit_collection")
  expect_identical(fits$methods, requested_methods)
  expect_identical(fits$rank, 1L)
  expect_false(fits$rank_selected)
  expect_null(fits$rank_selection)
  expect_identical(fits$data, data)

  for (method in requested_methods) {
    expect_s3_class(fits[[method]], "mr_rr_fit")
    expect_identical(fits[[method]]$method, method)
    expect_identical(fits[[method]]$rank, 1L)
    expect_identical(dim(fits[[method]]$AB), c(2L, 3L))
  }

  direct_corrected <- fit_mr_rr(
    data = data,
    method = "corrected",
    rank = 1
  )
  expect_equal(
    fits$corrected$AB,
    direct_corrected$AB,
    tolerance = 1e-12
  )
})


test_that("fit_all_mr_rr selects a common rank only once", {
  data <- fit_all_mr_rr_test_data()

  expected_selection <- select_mr_rr_rank(
    Y = data$Y,
    X = data$X,
    Sigma_X = data$Sigma_X,
    W = data$W,
    alpha = 0.05,
    min_rank = 1
  )

  fits <- fit_all_mr_rr(
    data = data,
    methods = c("naive", "corrected", "regularized"),
    rank = NULL,
    regularization_rate = 1e-6
  )

  expect_true(fits$rank_selected)
  expect_s3_class(fits$rank_selection, "mr_rr_rank_selection")
  expect_identical(fits$rank, expected_selection$selected_rank)

  for (method in fits$methods) {
    expect_identical(fits[[method]]$rank, fits$rank)
    expect_true(fits[[method]]$rank_selected)
    expect_identical(fits[[method]]$rank_selection, fits$rank_selection)
  }
})


test_that("fit_all_mr_rr can include sparse selection and refitting", {
  skip_if_fit_all_solver_unavailable()
  data <- fit_all_mr_rr_test_data()

  fits <- fit_all_mr_rr(
    data = data,
    methods = c("corrected", "sparse"),
    rank = 1,
    sparse_lambda = 0,
    sparse_refit = TRUE,
    sparse_threshold = 0,
    sparse_max_iter = 2,
    sparse_tol = 1e-10,
    sparse_solver = "OSQP"
  )

  expect_s3_class(fits$sparse, "mr_rr_fit")
  expect_true(fits$sparse$sparse_refitted)
  expect_true(is.list(fits$sparse$sparse_selection))
  expect_identical(
    fits$sparse$details$support,
    fits$sparse$sparse_selection$B != 0
  )
})


test_that("fit_all_mr_rr validates the method collection", {
  data <- fit_all_mr_rr_test_data()

  expect_error(
    fit_all_mr_rr(data, methods = character(0), rank = 1),
    "non-empty"
  )
  expect_error(
    fit_all_mr_rr(data, methods = c("corrected", "unknown"), rank = 1),
    "Unknown method"
  )
  expect_error(
    fit_all_mr_rr(data, methods = c("corrected", "corrected"), rank = 1),
    "duplicates"
  )
})
