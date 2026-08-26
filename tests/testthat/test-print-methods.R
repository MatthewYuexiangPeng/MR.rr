print_methods_test_data <- function() {
  set.seed(20260830)

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
  C <- matrix(
    c(
      0.40, -0.20, 0.10,
      0.15, 0.25, -0.10
    ),
    nrow = length(outcome_names),
    byrow = TRUE
  )
  Y <- X %*% t(C) +
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


test_that("print.mr_rr_fit displays a compact fit", {
  data <- print_methods_test_data()
  fit <- fit_mr_rr(data, method = "corrected", rank = 1)

  printed <- capture.output(returned <- withVisible(print(fit)))

  expect_false(returned$visible)
  expect_identical(returned$value, fit)
  expect_true(any(grepl("MR-rr fit", printed, fixed = TRUE)))
  expect_true(any(grepl("corrected", printed, fixed = TRUE)))
  expect_true(any(grepl("Rank", printed, fixed = TRUE)))
  expect_true(any(grepl("Estimated causal effect matrix", printed, fixed = TRUE)))
})


test_that("print.mr_rr_fit_collection displays a compact overview", {
  data <- print_methods_test_data()
  methods <- c("naive", "corrected", "regularized")
  fits <- fit_all_mr_rr(
    data = data,
    methods = methods,
    rank = 1,
    regularization_rate = 1e-6
  )

  printed <- capture.output(returned <- withVisible(print(fits)))

  expect_false(returned$visible)
  expect_identical(returned$value, fits)
  expect_true(any(grepl("MR-rr fit collection", printed, fixed = TRUE)))
  expect_true(any(grepl("Estimator overview", printed, fixed = TRUE)))
  expect_true(any(grepl("naive", printed, fixed = TRUE)))
  expect_true(any(grepl("corrected", printed, fixed = TRUE)))
  expect_true(any(grepl("regularized", printed, fixed = TRUE)))
})
