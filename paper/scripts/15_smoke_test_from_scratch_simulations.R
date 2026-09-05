#!/usr/bin/env Rscript

# From-scratch smoke test for the two manuscript simulation designs.
#
# This script intentionally does not read archived simulation-result RData.
# It reconstructs the data-generating mechanisms from the frozen raw inputs,
# generates one replicate per setting, runs every estimator used in Tables 3
# and 7, and writes compact diagnostic output outside the frozen snapshot.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    freeze_root <- file.path(
      current,
      "freeze",
      "current_analysis_20260823",
      "project"
    )

    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(freeze_root) &&
        dir.exists(file.path(current, "paper", "scripts"))) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      stop(
        "Could not locate the MR.rr repository root.",
        call. = FALSE
      )
    }

    current <- parent
  }
}


require_file <- function(path) {
  if (!file.exists(path)) {
    stop(
      "Required input file is missing: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}


parse_boolean <- function(value, variable_name) {
  normalized <- tolower(trimws(value))

  if (normalized %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }

  if (normalized %in% c("0", "false", "no", "n", "off")) {
    return(FALSE)
  }

  stop(
    variable_name,
    " must be one of true/false, yes/no, or 1/0.",
    call. = FALSE
  )
}


parse_integer <- function(value, variable_name, minimum = 0L) {
  parsed <- suppressWarnings(as.integer(value))

  if (length(parsed) != 1L || is.na(parsed) || parsed < minimum) {
    stop(
      variable_name,
      " must be an integer greater than or equal to ",
      minimum,
      ".",
      call. = FALSE
    )
  }

  parsed
}


assert_dimensions <- function(value, expected, label) {
  if (!identical(dim(value), expected)) {
    stop(
      label,
      " has dimensions ",
      paste(dim(value), collapse = " x "),
      "; expected ",
      paste(expected, collapse = " x "),
      ".",
      call. = FALSE
    )
  }

  invisible(value)
}


assert_finite <- function(value, label) {
  if (!all(is.finite(value))) {
    stop(
      label,
      " contains non-finite values.",
      call. = FALSE
    )
  }

  invisible(value)
}


assert_rounded_equal <- function(actual, expected, digits, label) {
  assert_dimensions(actual, dim(expected), label)

  rounded_actual <- round(actual, digits)

  if (!identical(
    as.numeric(rounded_actual),
    as.numeric(expected)
  )) {
    difference <- max(abs(rounded_actual - expected))

    stop(
      label,
      " does not reproduce the manuscript values; maximum rounded ",
      "difference = ",
      format(difference, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


second_moment <- function(value) {
  crossprod(value) / nrow(value)
}


compute_population_siv <- function(parameters, pz = 177L) {
  sigma_x <- (parameters$Sigma_X + t(parameters$Sigma_X)) / 2
  decomposition <- eigen(sigma_x, symmetric = TRUE)

  if (any(!is.finite(decomposition$values)) ||
      any(decomposition$values <= 0)) {
    stop(
      "Sigma_X is not positive definite.",
      call. = FALSE
    )
  }

  sigma_x_inverse_sqrt <-
    decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values)) %*%
    t(decomposition$vectors)

  standardized_strength <-
    sigma_x_inverse_sqrt %*%
    parameters$VX_tilde %*%
    sigma_x_inverse_sqrt
  standardized_strength <-
    (standardized_strength + t(standardized_strength)) / 2

  eigenvalues <- eigen(
    standardized_strength,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  sqrt(pz) * min(eigenvalues)
}


build_calibration <- function(lip_data, lip_correlation, pz = 177L) {
  px <- 9L
  py <- 3L

  if (!identical(nrow(lip_data), pz)) {
    stop(
      "The frozen simulation calibration requires exactly 177 instruments; ",
      "the raw input contains ",
      nrow(lip_data),
      ".",
      call. = FALSE
    )
  }

  required_columns <- c(
    "ImpMAF",
    paste0("gamma_exp", seq_len(px)),
    paste0("se_exp", seq_len(px)),
    paste0("se_out", 2:4)
  )
  missing_columns <- setdiff(required_columns, names(lip_data))

  if (length(missing_columns) > 0L) {
    stop(
      "The raw simulation input is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(lip_correlation) < 13L ||
      ncol(lip_correlation) < 13L) {
    stop(
      "The raw correlation matrix must be at least 13 x 13.",
      call. = FALSE
    )
  }

  variant_variance <-
    2 * lip_data$ImpMAF * (1 - lip_data$ImpMAF)
  sqrt_variant_variance <- sqrt(variant_variance)

  gamma_observed <- as.matrix(
    lip_data[, paste0("gamma_exp", seq_len(px))]
  )
  gamma_star_calibration <-
    gamma_observed * sqrt_variant_variance

  se_exposure <- as.matrix(
    lip_data[, paste0("se_exp", seq_len(px))]
  )
  se_exposure <- se_exposure * sqrt_variant_variance
  correlation_x <- as.matrix(
    lip_correlation[seq_len(px), seq_len(px)]
  )

  sigma_x_by_instrument <- lapply(
    seq_len(pz),
    function(instrument) {
      diagonal <- diag(se_exposure[instrument, ])
      diagonal %*% correlation_x %*% diagonal
    }
  )
  sigma_x_unweighted <-
    Reduce(`+`, sigma_x_by_instrument) / pz
  sigma_gg_unweighted <-
    second_moment(gamma_star_calibration) - sigma_x_unweighted

  se_outcome <- as.matrix(
    lip_data[, paste0("se_out", 2:4)]
  )
  se_outcome <- se_outcome * sqrt_variant_variance
  correlation_y <- as.matrix(
    lip_correlation[11:13, 11:13]
  )

  sigma_y_by_instrument <- lapply(
    seq_len(pz),
    function(instrument) {
      diagonal <- diag(se_outcome[instrument, ])
      diagonal %*% correlation_y %*% diagonal
    }
  )
  sigma_y <- Reduce(`+`, sigma_y_by_instrument) / pz

  values_to_check <- list(
    variant_variance = variant_variance,
    sigma_x_unweighted = sigma_x_unweighted,
    sigma_gg_unweighted = sigma_gg_unweighted,
    sigma_y = sigma_y
  )

  for (value_name in names(values_to_check)) {
    assert_finite(values_to_check[[value_name]], value_name)
  }

  list(
    pz = pz,
    px = px,
    py = py,
    variant_variance = variant_variance,
    sigma_x_unweighted = sigma_x_unweighted,
    sigma_gg_unweighted = sigma_gg_unweighted,
    sigma_y = sigma_y,
    weight_matrix = solve(sigma_y)
  )
}


generate_generic_effect <- function(
    px = 9L,
    py = 3L,
    rank = 2L,
    seed = 123L) {
  set.seed(seed)
  candidate <- matrix(
    stats::rnorm(py * px),
    nrow = py,
    ncol = px
  )
  decomposition <- svd(candidate)

  decomposition$u %*%
    diag(c(
      rep(1, rank),
      rep(0, min(px, py) - rank)
    )) %*%
    t(decomposition$v)
}


generate_sparse_loading_effect <- function(
    calibration,
    estimators,
    rank = 2L,
    seed = 123L) {
  set.seed(seed)

  Q <- qr.Q(
    qr(
      matrix(
        stats::rnorm(calibration$py * rank),
        nrow = calibration$py,
        ncol = rank
      )
    )
  )
  A <- solve(
    estimators$.sqrt_matrix(calibration$weight_matrix)
  ) %*% Q

  B <- matrix(
    0,
    nrow = rank,
    ncol = calibration$px
  )

  for (exposure in seq_len(calibration$px)) {
    nonzero_row <- sample(seq_len(rank), size = 1L)
    B[nonzero_row, exposure] <- stats::rnorm(
      1L,
      mean = 0,
      sd = 25
    )
  }

  list(
    A = A,
    B = B,
    C = A %*% B
  )
}


build_parameters <- function(
    causal_effect,
    exposure_error_weight,
    genetic_effect_weight,
    calibration,
    estimators,
    rank = 2L) {
  sigma_x <-
    exposure_error_weight * calibration$sigma_x_unweighted
  sigma_gg <-
    genetic_effect_weight * calibration$sigma_gg_unweighted
  sigma_y <- calibration$sigma_y
  weight_matrix <- solve(sigma_y)

  sigma_xx <- sigma_x + sigma_gg
  sigma_yx <- causal_effect %*% sigma_gg
  sigma_xy <- t(sigma_yx)
  weight_sqrt <- estimators$.sqrt_matrix(weight_matrix)
  weight_sqrt_inverse <- solve(weight_sqrt)

  naive_target <-
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_xx) %*%
    sigma_xy %*%
    weight_sqrt
  naive_vectors <- eigen(naive_target)$vectors[
    , seq_len(rank), drop = FALSE
  ]
  A <- weight_sqrt_inverse %*% naive_vectors
  B <-
    t(naive_vectors) %*%
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_xx)

  corrected_target <-
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_gg) %*%
    sigma_xy %*%
    weight_sqrt
  corrected_vectors <- eigen(corrected_target)$vectors[
    , seq_len(rank), drop = FALSE
  ]
  A_corrected <- weight_sqrt_inverse %*% corrected_vectors
  B_corrected <-
    t(corrected_vectors) %*%
    weight_sqrt %*%
    sigma_yx %*%
    solve(sigma_gg)

  parameters <- list(
    py = calibration$py,
    px = calibration$px,
    var_Z = calibration$variant_variance,
    VX_tilde = sigma_gg,
    Sigma_X = sigma_x,
    Sigma_Y = sigma_y,
    weight.matrix = weight_matrix,
    SigmaXX = sigma_xx,
    SigmaYX = sigma_yx,
    SigmaXY = sigma_xy,
    C = causal_effect,
    A = A,
    B = B,
    A_d = A_corrected,
    B_d = B_corrected,
    C_r = A_corrected %*% B_corrected,
    r_RR = rank,
    VY_tilde = causal_effect %*% sigma_gg %*% t(causal_effect)
  )

  parameters$iv_strength <- compute_population_siv(
    parameters,
    pz = calibration$pz
  )
  parameters
}


simulate_dataset <- function(parameters) {
  n <- length(parameters$var_Z)

  gamma_star <- MASS::mvrnorm(
    n = n,
    mu = rep(0, parameters$px),
    Sigma = parameters$VX_tilde,
    tol = 100
  )
  GAMMA_star <- gamma_star %*% t(parameters$C)

  gamma_hat <- matrix(
    0,
    nrow = n,
    ncol = parameters$px
  )
  GAMMA_hat <- matrix(
    0,
    nrow = n,
    ncol = parameters$py
  )

  for (instrument in seq_len(n)) {
    gamma_hat[instrument, ] <- MASS::mvrnorm(
      n = 1L,
      mu = gamma_star[instrument, ],
      Sigma = parameters$Sigma_X,
      tol = 100
    )
    GAMMA_hat[instrument, ] <- MASS::mvrnorm(
      n = 1L,
      mu = GAMMA_star[instrument, ],
      Sigma = parameters$Sigma_Y,
      tol = 100
    )
  }

  list(
    gamma_hat = gamma_hat,
    GAMMA_hat = GAMMA_hat
  )
}


fit_main_estimators <- function(
    simulated_data,
    parameters,
    regularization_rate,
    estimators,
    include_mrdag,
    mrdag_niter,
    mrdag_burnin) {
  GAMMA_hat <- simulated_data$GAMMA_hat
  gamma_hat <- simulated_data$gamma_hat

  fits <- list(
    naive_mr_rr = estimators$mr_rr_naive(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      W = parameters$weight.matrix
    )$AB,
    mr_rr = estimators$mr_rr(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      W = parameters$weight.matrix
    )$AB,
    regularized_mr_rr = estimators$mr_rr_regularized(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      regularization_rate = regularization_rate,
      W = parameters$weight.matrix
    )$AB,
    ivw = estimators$ivw_multiple_outcomes(
      Y = GAMMA_hat,
      X = gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    ),
    srivw = estimators$adivw_multiple_outcomes(
      Y = GAMMA_hat,
      X = gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    )
  )

  if (include_mrdag) {
    fits$mrdag <- estimators$Mr_DAG(
      Y = GAMMA_hat,
      X = gamma_hat,
      niter = mrdag_niter,
      burnin = mrdag_burnin
    )
  }

  fits
}


fit_sparse_estimator <- function(
    simulated_data,
    parameters,
    lambda,
    estimators) {
  estimators$mr_rr_sparse(
    GAMMA_hat = simulated_data$GAMMA_hat,
    gamma_hat = simulated_data$gamma_hat,
    W = parameters$weight.matrix,
    Sigma_X = parameters$Sigma_X,
    lambda = lambda,
    r = parameters$r_RR,
    max_iter = 100L,
    tol = 1e-2
  )
}


validate_effect_estimate <- function(
    estimate,
    parameters,
    method,
    require_low_rank = FALSE) {
  assert_dimensions(
    estimate,
    c(parameters$py, parameters$px),
    method
  )
  assert_finite(estimate, method)

  estimated_rank <- qr(estimate, tol = 1e-8)$rank

  if (require_low_rank && estimated_rank > parameters$r_RR) {
    stop(
      method,
      " returned rank ",
      estimated_rank,
      ", above working rank ",
      parameters$r_RR,
      ".",
      call. = FALSE
    )
  }

  invisible(estimated_rank)
}


make_summary_row <- function(
    design,
    scenario,
    phase,
    method,
    estimate,
    parameters,
    fixed_exposure,
    regularization_rate,
    sparse_lambda = NA_real_) {
  true_prediction <- as.vector(parameters$C %*% fixed_exposure)
  estimated_prediction <- as.vector(estimate %*% fixed_exposure)

  data.frame(
    design = design,
    scenario = scenario,
    phase = phase,
    method = method,
    me_weight = as.numeric(sub(
      "^me_([^_]+)_effect_.*$",
      "\\1",
      scenario
    )),
    effect_weight = as.numeric(sub(
      "^.*_effect_([^_]+)$",
      "\\1",
      scenario
    )),
    population_siv = parameters$iv_strength,
    regularization_rate = regularization_rate,
    sparse_lambda = sparse_lambda,
    estimated_rank = qr(estimate, tol = 1e-8)$rank,
    causal_rmse = sqrt(mean((estimate - parameters$C)^2)),
    prediction_rmse = sqrt(mean(
      (estimated_prediction - true_prediction)^2
    )),
    maximum_absolute_estimate = max(abs(estimate)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


run_design_smoke <- function(
    design_name,
    causal_effect,
    settings,
    regularization_rates,
    sparse_lambdas,
    calibration,
    estimators,
    fixed_exposure,
    include_mrdag,
    mrdag_niter,
    mrdag_burnin,
    seed = 123L) {
  parameters_by_setting <- setNames(
    vector("list", nrow(settings)),
    settings$scenario
  )

  for (setting_index in seq_len(nrow(settings))) {
    setting <- settings[setting_index, , drop = FALSE]
    parameters_by_setting[[setting$scenario]] <- build_parameters(
      causal_effect = causal_effect,
      exposure_error_weight = setting$me_weight,
      genetic_effect_weight = setting$effect_weight,
      calibration = calibration,
      estimators = estimators,
      rank = 2L
    )
  }

  main_results <- setNames(
    vector("list", nrow(settings)),
    settings$scenario
  )
  summary_rows <- list()
  summary_index <- 0L

  # The legacy main simulation resets the seed before entering its setting
  # loop. The full 1,000-replicate driver will use the same loop order.
  set.seed(seed)

  for (setting_index in seq_len(nrow(settings))) {
    setting <- settings[setting_index, , drop = FALSE]
    scenario <- setting$scenario
    parameters <- parameters_by_setting[[scenario]]

    cat(
      "  Main estimators:",
      design_name,
      scenario,
      "\n"
    )

    simulated_data <- simulate_dataset(parameters)
    fits <- fit_main_estimators(
      simulated_data = simulated_data,
      parameters = parameters,
      regularization_rate = regularization_rates[[scenario]],
      estimators = estimators,
      include_mrdag = include_mrdag,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin
    )

    for (method in names(fits)) {
      require_low_rank <- method %in% c(
        "naive_mr_rr",
        "mr_rr",
        "regularized_mr_rr"
      )
      validate_effect_estimate(
        fits[[method]],
        parameters,
        paste(design_name, scenario, method),
        require_low_rank = require_low_rank
      )

      summary_index <- summary_index + 1L
      summary_rows[[summary_index]] <- make_summary_row(
        design = design_name,
        scenario = scenario,
        phase = "main",
        method = method,
        estimate = fits[[method]],
        parameters = parameters,
        fixed_exposure = fixed_exposure,
        regularization_rate = regularization_rates[[scenario]]
      )
    }

    main_results[[scenario]] <- list(
      parameters = parameters,
      simulated_data = simulated_data,
      estimates = fits
    )
  }

  sparse_results <- setNames(
    vector("list", nrow(settings)),
    settings$scenario
  )

  # Sparse MR-rr was run as a separate simulation in the legacy workflow and
  # therefore has its own seed reset.
  set.seed(seed)

  for (setting_index in seq_len(nrow(settings))) {
    setting <- settings[setting_index, , drop = FALSE]
    scenario <- setting$scenario
    parameters <- parameters_by_setting[[scenario]]

    cat(
      "  Sparse estimator:",
      design_name,
      scenario,
      "\n"
    )

    simulated_data <- simulate_dataset(parameters)
    sparse_fit <- fit_sparse_estimator(
      simulated_data = simulated_data,
      parameters = parameters,
      lambda = rep(
        sparse_lambdas[[scenario]],
        parameters$px
      ),
      estimators = estimators
    )

    validate_effect_estimate(
      sparse_fit$AB,
      parameters,
      paste(design_name, scenario, "sparse_mr_rr"),
      require_low_rank = TRUE
    )
    assert_dimensions(
      sparse_fit$A,
      c(parameters$py, parameters$r_RR),
      paste(design_name, scenario, "sparse A")
    )
    assert_dimensions(
      sparse_fit$B,
      c(parameters$r_RR, parameters$px),
      paste(design_name, scenario, "sparse B")
    )

    reconstruction_difference <- max(
      abs(sparse_fit$AB - sparse_fit$A %*% sparse_fit$B)
    )

    if (!is.finite(reconstruction_difference) ||
        reconstruction_difference > 1e-8) {
      stop(
        design_name,
        " ",
        scenario,
        " sparse MR-rr failed AB = A %*% B.",
        call. = FALSE
      )
    }

    summary_index <- summary_index + 1L
    summary_rows[[summary_index]] <- make_summary_row(
      design = design_name,
      scenario = scenario,
      phase = "sparse",
      method = "sparse_mr_rr",
      estimate = sparse_fit$AB,
      parameters = parameters,
      fixed_exposure = fixed_exposure,
      regularization_rate = regularization_rates[[scenario]],
      sparse_lambda = sparse_lambdas[[scenario]]
    )

    sparse_results[[scenario]] <- list(
      parameters = parameters,
      simulated_data = simulated_data,
      estimate = sparse_fit
    )
  }

  list(
    parameters = parameters_by_setting,
    main = main_results,
    sparse = sparse_results,
    summary = do.call(rbind, summary_rows)
  )
}


run_from_scratch_simulation_smoke <- function() {
  repo_root <- locate_repo_root()
  freeze_root <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project"
  )
  estimator_file <- file.path(
    freeze_root,
    "scripts",
    "MR_rr_estimators.R"
  )
  data_file <- file.path(
    freeze_root,
    "data",
    "dat_1e-4.csv"
  )
  correlation_file <- file.path(
    freeze_root,
    "data",
    "rho_mat_1e-4.csv"
  )

  required_files <- c(
    estimator_file,
    data_file,
    correlation_file
  )
  invisible(lapply(required_files, require_file))

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package `MASS` is required.",
      call. = FALSE
    )
  }

  include_mrdag <- parse_boolean(
    Sys.getenv(
      "MRRR_SMOKE_INCLUDE_MRDAG",
      unset = "true"
    ),
    "MRRR_SMOKE_INCLUDE_MRDAG"
  )
  mrdag_niter <- parse_integer(
    Sys.getenv(
      "MRRR_SMOKE_MRDAG_NITER",
      unset = "1000"
    ),
    "MRRR_SMOKE_MRDAG_NITER",
    minimum = 1L
  )
  mrdag_burnin <- parse_integer(
    Sys.getenv(
      "MRRR_SMOKE_MRDAG_BURNIN",
      unset = "200"
    ),
    "MRRR_SMOKE_MRDAG_BURNIN",
    minimum = 0L
  )

  if (include_mrdag && mrdag_burnin >= mrdag_niter) {
    stop(
      "MRRR_SMOKE_MRDAG_BURNIN must be smaller than ",
      "MRRR_SMOKE_MRDAG_NITER.",
      call. = FALSE
    )
  }

  output_root_setting <- trimws(Sys.getenv(
    "MRRR_FULL_RUN_OUTPUT_DIR",
    unset = ""
  ))

  if (nzchar(output_root_setting)) {
    output_root <- normalizePath(
      output_root_setting,
      winslash = "/",
      mustWork = FALSE
    )
  } else {
    output_root <- file.path(
      repo_root,
      "paper",
      "output",
      "full_run",
      "smoke"
    )
  }
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(
    output_root,
    winslash = "/",
    mustWork = TRUE
  )

  cat("MR-rr from-scratch simulation smoke test\n")
  cat("Repository root:", repo_root, "\n")
  cat("Raw-data input:", data_file, "\n")
  cat("Estimator implementation:", estimator_file, "\n")
  cat("Archived result RData read: no\n")
  cat("Output directory:", output_root, "\n")
  cat("MrDAG included:", include_mrdag, "\n")

  estimators <- new.env(parent = globalenv())
  sys.source(estimator_file, envir = estimators)

  lip_data <- utils::read.csv(
    data_file,
    check.names = FALSE
  )
  lip_correlation <- utils::read.csv(
    correlation_file,
    check.names = FALSE
  )
  calibration <- build_calibration(
    lip_data = lip_data,
    lip_correlation = lip_correlation,
    pz = 177L
  )

  generic_effect <- generate_generic_effect(
    px = calibration$px,
    py = calibration$py,
    rank = 2L,
    seed = 123L
  )
  sparse_truth <- generate_sparse_loading_effect(
    calibration = calibration,
    estimators = estimators,
    rank = 2L,
    seed = 123L
  )

  set.seed(123L)
  fixed_exposure <- stats::rnorm(
    calibration$px,
    mean = 0,
    sd = 1
  )

  generic_effect_target <- matrix(
    c(
      -0.204, -0.168,  0.053, -0.020,  0.086,
       0.333,  0.137,  0.024, -0.187,
      -0.060,  0.095, -0.452,  0.429,  0.039,
       0.189, -0.172, -0.393, -0.601,
       0.354,  0.329, -0.215,  0.149, -0.146,
      -0.553, -0.293, -0.147,  0.181
    ),
    nrow = 3L,
    ncol = 9L,
    byrow = TRUE
  )
  sparse_loading_target <- matrix(
    c(
       0.000, -17.171, 42.255, 0.000, 13.727,
       0.000,  32.369, 0.000, -19.610,
       4.577,   0.000,  0.000, 8.995,  0.000,
     -13.896,   0.000, -49.165, 0.000
    ),
    nrow = 2L,
    ncol = 9L,
    byrow = TRUE
  )
  sparse_effect_target <- matrix(
    c(
      -0.063,  0.082, -0.202, -0.123, -0.065,
       0.191, -0.154,  0.675,  0.094,
      -0.030,  0.022, -0.054, -0.059, -0.018,
       0.091, -0.042,  0.323,  0.025,
      -0.029, -0.231,  0.569, -0.057,  0.185,
       0.088,  0.436,  0.310, -0.264
    ),
    nrow = 3L,
    ncol = 9L,
    byrow = TRUE
  )
  fixed_exposure_target <- c(
    -0.560, -0.230, 1.559,
     0.071,  0.129, 1.715,
     0.461, -1.265, -0.687
  )

  assert_rounded_equal(
    generic_effect,
    generic_effect_target,
    3L,
    "Generic true effect (Table 1)"
  )
  assert_rounded_equal(
    sparse_truth$B,
    sparse_loading_target,
    3L,
    "Sparse true loading (Table 4)"
  )
  assert_rounded_equal(
    sparse_truth$C,
    sparse_effect_target,
    3L,
    "Sparse true effect (Table 5)"
  )

  if (!identical(round(fixed_exposure, 3L), fixed_exposure_target)) {
    stop(
      "The fixed exposure vector does not match the manuscript simulation.",
      call. = FALSE
    )
  }

  orthonormality_error <- max(abs(
    t(sparse_truth$A) %*%
    calibration$weight_matrix %*%
    sparse_truth$A - diag(2L)
  ))

  if (!is.finite(orthonormality_error) ||
      orthonormality_error > 1e-8) {
    stop(
      "The sparse true A does not satisfy A^T W A = I.",
      call. = FALSE
    )
  }

  if (!identical(
    as.integer(colSums(sparse_truth$B != 0)),
    rep(1L, calibration$px)
  )) {
    stop(
      "The sparse true B must have one nonzero loading per exposure.",
      call. = FALSE
    )
  }

  settings <- data.frame(
    scenario = c(
      "me_2.5_effect_0.25",
      "me_2.5_effect_1",
      "me_1_effect_0.25",
      "me_1_effect_1"
    ),
    me_weight = c(2.5, 2.5, 1, 1),
    effect_weight = c(0.25, 1, 0.25, 1),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  generic_regularization_rates <- c(
    me_2.5_effect_0.25 = 1.007845e-10,
    me_2.5_effect_1 = 2.340371e-12,
    me_1_effect_0.25 = 1.578970e-12,
    me_1_effect_1 = 3.103420e-15
  )
  sparse_regularization_rates <- c(
    me_2.5_effect_0.25 = 9.109067e-11,
    me_2.5_effect_1 = 2.340371e-12,
    me_1_effect_0.25 = 1.578970e-12,
    me_1_effect_1 = 3.105502e-15
  )
  sparse_lambdas <- c(
    me_2.5_effect_0.25 = 1e-4,
    me_2.5_effect_1 = 1e-3,
    me_1_effect_0.25 = 1e-3,
    me_1_effect_1 = 1e-3
  )

  calibration_parameters <- lapply(
    seq_len(nrow(settings)),
    function(setting_index) {
      setting <- settings[setting_index, , drop = FALSE]
      build_parameters(
        causal_effect = generic_effect,
        exposure_error_weight = setting$me_weight,
        genetic_effect_weight = setting$effect_weight,
        calibration = calibration,
        estimators = estimators,
        rank = 2L
      )
    }
  )
  internal_siv <- vapply(
    calibration_parameters,
    function(parameters) parameters$iv_strength,
    numeric(1)
  )
  internal_siv_target <- c(3.60, 14.41, 9.00, 36.02)

  if (!identical(round(internal_siv, 2L), internal_siv_target)) {
    stop(
      "The reconstructed scaled instrument strengths do not match ",
      "the manuscript settings.",
      call. = FALSE
    )
  }

  cat("Raw-input calibration and true matrices: PASS\n")
  cat("Running generic low-rank design.\n")
  generic_result <- run_design_smoke(
    design_name = "generic_low_rank",
    causal_effect = generic_effect,
    settings = settings,
    regularization_rates = generic_regularization_rates,
    sparse_lambdas = sparse_lambdas,
    calibration = calibration,
    estimators = estimators,
    fixed_exposure = fixed_exposure,
    include_mrdag = include_mrdag,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin,
    seed = 123L
  )

  cat("Running sparse-loading design.\n")
  sparse_result <- run_design_smoke(
    design_name = "sparse_loading",
    causal_effect = sparse_truth$C,
    settings = settings,
    regularization_rates = sparse_regularization_rates,
    sparse_lambdas = sparse_lambdas,
    calibration = calibration,
    estimators = estimators,
    fixed_exposure = fixed_exposure,
    include_mrdag = include_mrdag,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin,
    seed = 123L
  )

  summary_table <- rbind(
    generic_result$summary,
    sparse_result$summary
  )
  rownames(summary_table) <- NULL

  expected_rows <- if (include_mrdag) 56L else 48L

  if (!identical(nrow(summary_table), expected_rows) ||
      any(!is.finite(summary_table$causal_rmse)) ||
      any(!is.finite(summary_table$prediction_rmse))) {
    stop(
      "The estimator summary has an invalid structure or values.",
      call. = FALSE
    )
  }

  input_manifest <- data.frame(
    role = c(
      "raw simulation data",
      "raw correlation matrix",
      "frozen legacy estimator"
    ),
    repository_path = c(
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "data",
        "dat_1e-4.csv"
      ),
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "data",
        "rho_mat_1e-4.csv"
      ),
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "scripts",
        "MR_rr_estimators.R"
      )
    ),
    md5 = unname(tools::md5sum(required_files)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  full_run_smoke <- list(
    metadata = list(
      schema_version = "1.0.0",
      generated_at_utc = format(
        Sys.time(),
        tz = "UTC",
        usetz = TRUE
      ),
      seed = 123L,
      replicates_per_setting = 1L,
      reads_archived_result_rdata = FALSE,
      mrdag_included = include_mrdag,
      mrdag_niter = if (include_mrdag) mrdag_niter else NA_integer_,
      mrdag_burnin = if (include_mrdag) mrdag_burnin else NA_integer_,
      note = paste(
        "Smoke mode uses one replicate per setting.",
        "A 1,000-replicate driver must retain the same setting-major loop",
        "and independent seed resets for main and sparse simulations."
      )
    ),
    input_manifest = input_manifest,
    settings = settings,
    fixed_exposure = fixed_exposure,
    truth = list(
      generic_C = generic_effect,
      sparse_A = sparse_truth$A,
      sparse_B = sparse_truth$B,
      sparse_C = sparse_truth$C
    ),
    generic = generic_result,
    sparse_loading = sparse_result,
    summary = summary_table
  )

  summary_file <- file.path(
    output_root,
    "simulation_smoke_summary.csv"
  )
  settings_file <- file.path(
    output_root,
    "simulation_settings.csv"
  )
  manifest_file <- file.path(
    output_root,
    "input_manifest.csv"
  )
  result_file <- file.path(
    output_root,
    "simulation_smoke_results.RData"
  )
  session_file <- file.path(
    output_root,
    "sessionInfo.txt"
  )

  utils::write.csv(
    summary_table,
    summary_file,
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    settings,
    settings_file,
    row.names = FALSE
  )
  utils::write.csv(
    input_manifest,
    manifest_file,
    row.names = FALSE
  )
  save(full_run_smoke, file = result_file)
  writeLines(
    capture.output(utils::sessionInfo()),
    session_file,
    useBytes = TRUE
  )

  output_files <- c(
    summary_file,
    settings_file,
    manifest_file,
    result_file,
    session_file
  )

  if (any(!file.exists(output_files))) {
    stop(
      "One or more smoke-test outputs were not written.",
      call. = FALSE
    )
  }

  cat("All requested estimator paths: PASS\n")
  cat("Estimator result rows:", nrow(summary_table), "\n")
  cat("Archived result RData independence: PASS\n")
  cat("Output files:", length(output_files), "\n")
  cat("MR-rr from-scratch simulation smoke test: PASS\n")

  invisible(full_run_smoke)
}


if (sys.nframe() == 0L) {
  run_from_scratch_simulation_smoke()
}
