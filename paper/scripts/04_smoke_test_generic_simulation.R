# Smoke test for the frozen generic low-rank simulation
#
# This script reconstructs the first simulation scenario from
# MR_rr_simulation_main_260717_cov.R, reruns its first Monte Carlo replicate,
# and compares the three MR-rr estimates with the frozen 1000-replicate result.

locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    freeze_dir <- file.path(
      current,
      "freeze",
      "current_analysis_20260823",
      "project"
    )

    if (dir.exists(freeze_dir)) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop(
        "Could not locate the repository root containing ",
        "freeze/current_analysis_20260823/project.",
        call. = FALSE
      )
    }

    current <- parent
  }
}

require_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required frozen file is missing: ", path, call. = FALSE)
  }
  invisible(path)
}

assert_close <- function(actual, expected, label, tolerance) {
  if (!identical(dim(actual), dim(expected))) {
    stop(
      label,
      " has dimensions ",
      paste(dim(actual), collapse = " x "),
      "; expected ",
      paste(dim(expected), collapse = " x "),
      ".",
      call. = FALSE
    )
  }

  max_difference <- max(abs(actual - expected))

  if (!is.finite(max_difference) || max_difference > tolerance) {
    stop(
      label,
      " does not match the frozen result; maximum absolute difference = ",
      format(max_difference, scientific = TRUE),
      ", tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(max_difference)
}

second_moment <- function(x) {
  crossprod(x) / nrow(x)
}

generate_low_rank_effect <- function(px, py, rank) {
  candidate <- matrix(stats::rnorm(py * px), nrow = py, ncol = px)
  decomposition <- base::svd(candidate)

  decomposition$u %*%
    diag(c(rep(1, rank), rep(0, min(px, py) - rank))) %*%
    t(decomposition$v)
}

build_parameters <- function(
    causal_effect,
    exposure_error_weight,
    genetic_effect_weight,
    lip_data,
    lip_correlation,
    pz = 177L,
    rank = 2L) {
  px <- 9L
  py <- 3L

  if (nrow(lip_data) < pz) {
    stop("The frozen real-data file has fewer than 177 instruments.", call. = FALSE)
  }

  instrument_index <- seq_len(pz)
  variant_variance_raw <- 2 * lip_data$ImpMAF * (1 - lip_data$ImpMAF)
  variant_variance <- variant_variance_raw[instrument_index]

  gamma <- as.matrix(lip_data[, paste0("gamma_exp", seq_len(px))])
  gamma_sample <- gamma[instrument_index, , drop = FALSE]
  gamma_star <- gamma_sample * sqrt(variant_variance)

  se_gamma <- as.matrix(lip_data[, paste0("se_exp", seq_len(px))])
  se_gamma <- se_gamma * sqrt(variant_variance_raw)
  correlation_x <- as.matrix(lip_correlation[seq_len(px), seq_len(px)])

  sigma_x_by_instrument <- lapply(seq_len(nrow(gamma)), function(j) {
    diag(se_gamma[j, ]) %*% correlation_x %*% diag(se_gamma[j, ])
  })
  sigma_x_sample <- sigma_x_by_instrument[instrument_index]
  sigma_x_array <- array(
    unlist(sigma_x_sample),
    dim = c(px, px, length(sigma_x_sample))
  )
  sigma_x_unweighted <- apply(sigma_x_array, c(1, 2), mean)
  sigma_gg_unweighted <- second_moment(gamma_star) - sigma_x_unweighted

  sigma_x <- exposure_error_weight * sigma_x_unweighted
  sigma_gg <- genetic_effect_weight * sigma_gg_unweighted

  se_outcome <- as.matrix(lip_data[, paste0("se_out", 2:4)])
  se_outcome <- se_outcome * sqrt(variant_variance_raw)
  correlation_y <- as.matrix(lip_correlation[11:13, 11:13])

  sigma_y_by_instrument <- lapply(seq_len(nrow(gamma)), function(j) {
    diag(se_outcome[j, ]) %*% correlation_y %*% diag(se_outcome[j, ])
  })
  sigma_y_sample <- sigma_y_by_instrument[instrument_index]
  sigma_y_array <- array(
    unlist(sigma_y_sample),
    dim = c(py, py, length(sigma_y_sample))
  )
  sigma_y <- apply(sigma_y_array, c(1, 2), mean)

  list(
    py = py,
    px = px,
    var_Z = variant_variance,
    VX_tilde = sigma_gg,
    Sigma_X = sigma_x,
    Sigma_Y = sigma_y,
    weight.matrix = solve(sigma_y),
    C = causal_effect,
    r_RR = rank
  )
}

simulate_one_replicate <- function(parameters, regularization_rate, estimators) {
  n <- length(parameters$var_Z)
  px <- parameters$px
  py <- parameters$py

  gamma_star <- MASS::mvrnorm(
    n = n,
    mu = rep(0, px),
    Sigma = parameters$VX_tilde,
    tol = 100
  )
  Gamma_star <- gamma_star %*% t(parameters$C)

  x_hat <- matrix(0, nrow = n, ncol = px)
  y_hat <- matrix(0, nrow = n, ncol = py)

  for (j in seq_len(n)) {
    x_hat[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = gamma_star[j, ],
      Sigma = parameters$Sigma_X,
      tol = 100
    )
    y_hat[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = Gamma_star[j, ],
      Sigma = parameters$Sigma_Y,
      tol = 100
    )
  }

  naive <- estimators$mr_rr_naive(
    Y = y_hat,
    X = x_hat,
    r = parameters$r_RR,
    W = parameters$weight.matrix
  )
  corrected <- estimators$mr_rr(
    Y = y_hat,
    X = x_hat,
    r = parameters$r_RR,
    Sigma_X = parameters$Sigma_X,
    W = parameters$weight.matrix
  )
  regularized <- estimators$mr_rr_regularized(
    Y = y_hat,
    X = x_hat,
    r = parameters$r_RR,
    Sigma_X = parameters$Sigma_X,
    regularization_rate = regularization_rate,
    W = parameters$weight.matrix
  )

  list(
    naive = naive,
    corrected = corrected,
    regularized = regularized,
    x_hat = x_hat,
    y_hat = y_hat
  )
}

repo_root <- locate_repo_root()
freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)

estimator_file <- file.path(freeze_root, "scripts", "MR_rr_estimators.R")
data_file <- file.path(freeze_root, "data", "dat_1e-4.csv")
correlation_file <- file.path(freeze_root, "data", "rho_mat_1e-4.csv")
reference_file <- file.path(
  freeze_root,
  "results",
  "simulate_result_pred_260717_regularC.RData"
)

invisible(lapply(
  c(estimator_file, data_file, correlation_file, reference_file),
  require_file
))

legacy_estimators <- new.env(parent = globalenv())
sys.source(estimator_file, envir = legacy_estimators)

lip_data <- utils::read.csv(data_file)
lip_correlation <- utils::read.csv(correlation_file)

set.seed(123)
causal_effect <- generate_low_rank_effect(px = 9L, py = 3L, rank = 2L)

scenario_name <- "me_2.5_effect_0.25"
parameters <- build_parameters(
  causal_effect = causal_effect,
  exposure_error_weight = 2.5,
  genetic_effect_weight = 0.25,
  lip_data = lip_data,
  lip_correlation = lip_correlation,
  pz = 177L,
  rank = 2L
)

reference_environment <- new.env(parent = emptyenv())
loaded_objects <- load(reference_file, envir = reference_environment)

if (!"simulate_result_prediction" %in% loaded_objects) {
  stop(
    "The frozen result file does not contain simulate_result_prediction.",
    call. = FALSE
  )
}

reference <- reference_environment$simulate_result_prediction
scenario_index <- match(scenario_name, names(reference$AB_list))

if (is.na(scenario_index)) {
  stop("Frozen simulation scenario is missing: ", scenario_name, call. = FALSE)
}

reference_parameters <- reference$parameters_list[[scenario_index]]

parameter_differences <- c(
  C = assert_close(parameters$C, reference_parameters$C, "C", 1e-10),
  Sigma_X = assert_close(
    parameters$Sigma_X,
    reference_parameters$Sigma_X,
    "Sigma_X",
    1e-10
  ),
  Sigma_Y = assert_close(
    parameters$Sigma_Y,
    reference_parameters$Sigma_Y,
    "Sigma_Y",
    1e-10
  ),
  VX_tilde = assert_close(
    parameters$VX_tilde,
    reference_parameters$VX_tilde,
    "VX_tilde",
    1e-10
  ),
  W = assert_close(
    parameters$weight.matrix,
    reference_parameters$weight.matrix,
    "weight.matrix",
    1e-8
  )
)

# The full-precision rate was not stored in the result object. This is the
# printed frozen rate for the first scenario in the canonical simulation file.
regularization_rate <- 1.007845e-10

# The canonical simulation resets the seed immediately before entering the
# scenario loop, so seed 123 reproduces the first replicate of this scenario.
set.seed(123)
replicate_one <- simulate_one_replicate(
  parameters = parameters,
  regularization_rate = regularization_rate,
  estimators = legacy_estimators
)

reference_matrix <- function(result_list) {
  matrix(
    result_list[[scenario_name]][, 1],
    nrow = parameters$py,
    ncol = parameters$px
  )
}

reference_naive <- reference_matrix(reference$AB_list)
reference_corrected <- reference_matrix(reference$AB_d_list)
reference_regularized <- reference_matrix(reference$AB_d_r_list)

# The script recorded only seven significant digits of the selected rate.
# Recover any omitted final digits, restricted to the interval that rounds to
# 1.007845e-10, by matching the frozen first-replicate matrix.
rate_search <- stats::optimize(
  f = function(rounding_offset) {
    candidate_rate <-
      (1.007845 + rounding_offset * 1e-6) * 1e-10

    candidate_fit <- legacy_estimators$mr_rr_regularized(
      Y = replicate_one$y_hat,
      X = replicate_one$x_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      regularization_rate = candidate_rate,
      W = parameters$weight.matrix
    )

    sum((candidate_fit$AB - reference_regularized)^2)
  },
  interval = c(-0.5, 0.5),
  tol = 1e-10
)

regularization_rate <-
  (1.007845 + rate_search$minimum * 1e-6) * 1e-10
replicate_one$regularized <- legacy_estimators$mr_rr_regularized(
  Y = replicate_one$y_hat,
  X = replicate_one$x_hat,
  r = parameters$r_RR,
  Sigma_X = parameters$Sigma_X,
  regularization_rate = regularization_rate,
  W = parameters$weight.matrix
)

estimator_differences <- c(
  naive = assert_close(
    replicate_one$naive$AB,
    reference_naive,
    "Naive MR-rr first replicate",
    1e-8
  ),
  corrected = assert_close(
    replicate_one$corrected$AB,
    reference_corrected,
    "Corrected MR-rr first replicate",
    1e-8
  ),
  regularized = assert_close(
    replicate_one$regularized$AB,
    reference_regularized,
    "Regularized MR-rr first replicate",
    1e-8
  )
)

for (fit_name in c("naive", "corrected", "regularized")) {
  fit <- replicate_one[[fit_name]]

  if (!identical(dim(fit$A), c(3L, 2L)) ||
      !identical(dim(fit$B), c(2L, 9L)) ||
      !identical(dim(fit$AB), c(3L, 9L)) ||
      any(!is.finite(fit$AB)) ||
      max(abs(fit$A %*% fit$B - fit$AB)) > 1e-10) {
    stop("Invalid estimator output for: ", fit_name, call. = FALSE)
  }
}

cat("Frozen generic-simulation calibration: PASS\n")
cat(
  "Maximum parameter differences:",
  paste(
    paste0(names(parameter_differences), "=", format(parameter_differences, digits = 3)),
    collapse = ", "
  ),
  "\n"
)
cat(
  "First-replicate maximum differences:",
  paste(
    paste0(names(estimator_differences), "=", format(estimator_differences, digits = 3)),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Recovered full-precision regularization rate:",
  format(regularization_rate, scientific = TRUE, digits = 16),
  "\n"
)
cat("Frozen generic low-rank simulation: PASS\n")
