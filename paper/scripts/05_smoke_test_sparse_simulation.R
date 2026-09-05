# Smoke test for the frozen sparse-loading simulation
#
# This script reconstructs the sparse rank-two causal effect used in
# Simulation_sparse_260719.R, reruns the first Monte Carlo replicate, and
# compares the four MR-rr estimates and the aligned sparse loading matrix with
# the frozen 1000-replicate result.

options(warn = 1)


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
    stop(
      "Required frozen file is missing: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}


assert_close <- function(actual, expected, label, tolerance) {
  if (!identical(dim(actual), dim(expected))) {
    actual_dimension <- if (is.null(dim(actual))) {
      paste0("length ", length(actual))
    } else {
      paste(dim(actual), collapse = " x ")
    }

    expected_dimension <- if (is.null(dim(expected))) {
      paste0("length ", length(expected))
    } else {
      paste(dim(expected), collapse = " x ")
    }

    stop(
      label,
      " has ",
      actual_dimension,
      "; expected ",
      expected_dimension,
      ".",
      call. = FALSE
    )
  }

  if (length(actual) != length(expected)) {
    stop(
      label,
      " has length ",
      length(actual),
      "; expected length ",
      length(expected),
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

  if (!identical(as.integer(pz), 177L)) {
    stop(
      "The frozen sparse simulation uses exactly 177 instruments.",
      call. = FALSE
    )
  }

  if (nrow(lip_data) < pz) {
    stop(
      "The frozen real-data file has fewer than 177 instruments.",
      call. = FALSE
    )
  }

  instrument_index <- seq_len(177L)
  variant_variance_raw <-
    2 * lip_data$ImpMAF * (1 - lip_data$ImpMAF)
  variant_variance <- variant_variance_raw[instrument_index]

  gamma <- as.matrix(
    lip_data[, paste0("gamma_exp", seq_len(px))]
  )
  gamma_sample <- gamma[instrument_index, , drop = FALSE]
  gamma_star <- gamma_sample * sqrt(variant_variance)

  se_gamma <- as.matrix(
    lip_data[, paste0("se_exp", seq_len(px))]
  )
  se_gamma <- se_gamma * sqrt(variant_variance_raw)
  correlation_x <- as.matrix(
    lip_correlation[seq_len(px), seq_len(px)]
  )

  sigma_x_by_instrument <- lapply(
    seq_len(nrow(gamma)),
    function(j) {
      diag(se_gamma[j, ]) %*%
        correlation_x %*%
        diag(se_gamma[j, ])
    }
  )
  sigma_x_sample <- sigma_x_by_instrument[instrument_index]
  sigma_x_array <- array(
    unlist(sigma_x_sample),
    dim = c(px, px, length(sigma_x_sample))
  )
  sigma_x_unweighted <- apply(
    sigma_x_array,
    c(1, 2),
    mean
  )
  sigma_gg_unweighted <-
    second_moment(gamma_star) - sigma_x_unweighted

  sigma_x <- exposure_error_weight * sigma_x_unweighted
  sigma_gg <- genetic_effect_weight * sigma_gg_unweighted

  se_outcome <- as.matrix(
    lip_data[, paste0("se_out", 2:4)]
  )
  se_outcome <- se_outcome * sqrt(variant_variance_raw)
  correlation_y <- as.matrix(
    lip_correlation[11:13, 11:13]
  )

  sigma_y_by_instrument <- lapply(
    seq_len(nrow(gamma)),
    function(j) {
      diag(se_outcome[j, ]) %*%
        correlation_y %*%
        diag(se_outcome[j, ])
    }
  )
  sigma_y_sample <- sigma_y_by_instrument[instrument_index]
  sigma_y_array <- array(
    unlist(sigma_y_sample),
    dim = c(py, py, length(sigma_y_sample))
  )
  sigma_y <- apply(
    sigma_y_array,
    c(1, 2),
    mean
  )

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


build_sparse_effect <- function(
    lip_data,
    lip_correlation,
    estimators,
    seed = 123L) {
  px <- 9L
  py <- 3L
  rank <- 2L

  variant_variance <-
    2 * lip_data$ImpMAF * (1 - lip_data$ImpMAF)

  se_outcome <- as.matrix(
    lip_data[, paste0("se_out", 2:(py + 1L))]
  )
  se_outcome <- se_outcome * sqrt(variant_variance)
  correlation_y <- as.matrix(
    lip_correlation[11:13, 11:13]
  )

  sigma_y_by_instrument <- lapply(
    seq_len(nrow(lip_data)),
    function(j) {
      diag(se_outcome[j, ]) %*%
        correlation_y %*%
        diag(se_outcome[j, ])
    }
  )
  sigma_y_array <- array(
    unlist(sigma_y_by_instrument),
    dim = c(py, py, length(sigma_y_by_instrument))
  )
  sigma_y <- apply(
    sigma_y_array,
    c(1, 2),
    mean
  )
  weight_matrix <- solve(sigma_y)

  # Reproduce the exact random-number order in Simulation_sparse_260719.R.
  set.seed(seed)
  Q <- qr.Q(
    qr(
      matrix(
        stats::rnorm(py * rank),
        nrow = py,
        ncol = rank
      )
    )
  )
  A <- solve(estimators$.sqrt_matrix(weight_matrix)) %*% Q

  B <- matrix(0, nrow = rank, ncol = px)

  for (j in seq_len(px)) {
    nonzero_row <- sample(seq_len(rank), size = 1L)
    B[nonzero_row, j] <- stats::rnorm(1L, mean = 0, sd = 25)
  }

  list(
    A = A,
    B = B,
    C = A %*% B,
    W = weight_matrix
  )
}


align_AB_to_reference <- function(A, B, B_reference) {
  if (!identical(dim(B), dim(B_reference))) {
    stop(
      "Sparse B and its reference have different dimensions.",
      call. = FALSE
    )
  }

  rank <- nrow(B)

  if (!identical(rank, 2L)) {
    stop(
      "The frozen sparse simulation alignment is defined here for rank two.",
      call. = FALSE
    )
  }

  permutations <- rbind(
    c(1L, 2L),
    c(2L, 1L)
  )
  sign_grid <- as.matrix(
    expand.grid(
      rep(list(c(-1, 1)), rank)
    )
  )

  best_error <- Inf
  best_A <- A
  best_B <- B

  for (i in seq_len(nrow(permutations))) {
    permutation <- permutations[i, ]
    B_permuted <- B[permutation, , drop = FALSE]
    A_permuted <- A[, permutation, drop = FALSE]

    for (j in seq_len(nrow(sign_grid))) {
      signs <- sign_grid[j, ]
      sign_matrix <- diag(as.numeric(signs), rank, rank)

      B_candidate <- sign_matrix %*% B_permuted
      A_candidate <- A_permuted %*% sign_matrix
      error <- sum((B_candidate - B_reference)^2)

      if (error < best_error) {
        best_error <- error
        best_A <- A_candidate
        best_B <- B_candidate
      }
    }
  }

  list(
    A = best_A,
    B = best_B,
    AB = best_A %*% best_B
  )
}


simulate_one_sparse_replicate <- function(
    parameters,
    B_true,
    lambda,
    regularization_rate,
    estimators) {
  n <- length(parameters$var_Z)
  px <- parameters$px
  py <- parameters$py

  gamma_star <- MASS::mvrnorm(
    n = n,
    mu = rep(0, px),
    Sigma = parameters$VX_tilde,
    tol = 100
  )
  GAMMA_star <- gamma_star %*% t(parameters$C)

  gamma_hat <- matrix(0, nrow = n, ncol = px)
  GAMMA_hat <- matrix(0, nrow = n, ncol = py)

  for (j in seq_len(n)) {
    gamma_hat[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = gamma_star[j, ],
      Sigma = parameters$Sigma_X,
      tol = 100
    )
    GAMMA_hat[j, ] <- MASS::mvrnorm(
      n = 1,
      mu = GAMMA_star[j, ],
      Sigma = parameters$Sigma_Y,
      tol = 100
    )
  }

  naive <- estimators$mr_rr_naive(
    Y = GAMMA_hat,
    X = gamma_hat,
    r = parameters$r_RR,
    W = parameters$weight.matrix
  )
  standard <- estimators$mr_rr(
    Y = GAMMA_hat,
    X = gamma_hat,
    r = parameters$r_RR,
    W = parameters$weight.matrix,
    Sigma_X = parameters$Sigma_X
  )
  regularized <- estimators$mr_rr_regularized(
    Y = GAMMA_hat,
    X = gamma_hat,
    r = parameters$r_RR,
    W = parameters$weight.matrix,
    Sigma_X = parameters$Sigma_X,
    regularization_rate = regularization_rate
  )
  sparse_unaligned <- estimators$mr_rr_sparse(
    GAMMA_hat = GAMMA_hat,
    gamma_hat = gamma_hat,
    W = parameters$weight.matrix,
    Sigma_X = parameters$Sigma_X,
    lambda = lambda,
    r = parameters$r_RR,
    max_iter = 100L
  )
  sparse <- align_AB_to_reference(
    A = sparse_unaligned$A,
    B = sparse_unaligned$B,
    B_reference = B_true
  )

  list(
    naive = naive,
    standard = standard,
    regularized = regularized,
    sparse = sparse,
    gamma_hat = gamma_hat,
    GAMMA_hat = GAMMA_hat
  )
}


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
reference_file <- file.path(
  freeze_root,
  "results",
  "simulation_sparseB_me1_effect1_eta_1e-3_260719.RData"
)

invisible(
  lapply(
    c(
      estimator_file,
      data_file,
      correlation_file,
      reference_file
    ),
    require_file
  )
)

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package `MASS` is required for this smoke test.", call. = FALSE)
}

legacy_estimators <- new.env(parent = globalenv())
sys.source(estimator_file, envir = legacy_estimators)

lip_data <- utils::read.csv(data_file)
lip_correlation <- utils::read.csv(correlation_file)

reference_environment <- new.env(parent = emptyenv())
loaded_objects <- load(reference_file, envir = reference_environment)
required_reference_objects <- c(
  "res",
  "C",
  "A",
  "B_sparse",
  "me_weight_sparse",
  "effect_weight_sparse",
  "sparse_lambda",
  "regularization_rate_sparse",
  "support_threshold"
)
missing_reference_objects <- setdiff(
  required_reference_objects,
  loaded_objects
)

if (length(missing_reference_objects) > 0L) {
  stop(
    "The frozen sparse result is missing objects: ",
    paste(missing_reference_objects, collapse = ", "),
    call. = FALSE
  )
}

reference <- reference_environment$res
reference_A <- reference_environment$A
reference_B <- reference_environment$B_sparse
reference_C <- reference_environment$C

required_result_fields <- c(
  "me_weight",
  "effect_weight",
  "lambda",
  "regularization_rate",
  "support_threshold",
  "C_true",
  "B_sparse",
  "result_AB_naive",
  "result_AB_standard",
  "result_AB_regularized",
  "result_AB_sparse",
  "result_B_sparse"
)
missing_result_fields <- setdiff(
  required_result_fields,
  names(reference)
)

if (length(missing_result_fields) > 0L) {
  stop(
    "The frozen `res` object is missing fields: ",
    paste(missing_result_fields, collapse = ", "),
    call. = FALSE
  )
}

if (!identical(dim(reference_A), c(3L, 2L)) ||
    !identical(dim(reference_B), c(2L, 9L)) ||
    !identical(dim(reference_C), c(3L, 9L))) {
  stop(
    "The frozen sparse effect has unexpected dimensions.",
    call. = FALSE
  )
}

if (any(vapply(
  reference[c(
    "result_AB_naive",
    "result_AB_standard",
    "result_AB_regularized",
    "result_AB_sparse",
    "result_B_sparse"
  )],
  ncol,
  integer(1)
) != 1000L)) {
  stop(
    "The frozen sparse result does not contain 1,000 replicates.",
    call. = FALSE
  )
}

sparse_effect <- build_sparse_effect(
  lip_data = lip_data,
  lip_correlation = lip_correlation,
  estimators = legacy_estimators,
  seed = 123L
)

effect_differences <- c(
  A = assert_close(
    sparse_effect$A,
    reference_A,
    "Sparse-effect A",
    1e-10
  ),
  B = assert_close(
    sparse_effect$B,
    reference_B,
    "Sparse-effect B",
    1e-10
  ),
  C = assert_close(
    sparse_effect$C,
    reference_C,
    "Sparse-effect C",
    1e-10
  ),
  C_in_res = assert_close(
    sparse_effect$C,
    reference$C_true,
    "Sparse-effect C stored in `res`",
    1e-10
  ),
  B_in_res = assert_close(
    sparse_effect$B,
    reference$B_sparse,
    "Sparse-effect B stored in `res`",
    1e-10
  )
)

orthonormality_error <- max(
  abs(
    crossprod(
      sparse_effect$A,
      sparse_effect$W %*% sparse_effect$A
    ) - diag(2L)
  )
)

if (!is.finite(orthonormality_error) ||
    orthonormality_error > 1e-8) {
  stop(
    "The regenerated A does not satisfy A^T W A = I.",
    call. = FALSE
  )
}

if (!identical(
  as.integer(colSums(sparse_effect$B != 0)),
  rep(1L, 9L)
)) {
  stop(
    "The regenerated true B does not have one nonzero per exposure.",
    call. = FALSE
  )
}

settings_differences <- c(
  me_weight = assert_close(
    reference_environment$me_weight_sparse,
    reference$me_weight,
    "Measurement-error weight",
    0
  ),
  effect_weight = assert_close(
    reference_environment$effect_weight_sparse,
    reference$effect_weight,
    "Genetic-effect weight",
    0
  ),
  lambda = assert_close(
    reference_environment$sparse_lambda,
    reference$lambda,
    "Sparse lambda",
    0
  ),
  regularization_rate = assert_close(
    reference_environment$regularization_rate_sparse,
    reference$regularization_rate,
    "Regularization rate",
    0
  ),
  support_threshold = assert_close(
    reference_environment$support_threshold,
    reference$support_threshold,
    "Support threshold",
    0
  )
)

parameters <- build_parameters(
  causal_effect = sparse_effect$C,
  exposure_error_weight = reference$me_weight,
  genetic_effect_weight = reference$effect_weight,
  lip_data = lip_data,
  lip_correlation = lip_correlation,
  pz = 177L,
  rank = 2L
)

if (max(abs(parameters$weight.matrix - sparse_effect$W)) > 1e-8) {
  stop(
    "The two regenerated outcome-weight matrices do not agree.",
    call. = FALSE
  )
}

# The canonical sparse script resets the seed immediately before the full
# Monte Carlo loop. Seed 123 therefore reproduces its first replicate.
set.seed(123)
replicate_one <- simulate_one_sparse_replicate(
  parameters = parameters,
  B_true = sparse_effect$B,
  lambda = reference$lambda,
  regularization_rate = reference$regularization_rate,
  estimators = legacy_estimators
)

reference_AB_matrix <- function(field) {
  matrix(
    reference[[field]][, 1L],
    nrow = parameters$py,
    ncol = parameters$px
  )
}

reference_sparse_B <- matrix(
  reference$result_B_sparse[, 1L],
  nrow = parameters$r_RR,
  ncol = parameters$px
)

closed_form_differences <- c(
  naive = assert_close(
    replicate_one$naive$AB,
    reference_AB_matrix("result_AB_naive"),
    "Naive MR-rr first sparse-simulation replicate",
    1e-8
  ),
  standard = assert_close(
    replicate_one$standard$AB,
    reference_AB_matrix("result_AB_standard"),
    "Standard MR-rr first sparse-simulation replicate",
    1e-8
  ),
  regularized = assert_close(
    replicate_one$regularized$AB,
    reference_AB_matrix("result_AB_regularized"),
    "Regularized MR-rr first sparse-simulation replicate",
    1e-8
  )
)

# CVXR/OSQP can differ slightly across compatible package builds. Exact support
# is required; a small numerical tolerance is allowed for nonzero coefficients.
sparse_tolerance <- 1e-3
sparse_differences <- c(
  AB = assert_close(
    replicate_one$sparse$AB,
    reference_AB_matrix("result_AB_sparse"),
    "Sparse MR-rr AB for the first replicate",
    sparse_tolerance
  ),
  B = assert_close(
    replicate_one$sparse$B,
    reference_sparse_B,
    "Aligned sparse MR-rr B for the first replicate",
    sparse_tolerance
  )
)

actual_support <-
  abs(replicate_one$sparse$B) >= reference$support_threshold
reference_support <-
  abs(reference_sparse_B) >= reference$support_threshold

if (!identical(actual_support, reference_support)) {
  stop(
    "The first-replicate sparse support differs from the frozen result.",
    call. = FALSE
  )
}

for (fit_name in c("naive", "standard", "regularized", "sparse")) {
  fit <- replicate_one[[fit_name]]

  if (!identical(dim(fit$A), c(3L, 2L)) ||
      !identical(dim(fit$B), c(2L, 9L)) ||
      !identical(dim(fit$AB), c(3L, 9L)) ||
      any(!is.finite(fit$AB)) ||
      max(abs(fit$A %*% fit$B - fit$AB)) > 1e-10) {
    stop(
      "Invalid estimator output for: ",
      fit_name,
      call. = FALSE
    )
  }
}

cat("Frozen sparse-simulation calibration: PASS\n")
cat(
  "Maximum effect-construction differences:",
  paste(
    paste0(
      names(effect_differences),
      "=",
      format(effect_differences, digits = 3)
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Frozen setting differences:",
  paste(
    paste0(
      names(settings_differences),
      "=",
      format(settings_differences, digits = 3)
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "First-replicate closed-form differences:",
  paste(
    paste0(
      names(closed_form_differences),
      "=",
      format(closed_form_differences, digits = 3)
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "First-replicate sparse differences:",
  paste(
    paste0(
      names(sparse_differences),
      "=",
      format(sparse_differences, digits = 3)
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "First-replicate selected coefficients:",
  sum(actual_support),
  "of",
  length(actual_support),
  "\n"
)
cat("Frozen sparse support for replicate one: PASS\n")
cat("Frozen sparse-loading simulation: PASS\n")
