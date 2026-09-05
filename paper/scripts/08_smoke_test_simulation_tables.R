# Smoke test for the frozen simulation summaries used in manuscript Tables 1-7
#
# This script validates the fixed true matrices, scaled instrument strengths,
# rank-selection distributions, point-estimation Bias/SD summaries, and pooled
# sparse-support recovery metrics. Bootstrap SE/CP columns are validated by
# 07_smoke_test_bootstrap_summaries.R. This script is read-only.

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


require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]

  if (length(missing) > 0L) {
    stop(
      paste(
        "Required frozen files are missing:",
        paste(missing, collapse = "\n  "),
        sep = "\n  "
      ),
      call. = FALSE
    )
  }

  invisible(paths)
}


load_required_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)

  if (!object_name %in% loaded) {
    stop(
      basename(path),
      " does not contain `",
      object_name,
      "`.",
      call. = FALSE
    )
  }

  environment[[object_name]]
}


assert_close <- function(
    actual,
    expected,
    label,
    tolerance = 1e-12) {
  if (!identical(dim(actual), dim(expected)) ||
      length(actual) != length(expected)) {
    stop(
      label,
      " has unexpected dimensions.",
      call. = FALSE
    )
  }

  difference <- max(abs(actual - expected))

  if (!is.finite(difference) || difference > tolerance) {
    stop(
      label,
      " differs by ",
      format(difference, scientific = TRUE),
      "; tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(difference)
}


summarize_vector <- function(x, digits = 3L) {
  c(
    med = round(stats::median(x, na.rm = TRUE), digits),
    q1 = round(
      unname(stats::quantile(x, 0.25, na.rm = TRUE)),
      digits
    ),
    q3 = round(
      unname(stats::quantile(x, 0.75, na.rm = TRUE)),
      digits
    )
  )
}


make_point_target <- function(
    setting_order,
    estimator_order,
    bias_med,
    bias_q1,
    bias_q3,
    sd_med,
    sd_q1,
    sd_q3) {
  expected_length <-
    length(setting_order) * length(estimator_order)
  supplied_lengths <- vapply(
    list(bias_med, bias_q1, bias_q3, sd_med, sd_q1, sd_q3),
    length,
    integer(1)
  )

  if (any(supplied_lengths != expected_length)) {
    stop(
      "A manuscript point-summary vector has incorrect length.",
      call. = FALSE
    )
  }

  data.frame(
    Setting = rep(
      setting_order,
      each = length(estimator_order)
    ),
    Estimator = rep(
      estimator_order,
      times = length(setting_order)
    ),
    Bias_med = bias_med,
    Bias_q1 = bias_q1,
    Bias_q3 = bias_q3,
    SD_med = sd_med,
    SD_q1 = sd_q1,
    SD_q3 = sd_q3,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


compare_point_tables <- function(
    actual,
    expected,
    label,
    tolerance = 1e-12) {
  expected_columns <- c(
    "Setting",
    "Estimator",
    "Bias_med",
    "Bias_q1",
    "Bias_q3",
    "SD_med",
    "SD_q1",
    "SD_q3"
  )

  if (!identical(names(actual), expected_columns) ||
      !identical(names(expected), expected_columns) ||
      !identical(dim(actual), dim(expected))) {
    stop(
      label,
      " has unexpected structure.",
      call. = FALSE
    )
  }

  if (!identical(actual$Setting, expected$Setting) ||
      !identical(actual$Estimator, expected$Estimator)) {
    stop(
      label,
      " has different labels or row ordering.",
      call. = FALSE
    )
  }

  numeric_columns <- setdiff(
    expected_columns,
    c("Setting", "Estimator")
  )

  difference <- max(
    abs(
      as.matrix(actual[, numeric_columns]) -
        as.matrix(expected[, numeric_columns])
    )
  )

  if (!is.finite(difference) || difference > tolerance) {
    stop(
      label,
      " differs by ",
      format(difference, scientific = TRUE),
      "; tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(difference)
}


compute_point_summary <- function(
    main_result,
    sparse_result,
    paper_scenario_names,
    setting_order,
    estimator_order,
    simulation_count = 1000L,
    effect_dimension = 27L) {
  result_sources <- list(
    "IVW" = main_result$result_C_ivw_list,
    "SRIVW" = main_result$result_C_adivw_list,
    "Naive MR-rr" = main_result$result_AB_list,
    "MR-rr" = main_result$result_AB_d_list,
    "Reg. MR-rr" = main_result$result_AB_d_r_list,
    "Sparse MR-rr" = sparse_result$result_C_sparse_list,
    "MrDAG" = main_result$result_MrDAG_list
  )

  if (!identical(names(result_sources), estimator_order)) {
    stop(
      "Internal estimator ordering is incorrect.",
      call. = FALSE
    )
  }

  rows <- list()
  row_index <- 1L

  for (setting_index in seq_along(paper_scenario_names)) {
    scenario_name <- paper_scenario_names[[setting_index]]

    for (estimator in estimator_order) {
      result_matrix <- result_sources[[estimator]][[scenario_name]]

      if (!identical(
        dim(result_matrix),
        c(effect_dimension, simulation_count)
      )) {
        stop(
          "Unexpected simulation matrix dimensions for ",
          estimator,
          " in ",
          scenario_name,
          ".",
          call. = FALSE
        )
      }

      entrywise_absolute_bias <- abs(
        rowMeans(result_matrix, na.rm = TRUE)
      )
      entrywise_sd <- apply(
        result_matrix,
        1L,
        stats::sd,
        na.rm = TRUE
      )

      bias_summary <- summarize_vector(entrywise_absolute_bias)
      sd_summary <- summarize_vector(entrywise_sd)

      if (any(!is.finite(c(bias_summary, sd_summary)))) {
        stop(
          "A point-estimation summary is non-finite for ",
          estimator,
          " in ",
          scenario_name,
          ".",
          call. = FALSE
        )
      }

      rows[[row_index]] <- data.frame(
        Setting = setting_order[[setting_index]],
        Estimator = estimator,
        Bias_med = unname(bias_summary[["med"]]),
        Bias_q1 = unname(bias_summary[["q1"]]),
        Bias_q3 = unname(bias_summary[["q3"]]),
        SD_med = unname(sd_summary[["med"]]),
        SD_q1 = unname(sd_summary[["q1"]]),
        SD_q3 = unname(sd_summary[["q3"]]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
      row_index <- row_index + 1L
    }
  }

  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}


compute_rank_counts <- function(
    rank_result,
    paper_scenario_names,
    simulation_count = 1000L) {
  count_matrix <- matrix(
    0L,
    nrow = length(paper_scenario_names),
    ncol = 4L
  )

  for (setting_index in seq_along(paper_scenario_names)) {
    scenario_name <- paper_scenario_names[[setting_index]]
    setting_result <- rank_result[[scenario_name]]

    if (is.null(setting_result) ||
        length(setting_result$selected_ranks) != simulation_count) {
      stop(
        "Rank-selection result is incomplete for ",
        scenario_name,
        ".",
        call. = FALSE
      )
    }

    selected_ranks <- as.integer(setting_result$selected_ranks)

    if (anyNA(selected_ranks) ||
        any(!selected_ranks %in% 1:4)) {
      stop(
        "Invalid selected rank for ",
        scenario_name,
        ".",
        call. = FALSE
      )
    }

    count_matrix[setting_index, ] <- tabulate(
      selected_ranks,
      nbins = 4L
    )

    reconstructed_correct_rate <-
      mean(selected_ranks == 2L)

    if (abs(
      reconstructed_correct_rate - setting_result$correct_rate
    ) > 1e-12) {
      stop(
        "Stored correct rank rate is inconsistent for ",
        scenario_name,
        ".",
        call. = FALSE
      )
    }
  }

  count_matrix
}


compute_population_siv <- function(parameters, pz = 177L) {
  Sigma_X <- parameters$Sigma_X
  Sigma_gg <- parameters$VX_tilde

  decomposition <- eigen(
    (Sigma_X + t(Sigma_X)) / 2,
    symmetric = TRUE
  )

  if (any(!is.finite(decomposition$values)) ||
      any(decomposition$values <= 0)) {
    stop(
      "Sigma_X is not positive definite.",
      call. = FALSE
    )
  }

  Sigma_X_inverse_sqrt <-
    decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values)) %*%
    t(decomposition$vectors)

  standardized_strength <-
    Sigma_X_inverse_sqrt %*%
    Sigma_gg %*%
    Sigma_X_inverse_sqrt
  standardized_strength <-
    (standardized_strength + t(standardized_strength)) / 2

  eigenvalues <- eigen(
    standardized_strength,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  sqrt(pz) * min(eigenvalues)
}


compute_support_metrics <- function(result) {
  true_loading <- result$B_sparse
  estimated_loading <- result$result_B_sparse
  threshold <- result$support_threshold

  if (!identical(dim(true_loading), c(2L, 9L)) ||
      !identical(dim(estimated_loading), c(18L, 1000L)) ||
      length(threshold) != 1L ||
      !is.finite(threshold)) {
    stop(
      "The support-recovery result has unexpected structure.",
      call. = FALSE
    )
  }

  true_support <- as.numeric(as.vector(true_loading) != 0)
  estimated_support <- apply(
    estimated_loading,
    2L,
    function(column) {
      as.numeric(abs(column) >= threshold)
    }
  )
  estimated_support <- matrix(
    estimated_support,
    nrow = length(true_support),
    ncol = ncol(estimated_loading)
  )

  true_support_matrix <- matrix(
    true_support,
    nrow = length(true_support),
    ncol = ncol(estimated_loading)
  )

  true_positive <- sum(
    estimated_support == 1 & true_support_matrix == 1
  )
  false_positive <- sum(
    estimated_support == 1 & true_support_matrix == 0
  )
  true_negative <- sum(
    estimated_support == 0 & true_support_matrix == 0
  )
  false_negative <- sum(
    estimated_support == 0 & true_support_matrix == 1
  )
  exact_recovery <- colSums(
    estimated_support == true_support_matrix
  ) == length(true_support)

  c(
    threshold = threshold,
    n_exact_recovery = sum(exact_recovery),
    exact_recovery_rate = mean(exact_recovery),
    sensitivity = true_positive /
      (true_positive + false_negative),
    specificity = true_negative /
      (true_negative + false_positive),
    precision = true_positive /
      (true_positive + false_positive),
    fdr = false_positive /
      (true_positive + false_positive)
  )
}


repo_root <- locate_repo_root()
freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)
results_root <- file.path(freeze_root, "results")

result_files <- list(
  generic_main = file.path(
    results_root,
    "simulate_result_pred_260717_regularC.RData"
  ),
  generic_sparse = file.path(
    results_root,
    "simulate_result_pred_sparse_260717_regularC.RData"
  ),
  generic_rank = file.path(
    results_root,
    "simulate_result_rank_260717_regularC.RData"
  ),
  sparse_main = file.path(
    results_root,
    "simulate_result_pred_260718_sparseC.RData"
  ),
  sparse_sparse = file.path(
    results_root,
    "simulate_result_pred_sparse_260718_sparseC.RData"
  ),
  sparse_rank = file.path(
    results_root,
    "simulate_result_rank_260718_sparseC.RData"
  ),
  support = file.path(
    results_root,
    "simulation_sparseB_me1_effect1_eta_1e-3_260719.RData"
  )
)

require_files(unlist(result_files, use.names = FALSE))
cat("Frozen simulation result inventory: PASS\n")

generic_main <- load_required_object(
  result_files$generic_main,
  "simulate_result_prediction"
)
generic_sparse <- load_required_object(
  result_files$generic_sparse,
  "sparse_results"
)
generic_rank <- load_required_object(
  result_files$generic_rank,
  "rank_results_all"
)
sparse_main <- load_required_object(
  result_files$sparse_main,
  "simulate_result_prediction"
)
sparse_sparse <- load_required_object(
  result_files$sparse_sparse,
  "sparse_results"
)
sparse_rank <- load_required_object(
  result_files$sparse_rank,
  "rank_results_all"
)
support_result <- load_required_object(
  result_files$support,
  "res"
)


internal_scenario_names <- c(
  "me_2.5_effect_0.25",
  "me_2.5_effect_1",
  "me_1_effect_0.25",
  "me_1_effect_1"
)
paper_scenario_names <- internal_scenario_names[c(1L, 3L, 2L, 4L)]
setting_order <- c(
  "me=2.5, effect=0.25",
  "me=1, effect=0.25",
  "me=2.5, effect=1",
  "me=1, effect=1"
)
estimator_order <- c(
  "IVW",
  "SRIVW",
  "Naive MR-rr",
  "MR-rr",
  "Reg. MR-rr",
  "Sparse MR-rr",
  "MrDAG"
)


# -------------------------------------------------------------------------
# Tables 1, 4, and 5: fixed true matrices
# -------------------------------------------------------------------------

generic_effect_target <- matrix(
  c(
    -0.204, -0.168, 0.053, -0.020, 0.086,
     0.333,  0.137, 0.024, -0.187,
    -0.060,  0.095, -0.452, 0.429, 0.039,
     0.189, -0.172, -0.393, -0.601,
     0.354,  0.329, -0.215, 0.149, -0.146,
    -0.553, -0.293, -0.147, 0.181
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

generic_true_effect <- generic_main$parameters_list[[1L]]$C
sparse_true_effect <- sparse_main$parameters_list[[1L]]$C

assert_close(
  round(generic_true_effect, 3L),
  generic_effect_target,
  "Manuscript Table 1"
)
assert_close(
  round(support_result$B_sparse, 3L),
  sparse_loading_target,
  "Manuscript Table 4"
)
assert_close(
  round(sparse_true_effect, 3L),
  sparse_effect_target,
  "Manuscript Table 5"
)

for (parameters in generic_main$parameters_list) {
  assert_close(
    parameters$C,
    generic_true_effect,
    "Generic C held fixed across settings",
    tolerance = 1e-12
  )
}
for (parameters in sparse_main$parameters_list) {
  assert_close(
    parameters$C,
    sparse_true_effect,
    "Sparse-loading C held fixed across settings",
    tolerance = 1e-12
  )
}

cat("Manuscript Tables 1, 4, and 5 true matrices: PASS\n")


# -------------------------------------------------------------------------
# Scaled instrument strengths and Tables 2 and 6 rank distributions
# -------------------------------------------------------------------------

paper_order_indices <- c(1L, 3L, 2L, 4L)
siv_target <- c(3.60, 9.00, 14.41, 36.02)

generic_population_siv <- vapply(
  generic_main$parameters_list,
  compute_population_siv,
  numeric(1)
)
sparse_population_siv <- vapply(
  sparse_main$parameters_list,
  compute_population_siv,
  numeric(1)
)

assert_close(
  round(generic_population_siv[paper_order_indices], 2L),
  siv_target,
  "Generic scaled instrument strengths"
)
assert_close(
  round(sparse_population_siv[paper_order_indices], 2L),
  siv_target,
  "Sparse-loading scaled instrument strengths"
)

generic_rank_target <- matrix(
  c(
    124L, 339L, 537L, 0L,
      1L, 690L, 309L, 0L,
      0L, 801L, 199L, 0L,
      0L, 928L,  72L, 0L
  ),
  nrow = 4L,
  ncol = 4L,
  byrow = TRUE
)
sparse_rank_target <- matrix(
  c(
    247L, 307L, 446L, 0L,
     63L, 667L, 270L, 0L,
      0L, 782L, 218L, 0L,
      0L, 905L,  95L, 0L
  ),
  nrow = 4L,
  ncol = 4L,
  byrow = TRUE
)

generic_rank_counts <- compute_rank_counts(
  generic_rank,
  paper_scenario_names
)
sparse_rank_counts <- compute_rank_counts(
  sparse_rank,
  paper_scenario_names
)

if (!identical(generic_rank_counts, generic_rank_target)) {
  stop(
    "Manuscript Table 2 rank counts differ from the frozen result.",
    call. = FALSE
  )
}
if (!identical(sparse_rank_counts, sparse_rank_target)) {
  stop(
    "Manuscript Table 6 rank counts differ from the frozen result.",
    call. = FALSE
  )
}

cat("Scaled instrument strengths: PASS\n")
cat("Manuscript Tables 2 and 6 rank distributions: PASS\n")


# -------------------------------------------------------------------------
# Tables 3 and 7: Bias and empirical SD columns
# -------------------------------------------------------------------------

point_targets <- list(
  generic = make_point_target(
    setting_order = setting_order,
    estimator_order = estimator_order,
    bias_med = c(
      0.100, 0.205, 0.105, 0.141, 0.026, 0.051, 0.160,
      0.047, 0.045, 0.062, 0.020, 0.020, 0.036, 0.144,
      0.037, 0.017, 0.039, 0.010, 0.010, 0.007, 0.110,
      0.017, 0.006, 0.016, 0.003, 0.003, 0.013, 0.077
    ),
    bias_q1 = c(
      0.041, 0.056, 0.041, 0.083, 0.020, 0.023, 0.108,
      0.023, 0.017, 0.025, 0.010, 0.010, 0.015, 0.090,
      0.017, 0.007, 0.017, 0.006, 0.006, 0.003, 0.053,
      0.006, 0.002, 0.008, 0.002, 0.002, 0.006, 0.042
    ),
    bias_q3 = c(
      0.148, 1.311, 0.144, 0.297, 0.058, 0.096, 0.264,
      0.098, 0.098, 0.084, 0.048, 0.048, 0.071, 0.194,
      0.064, 0.025, 0.060, 0.024, 0.024, 0.021, 0.137,
      0.033, 0.011, 0.033, 0.010, 0.010, 0.026, 0.121
    ),
    sd_med = c(
      0.089, 1.544, 0.077, 1.543, 0.210, 0.686, 0.035,
      0.111, 0.221, 0.079, 0.174, 0.173, 0.125, 0.043,
      0.068, 0.105, 0.050, 0.087, 0.087, 0.077, 0.061,
      0.072, 0.101, 0.052, 0.067, 0.067, 0.061, 0.065
    ),
    sd_q1 = c(
      0.084, 0.276, 0.064, 0.892, 0.167, 0.343, 0.026,
      0.092, 0.157, 0.074, 0.111, 0.111, 0.095, 0.025,
      0.052, 0.072, 0.042, 0.059, 0.058, 0.057, 0.037,
      0.053, 0.058, 0.045, 0.049, 0.049, 0.047, 0.039
    ),
    sd_q3 = c(
      0.110, 12.500, 0.085, 2.852, 0.288, 1.106, 0.063,
      0.142, 0.525, 0.112, 0.263, 0.262, 0.190, 0.074,
      0.083, 0.174, 0.066, 0.122, 0.122, 0.111, 0.077,
      0.092, 0.126, 0.071, 0.091, 0.091, 0.086, 0.082
    )
  ),
  sparse = make_point_target(
    setting_order = setting_order,
    estimator_order = estimator_order,
    bias_med = c(
      0.059, 0.074, 0.060, 0.141, 0.035, 0.042, 0.078,
      0.033, 0.032, 0.033, 0.014, 0.014, 0.032, 0.079,
      0.021, 0.013, 0.018, 0.008, 0.008, 0.010, 0.066,
      0.019, 0.018, 0.016, 0.016, 0.016, 0.013, 0.061
    ),
    bias_q1 = c(
      0.024, 0.016, 0.026, 0.049, 0.017, 0.013, 0.044,
      0.017, 0.018, 0.023, 0.005, 0.003, 0.020, 0.049,
      0.012, 0.006, 0.012, 0.002, 0.002, 0.004, 0.044,
      0.006, 0.011, 0.010, 0.008, 0.008, 0.006, 0.034
    ),
    bias_q3 = c(
      0.123, 0.191, 0.122, 0.273, 0.071, 0.124, 0.173,
      0.091, 0.067, 0.084, 0.032, 0.028, 0.073, 0.157,
      0.060, 0.039, 0.059, 0.022, 0.022, 0.020, 0.126,
      0.042, 0.036, 0.036, 0.026, 0.026, 0.026, 0.102
    ),
    sd_med = c(
      0.093, 1.103, 0.083, 2.431, 0.214, 0.631, 0.029,
      0.107, 0.267, 0.093, 0.193, 0.188, 0.121, 0.043,
      0.067, 0.097, 0.055, 0.084, 0.083, 0.060, 0.040,
      0.074, 0.088, 0.059, 0.069, 0.069, 0.051, 0.037
    ),
    sd_q1 = c(
      0.078, 0.147, 0.062, 1.525, 0.167, 0.380, 0.021,
      0.098, 0.204, 0.077, 0.161, 0.152, 0.086, 0.024,
      0.051, 0.075, 0.043, 0.062, 0.062, 0.053, 0.022,
      0.049, 0.057, 0.041, 0.049, 0.049, 0.043, 0.015
    ),
    sd_q3 = c(
      0.108, 3.952, 0.101, 6.084, 0.347, 1.055, 0.048,
      0.148, 0.772, 0.144, 0.357, 0.346, 0.169, 0.051,
      0.090, 0.174, 0.081, 0.142, 0.142, 0.102, 0.078,
      0.100, 0.142, 0.078, 0.106, 0.106, 0.085, 0.076
    )
  )
)

generic_point_summary <- compute_point_summary(
  main_result = generic_main,
  sparse_result = generic_sparse,
  paper_scenario_names = paper_scenario_names,
  setting_order = setting_order,
  estimator_order = estimator_order
)
sparse_point_summary <- compute_point_summary(
  main_result = sparse_main,
  sparse_result = sparse_sparse,
  paper_scenario_names = paper_scenario_names,
  setting_order = setting_order,
  estimator_order = estimator_order
)

generic_point_difference <- compare_point_tables(
  generic_point_summary,
  point_targets$generic,
  "Manuscript Table 3 Bias/SD"
)
sparse_point_difference <- compare_point_tables(
  sparse_point_summary,
  point_targets$sparse,
  "Manuscript Table 7 Bias/SD"
)

cat("Manuscript Table 3 Bias/empirical SD: PASS\n")
cat("Manuscript Table 7 Bias/empirical SD: PASS\n")


# -------------------------------------------------------------------------
# Section 7.4: pooled sparse-support recovery
# -------------------------------------------------------------------------

support_metrics <- compute_support_metrics(support_result)
support_target_rounded <- c(
  sensitivity = 0.951,
  specificity = 0.374,
  precision = 0.603,
  fdr = 0.397
)

if (!identical(
  as.integer(support_metrics[["n_exact_recovery"]]),
  0L
) ||
    support_metrics[["exact_recovery_rate"]] != 0) {
  stop(
    "Sparse exact-support recovery does not match Section 7.4.",
    call. = FALSE
  )
}

assert_close(
  round(
    support_metrics[names(support_target_rounded)],
    3L
  ),
  support_target_rounded,
  "Section 7.4 support-recovery metrics"
)

cat("Section 7.4 sparse-support recovery metrics: PASS\n")
cat(
  "Maximum point-summary differences: Table 3=",
  format(generic_point_difference, digits = 3),
  ", Table 7=",
  format(sparse_point_difference, digits = 3),
  "\n",
  sep = ""
)
cat("Frozen simulation summaries and manuscript Tables 1-7: PASS\n")
