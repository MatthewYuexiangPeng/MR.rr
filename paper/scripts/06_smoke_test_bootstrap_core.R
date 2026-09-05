# Smoke test for the frozen simulation-bootstrap core
#
# This script runs a deliberately small bootstrap job through the canonical
# frozen bootstrap_core.R implementation. It verifies input-setting mapping,
# parallel execution, stored bootstrap draws, percentile intervals, coverage,
# interval lengths, and aggregate summaries. It does not write result files.

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


assert_close <- function(actual, expected, label, tolerance = 1e-12) {
  if (!identical(dim(actual), dim(expected))) {
    stop(
      label,
      " has unexpected dimensions.",
      call. = FALSE
    )
  }

  if (length(actual) != length(expected)) {
    stop(
      label,
      " has unexpected length.",
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
bootstrap_core_file <- file.path(
  freeze_root,
  "scripts",
  "mian_sim_bootstrap_cluster_version",
  "scripts",
  "bootstrap_core.R"
)
bootstrap_input_file <- file.path(
  freeze_root,
  "results",
  "simulate_result_pred_260718_sparseC.RData"
)
sparse_effect_file <- file.path(
  freeze_root,
  "results",
  "simulation_sparseB_me1_effect1_eta_1e-3_260719.RData"
)

invisible(
  lapply(
    c(
      estimator_file,
      bootstrap_core_file,
      bootstrap_input_file,
      sparse_effect_file
    ),
    require_file
  )
)


# -------------------------------------------------------------------------
# Load and validate the canonical bootstrap input
# -------------------------------------------------------------------------

input_environment <- new.env(parent = emptyenv())
input_objects <- load(
  bootstrap_input_file,
  envir = input_environment
)

if (!"simulate_result_prediction" %in% input_objects) {
  stop(
    "The canonical bootstrap input lacks `simulate_result_prediction`.",
    call. = FALSE
  )
}

simulation_input <- input_environment$simulate_result_prediction

if (is.null(simulation_input$parameters_list) ||
    length(simulation_input$parameters_list) != 4L) {
  stop(
    "The bootstrap input does not contain four parameter settings.",
    call. = FALSE
  )
}

expected_scenario_names <- c(
  "me_2.5_effect_0.25",
  "me_2.5_effect_1",
  "me_1_effect_0.25",
  "me_1_effect_1"
)

if (!identical(
  names(simulation_input$AB_list),
  expected_scenario_names
)) {
  stop(
    "The canonical bootstrap scenarios have an unexpected order.",
    call. = FALSE
  )
}

sparse_effect_environment <- new.env(parent = emptyenv())
sparse_effect_objects <- load(
  sparse_effect_file,
  envir = sparse_effect_environment
)

if (!"C" %in% sparse_effect_objects) {
  stop(
    "The frozen sparse-effect file lacks `C`.",
    call. = FALSE
  )
}

true_sparse_effect <- sparse_effect_environment$C

input_effect_differences <- vapply(
  simulation_input$parameters_list,
  function(parameters) {
    assert_close(
      parameters$C,
      true_sparse_effect,
      "Bootstrap-input causal effect",
      tolerance = 1e-10
    )
  },
  numeric(1)
)

for (parameters in simulation_input$parameters_list) {
  required_parameter_names <- c(
    "py",
    "px",
    "VX_tilde",
    "Sigma_X",
    "Sigma_Y",
    "weight.matrix",
    "C",
    "r_RR"
  )

  if (!all(required_parameter_names %in% names(parameters))) {
    stop(
      "A bootstrap parameter setting is incomplete.",
      call. = FALSE
    )
  }

  if (!identical(as.integer(parameters$py), 3L) ||
      !identical(as.integer(parameters$px), 9L) ||
      !identical(as.integer(parameters$r_RR), 2L) ||
      !identical(dim(parameters$Sigma_X), c(9L, 9L)) ||
      !identical(dim(parameters$Sigma_Y), c(3L, 3L)) ||
      !identical(dim(parameters$C), c(3L, 9L))) {
    stop(
      "A bootstrap parameter setting has unexpected dimensions.",
      call. = FALSE
    )
  }
}

cat("Canonical sparse-B bootstrap input: PASS\n")


# -------------------------------------------------------------------------
# Load the frozen functions in an isolated environment
# -------------------------------------------------------------------------

bootstrap_environment <- new.env(parent = globalenv())

sys.source(
  estimator_file,
  envir = bootstrap_environment
)
sys.source(
  bootstrap_core_file,
  envir = bootstrap_environment
)

mapping_checks <- c(
  bootstrap_environment$.get_sim_index(1L, 1L, 2L),
  bootstrap_environment$.get_sim_index(1L, 2L, 2L),
  bootstrap_environment$.get_sim_index(2L, 1L, 2L),
  bootstrap_environment$.get_sim_index(2L, 2L, 2L)
)

if (!identical(as.integer(mapping_checks), 1:4)) {
  stop(
    "The frozen bootstrap setting-index mapping is incorrect.",
    call. = FALSE
  )
}


# -------------------------------------------------------------------------
# Run a minimal end-to-end parallel bootstrap
# -------------------------------------------------------------------------

smoke_iteration_count <- 2L
smoke_bootstrap_size <- 8L
expected_estimators <- c("Naive", "MR", "MR_r")

original_working_directory <- getwd()

# The frozen worker code sources scripts/MR_rr_estimators.R by a relative path.
# Running from freeze_root preserves that canonical behavior without editing it.
setwd(freeze_root)

bootstrap_result <- tryCatch(
  {
    set.seed(20260826)

    bootstrap_environment$nonpara_bootstrap_parallel(
      parameters_list = simulation_input$parameters_list,
      me_weight = "2.5",
      effect_weight = "0.25",
      regularization_rate = NULL,
      bootstrap_size = smoke_bootstrap_size,
      iteration = smoke_iteration_count,
      r_rank = 2L,
      n_cores = 1L,
      estimator_set = "mr_only",
      save_estimates = TRUE
    )
  },
  finally = {
    setwd(original_working_directory)
  }
)


# -------------------------------------------------------------------------
# Validate the returned object and automatic setting selection
# -------------------------------------------------------------------------

required_result_names <- c(
  "setting",
  "med_coverage",
  "ci_length_summary",
  "coverage_matrix",
  "ci_length_matrix",
  "estimates"
)

if (!all(required_result_names %in% names(bootstrap_result))) {
  stop(
    "The bootstrap result object is incomplete.",
    call. = FALSE
  )
}

if (!identical(
  names(bootstrap_result$coverage_matrix),
  expected_estimators
) ||
    !identical(
      names(bootstrap_result$ci_length_matrix),
      expected_estimators
    ) ||
    !identical(
      names(bootstrap_result$estimates),
      expected_estimators
    )) {
  stop(
    "The mr_only estimator ordering is incorrect.",
    call. = FALSE
  )
}

setting <- bootstrap_result$setting

if (!identical(setting$me_weight, "2.5") ||
    !identical(setting$effect_weight, "0.25") ||
    !identical(as.integer(setting$param_index), 1L) ||
    !identical(as.integer(setting$r_rank), 2L)) {
  stop(
    "The bootstrap result recorded the wrong setting.",
    call. = FALSE
  )
}

expected_regularization_rate <-
  bootstrap_environment$regularization_rate_list[[1L]]

assert_close(
  setting$regularization_rate,
  expected_regularization_rate,
  "Automatically selected regularization rate",
  tolerance = 0
)

true_effect_vector <- as.vector(
  simulation_input$parameters_list[[setting$param_index]]$C
)
effect_dimension <- length(true_effect_vector)

if (!identical(effect_dimension, 27L)) {
  stop(
    "The causal-effect vector does not have 27 entries.",
    call. = FALSE
  )
}

for (estimator in expected_estimators) {
  coverage_matrix <- bootstrap_result$coverage_matrix[[estimator]]
  length_matrix <- bootstrap_result$ci_length_matrix[[estimator]]
  estimate_list <- bootstrap_result$estimates[[estimator]]

  if (!identical(
    dim(coverage_matrix),
    c(effect_dimension, smoke_iteration_count)
  ) ||
      !identical(
        dim(length_matrix),
        c(effect_dimension, smoke_iteration_count)
      ) ||
      length(estimate_list) != smoke_iteration_count) {
    stop(
      "Unexpected bootstrap result dimensions for ",
      estimator,
      ".",
      call. = FALSE
    )
  }

  if (any(!is.finite(coverage_matrix)) ||
      any(!coverage_matrix %in% c(0, 1)) ||
      any(!is.finite(length_matrix)) ||
      any(length_matrix < 0)) {
    stop(
      "Invalid coverage or interval lengths for ",
      estimator,
      ".",
      call. = FALSE
    )
  }

  for (iteration_index in seq_len(smoke_iteration_count)) {
    estimates <- estimate_list[[iteration_index]]

    if (!identical(
      dim(estimates),
      c(smoke_bootstrap_size, effect_dimension)
    ) || any(!is.finite(estimates))) {
      stop(
        "Invalid saved bootstrap draws for ",
        estimator,
        ", iteration ",
        iteration_index,
        ".",
        call. = FALSE
      )
    }

    for (bootstrap_index in seq_len(smoke_bootstrap_size)) {
      estimated_effect <- matrix(
        estimates[bootstrap_index, ],
        nrow = 3L,
        ncol = 9L
      )
      singular_values <- svd(
        estimated_effect,
        nu = 0,
        nv = 0
      )$d

      rank_tolerance <-
        1e-10 + 1e-8 * max(singular_values)

      if (singular_values[[3L]] > rank_tolerance) {
        stop(
          estimator,
          " produced an estimate above rank two.",
          call. = FALSE
        )
      }
    }
  }
}

cat("Minimal parallel bootstrap execution: PASS\n")
cat("Bootstrap result structure: PASS\n")


# -------------------------------------------------------------------------
# Recompute every percentile interval from the stored bootstrap draws
# -------------------------------------------------------------------------

maximum_ci_difference <- setNames(
  numeric(length(expected_estimators)),
  expected_estimators
)
maximum_coverage_difference <- maximum_ci_difference
maximum_length_difference <- maximum_ci_difference

for (estimator in expected_estimators) {
  expected_coverage <- matrix(
    NA_real_,
    nrow = effect_dimension,
    ncol = smoke_iteration_count
  )
  expected_length <- matrix(
    NA_real_,
    nrow = effect_dimension,
    ncol = smoke_iteration_count
  )

  for (iteration_index in seq_len(smoke_iteration_count)) {
    estimates <-
      bootstrap_result$estimates[[estimator]][[iteration_index]]
    interval <- apply(
      estimates,
      2L,
      stats::quantile,
      probs = c(0.025, 0.975),
      na.rm = TRUE
    )

    expected_coverage[, iteration_index] <- as.numeric(
      true_effect_vector >= interval[1L, ] &
        true_effect_vector <= interval[2L, ]
    )
    expected_length[, iteration_index] <-
      interval[2L, ] - interval[1L, ]
  }

  actual_coverage <-
    bootstrap_result$coverage_matrix[[estimator]]
  actual_length <-
    bootstrap_result$ci_length_matrix[[estimator]]

  maximum_coverage_difference[[estimator]] <- assert_close(
    actual_coverage,
    expected_coverage,
    paste(estimator, "coverage matrix"),
    tolerance = 0
  )
  maximum_length_difference[[estimator]] <- assert_close(
    actual_length,
    expected_length,
    paste(estimator, "CI-length matrix"),
    tolerance = 1e-12
  )

  expected_median_coverage <- median(
    rowMeans(expected_coverage, na.rm = TRUE),
    na.rm = TRUE
  )
  assert_close(
    bootstrap_result$med_coverage[[estimator]],
    expected_median_coverage,
    paste(estimator, "median coverage summary"),
    tolerance = 0
  )

  entrywise_mean_length <- rowMeans(
    expected_length,
    na.rm = TRUE
  )
  expected_length_summary <- c(
    median = median(entrywise_mean_length, na.rm = TRUE),
    mean = mean(entrywise_mean_length, na.rm = TRUE),
    sd = stats::sd(entrywise_mean_length, na.rm = TRUE),
    min = min(entrywise_mean_length, na.rm = TRUE),
    max = max(entrywise_mean_length, na.rm = TRUE)
  )
  actual_length_summary <- unlist(
    bootstrap_result$ci_length_summary[[estimator]][
      names(expected_length_summary)
    ],
    use.names = FALSE
  )

  maximum_ci_difference[[estimator]] <- assert_close(
    actual_length_summary,
    unname(expected_length_summary),
    paste(estimator, "CI-length summary"),
    tolerance = 1e-12
  )
}

cat("Bootstrap percentile-CI reconstruction: PASS\n")
cat("Bootstrap coverage reconstruction: PASS\n")
cat(
  "Maximum input-effect differences:",
  paste(
    format(input_effect_differences, digits = 3),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Maximum reconstructed CI-length differences:",
  paste(
    paste0(
      names(maximum_length_difference),
      "=",
      format(maximum_length_difference, digits = 3)
    ),
    collapse = ", "
  ),
  "\n"
)
cat(
  "Legacy RNG note: worker streams are not explicitly seeded; ",
  "this smoke test validates the computation and result structure, ",
  "not equality to a historical random run.\n",
  sep = ""
)
cat("Frozen bootstrap core smoke test: PASS\n")
