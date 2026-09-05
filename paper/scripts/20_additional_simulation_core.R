#!/usr/bin/env Rscript

# Shared configuration, truth-construction, and seed utilities for the two
# additional MR-rr simulations. Running this file directly performs only
# deterministic design validation; it does not run Monte Carlo simulations.

options(warn = 1)


additional_required_config_columns <- c(
  "analysis_id",
  "setting_index",
  "scenario",
  "measurement_error_weight",
  "genetic_effect_weight",
  "truth_rank",
  "working_rank",
  "third_singular_value",
  "method_set",
  "bootstrap_rule",
  "monte_carlo_replicates",
  "bootstrap_size",
  "seed_base",
  "data_seed_group",
  "regularization_rate",
  "sparse_lambda",
  "status"
)


additional_setting_reference <- data.frame(
  setting_index = seq_len(4L),
  scenario = c(
    "me_2.5_effect_0.25",
    "me_1_effect_0.25",
    "me_2.5_effect_1",
    "me_1_effect_1"
  ),
  measurement_error_weight = c(2.5, 1, 2.5, 1),
  genetic_effect_weight = c(0.25, 0.25, 1, 1),
  regularization_rate = c(
    1.007845e-10,
    1.578970e-12,
    2.340371e-12,
    3.103420e-15
  ),
  sparse_lambda = c(1e-4, 1e-3, 1e-3, 1e-3),
  stringsAsFactors = FALSE,
  row.names = NULL
)


additional_find_repo_root <- function(start = getwd()) {
  file_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )
  starts <- start

  if (length(file_argument) > 0L) {
    script_path <- sub("^--file=", "", file_argument[[1L]])
    starts <- c(
      dirname(normalizePath(
        script_path,
        winslash = "/",
        mustWork = TRUE
      )),
      starts
    )
  }

  for (candidate_start in unique(starts)) {
    current <- normalizePath(
      candidate_start,
      winslash = "/",
      mustWork = TRUE
    )

    repeat {
      if (file.exists(file.path(current, "DESCRIPTION")) &&
          dir.exists(file.path(current, "paper", "scripts")) &&
          dir.exists(file.path(
            current,
            "freeze",
            "current_analysis_20260823",
            "project"
          ))) {
        return(current)
      }

      parent <- dirname(current)

      if (identical(parent, current)) {
        break
      }

      current <- parent
    }
  }

  stop(
    "Could not locate the MR.rr repository root.",
    call. = FALSE
  )
}


additional_require_file <- function(path, label = "Required file") {
  if (!file.exists(path)) {
    stop(
      label,
      " is missing: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}


additional_parse_numeric_column <- function(config, column) {
  parsed <- suppressWarnings(as.numeric(config[[column]]))

  if (any(!is.finite(parsed))) {
    stop(
      "Configuration column `",
      column,
      "` must contain only finite numeric values.",
      call. = FALSE
    )
  }

  parsed
}


additional_parse_integer_column <- function(
    config,
    column,
    minimum = 0L) {
  parsed_numeric <- additional_parse_numeric_column(config, column)

  if (any(parsed_numeric != round(parsed_numeric)) ||
      any(parsed_numeric < minimum) ||
      any(parsed_numeric > .Machine$integer.max)) {
    stop(
      "Configuration column `",
      column,
      "` must contain integers greater than or equal to ",
      minimum,
      ".",
      call. = FALSE
    )
  }

  as.integer(parsed_numeric)
}


additional_relative_equal <- function(
    actual,
    expected,
    tolerance = 1e-10) {
  difference <- abs(actual - expected)
  scale <- pmax(abs(expected), .Machine$double.xmin)
  difference <= tolerance * scale
}


additional_split_methods <- function(method_set) {
  methods <- strsplit(
    method_set,
    split = ";",
    fixed = TRUE
  )[[1L]]
  methods[nzchar(methods)]
}


additional_validate_shared_config <- function(config, path) {
  missing_columns <- setdiff(
    additional_required_config_columns,
    names(config)
  )
  extra_columns <- setdiff(
    names(config),
    additional_required_config_columns
  )

  if (length(missing_columns) > 0L) {
    stop(
      "Configuration is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(extra_columns) > 0L) {
    stop(
      "Configuration contains unexpected columns: ",
      paste(extra_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (nrow(config) == 0L) {
    stop(
      "Configuration has no data rows: ",
      path,
      call. = FALSE
    )
  }

  character_columns <- c(
    "analysis_id",
    "scenario",
    "method_set",
    "bootstrap_rule",
    "data_seed_group",
    "status"
  )

  for (column in character_columns) {
    config[[column]] <- trimws(as.character(config[[column]]))

    if (any(is.na(config[[column]]) | !nzchar(config[[column]]))) {
      stop(
        "Configuration column `",
        column,
        "` contains a missing or empty value.",
        call. = FALSE
      )
    }
  }

  integer_columns <- c(
    "setting_index",
    "truth_rank",
    "working_rank",
    "monte_carlo_replicates",
    "bootstrap_size",
    "seed_base"
  )
  integer_minima <- c(1L, 1L, 1L, 1L, 1L, 1L)

  for (column_index in seq_along(integer_columns)) {
    column <- integer_columns[[column_index]]
    config[[column]] <- additional_parse_integer_column(
      config,
      column,
      minimum = integer_minima[[column_index]]
    )
  }

  numeric_columns <- c(
    "measurement_error_weight",
    "genetic_effect_weight",
    "third_singular_value",
    "regularization_rate",
    "sparse_lambda"
  )

  for (column in numeric_columns) {
    config[[column]] <- additional_parse_numeric_column(config, column)
  }

  if (any(config$measurement_error_weight <= 0) ||
      any(config$genetic_effect_weight <= 0) ||
      any(config$third_singular_value < 0) ||
      any(config$regularization_rate < 0) ||
      any(config$sparse_lambda <= 0)) {
    stop(
      "Configuration contains an invalid nonpositive design value.",
      call. = FALSE
    )
  }

  reference_rows <- match(
    config$setting_index,
    additional_setting_reference$setting_index
  )

  if (any(is.na(reference_rows))) {
    stop(
      "Every setting index must be an integer from 1 through 4.",
      call. = FALSE
    )
  }

  reference <- additional_setting_reference[
    reference_rows,
    ,
    drop = FALSE
  ]

  if (!all(config$scenario == reference$scenario) ||
      !all(additional_relative_equal(
        config$measurement_error_weight,
        reference$measurement_error_weight
      )) ||
      !all(additional_relative_equal(
        config$genetic_effect_weight,
        reference$genetic_effect_weight
      )) ||
      !all(additional_relative_equal(
        config$regularization_rate,
        reference$regularization_rate
      )) ||
      !all(additional_relative_equal(
        config$sparse_lambda,
        reference$sparse_lambda
      ))) {
    stop(
      paste(
        "Configuration settings or tuning parameters differ from",
        "the locked generic-C reference settings."
      ),
      call. = FALSE
    )
  }

  expected_seed_group <- paste0(
    "generic_setting_",
    config$setting_index
  )

  if (!all(config$data_seed_group == expected_seed_group)) {
    stop(
      "Data seed groups do not match their setting indices.",
      call. = FALSE
    )
  }

  if (length(unique(config$monte_carlo_replicates)) != 1L ||
      unique(config$monte_carlo_replicates) != 1000L ||
      length(unique(config$bootstrap_size)) != 1L ||
      unique(config$bootstrap_size) != 300L ||
      length(unique(config$seed_base)) != 1L ||
      unique(config$seed_base) != 123L) {
    stop(
      paste(
        "Configuration must use 1,000 Monte Carlo replicates,",
        "bootstrap size 300, and seed base 123."
      ),
      call. = FALSE
    )
  }

  config
}


additional_validate_rank_config <- function(config) {
  expected_methods <- "regularized_mr_rr;sparse_mr_rr"
  expected_bootstrap_rule <- "regularized_mr_rr_only"

  if (nrow(config) != 12L ||
      !all(config$analysis_id == "rank_misspecification") ||
      !all(config$truth_rank == 2L) ||
      !all(config$third_singular_value == 0) ||
      !all(config$method_set == expected_methods) ||
      !all(config$bootstrap_rule == expected_bootstrap_rule) ||
      !all(config$status == "locked")) {
    stop(
      "Rank-misspecification configuration has an invalid design field.",
      call. = FALSE
    )
  }

  expected_pairs <- expand.grid(
    setting_index = seq_len(4L),
    working_rank = seq_len(3L),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  actual_pairs <- unique(
    config[, c("setting_index", "working_rank")]
  )
  expected_keys <- paste(
    expected_pairs$setting_index,
    expected_pairs$working_rank,
    sep = "/"
  )
  actual_keys <- paste(
    actual_pairs$setting_index,
    actual_pairs$working_rank,
    sep = "/"
  )

  if (anyDuplicated(paste(
    config$setting_index,
    config$working_rank,
    sep = "/"
  )) || !setequal(actual_keys, expected_keys)) {
    stop(
      "Rank configuration must contain each setting/rank pair exactly once.",
      call. = FALSE
    )
  }

  expected_method_vector <- c(
    "regularized_mr_rr",
    "sparse_mr_rr"
  )

  if (!all(vapply(
    config$method_set,
    function(value) identical(
      additional_split_methods(value),
      expected_method_vector
    ),
    logical(1)
  ))) {
    stop(
      "Rank configuration has an invalid method ordering.",
      call. = FALSE
    )
  }

  invisible(config)
}


additional_validate_approximate_config <- function(config) {
  expected_methods <- paste(
    c(
      "ivw",
      "srivw",
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr",
      "sparse_mr_rr",
      "mrdag"
    ),
    collapse = ";"
  )

  if (nrow(config) != 4L ||
      !all(config$analysis_id == "approximate_low_rank") ||
      !all(config$truth_rank == 3L) ||
      !all(config$working_rank == 2L) ||
      !all(config$method_set == expected_methods) ||
      !all(config$bootstrap_rule == "all_except_sparse_mr_rr") ||
      !all(config$status %in% c("provisional_delta", "locked"))) {
    stop(
      "Approximate-low-rank configuration has an invalid design field.",
      call. = FALSE
    )
  }

  if (!identical(sort(config$setting_index), seq_len(4L)) ||
      anyDuplicated(config$setting_index)) {
    stop(
      "Approximate-low-rank configuration must contain each setting once.",
      call. = FALSE
    )
  }

  delta_values <- unique(config$third_singular_value)

  if (length(delta_values) != 1L ||
      delta_values <= 0 ||
      delta_values >= 1) {
    stop(
      "Approximate-low-rank delta must be one common value between 0 and 1.",
      call. = FALSE
    )
  }

  invisible(config)
}


additional_read_config <- function(path, expected_analysis) {
  additional_require_file(path, "Configuration file")

  config <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    colClasses = "character",
    na.strings = character()
  )
  config <- additional_validate_shared_config(config, path)

  if (identical(expected_analysis, "rank_misspecification")) {
    additional_validate_rank_config(config)
  } else if (identical(expected_analysis, "approximate_low_rank")) {
    additional_validate_approximate_config(config)
  } else {
    stop(
      "Unknown expected analysis: ",
      expected_analysis,
      call. = FALSE
    )
  }

  config
}


additional_preserve_random_seed <- function(expression) {
  had_seed <- exists(
    ".Random.seed",
    envir = .GlobalEnv,
    inherits = FALSE
  )

  if (had_seed) {
    old_seed <- get(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )
  }

  on.exit({
    if (had_seed) {
      assign(
        ".Random.seed",
        old_seed,
        envir = .GlobalEnv
      )
    } else if (exists(
      ".Random.seed",
      envir = .GlobalEnv,
      inherits = FALSE
    )) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  force(expression)
}


additional_make_generic_truth <- function(
    third_singular_value = 0,
    px = 9L,
    py = 3L,
    seed = 123L) {
  if (length(third_singular_value) != 1L ||
      !is.finite(third_singular_value) ||
      third_singular_value < 0 ||
      third_singular_value >= 1) {
    stop(
      "`third_singular_value` must be one value in [0, 1).",
      call. = FALSE
    )
  }

  minimum_dimension <- min(px, py)

  if (minimum_dimension < 3L) {
    stop(
      "The generic design requires at least three singular values.",
      call. = FALSE
    )
  }

  decomposition <- additional_preserve_random_seed({
    set.seed(seed)
    candidate <- matrix(
      stats::rnorm(py * px),
      nrow = py,
      ncol = px
    )
    svd(candidate)
  })

  singular_values <- numeric(minimum_dimension)
  singular_values[seq_len(2L)] <- 1
  singular_values[[3L]] <- third_singular_value
  causal_effect <-
    decomposition$u %*%
    diag(
      singular_values,
      nrow = minimum_dimension,
      ncol = minimum_dimension
    ) %*%
    t(decomposition$v)

  list(
    C = causal_effect,
    U = decomposition$u,
    V = decomposition$v,
    singular_values = singular_values,
    numerical_rank = qr(causal_effect, tol = 1e-8)$rank,
    seed = as.integer(seed)
  )
}


additional_make_replicate_seed <- function(
    seed_base,
    setting_index,
    replicate_id) {
  values <- c(seed_base, setting_index, replicate_id)

  if (any(!is.finite(values)) ||
      any(values != round(values)) ||
      seed_base < 1 ||
      !setting_index %in% seq_len(4L) ||
      replicate_id < 1) {
    stop(
      "Invalid seed-base, setting-index, or replicate-id value.",
      call. = FALSE
    )
  }

  derived_seed <-
    as.double(seed_base) +
    (as.double(setting_index) - 1) * 100000 +
    as.double(replicate_id)

  if (derived_seed > .Machine$integer.max) {
    stop(
      "Derived replicate seed exceeds the R integer range.",
      call. = FALSE
    )
  }

  as.integer(derived_seed)
}


additional_load_baseline_core <- function(repo_root) {
  core_file <- file.path(
    repo_root,
    "paper",
    "scripts",
    "15_smoke_test_from_scratch_simulations.R"
  )
  additional_require_file(core_file, "Baseline simulation core")

  baseline <- new.env(parent = globalenv())
  sys.source(core_file, envir = baseline)

  required_functions <- c(
    "generate_generic_effect",
    "build_calibration",
    "build_parameters",
    "simulate_dataset"
  )
  missing_functions <- required_functions[
    !vapply(
      required_functions,
      exists,
      logical(1),
      envir = baseline,
      inherits = FALSE
    )
  ]

  if (length(missing_functions) > 0L) {
    stop(
      "Baseline core is missing functions: ",
      paste(missing_functions, collapse = ", "),
      call. = FALSE
    )
  }

  baseline
}


additional_load_estimators <- function(repo_root) {
  estimator_file <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project",
    "scripts",
    "MR_rr_estimators.R"
  )
  additional_require_file(
    estimator_file,
    "Frozen estimator implementation"
  )

  estimators <- new.env(parent = globalenv())
  sys.source(estimator_file, envir = estimators)
  estimators
}


additional_build_calibration <- function(
    repo_root,
    baseline,
    pz = 177L) {
  data_root <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project",
    "data"
  )
  data_file <- file.path(data_root, "dat_1e-4.csv")
  correlation_file <- file.path(data_root, "rho_mat_1e-4.csv")
  additional_require_file(data_file, "Frozen simulation data")
  additional_require_file(
    correlation_file,
    "Frozen correlation matrix"
  )

  lip_data <- utils::read.csv(
    data_file,
    check.names = FALSE
  )
  lip_correlation <- utils::read.csv(
    correlation_file,
    check.names = FALSE
  )

  baseline$build_calibration(
    lip_data = lip_data,
    lip_correlation = lip_correlation,
    pz = pz
  )
}


additional_build_parameters <- function(
    config_row,
    causal_effect,
    calibration,
    estimators,
    baseline) {
  if (nrow(config_row) != 1L) {
    stop(
      "`config_row` must contain exactly one row.",
      call. = FALSE
    )
  }

  baseline$build_parameters(
    causal_effect = causal_effect,
    exposure_error_weight = config_row$measurement_error_weight,
    genetic_effect_weight = config_row$genetic_effect_weight,
    calibration = calibration,
    estimators = estimators,
    rank = config_row$working_rank
  )
}


additional_simulate_dataset <- function(
    config_row,
    causal_effect,
    replicate_id,
    calibration,
    estimators,
    baseline) {
  parameters <- additional_build_parameters(
    config_row = config_row,
    causal_effect = causal_effect,
    calibration = calibration,
    estimators = estimators,
    baseline = baseline
  )
  replicate_seed <- additional_make_replicate_seed(
    seed_base = config_row$seed_base,
    setting_index = config_row$setting_index,
    replicate_id = replicate_id
  )

  set.seed(
    replicate_seed,
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )

  list(
    data = baseline$simulate_dataset(parameters),
    parameters = parameters,
    replicate_seed = replicate_seed
  )
}


additional_run_core_validation <- function() {
  repo_root <- additional_find_repo_root()
  rank_config_file <- file.path(
    repo_root,
    "paper",
    "config",
    "rank_misspecification_settings.csv"
  )
  approximate_config_file <- file.path(
    repo_root,
    "paper",
    "config",
    "approximate_low_rank_settings.csv"
  )
  rank_config <- additional_read_config(
    rank_config_file,
    "rank_misspecification"
  )
  approximate_config <- additional_read_config(
    approximate_config_file,
    "approximate_low_rank"
  )
  baseline <- additional_load_baseline_core(repo_root)

  exact_truth <- additional_make_generic_truth(
    third_singular_value = 0,
    seed = 123L
  )
  baseline_truth <- baseline$generate_generic_effect(
    px = 9L,
    py = 3L,
    rank = 2L,
    seed = 123L
  )
  exact_difference <- max(abs(exact_truth$C - baseline_truth))

  if (!is.finite(exact_difference) || exact_difference > 1e-12) {
    stop(
      "Exact-rank truth differs from the baseline generic truth.",
      call. = FALSE
    )
  }

  delta <- unique(approximate_config$third_singular_value)
  approximate_truth <- additional_make_generic_truth(
    third_singular_value = delta,
    seed = 123L
  )
  observed_singular_values <- svd(
    approximate_truth$C,
    nu = 0L,
    nv = 0L
  )$d
  expected_singular_values <- c(1, 1, delta)
  singular_value_difference <- max(abs(
    observed_singular_values - expected_singular_values
  ))

  if (exact_truth$numerical_rank != 2L ||
      approximate_truth$numerical_rank != 3L ||
      !is.finite(singular_value_difference) ||
      singular_value_difference > 1e-12 ||
      max(abs(exact_truth$U - approximate_truth$U)) > 0 ||
      max(abs(exact_truth$V - approximate_truth$V)) > 0) {
    stop(
      "Additional truth-matrix validation failed.",
      call. = FALSE
    )
  }

  seed_check <- vapply(
    seq_len(nrow(rank_config)),
    function(row_index) {
      additional_make_replicate_seed(
        seed_base = rank_config$seed_base[[row_index]],
        setting_index = rank_config$setting_index[[row_index]],
        replicate_id = 1L
      )
    },
    integer(1)
  )

  for (setting_index in seq_len(4L)) {
    selected_seeds <- seed_check[
      rank_config$setting_index == setting_index
    ]

    if (length(unique(selected_seeds)) != 1L) {
      stop(
        "Working ranks do not share a data seed in setting ",
        setting_index,
        ".",
        call. = FALSE
      )
    }
  }

  cat("MR-rr additional-simulation core validation\n")
  cat("Repository root:", repo_root, "\n")
  cat("Rank-misspecification rows:", nrow(rank_config), "\n")
  cat("Approximate-low-rank rows:", nrow(approximate_config), "\n")
  cat("Configured delta:", format(delta, digits = 6), "\n")
  cat(
    "Exact generic truth match: PASS (maximum difference =",
    format(exact_difference, scientific = TRUE),
    ")\n"
  )
  cat(
    "Approximate truth singular values: PASS (maximum difference =",
    format(singular_value_difference, scientific = TRUE),
    ")\n"
  )
  cat("Common seeds across working ranks: PASS\n")
  cat("Frozen files written: no\n")
  cat("Additional-simulation core validation: PASS\n")

  invisible(list(
    rank_config = rank_config,
    approximate_config = approximate_config,
    exact_truth = exact_truth,
    approximate_truth = approximate_truth
  ))
}


if (sys.nframe() == 0L) {
  additional_run_core_validation()
}
