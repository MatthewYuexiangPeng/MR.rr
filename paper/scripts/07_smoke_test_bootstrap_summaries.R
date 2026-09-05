# Smoke test for the frozen bootstrap summaries used in Tables 3 and 7
#
# This script reads all frozen no-MrDAG and MrDAG-only bootstrap result files,
# reconstructs the SE and CP summaries produced by make_summary.R, optionally
# compares them with summary CSV files, and checks the values printed in the
# manuscript. It is read-only and does not regenerate or overwrite any CSV.

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


rename_estimator <- function(estimator) {
  replacements <- c(
    adIVW = "SRIVW",
    Naive = "Naive MR-rr",
    MR = "MR-rr",
    MR_r = "Reg. MR-rr"
  )

  replace <- match(estimator, names(replacements))
  output <- estimator
  matched <- !is.na(replace)
  output[matched] <- unname(replacements[replace[matched]])
  output
}


summarize_median_iqr <- function(x, scale = 1, digits = 3L) {
  x <- scale * x

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


compute_average_se_by_entry <- function(
    estimates_list,
    effect_dimension,
    bootstrap_size,
    simulation_count) {
  if (length(estimates_list) != simulation_count) {
    stop(
      "An estimates list has an unexpected simulation count.",
      call. = FALSE
    )
  }

  dimensions_ok <- vapply(
    estimates_list,
    function(draws) {
      identical(
        dim(draws),
        c(bootstrap_size, effect_dimension)
      )
    },
    logical(1)
  )

  if (!all(dimensions_ok)) {
    stop(
      "A stored bootstrap-draw matrix has unexpected dimensions.",
      call. = FALSE
    )
  }

  se_matrix <- vapply(
    estimates_list,
    function(draws) {
      apply(
        draws,
        2L,
        stats::sd,
        na.rm = TRUE
      )
    },
    numeric(effect_dimension)
  )

  rowMeans(se_matrix, na.rm = TRUE)
}


compute_cp_by_entry <- function(
    coverage_matrix,
    effect_dimension,
    simulation_count) {
  if (!identical(
    dim(coverage_matrix),
    c(effect_dimension, simulation_count)
  )) {
    stop(
      "A coverage matrix has unexpected dimensions.",
      call. = FALSE
    )
  }

  rowMeans(coverage_matrix, na.rm = TRUE)
}


sort_summary_rows <- function(data, setting_order, estimator_order) {
  ordering <- order(
    match(data$Setting, setting_order),
    match(data$Estimator, estimator_order)
  )
  output <- data[ordering, , drop = FALSE]
  rownames(output) <- NULL
  output
}


compare_summary_tables <- function(
    actual,
    expected,
    label,
    tolerance = 1e-12) {
  expected_columns <- c(
    "Setting",
    "Estimator",
    "SE_med",
    "SE_q1",
    "SE_q3",
    "CP_med",
    "CP_q1",
    "CP_q3"
  )

  if (!identical(names(actual), expected_columns) ||
      !identical(names(expected), expected_columns)) {
    stop(
      label,
      " has unexpected columns.",
      call. = FALSE
    )
  }

  if (!identical(dim(actual), dim(expected))) {
    stop(
      label,
      " has unexpected dimensions.",
      call. = FALSE
    )
  }

  if (!identical(
    as.character(actual$Setting),
    as.character(expected$Setting)
  ) ||
  !identical(
    as.character(actual$Estimator),
    as.character(expected$Estimator)
  )) {
    stop(
      label,
      " has different setting or estimator labels.",
      call. = FALSE
    )
  }

  numeric_columns <- setdiff(
    expected_columns,
    c("Setting", "Estimator")
  )
  maximum_difference <- max(
    abs(
      as.matrix(actual[, numeric_columns]) -
        as.matrix(expected[, numeric_columns])
    )
  )

  if (!is.finite(maximum_difference) ||
      maximum_difference > tolerance) {
    stop(
      label,
      " differs by ",
      format(maximum_difference, scientific = TRUE),
      "; tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(maximum_difference)
}


make_manuscript_target <- function(
    setting_order,
    estimator_order,
    se_med,
    se_q1,
    se_q3,
    cp_med,
    cp_q1,
    cp_q3) {
  expected_length <-
    length(setting_order) * length(estimator_order)
  supplied_lengths <- vapply(
    list(se_med, se_q1, se_q3, cp_med, cp_q1, cp_q3),
    length,
    integer(1)
  )

  if (any(supplied_lengths != expected_length)) {
    stop(
      "A manuscript target vector has incorrect length.",
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
    SE_med = se_med,
    SE_q1 = se_q1,
    SE_q3 = se_q3,
    CP_med = cp_med,
    CP_q1 = cp_q1,
    CP_q3 = cp_q3,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


repo_root <- locate_repo_root()
freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)
bootstrap_results_override <- trimws(
  Sys.getenv(
    "MRRR_BOOTSTRAP_RESULTS_DIR",
    unset = ""
  )
)

if (nzchar(bootstrap_results_override)) {
  results_root <- normalizePath(
    bootstrap_results_override,
    winslash = "/",
    mustWork = TRUE
  )
} else {
  results_root <- file.path(freeze_root, "results")
}

cat("Bootstrap result directory:", results_root, "\n")

make_summary_file <- file.path(
  freeze_root,
  "scripts",
  "mian_sim_bootstrap_cluster_version",
  "scripts",
  "make_summary.R"
)

setting_specification <- data.frame(
  setting_index = 1:4,
  me_weight = c("2.5", "1", "2.5", "1"),
  effect_weight = c("0.25", "0.25", "1", "1"),
  parameter_index = c(1L, 3L, 2L, 4L),
  stringsAsFactors = FALSE
)

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
  "MrDAG"
)

designs <- list(
  regular_C = list(
    table_number = 3L,
    csv = "bootstrap_summary_SE_CP_med_IQR_regular_C.csv"
  ),
  sparseB_C = list(
    table_number = 7L,
    csv = "bootstrap_summary_SE_CP_med_IQR_sparseB_C.csv"
  )
)


# Values visually verified against manuscript Tables 3 and 7 (PDF dated
# 2026-08-23). Rows are ordered by setting and then estimator_order.
manuscript_targets <- list(
  regular_C = make_manuscript_target(
    setting_order = setting_order,
    estimator_order = estimator_order,
    se_med = c(
      0.095, 12.665, 0.081, 34.481, 0.289, 0.045,
      0.110, 7.331, 0.097, 8.513, 0.406, 0.054,
      0.066, 0.307, 0.050, 0.153, 0.113, 0.057,
      0.069, 0.091, 0.055, 0.071, 0.071, 0.060
    ),
    se_q1 = c(
      0.087, 0.761, 0.076, 13.373, 0.208, 0.039,
      0.099, 0.338, 0.085, 3.774, 0.207, 0.044,
      0.055, 0.200, 0.045, 0.096, 0.078, 0.042,
      0.055, 0.063, 0.044, 0.052, 0.052, 0.043
    ),
    se_q3 = c(
      0.111, 29.646, 0.095, 43.926, 0.356, 0.059,
      0.140, 15.133, 0.118, 13.353, 0.650, 0.065,
      0.081, 0.881, 0.068, 0.278, 0.173, 0.068,
      0.089, 0.128, 0.071, 0.093, 0.093, 0.069
    ),
    cp_med = c(
      85.9, 99.4, 78.6, 100.0, 99.6, 29.3,
      92.7, 98.1, 89.7, 99.2, 98.9, 47.0,
      90.9, 96.6, 86.8, 95.9, 95.9, 55.3,
      93.7, 95.3, 93.0, 95.0, 95.0, 68.2
    ),
    cp_q1 = c(
      60.2, 97.9, 53.0, 100.0, 99.2, 18.9,
      87.1, 97.0, 83.4, 99.0, 98.7, 31.9,
      81.5, 95.7, 77.8, 95.2, 95.1, 47.2,
      92.2, 94.6, 90.7, 94.3, 94.3, 56.5
    ),
    cp_q3 = c(
      91.3, 99.9, 89.3, 100.0, 99.8, 42.6,
      93.8, 99.3, 93.8, 99.3, 99.0, 53.2,
      94.2, 97.0, 92.8, 96.2, 96.2, 64.5,
      94.8, 95.7, 93.9, 95.3, 95.3, 74.5
    )
  ),
  sparseB_C = make_manuscript_target(
    setting_order = setting_order,
    estimator_order = estimator_order,
    se_med = c(
      0.092, 23.873, 0.086, 26.311, 0.256, 0.043,
      0.112, 13.422, 0.099, 7.912, 0.336, 0.048,
      0.061, 0.414, 0.054, 0.171, 0.103, 0.042,
      0.064, 0.084, 0.054, 0.068, 0.068, 0.041
    ),
    se_q1 = c(
      0.083, 0.161, 0.074, 15.434, 0.215, 0.031,
      0.099, 0.250, 0.090, 3.010, 0.224, 0.035,
      0.055, 0.125, 0.044, 0.121, 0.086, 0.026,
      0.055, 0.063, 0.047, 0.054, 0.054, 0.029
    ),
    se_q3 = c(
      0.112, 77.123, 0.106, 44.157, 0.394, 0.050,
      0.142, 36.761, 0.133, 13.921, 0.566, 0.061,
      0.080, 0.874, 0.075, 0.342, 0.181, 0.066,
      0.088, 0.122, 0.080, 0.105, 0.105, 0.068
    ),
    cp_med = c(
      90.7, 99.7, 89.7, 100.0, 99.7, 39.6,
      93.9, 98.1, 93.9, 99.3, 99.0, 48.7,
      92.7, 96.0, 91.9, 96.1, 96.0, 49.1,
      94.5, 95.0, 94.0, 94.9, 94.9, 59.1
    ),
    cp_q1 = c(
      75.7, 97.2, 74.2, 100.0, 99.3, 19.1,
      90.8, 97.1, 89.9, 99.1, 98.8, 26.2,
      87.9, 95.5, 86.1, 95.3, 95.2, 36.0,
      92.7, 94.6, 92.1, 94.4, 94.4, 48.4
    ),
    cp_q3 = c(
      94.1, 99.9, 95.6, 100.0, 99.9, 48.8,
      94.8, 99.6, 96.3, 99.5, 99.4, 59.2,
      94.3, 96.5, 93.7, 96.8, 96.8, 66.2,
      95.0, 95.7, 94.6, 95.4, 95.4, 70.4
    )
  )
)


effect_dimension <- 27L
simulation_count <- 1000L
bootstrap_size <- 300L
expected_no_mrdag_estimators <- c(
  "IVW",
  "adIVW",
  "Naive",
  "MR",
  "MR_r"
)
expected_mrdag_estimators <- "MrDAG"


# -------------------------------------------------------------------------
# Resolve every required RData file before loading large objects
# -------------------------------------------------------------------------

required_paths <- make_summary_file

for (design_name in names(designs)) {
  for (row_index in seq_len(nrow(setting_specification))) {
    setting <- setting_specification[row_index, ]

    for (estimator_set in c("no_mrdag", "mrdag_only")) {
      filename <- sprintf(
        "CI_%s_idx%d_me%s_eff%s_%s.RData",
        estimator_set,
        setting$setting_index,
        setting$me_weight,
        setting$effect_weight,
        design_name
      )
      required_paths <- c(
        required_paths,
        file.path(results_root, filename)
      )
    }
  }
}

require_files(required_paths)
cat("Frozen bootstrap result inventory: PASS\n")


# -------------------------------------------------------------------------
# Recompute both summary tables from the 16 frozen bootstrap result files
# -------------------------------------------------------------------------

design_differences <- list()

for (design_name in names(designs)) {
  rows <- list()
  row_counter <- 1L

  cat(
    "Recomputing bootstrap summary for ",
    design_name,
    "...\n",
    sep = ""
  )

  for (row_index in seq_len(nrow(setting_specification))) {
    setting_spec <- setting_specification[row_index, ]

    for (estimator_set in c("no_mrdag", "mrdag_only")) {
      filename <- sprintf(
        "CI_%s_idx%d_me%s_eff%s_%s.RData",
        estimator_set,
        setting_spec$setting_index,
        setting_spec$me_weight,
        setting_spec$effect_weight,
        design_name
      )
      result_path <- file.path(results_root, filename)
      result_environment <- new.env(parent = emptyenv())
      loaded_objects <- load(
        result_path,
        envir = result_environment
      )

      if (!"result" %in% loaded_objects) {
        stop(
          "Bootstrap file lacks `result`: ",
          result_path,
          call. = FALSE
        )
      }

      result <- result_environment$result
      result_setting <- result$setting

      if (!identical(
        as.character(result_setting$me_weight),
        setting_spec$me_weight
      ) ||
      !identical(
        as.character(result_setting$effect_weight),
        setting_spec$effect_weight
      ) ||
      !identical(
        as.integer(result_setting$param_index),
        setting_spec$parameter_index
      )) {
        stop(
          "Setting metadata mismatch in: ",
          result_path,
          call. = FALSE
        )
      }

      expected_estimators <- if (estimator_set == "no_mrdag") {
        expected_no_mrdag_estimators
      } else {
        expected_mrdag_estimators
      }

      if (!identical(
        names(result$coverage_matrix),
        expected_estimators
      ) ||
      !identical(
        names(result$estimates),
        expected_estimators
      )) {
        stop(
          "Estimator ordering mismatch in: ",
          result_path,
          call. = FALSE
        )
      }

      setting_label <- sprintf(
        "me=%s, effect=%s",
        result_setting$me_weight,
        result_setting$effect_weight
      )

      for (estimator in expected_estimators) {
        cp_by_entry <- compute_cp_by_entry(
          coverage_matrix = result$coverage_matrix[[estimator]],
          effect_dimension = effect_dimension,
          simulation_count = simulation_count
        )
        cp_summary <- summarize_median_iqr(
          cp_by_entry,
          scale = 100,
          digits = 1L
        )

        if (is.null(result$estimates[[estimator]])) {
          stop(
            "Saved bootstrap estimates are missing for ",
            estimator,
            " in ",
            result_path,
            ".",
            call. = FALSE
          )
        }

        average_se_by_entry <- compute_average_se_by_entry(
          estimates_list = result$estimates[[estimator]],
          effect_dimension = effect_dimension,
          bootstrap_size = bootstrap_size,
          simulation_count = simulation_count
        )
        se_summary <- summarize_median_iqr(
          average_se_by_entry,
          scale = 1,
          digits = 3L
        )

        if (any(!is.finite(c(cp_summary, se_summary)))) {
          stop(
            "A recomputed summary is non-finite for ",
            estimator,
            " in ",
            result_path,
            ".",
            call. = FALSE
          )
        }

        rows[[row_counter]] <- data.frame(
          Setting = setting_label,
          Estimator = rename_estimator(estimator),
          SE_med = unname(se_summary[["med"]]),
          SE_q1 = unname(se_summary[["q1"]]),
          SE_q3 = unname(se_summary[["q3"]]),
          CP_med = unname(cp_summary[["med"]]),
          CP_q1 = unname(cp_summary[["q1"]]),
          CP_q3 = unname(cp_summary[["q3"]]),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
        row_counter <- row_counter + 1L
      }

      rm(result, result_environment)
      invisible(gc(verbose = FALSE))
    }
  }

  recomputed_summary <- do.call(rbind, rows)
  recomputed_summary <- sort_summary_rows(
    recomputed_summary,
    setting_order = setting_order,
    estimator_order = estimator_order
  )

  frozen_summary_path <- file.path(
    results_root,
    designs[[design_name]]$csv
  )

  if (file.exists(frozen_summary_path)) {
    frozen_summary <- utils::read.csv(
      frozen_summary_path,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    frozen_summary <- sort_summary_rows(
      frozen_summary,
      setting_order = setting_order,
      estimator_order = estimator_order
    )

    frozen_difference <- compare_summary_tables(
      actual = recomputed_summary,
      expected = frozen_summary,
      label = paste(design_name, "recomputed versus frozen CSV")
    )
  } else {
    frozen_difference <- NA_real_
    cat(
      "Optional frozen summary CSV not supplied: ",
      basename(frozen_summary_path),
      "\n",
      sep = ""
    )
  }

  manuscript_difference <- compare_summary_tables(
    actual = recomputed_summary,
    expected = manuscript_targets[[design_name]],
    label = paste(
      design_name,
      "recomputed versus manuscript Table",
      designs[[design_name]]$table_number
    )
  )

  design_differences[[design_name]] <- c(
    frozen_csv = frozen_difference,
    manuscript = manuscript_difference
  )

  cat(
    "Table ",
    designs[[design_name]]$table_number,
    " bootstrap SE/CP summary: PASS\n",
    sep = ""
  )
}

cat(
  "Maximum summary differences: ",
  paste(
    vapply(
      names(design_differences),
      function(design_name) {
        differences <- design_differences[[design_name]]
        paste0(
          design_name,
          "[CSV=",
          format(differences[["frozen_csv"]], digits = 3),
          ", manuscript=",
          format(differences[["manuscript"]], digits = 3),
          "]"
        )
      },
      character(1)
    ),
    collapse = "; "
  ),
  "\n",
  sep = ""
)
cat("Frozen bootstrap RData and manuscript tables: PASS\n")
