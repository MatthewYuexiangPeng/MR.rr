#!/usr/bin/env Rscript

# Reproduce the numerical manuscript tables from the archived result objects.
# Generated CSV files are written to paper/output/tables/. The frozen project
# and any externally supplied bootstrap results are read-only inputs.

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


require_files <- function(paths) {
  missing <- paths[!file.exists(paths)]

  if (length(missing) > 0L) {
    stop(
      paste(
        "Required files are missing:",
        paste(missing, collapse = "\n  "),
        sep = "\n  "
      ),
      call. = FALSE
    )
  }

  invisible(paths)
}


summarize_median_iqr <- function(
    values,
    scale = 1,
    digits = 3L) {
  values <- scale * values

  c(
    med = round(stats::median(values, na.rm = TRUE), digits),
    q1 = round(
      unname(stats::quantile(
        values,
        0.25,
        na.rm = TRUE
      )),
      digits
    ),
    q3 = round(
      unname(stats::quantile(
        values,
        0.75,
        na.rm = TRUE
      )),
      digits
    )
  )
}


compute_average_se_by_entry <- function(
    estimates_list,
    effect_dimension = 27L,
    bootstrap_size = 300L,
    simulation_count = 1000L) {
  if (length(estimates_list) != simulation_count) {
    stop(
      "A bootstrap estimate list has an unexpected simulation count.",
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
      "A bootstrap draw matrix has unexpected dimensions.",
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
    effect_dimension = 27L,
    simulation_count = 1000L) {
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


rename_estimator <- function(estimator) {
  replacements <- c(
    IVW = "IVW",
    adIVW = "SRIVW",
    Naive = "Naive MR-rr",
    MR = "MR-rr",
    MR_r = "Reg. MR-rr",
    MrDAG = "MrDAG"
  )

  if (!estimator %in% names(replacements)) {
    stop(
      "Unknown bootstrap estimator: ",
      estimator,
      call. = FALSE
    )
  }

  unname(replacements[[estimator]])
}


recompute_bootstrap_summaries <- function(results_root) {
  setting_specification <- data.frame(
    setting_index = 1:4,
    me_weight = c("2.5", "1", "2.5", "1"),
    effect_weight = c("0.25", "0.25", "1", "1"),
    parameter_index = c(1L, 3L, 2L, 4L),
    setting_label = c(
      "me=2.5, effect=0.25",
      "me=1, effect=0.25",
      "me=2.5, effect=1",
      "me=1, effect=1"
    ),
    stringsAsFactors = FALSE
  )

  estimator_sets <- list(
    no_mrdag = c(
      "IVW",
      "adIVW",
      "Naive",
      "MR",
      "MR_r"
    ),
    mrdag_only = "MrDAG"
  )
  design_names <- c("regular_C", "sparseB_C")

  required_paths <- character(0)

  for (design_name in design_names) {
    for (setting_index in seq_len(nrow(setting_specification))) {
      setting <- setting_specification[setting_index, ]

      for (estimator_set in names(estimator_sets)) {
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

  summaries <- list()

  for (design_name in design_names) {
    rows <- list()
    row_index <- 1L

    cat(
      "Recomputing bootstrap summaries for ",
      design_name,
      "...\n",
      sep = ""
    )

    for (setting_index in seq_len(nrow(setting_specification))) {
      setting <- setting_specification[setting_index, ]

      for (estimator_set in names(estimator_sets)) {
        filename <- sprintf(
          "CI_%s_idx%d_me%s_eff%s_%s.RData",
          estimator_set,
          setting$setting_index,
          setting$me_weight,
          setting$effect_weight,
          design_name
        )
        result_path <- file.path(results_root, filename)
        result_environment <- new.env(parent = emptyenv())
        loaded <- load(
          result_path,
          envir = result_environment
        )

        if (!"result" %in% loaded) {
          stop(
            "Bootstrap file does not contain `result`: ",
            result_path,
            call. = FALSE
          )
        }

        result <- result_environment$result

        if (!identical(
          as.character(result$setting$me_weight),
          setting$me_weight
        ) ||
            !identical(
              as.character(result$setting$effect_weight),
              setting$effect_weight
            ) ||
            !identical(
              as.integer(result$setting$param_index),
              setting$parameter_index
            )) {
          stop(
            "Bootstrap setting metadata mismatch in: ",
            result_path,
            call. = FALSE
          )
        }

        expected_estimators <- estimator_sets[[estimator_set]]

        if (!identical(
          names(result$coverage_matrix),
          expected_estimators
        ) ||
            !identical(
              names(result$estimates),
              expected_estimators
            )) {
          stop(
            "Bootstrap estimator ordering mismatch in: ",
            result_path,
            call. = FALSE
          )
        }

        for (estimator in expected_estimators) {
          cp_by_entry <- compute_cp_by_entry(
            result$coverage_matrix[[estimator]]
          )
          cp_summary <- summarize_median_iqr(
            cp_by_entry,
            scale = 100,
            digits = 1L
          )

          if (is.null(result$estimates[[estimator]])) {
            stop(
              "Bootstrap estimates are missing for ",
              estimator,
              " in ",
              result_path,
              ".",
              call. = FALSE
            )
          }

          average_se_by_entry <- compute_average_se_by_entry(
            result$estimates[[estimator]]
          )
          se_summary <- summarize_median_iqr(
            average_se_by_entry,
            digits = 3L
          )

          if (any(!is.finite(c(cp_summary, se_summary)))) {
            stop(
              "A bootstrap summary is non-finite for ",
              estimator,
              " in ",
              result_path,
              ".",
              call. = FALSE
            )
          }

          rows[[row_index]] <- data.frame(
            Setting = setting$setting_label,
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
          row_index <- row_index + 1L
        }

        rm(result, result_environment)
        invisible(gc(verbose = FALSE))
      }
    }

    summary <- do.call(rbind, rows)
    setting_order <- setting_specification$setting_label
    estimator_order <- c(
      "IVW",
      "SRIVW",
      "Naive MR-rr",
      "MR-rr",
      "Reg. MR-rr",
      "MrDAG"
    )
    ordering <- order(
      match(summary$Setting, setting_order),
      match(summary$Estimator, estimator_order)
    )
    summary <- summary[ordering, , drop = FALSE]
    rownames(summary) <- NULL
    summaries[[design_name]] <- summary
  }

  summaries
}


combine_point_and_bootstrap <- function(
    point_summary,
    bootstrap_summary) {
  point_key <- paste(
    point_summary$Setting,
    point_summary$Estimator,
    sep = "|"
  )
  bootstrap_key <- paste(
    bootstrap_summary$Setting,
    bootstrap_summary$Estimator,
    sep = "|"
  )
  bootstrap_match <- match(point_key, bootstrap_key)
  sparse_rows <- point_summary$Estimator == "Sparse MR-rr"

  if (any(is.na(bootstrap_match[!sparse_rows])) ||
      any(!is.na(bootstrap_match[sparse_rows]))) {
    stop(
      "Point and bootstrap summary rows could not be aligned.",
      call. = FALSE
    )
  }

  bootstrap_columns <- c(
    "SE_med",
    "SE_q1",
    "SE_q3",
    "CP_med",
    "CP_q1",
    "CP_q3"
  )
  bootstrap_values <- bootstrap_summary[
    bootstrap_match,
    bootstrap_columns,
    drop = FALSE
  ]

  output <- cbind(
    point_summary,
    bootstrap_values
  )
  rownames(output) <- NULL
  output
}


matrix_table <- function(
    matrix,
    row_name,
    row_labels,
    column_labels,
    digits) {
  if (!identical(
    dim(matrix),
    c(length(row_labels), length(column_labels))
  )) {
    stop(
      "A manuscript matrix has unexpected dimensions.",
      call. = FALSE
    )
  }

  values <- as.data.frame(
    round(matrix, digits),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(values) <- column_labels

  output <- data.frame(
    row_labels,
    values,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(output)[[1L]] <- row_name
  output
}


make_rank_table <- function(rank_counts, siv) {
  if (!identical(dim(rank_counts), c(4L, 4L)) ||
      length(siv) != 4L) {
    stop(
      "A rank-selection result has unexpected dimensions.",
      call. = FALSE
    )
  }

  data.frame(
    Index = 1:4,
    Sigma_X_multiplier = c(2.5, 1.0, 2.5, 1.0),
    Sigma_gamma_gamma_multiplier = c(0.25, 0.25, 1.0, 1.0),
    SIV = round(siv, 2L),
    Rank_1 = rank_counts[, 1L],
    Rank_2 = rank_counts[, 2L],
    Rank_3 = rank_counts[, 3L],
    Rank_4 = rank_counts[, 4L],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


write_output_csv <- function(data, path) {
  utils::write.csv(
    data,
    path,
    row.names = FALSE,
    na = ""
  )
  cat("Wrote:", path, "\n")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}


repo_root <- locate_repo_root()
setwd(repo_root)

scripts_root <- file.path(repo_root, "paper", "scripts")
simulation_check_path <- file.path(
  scripts_root,
  "08_smoke_test_simulation_tables.R"
)
real_data_check_path <- file.path(
  scripts_root,
  "09_smoke_test_real_data_manuscript.R"
)
require_files(c(
  simulation_check_path,
  real_data_check_path
))

cat("Loading verified simulation summaries.\n")
simulation <- new.env(parent = globalenv())
sys.source(
  simulation_check_path,
  envir = simulation
)

cat("\nLoading verified real-data summaries.\n")
real_data <- new.env(parent = globalenv())
sys.source(
  real_data_check_path,
  envir = real_data
)

freeze_results_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project",
  "results"
)
bootstrap_results_root <- trimws(Sys.getenv(
  "MRRR_BOOTSTRAP_RESULTS_DIR",
  unset = ""
))

if (nzchar(bootstrap_results_root)) {
  bootstrap_results_root <- normalizePath(
    bootstrap_results_root,
    winslash = "/",
    mustWork = TRUE
  )
} else {
  bootstrap_results_root <- freeze_results_root
}

cat(
  "\nSimulation-bootstrap result directory:",
  bootstrap_results_root,
  "\n"
)

bootstrap_summaries <- recompute_bootstrap_summaries(
  bootstrap_results_root
)

table_01 <- matrix_table(
  simulation$generic_true_effect,
  row_name = "Outcome",
  row_labels = paste0("Y", 1:3),
  column_labels = paste0("X", 1:9),
  digits = 3L
)
table_02 <- make_rank_table(
  simulation$generic_rank_counts,
  simulation$generic_population_siv[
    simulation$paper_order_indices
  ]
)
table_03 <- combine_point_and_bootstrap(
  simulation$generic_point_summary,
  bootstrap_summaries$regular_C
)
table_04 <- matrix_table(
  simulation$support_result$B_sparse,
  row_name = "Pathway",
  row_labels = paste0("Pathway ", 1:2),
  column_labels = paste0("X", 1:9),
  digits = 3L
)
table_05 <- matrix_table(
  simulation$sparse_true_effect,
  row_name = "Outcome",
  row_labels = paste0("Y", 1:3),
  column_labels = paste0("X", 1:9),
  digits = 3L
)
table_06 <- make_rank_table(
  simulation$sparse_rank_counts,
  simulation$sparse_population_siv[
    simulation$paper_order_indices
  ]
)
table_07 <- combine_point_and_bootstrap(
  simulation$sparse_point_summary,
  bootstrap_summaries$sparseB_C
)

outcome_names <- c("LAS", "CES", "SVS")
protein_names <- c(
  "MMP12",
  "CNTN1",
  "FGL1",
  "MXRA8",
  "CNTFR",
  "SCG3",
  "HTRA1",
  "CLEC3B",
  "ANTXR2"
)

table_08_A <- matrix_table(
  real_data$rank1_A_actual,
  row_name = "Outcome",
  row_labels = outcome_names,
  column_labels = c(
    "MR-rr",
    "Reg. MR-rr",
    "Sparse MR-rr"
  ),
  digits = 4L
)
table_08_B <- matrix_table(
  real_data$rank1_B_actual,
  row_name = "Protein",
  row_labels = protein_names,
  column_labels = c(
    "MR-rr",
    "Reg. MR-rr",
    "Sparse MR-rr"
  ),
  digits = 3L
)

table_S3_A <- matrix_table(
  real_data$rank2_A_actual,
  row_name = "Outcome",
  row_labels = outcome_names,
  column_labels = c(
    "MR-rr Pathway 1",
    "MR-rr Pathway 2",
    "Reg. MR-rr Pathway 1",
    "Reg. MR-rr Pathway 2",
    "Sparse MR-rr Pathway 1",
    "Sparse MR-rr Pathway 2"
  ),
  digits = 4L
)
table_S3_B <- matrix_table(
  real_data$rank2_B_actual,
  row_name = "Protein",
  row_labels = protein_names,
  column_labels = c(
    "MR-rr Pathway 1",
    "MR-rr Pathway 2",
    "Reg. MR-rr Pathway 1",
    "Reg. MR-rr Pathway 2",
    "Sparse MR-rr Pathway 1",
    "Sparse MR-rr Pathway 2"
  ),
  digits = 3L
)

support_recovery <- data.frame(
  Metric = names(simulation$support_metrics),
  Value = unname(simulation$support_metrics),
  stringsAsFactors = FALSE,
  row.names = NULL
)

real_data_summary <- data.frame(
  Instruments = 177L,
  Exposures = 9L,
  Outcomes = 3L,
  Estimated_SIV = round(real_data$estimated_siv, 2L),
  Primary_rank = as.integer(real_data$rank1_result$r_RR),
  Primary_sparse_eta = real_data$rank1_result$sparse_eta,
  Sensitivity_rank = as.integer(real_data$rank2_result$r_RR),
  Sensitivity_sparse_eta = real_data$rank2_result$sparse_eta,
  stringsAsFactors = FALSE,
  row.names = NULL
)

output_root <- file.path(
  repo_root,
  "paper",
  "output",
  "tables"
)
dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

output_specification <- list(
  "table_01_true_effect_matrix.csv" = table_01,
  "table_02_rank_selection_generic.csv" = table_02,
  "table_03_simulation_generic.csv" = table_03,
  "table_04_true_sparse_loading_matrix.csv" = table_04,
  "table_05_true_sparse_effect_matrix.csv" = table_05,
  "table_06_rank_selection_sparse.csv" = table_06,
  "table_07_simulation_sparse.csv" = table_07,
  "table_08_panel_A_rank1_latent_effects.csv" = table_08_A,
  "table_08_panel_B_rank1_protein_loadings.csv" = table_08_B,
  "supp_table_S3_panel_A_rank2_latent_effects.csv" = table_S3_A,
  "supp_table_S3_panel_B_rank2_protein_loadings.csv" = table_S3_B,
  "section_7_4_sparse_support_recovery.csv" = support_recovery,
  "section_8_real_data_summary.csv" = real_data_summary,
  "figure_02_interval_exclusions.csv" = real_data$rank1_inventory,
  "supp_figure_S12_interval_exclusions.csv" = real_data$rank2_inventory
)

generated_paths <- vapply(
  names(output_specification),
  function(filename) {
    write_output_csv(
      output_specification[[filename]],
      file.path(output_root, filename)
    )
  },
  character(1)
)

manifest <- data.frame(
  File = basename(generated_paths),
  MD5 = unname(tools::md5sum(generated_paths)),
  stringsAsFactors = FALSE,
  row.names = NULL
)
manifest_path <- write_output_csv(
  manifest,
  file.path(output_root, "manifest.csv")
)

cat("\nGenerated table directory:", output_root, "\n")
cat("Generated numerical artifacts:", length(generated_paths), "\n")
cat("Manifest:", manifest_path, "\n")
cat("Manuscript table reproduction: PASS\n")
