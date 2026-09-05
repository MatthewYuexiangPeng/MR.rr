#!/usr/bin/env Rscript

# One-command reproduction of all manuscript tables and figures from archived
# result objects. This is the fast artifact-reproduction workflow; it does not
# rerun the full Monte Carlo or bootstrap computations.

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


find_rscript <- function() {
  executable_name <- if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else {
    "Rscript"
  }
  candidate <- file.path(
    R.home("bin"),
    executable_name
  )

  if (file.exists(candidate)) {
    return(normalizePath(
      candidate,
      winslash = "/",
      mustWork = TRUE
    ))
  }

  candidate <- Sys.which("Rscript")

  if (length(candidate) == 1L && nzchar(candidate)) {
    return(normalizePath(
      candidate,
      winslash = "/",
      mustWork = TRUE
    ))
  }

  stop(
    "Could not locate the Rscript executable.",
    call. = FALSE
  )
}


format_elapsed <- function(seconds) {
  if (seconds < 60) {
    return(sprintf("%.1f s", seconds))
  }

  sprintf(
    "%d min %.1f s",
    floor(seconds / 60),
    seconds %% 60
  )
}


repo_root <- locate_repo_root()
scripts_root <- file.path(repo_root, "paper", "scripts")
output_root <- file.path(repo_root, "paper", "output")
rscript <- find_rscript()
setwd(repo_root)

workflow <- data.frame(
  script = c(
    "11_reproduce_manuscript_tables.R",
    "12_reproduce_real_data_figures.R",
    "13_reproduce_simulation_figures.R"
  ),
  description = c(
    "Tables 1-8, Table S3, and numerical inventories",
    "Figure 2 and supplementary Figures S2 and S12",
    "Supplementary Figures S1 and S3-S11"
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)
script_paths <- file.path(scripts_root, workflow$script)
missing_scripts <- script_paths[!file.exists(script_paths)]

if (length(missing_scripts) > 0L) {
  stop(
    paste(
      "Required reproduction scripts are missing:",
      paste(missing_scripts, collapse = "\n  "),
      sep = "\n  "
    ),
    call. = FALSE
  )
}

cat("MR-rr archived-result artifact reproduction\n")
cat("Repository root:", repo_root, "\n")
cat("Rscript:", rscript, "\n")

bootstrap_results_directory <- trimws(Sys.getenv(
  "MRRR_BOOTSTRAP_RESULTS_DIR",
  unset = ""
))

if (nzchar(bootstrap_results_directory)) {
  bootstrap_results_directory <- normalizePath(
    bootstrap_results_directory,
    winslash = "/",
    mustWork = TRUE
  )
  Sys.setenv(
    MRRR_BOOTSTRAP_RESULTS_DIR = bootstrap_results_directory
  )
  cat(
    "Simulation-bootstrap results:",
    bootstrap_results_directory,
    "\n"
  )
} else {
  cat(
    "Simulation-bootstrap results:",
    "frozen results directory fallback\n"
  )
}

real_data_results_directory <- trimws(Sys.getenv(
  "MRRR_REAL_DATA_RESULTS_DIR",
  unset = ""
))

if (nzchar(real_data_results_directory)) {
  real_data_results_directory <- normalizePath(
    real_data_results_directory,
    winslash = "/",
    mustWork = TRUE
  )
  Sys.setenv(
    MRRR_REAL_DATA_RESULTS_DIR = real_data_results_directory
  )
  cat(
    "Real-data results:",
    real_data_results_directory,
    "\n"
  )
} else {
  cat("Real-data results: frozen results directory\n")
}

cat("Artifact-generation scripts:", nrow(workflow), "\n\n")

elapsed_seconds <- numeric(nrow(workflow))
workflow_start <- proc.time()[["elapsed"]]

for (script_index in seq_len(nrow(workflow))) {
  script <- workflow[script_index, , drop = FALSE]
  script_path <- script_paths[[script_index]]

  cat(strrep("=", 78L), "\n", sep = "")
  cat(
    sprintf(
      "[%d/%d] %s\n%s\n",
      script_index,
      nrow(workflow),
      script$script,
      script$description
    )
  )
  cat(strrep("=", 78L), "\n", sep = "")

  script_start <- proc.time()[["elapsed"]]
  status <- system2(
    command = rscript,
    args = c(
      "--vanilla",
      shQuote(script_path)
    ),
    stdout = "",
    stderr = ""
  )
  elapsed_seconds[[script_index]] <-
    proc.time()[["elapsed"]] - script_start

  if (!identical(as.integer(status), 0L)) {
    stop(
      script$script,
      " failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }

  cat(
    "Completed ",
    script$script,
    " in ",
    format_elapsed(elapsed_seconds[[script_index]]),
    "\n\n",
    sep = ""
  )
}

table_files <- c(
  "table_01_true_effect_matrix.csv",
  "table_02_rank_selection_generic.csv",
  "table_03_simulation_generic.csv",
  "table_04_true_sparse_loading_matrix.csv",
  "table_05_true_sparse_effect_matrix.csv",
  "table_06_rank_selection_sparse.csv",
  "table_07_simulation_sparse.csv",
  "table_08_panel_A_rank1_latent_effects.csv",
  "table_08_panel_B_rank1_protein_loadings.csv",
  "supp_table_S3_panel_A_rank2_latent_effects.csv",
  "supp_table_S3_panel_B_rank2_protein_loadings.csv",
  "section_7_4_sparse_support_recovery.csv",
  "section_8_real_data_summary.csv",
  "figure_02_interval_exclusions.csv",
  "supp_figure_S12_interval_exclusions.csv"
)
figure_files <- c(
  "Figure_02_real_data_rank1.png",
  "Figure_S01_simulation_eta_paths.png",
  "Figure_S02_sparse_eta_paths.png",
  "Figure_S03_generic_strong_entrywise.png",
  "Figure_S04_generic_weak_entrywise.png",
  "Figure_S05_sparse_strong_entrywise.png",
  "Figure_S06_sparse_weak_entrywise.png",
  "Figure_S07_sparse_loading_B.png",
  "Figure_S08_generic_strong_prediction.png",
  "Figure_S09_generic_weak_prediction.png",
  "Figure_S10_sparse_strong_prediction.png",
  "Figure_S11_sparse_weak_prediction.png",
  "Figure_S12_real_data_rank2.png"
)
figure_source_files <- c(
  "Figure_S01_simulation_eta_paths.csv",
  "Figure_S02_sparse_eta_paths.csv"
)
artifact_paths <- c(
  file.path(output_root, "tables", table_files),
  file.path(output_root, "figures", figure_files),
  file.path(output_root, "figures", figure_source_files)
)
missing_artifacts <- artifact_paths[!file.exists(artifact_paths)]

if (length(missing_artifacts) > 0L) {
  stop(
    paste(
      "Expected manuscript artifacts are missing:",
      paste(missing_artifacts, collapse = "\n  "),
      sep = "\n  "
    ),
    call. = FALSE
  )
}

artifact_manifest <- data.frame(
  Type = c(
    rep("table", length(table_files)),
    rep("figure", length(figure_files)),
    rep("figure_source_data", length(figure_source_files))
  ),
  File = c(
    file.path("tables", table_files),
    file.path("figures", figure_files),
    file.path("figures", figure_source_files)
  ),
  Bytes = as.numeric(file.info(artifact_paths)$size),
  MD5 = unname(tools::md5sum(artifact_paths)),
  stringsAsFactors = FALSE,
  row.names = NULL
)

if (nrow(artifact_manifest) != 30L ||
    any(artifact_manifest$Bytes <= 0) ||
    any(!nzchar(artifact_manifest$MD5))) {
  stop(
    "The archived-result artifact manifest is incomplete.",
    call. = FALSE
  )
}

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)
artifact_manifest_path <- file.path(
  output_root,
  "archived_artifact_manifest.csv"
)
utils::write.csv(
  artifact_manifest,
  artifact_manifest_path,
  row.names = FALSE
)

workflow_summary <- data.frame(
  Script = workflow$script,
  Description = workflow$description,
  Elapsed_seconds = round(elapsed_seconds, 3L),
  Status = "PASS",
  stringsAsFactors = FALSE,
  row.names = NULL
)
workflow_summary_path <- file.path(
  output_root,
  "archived_reproduction_run.csv"
)
utils::write.csv(
  workflow_summary,
  workflow_summary_path,
  row.names = FALSE
)

total_elapsed <- proc.time()[["elapsed"]] - workflow_start

cat(strrep("=", 78L), "\n", sep = "")
print(workflow_summary, row.names = FALSE)
cat(strrep("=", 78L), "\n", sep = "")
cat("Verified manuscript artifacts: 30\n")
cat("Total elapsed time:", format_elapsed(total_elapsed), "\n")
cat("Archived-result manuscript reproduction: PASS\n")
