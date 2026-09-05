#!/usr/bin/env Rscript

options(warn = 1)

.find_repo_root <- function() {
  file_arg <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  starts <- getwd()

  if (length(file_arg) > 0L) {
    script_path <- sub("^--file=", "", file_arg[[1L]])
    starts <- c(dirname(normalizePath(
      script_path,
      winslash = "/",
      mustWork = TRUE
    )), starts)
  }

  for (start in unique(starts)) {
    current <- normalizePath(
      start,
      winslash = "/",
      mustWork = TRUE
    )

    repeat {
      is_repo_root <-
        file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(
          current,
          "freeze",
          "current_analysis_20260823",
          "project"
        ))

      if (is_repo_root) {
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


repo_root <- .find_repo_root()

freeze_relative <- file.path(
  "freeze",
  "current_analysis_20260823",
  "project"
)

freeze_root <- file.path(repo_root, freeze_relative)

required_files <- c(
  "FREEZE_README.md",
  "SHA256SUMS.txt",
  "metadata/package_versions.csv",
  "metadata/sessionInfo.txt",
  "data/Kennetu_2016_download_links_updated.xlsx",
  "data/lipids_total24_5e-08.csv",
  "data/lipids_total24_5e-08_cor_mat.csv",
  "data/dat_1e-4.csv",
  "data/rho_mat_1e-4.csv",
  "data/traits_1e-4.csv",
  "scripts/MR_rr_estimators.R",
  "scripts/MR_rr_simulation_main_260717_cov.R",
  "scripts/Simulation_sparse_260719.R",
  "scripts/real_data_260820.R",
  paste0(
    "scripts/mian_sim_bootstrap_cluster_version/",
    "scripts/run_bootstrap_one.R"
  ),
  paste0(
    "scripts/mian_sim_bootstrap_cluster_version/",
    "scripts/bootstrap_core.R"
  ),
  paste0(
    "scripts/mian_sim_bootstrap_cluster_version/",
    "scripts/make_summary.R"
  ),
  paste0(
    "scripts/mian_sim_bootstrap_cluster_version/",
    "run_no_mrdag.sbatch"
  ),
  paste0(
    "scripts/mian_sim_bootstrap_cluster_version/",
    "run_mrdag.sbatch"
  )
)

missing_files <- required_files[
  !file.exists(file.path(freeze_root, required_files))
]

cat("MR-rr reproducibility preflight\n")
cat("Repository root:", repo_root, "\n")
cat("Frozen project:", freeze_root, "\n")
cat("R version:", R.version.string, "\n\n")

if (length(missing_files) > 0L) {
  stop(
    paste(
      "Required frozen files are missing:",
      paste(missing_files, collapse = "\n  "),
      sep = "\n  "
    ),
    call. = FALSE
  )
}

cat("Required frozen files: PASS\n")


setwd(repo_root)

freeze_status <- system2(
  "git",
  c(
    "status",
    "--porcelain",
    "--",
    "freeze/current_analysis_20260823"
  ),
  stdout = TRUE,
  stderr = TRUE
)

if (length(freeze_status) > 0L) {
  stop(
    paste(
      "The frozen directory has working-tree changes:",
      paste(freeze_status, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

freeze_diff <- system2(
  "git",
  c(
    "diff",
    "--name-only",
    "analysis-freeze-2026-08-23",
    "--",
    "freeze/current_analysis_20260823"
  ),
  stdout = TRUE,
  stderr = TRUE
)

if (length(freeze_diff) > 0L) {
  stop(
    paste(
      "The frozen directory differs from its reference tag:",
      paste(freeze_diff, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

cat("Frozen snapshot integrity: PASS\n")


required_packages <- c(
  "mr.divw",
  "GRAPPLE",
  "MrDAG",
  "CVXR",
  "ADMM",
  "foreach",
  "doParallel",
  "ggplot2",
  "dplyr",
  "tidyr",
  "purrr",
  "latex2exp",
  "readxl",
  "tibble",
  "readr",
  "patchwork",
  "reshape2",
  "ggdist",
  "pheatmap",
  "scales",
  "MASS"
)

installed <- vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)

installed_versions <- vapply(
  required_packages,
  function(package) {
    if (requireNamespace(package, quietly = TRUE)) {
      as.character(utils::packageVersion(package))
    } else {
      NA_character_
    }
  },
  character(1)
)

version_file <- file.path(
  freeze_root,
  "metadata",
  "package_versions.csv"
)

frozen_versions_data <- utils::read.csv(
  version_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

frozen_package_names <- as.character(
  frozen_versions_data[[1L]]
)

frozen_package_versions <- as.character(
  frozen_versions_data[[2L]]
)

frozen_versions <- frozen_package_versions[
  match(required_packages, frozen_package_names)
]

package_check <- data.frame(
  package = required_packages,
  installed = installed,
  installed_version = installed_versions,
  frozen_version = frozen_versions,
  stringsAsFactors = FALSE
)

print(package_check, row.names = FALSE)

if (any(!installed)) {
  stop(
    paste(
      "Missing packages:",
      paste(required_packages[!installed], collapse = ", ")
    ),
    call. = FALSE
  )
}

version_mismatch <-
  !is.na(frozen_versions) &
  installed_versions != frozen_versions

if (any(version_mismatch)) {
  cat(
    "\nVersion differences detected for:",
    paste(
      required_packages[version_mismatch],
      collapse = ", "
    ),
    "\n"
  )
}

cat("\nReproducibility environment preflight: PASS\n")
