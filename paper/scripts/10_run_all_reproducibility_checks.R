#!/usr/bin/env Rscript

# Run every MR-rr paper reproducibility check in an independent R process.
# The checks are read-only. Any nonzero exit status stops the sequence.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
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
      )) &&
      dir.exists(file.path(current, "paper", "scripts"))

    if (is_repo_root) {
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
rscript <- find_rscript()

checks <- data.frame(
  id = sprintf("%02d", 0:9),
  script = c(
    "00_check_environment.R",
    "01_smoke_test_estimators.R",
    "02_smoke_test_real_data_core.R",
    "03_smoke_test_real_data_sparse.R",
    "04_smoke_test_generic_simulation.R",
    "05_smoke_test_sparse_simulation.R",
    "06_smoke_test_bootstrap_core.R",
    "07_smoke_test_bootstrap_summaries.R",
    "08_smoke_test_simulation_tables.R",
    "09_smoke_test_real_data_manuscript.R"
  ),
  description = c(
    "environment and frozen-snapshot integrity",
    "legacy estimator implementation",
    "primary real-data core",
    "primary real-data sparse refit",
    "generic simulation first replicate",
    "sparse-loading simulation first replicate",
    "minimal simulation-bootstrap core",
    "Tables 3 and 7 bootstrap summaries",
    "simulation manuscript Tables 1-7",
    "real-data manuscript results"
  ),
  stringsAsFactors = FALSE
)

script_paths <- file.path(scripts_root, checks$script)
missing_scripts <- script_paths[!file.exists(script_paths)]

if (length(missing_scripts) > 0L) {
  stop(
    paste(
      "Required reproducibility scripts are missing:",
      paste(missing_scripts, collapse = "\n  "),
      sep = "\n  "
    ),
    call. = FALSE
  )
}

setwd(repo_root)

cat("MR-rr complete reproducibility check\n")
cat("Repository root:", repo_root, "\n")
cat("Rscript:", rscript, "\n")

bootstrap_results_directory <- Sys.getenv(
  "MRRR_BOOTSTRAP_RESULTS_DIR",
  unset = ""
)

if (nzchar(bootstrap_results_directory)) {
  cat(
    "External simulation-bootstrap results:",
    bootstrap_results_directory,
    "\n"
  )
} else {
  cat(
    "Simulation-bootstrap results: frozen results directory fallback\n"
  )
}

real_data_results_directory <- Sys.getenv(
  "MRRR_REAL_DATA_RESULTS_DIR",
  unset = ""
)

if (nzchar(real_data_results_directory)) {
  cat(
    "External real-data results:",
    real_data_results_directory,
    "\n"
  )
} else {
  cat("Real-data results: frozen results directory\n")
}

cat("Checks scheduled:", nrow(checks), "\n\n")

elapsed_seconds <- numeric(nrow(checks))
suite_start <- proc.time()[["elapsed"]]

for (check_index in seq_len(nrow(checks))) {
  check <- checks[check_index, , drop = FALSE]
  script_path <- script_paths[[check_index]]

  cat(strrep("=", 78L), "\n", sep = "")
  cat(
    sprintf(
      "[%d/%d] %s - %s\n",
      check_index,
      nrow(checks),
      check$script,
      check$description
    )
  )
  cat(strrep("=", 78L), "\n", sep = "")

  check_start <- proc.time()[["elapsed"]]

  status <- system2(
    command = rscript,
    args = c(
      "--vanilla",
      shQuote(script_path)
    ),
    stdout = "",
    stderr = ""
  )

  elapsed_seconds[[check_index]] <-
    proc.time()[["elapsed"]] - check_start

  if (!identical(as.integer(status), 0L)) {
    stop(
      check$script,
      " failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }

  cat(
    "Completed ",
    check$script,
    " in ",
    format_elapsed(elapsed_seconds[[check_index]]),
    "\n\n",
    sep = ""
  )
}

total_elapsed <- proc.time()[["elapsed"]] - suite_start

summary_table <- data.frame(
  check = checks$id,
  script = checks$script,
  elapsed = vapply(
    elapsed_seconds,
    format_elapsed,
    character(1)
  ),
  status = "PASS",
  stringsAsFactors = FALSE
)

cat(strrep("=", 78L), "\n", sep = "")
print(summary_table, row.names = FALSE)
cat(strrep("=", 78L), "\n", sep = "")
cat(
  "Total elapsed time: ",
  format_elapsed(total_elapsed),
  "\n",
  sep = ""
)
cat("All MR-rr reproducibility checks: PASS\n")
