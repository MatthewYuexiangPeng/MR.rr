#!/usr/bin/env Rscript

# Dispatch one row of the full simulation task manifest.
#
# This script is the boundary between a scheduler and the portable chunk
# worker. Slurm supplies a resource class and array index; this dispatcher
# validates the corresponding CSV row and calls 16_run_simulation_chunk.R.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        file.exists(file.path(
          current,
          "paper",
          "scripts",
          "16_run_simulation_chunk.R"
        )) &&
        file.exists(file.path(
          current,
          "paper",
          "config",
          "full_simulation_tasks.csv"
        ))) {
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


parse_named_arguments <- function(arguments) {
  if (length(arguments) == 0L) {
    return(list())
  }

  if (any(!grepl("^--[^=]+=.*$", arguments))) {
    stop(
      "Every argument must use --name=value syntax.",
      call. = FALSE
    )
  }

  keys <- sub("^--([^=]+)=.*$", "\\1", arguments)
  values <- sub("^--[^=]+=(.*)$", "\\1", arguments)

  if (anyDuplicated(keys)) {
    stop(
      "Duplicate command-line arguments: ",
      paste(unique(keys[duplicated(keys)]), collapse = ", "),
      call. = FALSE
    )
  }

  stats::setNames(as.list(values), keys)
}


parse_integer_argument <- function(arguments, name, minimum = 1L) {
  value <- arguments[[name]]
  parsed <- suppressWarnings(as.integer(value))

  if (is.null(value) || length(parsed) != 1L ||
      is.na(parsed) || parsed < minimum) {
    stop(
      "--",
      name,
      " must be an integer greater than or equal to ",
      minimum,
      ".",
      call. = FALSE
    )
  }

  parsed
}


parse_boolean_argument <- function(
    arguments,
    name,
    default = FALSE) {
  value <- arguments[[name]]

  if (is.null(value)) {
    return(default)
  }

  normalized <- tolower(trimws(value))

  if (normalized %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }

  if (normalized %in% c("0", "false", "no", "n", "off")) {
    return(FALSE)
  }

  stop(
    "--",
    name,
    " must be one of true/false, yes/no, or 1/0.",
    call. = FALSE
  )
}


validate_argument_names <- function(arguments, allowed_names) {
  unknown_names <- setdiff(names(arguments), allowed_names)

  if (length(unknown_names) > 0L) {
    stop(
      "Unknown arguments: ",
      paste(paste0("--", unknown_names), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


resolve_repo_path <- function(path, repo_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(
      path,
      winslash = "/",
      mustWork = FALSE
    ))
  }

  normalizePath(
    file.path(repo_root, path),
    winslash = "/",
    mustWork = FALSE
  )
}


find_rscript <- function() {
  executable_name <- if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else {
    "Rscript"
  }
  candidate <- file.path(R.home("bin"), executable_name)

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


format_command <- function(command, arguments) {
  paste(
    c(shQuote(command), vapply(arguments, shQuote, character(1))),
    collapse = " "
  )
}


dispatch_simulation_task <- function() {
  arguments <- parse_named_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "resource-class",
    "array-index",
    "manifest",
    "dry-run",
    "overwrite"
  )
  validate_argument_names(arguments, allowed_arguments)

  resource_class <- arguments[["resource-class"]]

  if (is.null(resource_class) ||
      !resource_class %in% c("standard", "mrdag")) {
    stop(
      "--resource-class must be standard or mrdag.",
      call. = FALSE
    )
  }

  array_index <- parse_integer_argument(
    arguments,
    "array-index",
    minimum = 1L
  )
  dry_run <- parse_boolean_argument(
    arguments,
    "dry-run",
    default = FALSE
  )
  overwrite <- parse_boolean_argument(
    arguments,
    "overwrite",
    default = FALSE
  )

  repo_root <- locate_repo_root()
  manifest_argument <- arguments[["manifest"]]
  manifest_file <- if (is.null(manifest_argument) ||
      !nzchar(trimws(manifest_argument))) {
    file.path(
      repo_root,
      "paper",
      "config",
      "full_simulation_tasks.csv"
    )
  } else {
    resolve_repo_path(manifest_argument, repo_root)
  }

  if (!file.exists(manifest_file)) {
    stop(
      "Task manifest is missing: ",
      manifest_file,
      call. = FALSE
    )
  }

  manifest <- utils::read.csv(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(
    "task_id",
    "resource_class",
    "array_index",
    "design",
    "phase",
    "setting",
    "replicate_start",
    "replicate_end",
    "seed_base",
    "mrdag_niter",
    "mrdag_burnin",
    "worker_script",
    "output_directory",
    "expected_output_file"
  )
  missing_columns <- setdiff(required_columns, names(manifest))

  if (length(missing_columns) > 0L) {
    stop(
      "Task manifest is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  selected <- manifest[
    manifest$resource_class == resource_class &
    manifest$array_index == array_index,
    ,
    drop = FALSE
  ]

  if (nrow(selected) != 1L) {
    stop(
      sprintf(
        paste(
          "Expected exactly one manifest row for",
          "resource_class=%s and array_index=%d; found %d."
        ),
        resource_class,
        array_index,
        nrow(selected)
      ),
      call. = FALSE
    )
  }

  worker_file <- resolve_repo_path(
    selected$worker_script,
    repo_root
  )
  output_directory <- resolve_repo_path(
    selected$output_directory,
    repo_root
  )
  expected_output_file <- file.path(
    output_directory,
    selected$expected_output_file
  )

  if (!file.exists(worker_file)) {
    stop(
      "Chunk worker is missing: ",
      worker_file,
      call. = FALSE
    )
  }

  worker_arguments <- c(
    "--vanilla",
    shQuote(worker_file),
    paste0("--design=", selected$design),
    paste0("--phase=", selected$phase),
    paste0("--setting=", selected$setting),
    paste0("--start=", selected$replicate_start),
    paste0("--end=", selected$replicate_end),
    paste0("--seed-base=", selected$seed_base),
    paste0("--output-dir=", shQuote(output_directory)),
    paste0("--overwrite=", tolower(as.character(overwrite)))
  )

  if (selected$phase == "main_mrdag") {
    if (is.na(selected$mrdag_niter) ||
        is.na(selected$mrdag_burnin)) {
      stop(
        "MrDAG task is missing niter or burnin in the manifest.",
        call. = FALSE
      )
    }

    worker_arguments <- c(
      worker_arguments,
      paste0("--mrdag-niter=", selected$mrdag_niter),
      paste0("--mrdag-burnin=", selected$mrdag_burnin)
    )
  }

  rscript <- find_rscript()

  cat("MR-rr simulation manifest dispatcher\n")
  cat("Repository root:", repo_root, "\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Resource class:", resource_class, "\n")
  cat("Array index:", array_index, "\n")
  cat("Global task ID:", selected$task_id, "\n")
  cat("Design:", selected$design, "\n")
  cat("Phase:", selected$phase, "\n")
  cat("Setting:", selected$setting, "\n")
  cat(
    "Replicates:",
    selected$replicate_start,
    "through",
    selected$replicate_end,
    "\n"
  )
  cat("Expected output:", expected_output_file, "\n")
  cat(
    "Command:",
    format_command(rscript, worker_arguments),
    "\n"
  )

  if (dry_run) {
    cat("Manifest task dry run: PASS\n")
    return(invisible(selected))
  }

  old_working_directory <- getwd()
  on.exit(setwd(old_working_directory), add = TRUE)
  setwd(repo_root)

  status <- system2(
    command = rscript,
    args = worker_arguments,
    stdout = "",
    stderr = ""
  )

  if (!identical(as.integer(status), 0L)) {
    stop(
      "Chunk worker failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }

  if (!file.exists(expected_output_file)) {
    stop(
      "Chunk worker succeeded but the expected output is missing: ",
      expected_output_file,
      call. = FALSE
    )
  }

  chunk <- readRDS(expected_output_file)
  metadata <- chunk$metadata
  metadata_matches <-
    identical(as.character(metadata$design), selected$design) &&
    identical(as.character(metadata$phase), selected$phase) &&
    identical(as.integer(metadata$setting_index), as.integer(selected$setting)) &&
    identical(
      as.integer(metadata$replicate_start),
      as.integer(selected$replicate_start)
    ) &&
    identical(
      as.integer(metadata$replicate_end),
      as.integer(selected$replicate_end)
    )

  if (!metadata_matches ||
      !isTRUE(metadata$complete) ||
      !identical(metadata$archived_result_rdata_read, FALSE) ||
      nrow(chunk$errors) != 0L) {
    stop(
      "The completed chunk failed dispatcher validation.",
      call. = FALSE
    )
  }

  cat("Manifest task output validation: PASS\n")
  cat("Manifest task: PASS\n")
  invisible(selected)
}


if (sys.nframe() == 0L) {
  dispatch_simulation_task()
}
