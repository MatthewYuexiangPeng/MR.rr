#!/usr/bin/env Rscript

# Dispatch one row of the additional-simulation point-task manifest.
#
# Slurm supplies a resource class and an array index. This dispatcher checks
# the manifest row, verifies the locked input hashes, calls the portable point
# worker, and validates the completed chunk before returning success.

options(warn = 1)


additional_dispatch_parse_arguments <- function(arguments) {
  if (length(arguments) == 0L) {
    return(list())
  }

  if (any(!grepl("^--[^=]+=.*$", arguments))) {
    stop("Every argument must use --name=value syntax.", call. = FALSE)
  }

  keys <- sub("^--([^=]+)=.*$", "\\1", arguments)
  values <- sub("^--[^=]+=(.*)$", "\\1", arguments)

  if (anyDuplicated(keys)) {
    stop(
      "Duplicate arguments: ",
      paste(unique(keys[duplicated(keys)]), collapse = ", "),
      call. = FALSE
    )
  }

  stats::setNames(as.list(values), keys)
}


additional_dispatch_validate_names <- function(arguments, allowed) {
  unknown <- setdiff(names(arguments), allowed)

  if (length(unknown) > 0L) {
    stop(
      "Unknown arguments: ",
      paste(paste0("--", unknown), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


additional_dispatch_parse_integer <- function(
    arguments,
    name,
    minimum = 1L) {
  value <- arguments[[name]]
  numeric_value <- suppressWarnings(as.numeric(value))

  if (is.null(value) ||
      length(numeric_value) != 1L ||
      !is.finite(numeric_value) ||
      numeric_value != round(numeric_value) ||
      numeric_value < minimum ||
      numeric_value > .Machine$integer.max) {
    stop(
      "--",
      name,
      " must be an integer greater than or equal to ",
      minimum,
      ".",
      call. = FALSE
    )
  }

  as.integer(numeric_value)
}


additional_dispatch_parse_boolean <- function(
    arguments,
    name,
    default = FALSE) {
  value <- arguments[[name]]

  if (is.null(value)) {
    return(default)
  }

  additional_dispatch_boolean_value(value, paste0("--", name))
}


additional_dispatch_boolean_value <- function(value, label) {
  if (length(value) != 1L || is.na(value)) {
    stop(label, " must contain exactly one Boolean value.", call. = FALSE)
  }

  if (is.logical(value)) {
    return(value)
  }

  normalized <- tolower(trimws(as.character(value)))

  if (normalized %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }

  if (normalized %in% c("0", "false", "no", "n", "off")) {
    return(FALSE)
  }

  stop(
    label,
    " must be one of true/false, yes/no, or 1/0.",
    call. = FALSE
  )
}


additional_dispatch_integer_value <- function(
    selected,
    name,
    minimum = 0L,
    allow_na = FALSE) {
  value <- selected[[name]]

  if (length(value) != 1L) {
    stop("Manifest field must be scalar: ", name, call. = FALSE)
  }

  numeric_value <- suppressWarnings(as.numeric(value))

  if (allow_na && (is.na(value) || is.na(numeric_value))) {
    return(NA_integer_)
  }

  if (!is.finite(numeric_value) ||
      numeric_value != round(numeric_value) ||
      numeric_value < minimum ||
      numeric_value > .Machine$integer.max) {
    stop("Manifest field is not a valid integer: ", name, call. = FALSE)
  }

  as.integer(numeric_value)
}


additional_dispatch_text_value <- function(selected, name) {
  value <- selected[[name]]

  if (length(value) != 1L || is.na(value) || !nzchar(trimws(value))) {
    stop("Manifest field is empty or invalid: ", name, call. = FALSE)
  }

  trimws(as.character(value))
}


additional_dispatch_resolve_path <- function(path, repo_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(
    file.path(repo_root, path),
    winslash = "/",
    mustWork = FALSE
  )
}


additional_dispatch_find_rscript <- function() {
  executable_name <- if (.Platform$OS.type == "windows") {
    "Rscript.exe"
  } else {
    "Rscript"
  }
  candidate <- file.path(R.home("bin"), executable_name)

  if (file.exists(candidate)) {
    return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }

  candidate <- Sys.which("Rscript")

  if (length(candidate) == 1L && nzchar(candidate)) {
    return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate the Rscript executable.", call. = FALSE)
}


additional_dispatch_format_command <- function(command, arguments) {
  paste(c(shQuote(command), arguments), collapse = " ")
}


additional_dispatch_check_md5 <- function(path, expected, label) {
  if (!file.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }

  if (length(expected) != 1L ||
      is.na(expected) ||
      !grepl("^[0-9a-fA-F]{32}$", expected)) {
    stop("Manifest has an invalid ", label, " MD5.", call. = FALSE)
  }

  observed <- unname(tools::md5sum(path))

  if (is.na(observed) ||
      !identical(tolower(observed), tolower(as.character(expected)))) {
    stop(
      label,
      " MD5 differs from the manifest. Regenerate the manifest after ",
      "confirming the intended code/configuration version.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


additional_dispatch_validate_completed_chunk <- function(
    output_file,
    selected,
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end,
    seed_base,
    allow_provisional,
    mrdag_niter,
    mrdag_burnin) {
  if (!file.exists(output_file)) {
    stop(
      "Point worker succeeded but the expected output is missing: ",
      output_file,
      call. = FALSE
    )
  }

  chunk <- tryCatch(
    readRDS(output_file),
    error = function(error) error
  )

  if (inherits(chunk, "error") || !is.list(chunk)) {
    stop("The completed point chunk cannot be read.", call. = FALSE)
  }

  metadata <- chunk$metadata
  required_metadata <- c(
    "analysis",
    "phase",
    "setting_index",
    "replicate_start",
    "replicate_end",
    "replicate_count",
    "seed_base",
    "mrdag_niter",
    "mrdag_burnin",
    "allow_provisional",
    "complete",
    "archived_simulation_result_rdata_read",
    "frozen_snapshot_written"
  )

  if (!is.list(metadata) ||
      length(setdiff(required_metadata, names(metadata))) > 0L) {
    stop("The completed point chunk has incomplete metadata.", call. = FALSE)
  }

  expected_count <- replicate_end - replicate_start + 1L
  metadata_matches <-
    identical(as.character(metadata$analysis), analysis) &&
    identical(as.character(metadata$phase), phase) &&
    identical(as.integer(metadata$setting_index), setting_index) &&
    identical(as.integer(metadata$replicate_start), replicate_start) &&
    identical(as.integer(metadata$replicate_end), replicate_end) &&
    identical(as.integer(metadata$replicate_count), expected_count) &&
    identical(as.integer(metadata$seed_base), seed_base) &&
    identical(isTRUE(metadata$allow_provisional), allow_provisional)

  if (phase == "mrdag") {
    metadata_matches <- metadata_matches &&
      identical(as.integer(metadata$mrdag_niter), mrdag_niter) &&
      identical(as.integer(metadata$mrdag_burnin), mrdag_burnin)
  }

  errors_are_empty <- is.data.frame(chunk$errors) && nrow(chunk$errors) == 0L
  success_counts <- chunk$successful_replicates
  success_counts_are_complete <-
    is.numeric(success_counts) &&
    length(success_counts) > 0L &&
    all(is.finite(success_counts)) &&
    all(as.integer(success_counts) == expected_count)

  if (!metadata_matches ||
      !isTRUE(metadata$complete) ||
      !identical(metadata$archived_simulation_result_rdata_read, FALSE) ||
      !identical(metadata$frozen_snapshot_written, FALSE) ||
      !errors_are_empty ||
      !success_counts_are_complete ||
      !is.list(chunk$estimates) ||
      length(chunk$estimates) == 0L) {
    stop("The completed point chunk failed dispatcher validation.", call. = FALSE)
  }

  expected_output_name <- additional_dispatch_text_value(
    selected,
    "expected_output_file"
  )

  if (!identical(basename(output_file), expected_output_name)) {
    stop("The completed chunk filename differs from the manifest.", call. = FALSE)
  }

  invisible(chunk)
}


additional_dispatch_point_task <- function() {
  arguments <- additional_dispatch_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "resource-class",
    "array-index",
    "manifest",
    "dry-run",
    "overwrite"
  )
  additional_dispatch_validate_names(arguments, allowed_arguments)
  resource_class <- arguments[["resource-class"]]

  if (is.null(resource_class) ||
      !resource_class %in% c("standard", "mrdag")) {
    stop("--resource-class must be standard or mrdag.", call. = FALSE)
  }

  array_index <- additional_dispatch_parse_integer(
    arguments,
    "array-index",
    minimum = 1L
  )
  dry_run <- additional_dispatch_parse_boolean(
    arguments,
    "dry-run",
    default = FALSE
  )
  overwrite <- additional_dispatch_parse_boolean(
    arguments,
    "overwrite",
    default = FALSE
  )
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the dispatcher script path.", call. = FALSE)
  }

  dispatcher_file <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  scripts_directory <- dirname(dispatcher_file)
  repo_root <- normalizePath(
    file.path(scripts_directory, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
    stop("Could not locate the MR.rr repository root.", call. = FALSE)
  }

  manifest_argument <- arguments[["manifest"]]
  manifest_file <- if (is.null(manifest_argument) ||
      !nzchar(trimws(manifest_argument))) {
    file.path(repo_root, "paper", "config", "additional_point_tasks.csv")
  } else {
    additional_dispatch_resolve_path(trimws(manifest_argument), repo_root)
  }

  if (!file.exists(manifest_file)) {
    stop("Point-task manifest is missing: ", manifest_file, call. = FALSE)
  }

  manifest <- utils::read.csv(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_columns <- c(
    "manifest_schema_version",
    "task_id",
    "resource_class",
    "array_index",
    "analysis",
    "phase",
    "setting_index",
    "configuration_status",
    "replicate_start",
    "replicate_end",
    "replicate_count",
    "seed_base",
    "mrdag_niter",
    "mrdag_burnin",
    "allow_provisional",
    "worker_script",
    "configuration_file",
    "output_directory",
    "expected_output_file",
    "configuration_md5",
    "core_md5",
    "worker_md5"
  )
  missing_columns <- setdiff(required_columns, names(manifest))

  if (length(missing_columns) > 0L) {
    stop(
      "Point-task manifest is missing columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  manifest_keys <- paste(
    manifest$resource_class,
    manifest$array_index,
    sep = "/"
  )

  if (nrow(manifest) == 0L ||
      anyNA(manifest_keys) ||
      anyDuplicated(manifest_keys) ||
      anyDuplicated(manifest$task_id) ||
      anyDuplicated(manifest$expected_output_file)) {
    stop("Point-task manifest keys are invalid or duplicated.", call. = FALSE)
  }

  selected <- manifest[
    manifest$resource_class == resource_class &
      suppressWarnings(as.integer(manifest$array_index)) == array_index,
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

  schema_version <- additional_dispatch_text_value(
    selected,
    "manifest_schema_version"
  )

  if (!identical(schema_version, "1.0.0")) {
    stop("Unsupported point-task manifest schema: ", schema_version, call. = FALSE)
  }

  task_id <- additional_dispatch_integer_value(selected, "task_id", 1L)
  selected_array_index <- additional_dispatch_integer_value(
    selected,
    "array_index",
    1L
  )
  analysis <- additional_dispatch_text_value(selected, "analysis")
  phase <- additional_dispatch_text_value(selected, "phase")
  setting_index <- additional_dispatch_integer_value(
    selected,
    "setting_index",
    1L
  )
  configuration_status <- additional_dispatch_text_value(
    selected,
    "configuration_status"
  )
  replicate_start <- additional_dispatch_integer_value(
    selected,
    "replicate_start",
    1L
  )
  replicate_end <- additional_dispatch_integer_value(
    selected,
    "replicate_end",
    1L
  )
  replicate_count <- additional_dispatch_integer_value(
    selected,
    "replicate_count",
    1L
  )
  seed_base <- additional_dispatch_integer_value(selected, "seed_base", 1L)
  allow_provisional <- additional_dispatch_boolean_value(
    selected$allow_provisional,
    "Manifest field allow_provisional"
  )

  if (!identical(selected_array_index, array_index)) {
    stop("Selected manifest array index is inconsistent.", call. = FALSE)
  }

  if (!analysis %in% c("rank_misspecification", "approximate_low_rank") ||
      !phase %in% c("standard", "mrdag") ||
      !setting_index %in% seq_len(4L) ||
      replicate_end < replicate_start ||
      replicate_count != replicate_end - replicate_start + 1L) {
    stop("Selected point-task manifest row is internally invalid.", call. = FALSE)
  }

  if ((resource_class == "standard" && phase != "standard") ||
      (resource_class == "mrdag" && phase != "mrdag") ||
      (analysis == "rank_misspecification" && phase != "standard") ||
      (phase == "mrdag" && analysis != "approximate_low_rank")) {
    stop("Analysis, phase, and resource class are inconsistent.", call. = FALSE)
  }

  if (analysis == "rank_misspecification" && allow_provisional) {
    stop("Rank-misspecification tasks cannot be provisional.", call. = FALSE)
  }

  if (configuration_status != "locked") {
    if (analysis != "approximate_low_rank" ||
        !allow_provisional ||
        replicate_count > 5L) {
      stop(
        "A non-locked configuration is allowed only for a small ",
        "approximate-low-rank validation task.",
        call. = FALSE
      )
    }
  }

  mrdag_niter <- additional_dispatch_integer_value(
    selected,
    "mrdag_niter",
    minimum = 1L,
    allow_na = phase != "mrdag"
  )
  mrdag_burnin <- additional_dispatch_integer_value(
    selected,
    "mrdag_burnin",
    minimum = 0L,
    allow_na = phase != "mrdag"
  )

  if (phase == "mrdag" && mrdag_burnin >= mrdag_niter) {
    stop("MrDAG burn-in must be smaller than niter.", call. = FALSE)
  }

  worker_file <- additional_dispatch_resolve_path(
    additional_dispatch_text_value(selected, "worker_script"),
    repo_root
  )
  expected_worker_file <- normalizePath(
    file.path(scripts_directory, "22_run_additional_simulation_chunk.R"),
    winslash = "/",
    mustWork = FALSE
  )
  core_file <- normalizePath(
    file.path(scripts_directory, "20_additional_simulation_core.R"),
    winslash = "/",
    mustWork = FALSE
  )
  configuration_file <- additional_dispatch_resolve_path(
    additional_dispatch_text_value(selected, "configuration_file"),
    repo_root
  )
  expected_configuration_file <- normalizePath(
    file.path(
      repo_root,
      "paper",
      "config",
      paste0(analysis, "_settings.csv")
    ),
    winslash = "/",
    mustWork = FALSE
  )

  if (!identical(worker_file, expected_worker_file)) {
    stop("Manifest worker path is not the expected point worker.", call. = FALSE)
  }

  if (!identical(configuration_file, expected_configuration_file)) {
    stop("Manifest configuration path does not match the analysis.", call. = FALSE)
  }

  additional_dispatch_check_md5(
    worker_file,
    selected$worker_md5,
    "Point worker"
  )
  additional_dispatch_check_md5(
    core_file,
    selected$core_md5,
    "Additional-simulation core"
  )
  additional_dispatch_check_md5(
    configuration_file,
    selected$configuration_md5,
    "Analysis configuration"
  )

  output_directory <- additional_dispatch_resolve_path(
    additional_dispatch_text_value(selected, "output_directory"),
    repo_root
  )
  expected_output_name <- additional_dispatch_text_value(
    selected,
    "expected_output_file"
  )
  canonical_output_name <- sprintf(
    "%s__%s__setting-%d__rep-%04d-%04d.rds",
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end
  )

  if (!identical(expected_output_name, basename(expected_output_name)) ||
      !identical(expected_output_name, canonical_output_name)) {
    stop("Manifest expected output filename is invalid.", call. = FALSE)
  }

  expected_output_file <- file.path(output_directory, expected_output_name)
  worker_arguments <- c(
    "--vanilla",
    shQuote(worker_file),
    paste0("--analysis=", analysis),
    paste0("--phase=", phase),
    paste0("--setting=", setting_index),
    paste0("--start=", replicate_start),
    paste0("--end=", replicate_end),
    paste0("--seed-base=", seed_base),
    paste0("--output-dir=", shQuote(output_directory)),
    paste0("--allow-provisional=", tolower(as.character(allow_provisional))),
    paste0("--overwrite=", tolower(as.character(overwrite)))
  )

  if (phase == "mrdag") {
    worker_arguments <- c(
      worker_arguments,
      paste0("--mrdag-niter=", mrdag_niter),
      paste0("--mrdag-burnin=", mrdag_burnin)
    )
  }

  rscript <- additional_dispatch_find_rscript()

  cat("MR-rr additional point-task dispatcher\n")
  cat("Repository root:", repo_root, "\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Resource class:", resource_class, "\n")
  cat("Array index:", array_index, "\n")
  cat("Global task ID:", task_id, "\n")
  cat("Analysis:", analysis, "\n")
  cat("Phase:", phase, "\n")
  cat("Setting:", setting_index, "\n")
  cat("Replicates:", replicate_start, "through", replicate_end, "\n")
  cat("Configuration status:", configuration_status, "\n")
  cat("Expected output:", expected_output_file, "\n")
  cat(
    "Command:",
    additional_dispatch_format_command(rscript, worker_arguments),
    "\n"
  )

  if (dry_run) {
    cat("Additional point-task dry run: PASS\n")
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
      "Point worker failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }

  additional_dispatch_validate_completed_chunk(
    output_file = expected_output_file,
    selected = selected,
    analysis = analysis,
    phase = phase,
    setting_index = setting_index,
    replicate_start = replicate_start,
    replicate_end = replicate_end,
    seed_base = seed_base,
    allow_provisional = allow_provisional,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin
  )
  cat("Additional point-task output validation: PASS\n")
  cat("Additional point task: PASS\n")
  invisible(selected)
}


if (sys.nframe() == 0L) {
  additional_dispatch_point_task()
}
