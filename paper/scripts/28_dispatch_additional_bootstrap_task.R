#!/usr/bin/env Rscript

# Dispatch one row of the additional-simulation bootstrap-task manifest.
#
# Slurm supplies a resource class and array index. This dispatcher validates
# the manifest row and all locked hashes, invokes the bootstrap chunk worker,
# and validates the completed output before returning success.

options(warn = 1)


additional_bootstrap_dispatch_expected_keys <- function(analysis, phase) {
  if (analysis == "rank_misspecification") {
    return(paste0("r", 1:3, "__regularized_mr_rr"))
  }

  if (phase == "standard") {
    return(c(
      "ivw",
      "srivw",
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr"
    ))
  }

  "mrdag"
}


additional_bootstrap_dispatch_check_input_manifest <- function(
    input_manifest,
    repo_root) {
  required_columns <- c("repository_path", "md5")

  if (!is.data.frame(input_manifest) ||
      !all(required_columns %in% names(input_manifest)) ||
      nrow(input_manifest) == 0L ||
      anyNA(input_manifest[, required_columns]) ||
      any(!nzchar(input_manifest$repository_path)) ||
      any(!grepl("^[0-9a-fA-F]{32}$", input_manifest$md5))) {
    stop("Bootstrap output input provenance is invalid.", call. = FALSE)
  }

  for (row_index in seq_len(nrow(input_manifest))) {
    path <- normalizePath(
      file.path(
        repo_root,
        gsub("\\\\", "/", input_manifest$repository_path[[row_index]])
      ),
      winslash = "/",
      mustWork = FALSE
    )

    if (!file.exists(path)) {
      stop("Bootstrap provenance input is missing: ", path, call. = FALSE)
    }

    observed <- unname(tools::md5sum(path))
    expected <- tolower(input_manifest$md5[[row_index]])

    if (is.na(observed) || tolower(observed) != expected) {
      stop(
        "Bootstrap provenance input has changed: ",
        input_manifest$repository_path[[row_index]],
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}


additional_bootstrap_dispatch_validate_matrix_list <- function(
    values,
    expected_keys,
    effect_dimension,
    replicate_count,
    label,
    predicate = NULL) {
  if (!is.list(values) || !identical(names(values), expected_keys)) {
    stop(label, " result keys are invalid.", call. = FALSE)
  }

  for (key in expected_keys) {
    value <- values[[key]]

    if (!is.matrix(value) ||
        !identical(dim(value), c(effect_dimension, replicate_count)) ||
        any(!is.finite(value))) {
      stop(label, " matrix is invalid for ", key, ".", call. = FALSE)
    }

    if (!is.null(predicate) && !all(predicate(value))) {
      stop(label, " values are invalid for ", key, ".", call. = FALSE)
    }
  }

  invisible(TRUE)
}


additional_bootstrap_dispatch_validate_completed_chunk <- function(
    output_file,
    selected,
    repo_root,
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end,
    seed_base,
    bootstrap_size,
    configured_bootstrap_size,
    bootstrap_size_override,
    requested_cores,
    allow_provisional,
    mrdag_niter,
    mrdag_burnin,
    core,
    worker,
    point_dispatch) {
  if (!file.exists(output_file)) {
    stop(
      "Bootstrap worker succeeded but the expected output is missing: ",
      output_file,
      call. = FALSE
    )
  }

  chunk <- tryCatch(
    readRDS(output_file),
    error = function(error) error
  )

  if (inherits(chunk, "error") || !is.list(chunk)) {
    stop("The completed bootstrap chunk cannot be read.", call. = FALSE)
  }

  metadata <- chunk$metadata
  required_metadata <- c(
    "schema_version",
    "analysis",
    "phase",
    "setting_index",
    "replicate_start",
    "replicate_end",
    "replicate_ids",
    "replicate_count",
    "seed_base",
    "data_seeds",
    "resample_seeds",
    "rng_kind",
    "bootstrap_size",
    "configured_bootstrap_size",
    "bootstrap_size_override",
    "minimum_successful_draws",
    "requested_cores",
    "effective_cores",
    "mrdag_niter",
    "mrdag_burnin",
    "allow_provisional",
    "complete",
    "all_bootstrap_draws_successful",
    "point_result_files_read",
    "archived_simulation_result_rdata_read",
    "frozen_snapshot_written"
  )

  if (!is.list(metadata) ||
      length(setdiff(required_metadata, names(metadata))) > 0L) {
    stop("The completed bootstrap chunk has incomplete metadata.", call. = FALSE)
  }

  expected_ids <- seq.int(replicate_start, replicate_end)
  replicate_count <- length(expected_ids)
  metadata_matches <-
    identical(as.character(metadata$schema_version), "1.0.0") &&
    identical(as.character(metadata$analysis), analysis) &&
    identical(as.character(metadata$phase), phase) &&
    identical(as.integer(metadata$setting_index), setting_index) &&
    identical(as.integer(metadata$replicate_start), replicate_start) &&
    identical(as.integer(metadata$replicate_end), replicate_end) &&
    identical(as.integer(metadata$replicate_ids), expected_ids) &&
    identical(as.integer(metadata$replicate_count), replicate_count) &&
    identical(as.integer(metadata$seed_base), seed_base) &&
    identical(as.integer(metadata$bootstrap_size), bootstrap_size) &&
    identical(
      as.integer(metadata$configured_bootstrap_size),
      configured_bootstrap_size
    ) &&
    identical(
      isTRUE(metadata$bootstrap_size_override),
      bootstrap_size_override
    ) &&
    identical(as.integer(metadata$requested_cores), requested_cores) &&
    identical(isTRUE(metadata$allow_provisional), allow_provisional)

  if (phase == "mrdag") {
    metadata_matches <- metadata_matches &&
      identical(as.integer(metadata$mrdag_niter), mrdag_niter) &&
      identical(as.integer(metadata$mrdag_burnin), mrdag_burnin)
  }

  expected_data_seeds <- vapply(
    expected_ids,
    function(replicate_id) {
      core$additional_make_replicate_seed(
        seed_base = seed_base,
        setting_index = setting_index,
        replicate_id = replicate_id
      )
    },
    integer(1)
  )
  expected_resample_seeds <- vapply(
    expected_data_seeds,
    worker$additional_bootstrap_make_resample_seed,
    integer(1)
  )

  if (!metadata_matches ||
      !identical(as.integer(metadata$data_seeds), expected_data_seeds) ||
      !identical(
        as.integer(metadata$resample_seeds),
        expected_resample_seeds
      ) ||
      !identical(
        as.character(metadata$rng_kind),
        c("L'Ecuyer-CMRG", "Inversion", "Rejection")
      ) ||
      !isTRUE(metadata$complete) ||
      !identical(metadata$point_result_files_read, FALSE) ||
      !identical(metadata$archived_simulation_result_rdata_read, FALSE) ||
      !identical(metadata$frozen_snapshot_written, FALSE)) {
    stop(
      "The completed bootstrap chunk failed metadata validation.",
      call. = FALSE
    )
  }

  expected_keys <- additional_bootstrap_dispatch_expected_keys(
    analysis,
    phase
  )
  specification <- chunk$result_specification

  if (!is.data.frame(specification) ||
      !all(c(
        "result_key",
        "config_row_index",
        "working_rank",
        "method"
      ) %in% names(specification)) ||
      !identical(as.character(specification$result_key), expected_keys)) {
    stop("Bootstrap result specification is invalid.", call. = FALSE)
  }

  truth <- chunk$truth
  effect_dimension <- 27L

  if (!is.list(truth) ||
      !is.matrix(truth$C) ||
      !identical(dim(truth$C), c(3L, 9L)) ||
      any(!is.finite(truth$C)) ||
      length(truth$singular_values) != 3L ||
      any(!is.finite(truth$singular_values)) ||
      !is.numeric(truth$numerical_rank)) {
    stop("Bootstrap truth object is invalid.", call. = FALSE)
  }

  additional_bootstrap_dispatch_validate_matrix_list(
    chunk$standard_errors,
    expected_keys,
    effect_dimension,
    replicate_count,
    "Bootstrap SE",
    predicate = function(value) value >= 0
  )
  additional_bootstrap_dispatch_validate_matrix_list(
    chunk$ci_lower,
    expected_keys,
    effect_dimension,
    replicate_count,
    "Bootstrap lower interval"
  )
  additional_bootstrap_dispatch_validate_matrix_list(
    chunk$ci_upper,
    expected_keys,
    effect_dimension,
    replicate_count,
    "Bootstrap upper interval"
  )
  additional_bootstrap_dispatch_validate_matrix_list(
    chunk$coverage,
    expected_keys,
    effect_dimension,
    replicate_count,
    "Bootstrap coverage",
    predicate = function(value) value %in% c(0, 1)
  )

  for (key in expected_keys) {
    if (any(chunk$ci_lower[[key]] > chunk$ci_upper[[key]])) {
      stop("Bootstrap interval bounds are reversed for ", key, ".", call. = FALSE)
    }

    reconstructed_coverage <- (
      as.vector(truth$C) >= chunk$ci_lower[[key]] &
        as.vector(truth$C) <= chunk$ci_upper[[key]]
    )

    if (!identical(
      as.integer(chunk$coverage[[key]]),
      as.integer(reconstructed_coverage)
    )) {
      stop("Bootstrap coverage cannot be reconstructed for ", key, ".", call. = FALSE)
    }
  }

  successful_draws <- chunk$successful_draws
  minimum_successful <- as.integer(metadata$minimum_successful_draws)

  if (!is.matrix(successful_draws) ||
      !identical(dim(successful_draws), c(length(expected_keys), replicate_count)) ||
      !identical(rownames(successful_draws), expected_keys) ||
      any(!is.finite(successful_draws)) ||
      any(successful_draws != round(successful_draws)) ||
      any(successful_draws < minimum_successful) ||
      any(successful_draws > bootstrap_size)) {
    stop("Bootstrap successful-draw counts are invalid.", call. = FALSE)
  }

  expected_all_successful <- all(successful_draws == bootstrap_size)

  if (!identical(
    isTRUE(metadata$all_bootstrap_draws_successful),
    expected_all_successful
  )) {
    stop("Bootstrap all-successful flag is invalid.", call. = FALSE)
  }

  errors <- chunk$errors

  if (!is.data.frame(errors) ||
      !all(c(
        "replicate_id",
        "bootstrap_id",
        "result_key",
        "message"
      ) %in% names(errors)) ||
      nrow(errors) != sum(bootstrap_size - successful_draws)) {
    stop("Bootstrap estimator-error records are invalid.", call. = FALSE)
  }

  fit_seeds <- chunk$fit_seeds

  if (!is.list(fit_seeds) || !identical(names(fit_seeds), expected_keys)) {
    stop("Bootstrap fit-seed keys are invalid.", call. = FALSE)
  }

  for (specification_index in seq_len(nrow(specification))) {
    key <- specification$result_key[[specification_index]]
    observed <- fit_seeds[[key]]
    expected <- matrix(
      NA_integer_,
      nrow = bootstrap_size,
      ncol = replicate_count
    )

    for (replicate_position in seq_along(expected_ids)) {
      expected[, replicate_position] <- vapply(
        seq_len(bootstrap_size),
        function(bootstrap_id) {
          worker$additional_bootstrap_make_fit_seed(
            data_seed = expected_data_seeds[[replicate_position]],
            method = specification$method[[specification_index]],
            working_rank = specification$working_rank[[specification_index]],
            bootstrap_id = bootstrap_id
          )
        },
        integer(1)
      )
    }

    if (!is.matrix(observed) ||
        !identical(dim(observed), dim(expected)) ||
        !identical(as.integer(observed), as.integer(expected))) {
      stop("Bootstrap fit seeds are invalid for ", key, ".", call. = FALSE)
    }
  }

  additional_bootstrap_dispatch_check_input_manifest(
    chunk$input_manifest,
    repo_root
  )
  expected_output_name <- point_dispatch$additional_dispatch_text_value(
    selected,
    "expected_output_file"
  )

  if (!identical(basename(output_file), expected_output_name)) {
    stop("The completed bootstrap filename differs from the manifest.", call. = FALSE)
  }

  invisible(chunk)
}


additional_dispatch_bootstrap_task <- function() {
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
  point_dispatch_file <- file.path(
    scripts_directory,
    "24_dispatch_additional_point_task.R"
  )

  if (!file.exists(point_dispatch_file)) {
    stop("Point dispatcher dependency is missing.", call. = FALSE)
  }

  point_dispatch <- new.env(parent = globalenv())
  sys.source(point_dispatch_file, envir = point_dispatch)
  arguments <- point_dispatch$additional_dispatch_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "resource-class",
    "array-index",
    "manifest",
    "dry-run",
    "overwrite"
  )
  point_dispatch$additional_dispatch_validate_names(
    arguments,
    allowed_arguments
  )
  resource_class <- arguments[["resource-class"]]

  if (is.null(resource_class) ||
      !resource_class %in% c("standard", "mrdag")) {
    stop("--resource-class must be standard or mrdag.", call. = FALSE)
  }

  array_index <- point_dispatch$additional_dispatch_parse_integer(
    arguments,
    "array-index",
    minimum = 1L
  )
  dry_run <- point_dispatch$additional_dispatch_parse_boolean(
    arguments,
    "dry-run",
    default = FALSE
  )
  overwrite <- point_dispatch$additional_dispatch_parse_boolean(
    arguments,
    "overwrite",
    default = FALSE
  )
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
    file.path(
      repo_root,
      "paper",
      "config",
      "additional_bootstrap_tasks.csv"
    )
  } else {
    point_dispatch$additional_dispatch_resolve_path(
      trimws(manifest_argument),
      repo_root
    )
  }

  if (!file.exists(manifest_file)) {
    stop("Bootstrap-task manifest is missing: ", manifest_file, call. = FALSE)
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
    "bootstrap_size",
    "configured_bootstrap_size",
    "bootstrap_size_override",
    "requested_cores",
    "mrdag_niter",
    "mrdag_burnin",
    "allow_provisional",
    "worker_script",
    "point_worker_script",
    "configuration_file",
    "output_directory",
    "expected_output_file",
    "configuration_md5",
    "core_md5",
    "point_worker_md5",
    "worker_md5"
  )
  missing_columns <- setdiff(required_columns, names(manifest))

  if (length(missing_columns) > 0L) {
    stop(
      "Bootstrap-task manifest is missing columns: ",
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
    stop("Bootstrap-task manifest keys are invalid.", call. = FALSE)
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
          "Expected exactly one bootstrap manifest row for",
          "resource_class=%s and array_index=%d; found %d."
        ),
        resource_class,
        array_index,
        nrow(selected)
      ),
      call. = FALSE
    )
  }

  text_value <- point_dispatch$additional_dispatch_text_value
  integer_value <- point_dispatch$additional_dispatch_integer_value
  boolean_value <- point_dispatch$additional_dispatch_boolean_value
  schema_version <- text_value(selected, "manifest_schema_version")

  if (!identical(schema_version, "1.0.0")) {
    stop(
      "Unsupported bootstrap-task manifest schema: ",
      schema_version,
      call. = FALSE
    )
  }

  task_id <- integer_value(selected, "task_id", 1L)
  selected_array_index <- integer_value(selected, "array_index", 1L)
  analysis <- text_value(selected, "analysis")
  phase <- text_value(selected, "phase")
  setting_index <- integer_value(selected, "setting_index", 1L)
  configuration_status <- text_value(selected, "configuration_status")
  replicate_start <- integer_value(selected, "replicate_start", 1L)
  replicate_end <- integer_value(selected, "replicate_end", 1L)
  replicate_count <- integer_value(selected, "replicate_count", 1L)
  seed_base <- integer_value(selected, "seed_base", 1L)
  bootstrap_size <- integer_value(selected, "bootstrap_size", 2L)
  configured_bootstrap_size <- integer_value(
    selected,
    "configured_bootstrap_size",
    2L
  )
  bootstrap_size_override <- boolean_value(
    selected$bootstrap_size_override,
    "Manifest field bootstrap_size_override"
  )
  requested_cores <- integer_value(selected, "requested_cores", 1L)
  allow_provisional <- boolean_value(
    selected$allow_provisional,
    "Manifest field allow_provisional"
  )

  if (!identical(selected_array_index, array_index) ||
      !analysis %in% c("rank_misspecification", "approximate_low_rank") ||
      !phase %in% c("standard", "mrdag") ||
      !setting_index %in% seq_len(4L) ||
      replicate_end < replicate_start ||
      replicate_count != replicate_end - replicate_start + 1L ||
      bootstrap_size_override !=
        (bootstrap_size != configured_bootstrap_size)) {
    stop("Selected bootstrap manifest row is invalid.", call. = FALSE)
  }

  if ((resource_class == "standard" && phase != "standard") ||
      (resource_class == "mrdag" && phase != "mrdag") ||
      (analysis == "rank_misspecification" && phase != "standard") ||
      (phase == "mrdag" && analysis != "approximate_low_rank")) {
    stop(
      "Analysis, phase, and bootstrap resource class are inconsistent.",
      call. = FALSE
    )
  }

  if (bootstrap_size_override &&
      (bootstrap_size > 20L || replicate_count > 5L)) {
    stop("Bootstrap override is not a small validation task.", call. = FALSE)
  }

  if (analysis == "rank_misspecification" && allow_provisional) {
    stop("Rank-misspecification tasks cannot be provisional.", call. = FALSE)
  }

  if (configuration_status != "locked" &&
      (analysis != "approximate_low_rank" ||
        !allow_provisional ||
        replicate_count > 5L)) {
    stop(
      "Non-locked bootstrap configuration is not a small validation task.",
      call. = FALSE
    )
  }

  mrdag_niter <- integer_value(
    selected,
    "mrdag_niter",
    minimum = 1L,
    allow_na = phase != "mrdag"
  )
  mrdag_burnin <- integer_value(
    selected,
    "mrdag_burnin",
    minimum = 0L,
    allow_na = phase != "mrdag"
  )

  if (phase == "mrdag" && mrdag_burnin >= mrdag_niter) {
    stop("MrDAG burn-in must be smaller than niter.", call. = FALSE)
  }

  resolve_path <- point_dispatch$additional_dispatch_resolve_path
  worker_file <- resolve_path(
    text_value(selected, "worker_script"),
    repo_root
  )
  point_worker_file <- resolve_path(
    text_value(selected, "point_worker_script"),
    repo_root
  )
  core_file <- normalizePath(
    file.path(scripts_directory, "20_additional_simulation_core.R"),
    winslash = "/",
    mustWork = FALSE
  )
  expected_worker_file <- normalizePath(
    file.path(scripts_directory, "26_run_additional_bootstrap_chunk.R"),
    winslash = "/",
    mustWork = FALSE
  )
  expected_point_worker_file <- normalizePath(
    file.path(scripts_directory, "22_run_additional_simulation_chunk.R"),
    winslash = "/",
    mustWork = FALSE
  )
  configuration_file <- resolve_path(
    text_value(selected, "configuration_file"),
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

  if (!identical(worker_file, expected_worker_file) ||
      !identical(point_worker_file, expected_point_worker_file) ||
      !identical(configuration_file, expected_configuration_file)) {
    stop("Bootstrap manifest contains an unexpected input path.", call. = FALSE)
  }

  check_md5 <- point_dispatch$additional_dispatch_check_md5
  check_md5(worker_file, selected$worker_md5, "Bootstrap worker")
  check_md5(
    point_worker_file,
    selected$point_worker_md5,
    "Point worker dependency"
  )
  check_md5(core_file, selected$core_md5, "Additional-simulation core")
  check_md5(
    configuration_file,
    selected$configuration_md5,
    "Analysis configuration"
  )
  core <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  worker <- new.env(parent = globalenv())
  sys.source(worker_file, envir = worker)
  output_directory <- resolve_path(
    text_value(selected, "output_directory"),
    repo_root
  )
  expected_output_name <- text_value(selected, "expected_output_file")
  canonical_output_name <- sprintf(
    "%s__%s__setting-%d__rep-%04d-%04d__bt-%04d.rds",
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end,
    bootstrap_size
  )

  if (!identical(expected_output_name, basename(expected_output_name)) ||
      !identical(expected_output_name, canonical_output_name)) {
    stop("Bootstrap manifest output filename is invalid.", call. = FALSE)
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
    paste0("--bootstrap-size=", bootstrap_size),
    paste0(
      "--allow-bootstrap-override=",
      tolower(as.character(bootstrap_size_override))
    ),
    paste0("--n-cores=", requested_cores),
    paste0("--output-dir=", shQuote(output_directory)),
    paste0(
      "--allow-provisional=",
      tolower(as.character(allow_provisional))
    ),
    paste0("--overwrite=", tolower(as.character(overwrite)))
  )

  if (phase == "mrdag") {
    worker_arguments <- c(
      worker_arguments,
      paste0("--mrdag-niter=", mrdag_niter),
      paste0("--mrdag-burnin=", mrdag_burnin)
    )
  }

  rscript <- point_dispatch$additional_dispatch_find_rscript()

  cat("MR-rr additional bootstrap-task dispatcher\n")
  cat("Repository root:", repo_root, "\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Resource class:", resource_class, "\n")
  cat("Array index:", array_index, "\n")
  cat("Global task ID:", task_id, "\n")
  cat("Analysis:", analysis, "\n")
  cat("Phase:", phase, "\n")
  cat("Setting:", setting_index, "\n")
  cat("Replicates:", replicate_start, "through", replicate_end, "\n")
  cat("Bootstrap draws:", bootstrap_size, "\n")
  cat("Requested cores:", requested_cores, "\n")
  cat("Configuration status:", configuration_status, "\n")
  cat("Expected output:", expected_output_file, "\n")
  cat(
    "Command:",
    point_dispatch$additional_dispatch_format_command(
      rscript,
      worker_arguments
    ),
    "\n"
  )

  if (dry_run) {
    cat("Additional bootstrap-task dry run: PASS\n")
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
      "Bootstrap worker failed with exit status ",
      status,
      ".",
      call. = FALSE
    )
  }

  additional_bootstrap_dispatch_validate_completed_chunk(
    output_file = expected_output_file,
    selected = selected,
    repo_root = repo_root,
    analysis = analysis,
    phase = phase,
    setting_index = setting_index,
    replicate_start = replicate_start,
    replicate_end = replicate_end,
    seed_base = seed_base,
    bootstrap_size = bootstrap_size,
    configured_bootstrap_size = configured_bootstrap_size,
    bootstrap_size_override = bootstrap_size_override,
    requested_cores = requested_cores,
    allow_provisional = allow_provisional,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin,
    core = core,
    worker = worker,
    point_dispatch = point_dispatch
  )
  cat("Additional bootstrap-task output validation: PASS\n")
  cat("Additional bootstrap task: PASS\n")
  invisible(selected)
}


if (sys.nframe() == 0L) {
  additional_dispatch_bootstrap_task()
}
