#!/usr/bin/env Rscript

# Build the task manifest for additional-simulation bootstrap inference.
#
# The standard resource class contains rank misspecification and approximate
# low-rank non-MrDAG tasks. The MrDAG class contains only approximate-low-rank
# MrDAG tasks. A formal manifest requires the approximate delta to be locked.

options(warn = 1)


additional_bootstrap_manifest_make_block <- function(
    analysis,
    phase,
    resource_class,
    config,
    ranges,
    bootstrap_size,
    bootstrap_size_override,
    requested_cores,
    output_directory,
    worker_path,
    point_worker_path,
    config_path,
    config_md5,
    core_md5,
    point_worker_md5,
    worker_md5,
    allow_provisional,
    mrdag_niter,
    mrdag_burnin,
    point_manifest) {
  rows <- list()
  row_index <- 0L

  for (setting_index in seq_len(4L)) {
    setting <- point_manifest$additional_manifest_setting_row(
      config,
      setting_index
    )

    for (range_index in seq_len(nrow(ranges))) {
      range <- ranges[range_index, , drop = FALSE]
      row_index <- row_index + 1L
      expected_output_file <- sprintf(
        "%s__%s__setting-%d__rep-%04d-%04d__bt-%04d.rds",
        analysis,
        phase,
        setting_index,
        range$replicate_start,
        range$replicate_end,
        bootstrap_size
      )
      rows[[row_index]] <- data.frame(
        manifest_schema_version = "1.0.0",
        task_id = NA_integer_,
        resource_class = resource_class,
        array_index = NA_integer_,
        analysis = analysis,
        phase = phase,
        setting_index = setting_index,
        scenario = setting$scenario,
        measurement_error_weight = setting$measurement_error_weight,
        genetic_effect_weight = setting$genetic_effect_weight,
        third_singular_value = setting$third_singular_value,
        configuration_status = setting$status,
        chunk_index = range$chunk_index,
        replicate_start = range$replicate_start,
        replicate_end = range$replicate_end,
        replicate_count = range$replicate_count,
        seed_base = setting$seed_base,
        bootstrap_size = bootstrap_size,
        configured_bootstrap_size = setting$bootstrap_size,
        bootstrap_size_override = bootstrap_size_override,
        requested_cores = requested_cores,
        mrdag_niter = if (phase == "mrdag") {
          mrdag_niter
        } else {
          NA_integer_
        },
        mrdag_burnin = if (phase == "mrdag") {
          mrdag_burnin
        } else {
          NA_integer_
        },
        allow_provisional = allow_provisional,
        worker_script = worker_path,
        point_worker_script = point_worker_path,
        configuration_file = config_path,
        output_directory = output_directory,
        expected_output_file = expected_output_file,
        configuration_md5 = config_md5,
        core_md5 = core_md5,
        point_worker_md5 = point_worker_md5,
        worker_md5 = worker_md5,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  do.call(rbind, rows)
}


additional_bootstrap_manifest_validate_coverage <- function(
    manifest,
    total_replicates) {
  group_keys <- unique(
    manifest[, c("analysis", "phase", "setting_index")]
  )
  rownames(group_keys) <- NULL

  if (nrow(group_keys) != 12L) {
    stop(
      "Bootstrap manifest does not contain the expected 12 groups.",
      call. = FALSE
    )
  }

  for (group_index in seq_len(nrow(group_keys))) {
    key <- group_keys[group_index, , drop = FALSE]
    selected <- manifest[
      manifest$analysis == key$analysis &
        manifest$phase == key$phase &
        manifest$setting_index == key$setting_index,
      ,
      drop = FALSE
    ]
    replicate_ids <- unlist(Map(
      seq.int,
      selected$replicate_start,
      selected$replicate_end
    ))

    if (!identical(
      sort(as.integer(replicate_ids)),
      seq_len(total_replicates)
    )) {
      stop(
        "Bootstrap manifest replicate coverage failed for ",
        key$analysis,
        "/",
        key$phase,
        "/setting=",
        key$setting_index,
        ".",
        call. = FALSE
      )
    }
  }

  invisible(group_keys)
}


additional_make_bootstrap_manifest <- function() {
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the manifest-generator path.", call. = FALSE)
  }

  script_path <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  scripts_directory <- dirname(script_path)
  point_manifest_file <- file.path(
    scripts_directory,
    "23_make_additional_point_manifest.R"
  )
  core_file <- file.path(
    scripts_directory,
    "20_additional_simulation_core.R"
  )
  point_worker_file <- file.path(
    scripts_directory,
    "22_run_additional_simulation_chunk.R"
  )
  worker_file <- file.path(
    scripts_directory,
    "26_run_additional_bootstrap_chunk.R"
  )
  required_scripts <- c(
    point_manifest_file,
    core_file,
    point_worker_file,
    worker_file
  )
  missing_scripts <- required_scripts[!file.exists(required_scripts)]

  if (length(missing_scripts) > 0L) {
    stop(
      "Bootstrap manifest dependencies are missing: ",
      paste(missing_scripts, collapse = ", "),
      call. = FALSE
    )
  }

  point_manifest <- new.env(parent = globalenv())
  sys.source(point_manifest_file, envir = point_manifest)
  arguments <- point_manifest$additional_manifest_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "total-replicates",
    "standard-chunk-size",
    "mrdag-chunk-size",
    "bootstrap-size",
    "allow-bootstrap-override",
    "standard-cores",
    "mrdag-cores",
    "mrdag-niter",
    "mrdag-burnin",
    "output",
    "output-directory",
    "allow-provisional",
    "overwrite"
  )
  point_manifest$additional_manifest_validate_names(
    arguments,
    allowed_arguments
  )
  total_replicates <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "total-replicates",
    default = 1000L
  )
  standard_chunk_size <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "standard-chunk-size",
    default = 10L
  )
  mrdag_chunk_size <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "mrdag-chunk-size",
    default = 10L
  )
  requested_bootstrap_size <- if (is.null(arguments[["bootstrap-size"]])) {
    NULL
  } else {
    point_manifest$additional_manifest_parse_integer(
      arguments,
      "bootstrap-size",
      minimum = 2L
    )
  }
  allow_bootstrap_override <-
    point_manifest$additional_manifest_parse_boolean(
      arguments,
      "allow-bootstrap-override",
      default = FALSE
    )
  standard_cores <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "standard-cores",
    default = 16L
  )
  mrdag_cores <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "mrdag-cores",
    default = 16L
  )
  mrdag_niter <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "mrdag-niter",
    default = 1000L
  )
  mrdag_burnin <- point_manifest$additional_manifest_parse_integer(
    arguments,
    "mrdag-burnin",
    default = 200L,
    minimum = 0L
  )
  allow_provisional <- point_manifest$additional_manifest_parse_boolean(
    arguments,
    "allow-provisional",
    default = FALSE
  )
  overwrite <- point_manifest$additional_manifest_parse_boolean(
    arguments,
    "overwrite",
    default = FALSE
  )

  if (mrdag_burnin >= mrdag_niter) {
    stop("--mrdag-burnin must be smaller than --mrdag-niter.", call. = FALSE)
  }

  core <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  repo_root <- core$additional_find_repo_root(dirname(script_path))
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
  rank_config <- core$additional_read_config(
    rank_config_file,
    "rank_misspecification"
  )
  approximate_config <- core$additional_read_config(
    approximate_config_file,
    "approximate_low_rank"
  )
  configured_total <- unique(c(
    rank_config$monte_carlo_replicates,
    approximate_config$monte_carlo_replicates
  ))

  if (length(configured_total) != 1L ||
      total_replicates > configured_total) {
    stop(
      "--total-replicates exceeds the common configured total.",
      call. = FALSE
    )
  }

  configured_bootstrap_size <- unique(c(
    rank_config$bootstrap_size,
    approximate_config$bootstrap_size
  ))

  if (length(configured_bootstrap_size) != 1L) {
    stop("Configured bootstrap sizes are inconsistent.", call. = FALSE)
  }

  bootstrap_size <- if (is.null(requested_bootstrap_size)) {
    configured_bootstrap_size
  } else {
    requested_bootstrap_size
  }
  bootstrap_size_override <- bootstrap_size != configured_bootstrap_size

  if (bootstrap_size_override && !allow_bootstrap_override) {
    stop(
      paste0(
        "--bootstrap-size differs from the configured value ",
        configured_bootstrap_size,
        "; use --allow-bootstrap-override=true only for a smoke test."
      ),
      call. = FALSE
    )
  }

  if (allow_bootstrap_override && bootstrap_size > 20L) {
    stop(
      "A bootstrap-size override is limited to at most 20 draws.",
      call. = FALSE
    )
  }

  approximate_is_provisional <- any(approximate_config$status != "locked")

  if (approximate_is_provisional) {
    if (!allow_provisional) {
      stop(
        paste(
          "Approximate-low-rank delta is still provisional.",
          "A formal bootstrap manifest cannot be generated."
        ),
        call. = FALSE
      )
    }

    if (total_replicates > 5L) {
      stop(
        "A provisional bootstrap manifest is limited to 5 replicates.",
        call. = FALSE
      )
    }
  }

  if ((bootstrap_size_override || approximate_is_provisional) &&
      total_replicates > 5L) {
    stop(
      "Validation manifests are limited to at most 5 replicates.",
      call. = FALSE
    )
  }

  output_argument <- arguments[["output"]]
  output_file <- if (is.null(output_argument) ||
      !nzchar(trimws(output_argument))) {
    file.path(
      repo_root,
      "paper",
      "config",
      "additional_bootstrap_tasks.csv"
    )
  } else {
    point_manifest$additional_manifest_resolve_path(
      trimws(output_argument),
      repo_root
    )
  }
  task_output_argument <- arguments[["output-directory"]]
  task_output_directory <- if (is.null(task_output_argument) ||
      !nzchar(trimws(task_output_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "bootstrap_chunks"
    )
  } else {
    point_manifest$additional_manifest_resolve_path(
      trimws(task_output_argument),
      repo_root
    )
  }

  if (file.exists(output_file) && !overwrite) {
    stop(
      "Manifest already exists; use --overwrite=true to replace it: ",
      output_file,
      call. = FALSE
    )
  }

  standard_ranges <- point_manifest$additional_manifest_make_ranges(
    total_replicates,
    standard_chunk_size
  )
  mrdag_ranges <- point_manifest$additional_manifest_make_ranges(
    total_replicates,
    mrdag_chunk_size
  )
  worker_path <- point_manifest$additional_manifest_repo_path(
    worker_file,
    repo_root
  )
  point_worker_path <- point_manifest$additional_manifest_repo_path(
    point_worker_file,
    repo_root
  )
  task_output_path <- point_manifest$additional_manifest_repo_path(
    task_output_directory,
    repo_root
  )
  core_md5 <- unname(tools::md5sum(core_file))
  point_worker_md5 <- unname(tools::md5sum(point_worker_file))
  worker_md5 <- unname(tools::md5sum(worker_file))
  blocks <- list(
    additional_bootstrap_manifest_make_block(
      analysis = "rank_misspecification",
      phase = "standard",
      resource_class = "standard",
      config = rank_config,
      ranges = standard_ranges,
      bootstrap_size = bootstrap_size,
      bootstrap_size_override = bootstrap_size_override,
      requested_cores = standard_cores,
      output_directory = task_output_path,
      worker_path = worker_path,
      point_worker_path = point_worker_path,
      config_path = point_manifest$additional_manifest_repo_path(
        rank_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(rank_config_file)),
      core_md5 = core_md5,
      point_worker_md5 = point_worker_md5,
      worker_md5 = worker_md5,
      allow_provisional = FALSE,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin,
      point_manifest = point_manifest
    ),
    additional_bootstrap_manifest_make_block(
      analysis = "approximate_low_rank",
      phase = "standard",
      resource_class = "standard",
      config = approximate_config,
      ranges = standard_ranges,
      bootstrap_size = bootstrap_size,
      bootstrap_size_override = bootstrap_size_override,
      requested_cores = standard_cores,
      output_directory = task_output_path,
      worker_path = worker_path,
      point_worker_path = point_worker_path,
      config_path = point_manifest$additional_manifest_repo_path(
        approximate_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(approximate_config_file)),
      core_md5 = core_md5,
      point_worker_md5 = point_worker_md5,
      worker_md5 = worker_md5,
      allow_provisional = allow_provisional,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin,
      point_manifest = point_manifest
    ),
    additional_bootstrap_manifest_make_block(
      analysis = "approximate_low_rank",
      phase = "mrdag",
      resource_class = "mrdag",
      config = approximate_config,
      ranges = mrdag_ranges,
      bootstrap_size = bootstrap_size,
      bootstrap_size_override = bootstrap_size_override,
      requested_cores = mrdag_cores,
      output_directory = task_output_path,
      worker_path = worker_path,
      point_worker_path = point_worker_path,
      config_path = point_manifest$additional_manifest_repo_path(
        approximate_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(approximate_config_file)),
      core_md5 = core_md5,
      point_worker_md5 = point_worker_md5,
      worker_md5 = worker_md5,
      allow_provisional = allow_provisional,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin,
      point_manifest = point_manifest
    )
  )
  manifest <- do.call(rbind, blocks)
  rownames(manifest) <- NULL
  manifest$task_id <- seq_len(nrow(manifest))
  manifest$array_index <- ave(
    manifest$task_id,
    manifest$resource_class,
    FUN = seq_along
  )
  manifest$array_index <- as.integer(manifest$array_index)
  manifest <- manifest[, c(
    "manifest_schema_version",
    "task_id",
    "resource_class",
    "array_index",
    setdiff(
      names(manifest),
      c(
        "manifest_schema_version",
        "task_id",
        "resource_class",
        "array_index"
      )
    )
  )]
  task_keys <- paste(
    manifest$resource_class,
    manifest$array_index,
    sep = "/"
  )

  if (anyDuplicated(task_keys) ||
      anyDuplicated(manifest$expected_output_file) ||
      any(manifest$replicate_count < 1L) ||
      any(manifest$replicate_start > manifest$replicate_end) ||
      any(manifest$bootstrap_size != bootstrap_size) ||
      any(manifest$requested_cores < 1L)) {
    stop(
      "Generated bootstrap-task manifest is internally invalid.",
      call. = FALSE
    )
  }

  additional_bootstrap_manifest_validate_coverage(
    manifest,
    total_replicates
  )
  point_manifest$additional_manifest_write_csv(manifest, output_file)
  resource_counts <- table(manifest$resource_class)

  cat("MR-rr additional bootstrap-task manifest\n")
  cat("Repository root:", repo_root, "\n")
  cat("Total replicates per setting:", total_replicates, "\n")
  cat("Bootstrap draws per replicate:", bootstrap_size, "\n")
  cat("Configured bootstrap draws:", configured_bootstrap_size, "\n")
  cat("Standard chunk size:", standard_chunk_size, "\n")
  cat("MrDAG chunk size:", mrdag_chunk_size, "\n")
  cat("Standard cores per task:", standard_cores, "\n")
  cat("MrDAG cores per task:", mrdag_cores, "\n")
  cat(
    "Approximate delta:",
    unique(approximate_config$third_singular_value),
    "\n"
  )
  cat(
    "Approximate configuration status:",
    paste(unique(approximate_config$status), collapse = ", "),
    "\n"
  )
  cat("Total tasks:", nrow(manifest), "\n")
  cat("Standard tasks:", unname(resource_counts[["standard"]]), "\n")
  cat("MrDAG tasks:", unname(resource_counts[["mrdag"]]), "\n")
  cat("Manifest:", output_file, "\n")
  cat("Bootstrap-task manifest validation: PASS\n")
  invisible(manifest)
}


if (sys.nframe() == 0L) {
  additional_make_bootstrap_manifest()
}
