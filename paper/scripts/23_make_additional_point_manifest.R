#!/usr/bin/env Rscript

# Build the task manifest for additional-simulation point estimates.
#
# Standard and MrDAG tasks receive independent Slurm array indices. A full
# manifest cannot be generated while approximate-low-rank delta is provisional.

options(warn = 1)


additional_manifest_parse_arguments <- function(arguments) {
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


additional_manifest_validate_names <- function(arguments, allowed) {
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


additional_manifest_parse_integer <- function(
    arguments,
    name,
    default,
    minimum = 1L) {
  value <- arguments[[name]]

  if (is.null(value)) {
    value <- default
  }

  numeric_value <- suppressWarnings(as.numeric(value))

  if (length(numeric_value) != 1L ||
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


additional_manifest_parse_boolean <- function(
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


additional_manifest_resolve_path <- function(path, repo_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(
    file.path(repo_root, path),
    winslash = "/",
    mustWork = FALSE
  )
}


additional_manifest_repo_path <- function(path, repo_root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  normalized_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(normalized_root, "/")

  if (startsWith(normalized_path, prefix)) {
    return(substring(normalized_path, nchar(prefix) + 1L))
  }

  normalized_path
}


additional_manifest_make_ranges <- function(total, chunk_size) {
  starts <- seq.int(1L, total, by = chunk_size)
  ends <- pmin(starts + chunk_size - 1L, total)

  data.frame(
    chunk_index = seq_along(starts),
    replicate_start = starts,
    replicate_end = ends,
    replicate_count = ends - starts + 1L,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


additional_manifest_setting_row <- function(config, setting_index) {
  selected <- config[
    config$setting_index == setting_index,
    ,
    drop = FALSE
  ]

  if (nrow(selected) == 0L) {
    stop("Setting is absent from configuration: ", setting_index, call. = FALSE)
  }

  fields <- c(
    "scenario",
    "measurement_error_weight",
    "genetic_effect_weight",
    "third_singular_value",
    "monte_carlo_replicates",
    "bootstrap_size",
    "seed_base",
    "data_seed_group",
    "status"
  )

  for (field in fields) {
    if (length(unique(selected[[field]])) != 1L) {
      stop(
        "Configuration field varies within setting ",
        setting_index,
        ": ",
        field,
        call. = FALSE
      )
    }
  }

  selected[1L, , drop = FALSE]
}


additional_manifest_make_block <- function(
    analysis,
    phase,
    resource_class,
    config,
    ranges,
    output_directory,
    worker_path,
    config_path,
    config_md5,
    core_md5,
    worker_md5,
    allow_provisional,
    mrdag_niter,
    mrdag_burnin) {
  rows <- list()
  row_index <- 0L

  for (setting_index in seq_len(4L)) {
    setting <- additional_manifest_setting_row(config, setting_index)

    for (range_index in seq_len(nrow(ranges))) {
      range <- ranges[range_index, , drop = FALSE]
      row_index <- row_index + 1L
      expected_output_file <- sprintf(
        "%s__%s__setting-%d__rep-%04d-%04d.rds",
        analysis,
        phase,
        setting_index,
        range$replicate_start,
        range$replicate_end
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
        configuration_file = config_path,
        output_directory = output_directory,
        expected_output_file = expected_output_file,
        configuration_md5 = config_md5,
        core_md5 = core_md5,
        worker_md5 = worker_md5,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  do.call(rbind, rows)
}


additional_manifest_write_csv <- function(value, path) {
  directory <- dirname(path)
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = directory,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  utils::write.csv(
    value,
    temporary_file,
    row.names = FALSE,
    na = ""
  )

  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = directory,
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop("Could not preserve existing manifest: ", path, call. = FALSE)
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop("Could not install completed manifest: ", path, call. = FALSE)
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


additional_make_point_manifest <- function() {
  arguments <- additional_manifest_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "total-replicates",
    "standard-chunk-size",
    "mrdag-chunk-size",
    "mrdag-niter",
    "mrdag-burnin",
    "output",
    "output-directory",
    "allow-provisional",
    "overwrite"
  )
  additional_manifest_validate_names(arguments, allowed_arguments)
  total_replicates <- additional_manifest_parse_integer(
    arguments,
    "total-replicates",
    default = 1000L
  )
  standard_chunk_size <- additional_manifest_parse_integer(
    arguments,
    "standard-chunk-size",
    default = 100L
  )
  mrdag_chunk_size <- additional_manifest_parse_integer(
    arguments,
    "mrdag-chunk-size",
    default = 25L
  )
  mrdag_niter <- additional_manifest_parse_integer(
    arguments,
    "mrdag-niter",
    default = 1000L
  )
  mrdag_burnin <- additional_manifest_parse_integer(
    arguments,
    "mrdag-burnin",
    default = 200L,
    minimum = 0L
  )
  allow_provisional <- additional_manifest_parse_boolean(
    arguments,
    "allow-provisional",
    default = FALSE
  )
  overwrite <- additional_manifest_parse_boolean(
    arguments,
    "overwrite",
    default = FALSE
  )

  if (mrdag_burnin >= mrdag_niter) {
    stop("--mrdag-burnin must be smaller than --mrdag-niter.", call. = FALSE)
  }

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
  core_file <- file.path(dirname(script_path), "20_additional_simulation_core.R")
  worker_file <- file.path(
    dirname(script_path),
    "22_run_additional_simulation_chunk.R"
  )

  if (!file.exists(core_file) || !file.exists(worker_file)) {
    stop("Additional simulation core or point worker is missing.", call. = FALSE)
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

  if (length(configured_total) != 1L || total_replicates > configured_total) {
    stop(
      "--total-replicates exceeds the common configured total.",
      call. = FALSE
    )
  }

  approximate_is_provisional <- any(approximate_config$status != "locked")

  if (approximate_is_provisional) {
    if (!allow_provisional) {
      stop(
        paste(
          "Approximate-low-rank delta is still provisional.",
          "A full point manifest cannot be generated."
        ),
        call. = FALSE
      )
    }

    if (total_replicates > 5L) {
      stop(
        "A provisional manifest is limited to at most 5 replicates.",
        call. = FALSE
      )
    }
  }

  output_argument <- arguments[["output"]]
  output_file <- if (is.null(output_argument) ||
      !nzchar(trimws(output_argument))) {
    file.path(
      repo_root,
      "paper",
      "config",
      "additional_point_tasks.csv"
    )
  } else {
    additional_manifest_resolve_path(trimws(output_argument), repo_root)
  }
  task_output_argument <- arguments[["output-directory"]]
  task_output_directory <- if (is.null(task_output_argument) ||
      !nzchar(trimws(task_output_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "point_chunks"
    )
  } else {
    additional_manifest_resolve_path(
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

  standard_ranges <- additional_manifest_make_ranges(
    total_replicates,
    standard_chunk_size
  )
  mrdag_ranges <- additional_manifest_make_ranges(
    total_replicates,
    mrdag_chunk_size
  )
  worker_path <- additional_manifest_repo_path(worker_file, repo_root)
  task_output_path <- additional_manifest_repo_path(
    task_output_directory,
    repo_root
  )
  core_md5 <- unname(tools::md5sum(core_file))
  worker_md5 <- unname(tools::md5sum(worker_file))
  blocks <- list(
    additional_manifest_make_block(
      analysis = "rank_misspecification",
      phase = "standard",
      resource_class = "standard",
      config = rank_config,
      ranges = standard_ranges,
      output_directory = task_output_path,
      worker_path = worker_path,
      config_path = additional_manifest_repo_path(
        rank_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(rank_config_file)),
      core_md5 = core_md5,
      worker_md5 = worker_md5,
      allow_provisional = FALSE,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin
    ),
    additional_manifest_make_block(
      analysis = "approximate_low_rank",
      phase = "standard",
      resource_class = "standard",
      config = approximate_config,
      ranges = standard_ranges,
      output_directory = task_output_path,
      worker_path = worker_path,
      config_path = additional_manifest_repo_path(
        approximate_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(approximate_config_file)),
      core_md5 = core_md5,
      worker_md5 = worker_md5,
      allow_provisional = allow_provisional,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin
    ),
    additional_manifest_make_block(
      analysis = "approximate_low_rank",
      phase = "mrdag",
      resource_class = "mrdag",
      config = approximate_config,
      ranges = mrdag_ranges,
      output_directory = task_output_path,
      worker_path = worker_path,
      config_path = additional_manifest_repo_path(
        approximate_config_file,
        repo_root
      ),
      config_md5 = unname(tools::md5sum(approximate_config_file)),
      core_md5 = core_md5,
      worker_md5 = worker_md5,
      allow_provisional = allow_provisional,
      mrdag_niter = mrdag_niter,
      mrdag_burnin = mrdag_burnin
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
      any(manifest$replicate_start > manifest$replicate_end)) {
    stop("Generated point-task manifest is internally invalid.", call. = FALSE)
  }

  for (analysis in unique(manifest$analysis)) {
    for (phase in unique(manifest$phase[manifest$analysis == analysis])) {
      selected <- manifest[
        manifest$analysis == analysis & manifest$phase == phase,
        ,
        drop = FALSE
      ]

      for (setting_index in seq_len(4L)) {
        replicate_ids <- unlist(Map(
          seq.int,
          selected$replicate_start[selected$setting_index == setting_index],
          selected$replicate_end[selected$setting_index == setting_index]
        ))

        if (!identical(sort(as.integer(replicate_ids)), seq_len(total_replicates))) {
          stop(
            "Manifest replicate coverage failed for ",
            analysis,
            "/",
            phase,
            "/setting=",
            setting_index,
            ".",
            call. = FALSE
          )
        }
      }
    }
  }

  additional_manifest_write_csv(manifest, output_file)
  resource_counts <- table(manifest$resource_class)

  cat("MR-rr additional point-task manifest\n")
  cat("Repository root:", repo_root, "\n")
  cat("Total replicates per setting:", total_replicates, "\n")
  cat("Standard chunk size:", standard_chunk_size, "\n")
  cat("MrDAG chunk size:", mrdag_chunk_size, "\n")
  cat("Approximate delta:", unique(approximate_config$third_singular_value), "\n")
  cat(
    "Approximate configuration status:",
    paste(unique(approximate_config$status), collapse = ", "),
    "\n"
  )
  cat("Total tasks:", nrow(manifest), "\n")
  cat("Standard tasks:", unname(resource_counts[["standard"]]), "\n")
  cat("MrDAG tasks:", unname(resource_counts[["mrdag"]]), "\n")
  cat("Manifest:", output_file, "\n")
  cat("Point-task manifest validation: PASS\n")
  invisible(manifest)
}


if (sys.nframe() == 0L) {
  additional_make_point_manifest()
}
