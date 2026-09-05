#!/usr/bin/env Rscript

# Build the complete task manifest for the from-scratch MR-rr simulations.
#
# The manifest separates standard tasks from MrDAG tasks so that Slurm can use
# different time limits and concurrency caps. It contains no cluster-specific
# paths and can also be consumed by a non-Slurm scheduler.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "paper", "scripts")) &&
        file.exists(file.path(
          current,
          "paper",
          "scripts",
          "16_run_simulation_chunk.R"
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


parse_integer_argument <- function(
    arguments,
    name,
    default,
    minimum = 1L) {
  value <- arguments[[name]]

  if (is.null(value)) {
    value <- default
  }

  parsed <- suppressWarnings(as.integer(value))

  if (length(parsed) != 1L || is.na(parsed) || parsed < minimum) {
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


make_chunk_ranges <- function(total_replicates, chunk_size) {
  starts <- seq.int(
    from = 1L,
    to = total_replicates,
    by = chunk_size
  )
  ends <- pmin(
    starts + chunk_size - 1L,
    total_replicates
  )

  data.frame(
    chunk_index = seq_along(starts),
    replicate_start = starts,
    replicate_end = ends,
    replicate_count = ends - starts + 1L,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


make_task_block <- function(
    resource_class,
    designs,
    phases,
    settings,
    chunk_ranges,
    seed_base,
    mrdag_niter,
    mrdag_burnin) {
  rows <- list()
  row_index <- 0L

  for (design in designs) {
    for (phase in phases) {
      for (setting_index in settings$setting_index) {
        setting <- settings[
          settings$setting_index == setting_index,
          ,
          drop = FALSE
        ]

        for (chunk_row in seq_len(nrow(chunk_ranges))) {
          chunk <- chunk_ranges[chunk_row, , drop = FALSE]
          row_index <- row_index + 1L

          expected_output_file <- sprintf(
            "%s__%s__setting-%d__rep-%04d-%04d.rds",
            design,
            phase,
            setting_index,
            chunk$replicate_start,
            chunk$replicate_end
          )

          rows[[row_index]] <- data.frame(
            task_id = NA_integer_,
            resource_class = resource_class,
            array_index = NA_integer_,
            design = design,
            phase = phase,
            setting = setting_index,
            scenario = setting$scenario,
            me_weight = setting$me_weight,
            effect_weight = setting$effect_weight,
            chunk_index = chunk$chunk_index,
            replicate_start = chunk$replicate_start,
            replicate_end = chunk$replicate_end,
            replicate_count = chunk$replicate_count,
            seed_base = seed_base,
            mrdag_niter = if (phase == "main_mrdag") {
              mrdag_niter
            } else {
              NA_integer_
            },
            mrdag_burnin = if (phase == "main_mrdag") {
              mrdag_burnin
            } else {
              NA_integer_
            },
            worker_script = file.path(
              "paper",
              "scripts",
              "16_run_simulation_chunk.R"
            ),
            output_directory = file.path(
              "paper",
              "output",
              "full_run",
              "chunks"
            ),
            expected_output_file = expected_output_file,
            stringsAsFactors = FALSE,
            row.names = NULL
          )
        }
      }
    }
  }

  task_block <- do.call(rbind, rows)
  rownames(task_block) <- NULL
  task_block$array_index <- seq_len(nrow(task_block))
  task_block
}


validate_task_coverage <- function(
    manifest,
    total_replicates,
    designs,
    phases,
    settings) {
  expected_groups <- expand.grid(
    design = designs,
    phase = phases,
    setting = settings$setting_index,
    stringsAsFactors = FALSE
  )

  for (group_index in seq_len(nrow(expected_groups))) {
    group <- expected_groups[group_index, , drop = FALSE]
    selected <- manifest[
      manifest$design == group$design &
      manifest$phase == group$phase &
      manifest$setting == group$setting,
      ,
      drop = FALSE
    ]

    if (nrow(selected) == 0L) {
      stop(
        "Task manifest is missing group: ",
        paste(
          group$design,
          group$phase,
          group$setting,
          sep = "/"
        ),
        call. = FALSE
      )
    }

    covered_replicates <- unlist(
      Map(
        seq.int,
        selected$replicate_start,
        selected$replicate_end
      ),
      use.names = FALSE
    )

    if (!identical(
      sort(as.integer(covered_replicates)),
      seq_len(total_replicates)
    )) {
      stop(
        "Task manifest has a gap, overlap, or out-of-range replicate in ",
        paste(
          group$design,
          group$phase,
          group$setting,
          sep = "/"
        ),
        ".",
        call. = FALSE
      )
    }
  }

  actual_groups <- unique(
    manifest[, c("design", "phase", "setting")]
  )

  if (nrow(actual_groups) != nrow(expected_groups)) {
    stop(
      "Task manifest contains unexpected design/phase/setting groups.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


write_csv_atomically <- function(data, path) {
  output_directory <- dirname(path)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = output_directory,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)

  utils::write.csv(
    data,
    temporary_file,
    row.names = FALSE,
    na = ""
  )

  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = output_directory,
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop(
        "Could not preserve the existing manifest before replacement.",
        call. = FALSE
      )
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop(
      "Could not move the completed manifest to: ",
      path,
      call. = FALSE
    )
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


make_full_simulation_manifest <- function() {
  arguments <- parse_named_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "total-replicates",
    "standard-chunk-size",
    "mrdag-chunk-size",
    "seed-base",
    "mrdag-niter",
    "mrdag-burnin",
    "output",
    "overwrite"
  )
  validate_argument_names(arguments, allowed_arguments)

  total_replicates <- parse_integer_argument(
    arguments,
    "total-replicates",
    default = 1000L
  )
  standard_chunk_size <- parse_integer_argument(
    arguments,
    "standard-chunk-size",
    default = 100L
  )
  mrdag_chunk_size <- parse_integer_argument(
    arguments,
    "mrdag-chunk-size",
    default = 25L
  )
  seed_base <- parse_integer_argument(
    arguments,
    "seed-base",
    default = 123L
  )
  mrdag_niter <- parse_integer_argument(
    arguments,
    "mrdag-niter",
    default = 1000L
  )
  mrdag_burnin <- parse_integer_argument(
    arguments,
    "mrdag-burnin",
    default = 200L,
    minimum = 0L
  )
  overwrite <- parse_boolean_argument(
    arguments,
    "overwrite",
    default = FALSE
  )

  if (mrdag_burnin >= mrdag_niter) {
    stop(
      "--mrdag-burnin must be smaller than --mrdag-niter.",
      call. = FALSE
    )
  }

  repo_root <- locate_repo_root()
  output_argument <- arguments[["output"]]

  if (is.null(output_argument) || !nzchar(trimws(output_argument))) {
    output_file <- file.path(
      repo_root,
      "paper",
      "config",
      "full_simulation_tasks.csv"
    )
  } else if (grepl("^(/|[A-Za-z]:[/\\\\])", output_argument)) {
    output_file <- output_argument
  } else {
    output_file <- file.path(repo_root, output_argument)
  }

  output_directory <- dirname(output_file)
  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  output_file <- normalizePath(
    output_file,
    winslash = "/",
    mustWork = FALSE
  )

  if (file.exists(output_file) && !overwrite) {
    stop(
      "Manifest already exists: ",
      output_file,
      "\nUse --overwrite=true only when deliberately changing it.",
      call. = FALSE
    )
  }

  designs <- c(
    "generic_low_rank",
    "sparse_loading"
  )
  settings <- data.frame(
    setting_index = seq_len(4L),
    scenario = c(
      "me_2.5_effect_0.25",
      "me_1_effect_0.25",
      "me_2.5_effect_1",
      "me_1_effect_1"
    ),
    me_weight = c(2.5, 1, 2.5, 1),
    effect_weight = c(0.25, 0.25, 1, 1),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  standard_phases <- c(
    "main_no_mrdag",
    "sparse"
  )
  mrdag_phases <- "main_mrdag"
  all_phases <- c(standard_phases, mrdag_phases)

  standard_ranges <- make_chunk_ranges(
    total_replicates = total_replicates,
    chunk_size = standard_chunk_size
  )
  mrdag_ranges <- make_chunk_ranges(
    total_replicates = total_replicates,
    chunk_size = mrdag_chunk_size
  )

  standard_tasks <- make_task_block(
    resource_class = "standard",
    designs = designs,
    phases = standard_phases,
    settings = settings,
    chunk_ranges = standard_ranges,
    seed_base = seed_base,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin
  )
  mrdag_tasks <- make_task_block(
    resource_class = "mrdag",
    designs = designs,
    phases = mrdag_phases,
    settings = settings,
    chunk_ranges = mrdag_ranges,
    seed_base = seed_base,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin
  )

  manifest <- rbind(
    standard_tasks,
    mrdag_tasks
  )
  rownames(manifest) <- NULL
  manifest$task_id <- seq_len(nrow(manifest))

  expected_standard_tasks <-
    length(designs) *
    length(standard_phases) *
    nrow(settings) *
    nrow(standard_ranges)
  expected_mrdag_tasks <-
    length(designs) *
    length(mrdag_phases) *
    nrow(settings) *
    nrow(mrdag_ranges)

  if (!identical(nrow(standard_tasks), expected_standard_tasks) ||
      !identical(nrow(mrdag_tasks), expected_mrdag_tasks)) {
    stop(
      "The generated task counts are inconsistent.",
      call. = FALSE
    )
  }

  if (!identical(
    standard_tasks$array_index,
    seq_len(nrow(standard_tasks))
  ) || !identical(
    mrdag_tasks$array_index,
    seq_len(nrow(mrdag_tasks))
  )) {
    stop(
      "Array indices are not consecutive within resource class.",
      call. = FALSE
    )
  }

  if (anyDuplicated(manifest$expected_output_file) ||
      anyDuplicated(manifest$task_id)) {
    stop(
      "The task manifest contains duplicate task IDs or output files.",
      call. = FALSE
    )
  }

  validate_task_coverage(
    manifest = manifest,
    total_replicates = total_replicates,
    designs = designs,
    phases = all_phases,
    settings = settings
  )

  write_csv_atomically(manifest, output_file)

  regenerated <- utils::read.csv(
    output_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (nrow(regenerated) != nrow(manifest) ||
      !identical(regenerated$task_id, manifest$task_id) ||
      !identical(
        regenerated$expected_output_file,
        manifest$expected_output_file
      )) {
    stop(
      "The written manifest failed its read-back check.",
      call. = FALSE
    )
  }

  cat("MR-rr full simulation task manifest\n")
  cat("Total simulation replicates per setting:", total_replicates, "\n")
  cat(
    "Standard chunk size:",
    standard_chunk_size,
    "| tasks:",
    nrow(standard_tasks),
    "\n"
  )
  cat(
    "MrDAG chunk size:",
    mrdag_chunk_size,
    "| tasks:",
    nrow(mrdag_tasks),
    "\n"
  )
  cat("Total tasks:", nrow(manifest), "\n")
  cat("Manifest:", output_file, "\n")
  cat("Task manifest coverage: PASS\n")

  invisible(manifest)
}


if (sys.nframe() == 0L) {
  make_full_simulation_manifest()
}
