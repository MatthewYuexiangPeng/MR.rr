#!/usr/bin/env Rscript

# Validate and merge point-estimation chunks for the additional simulations.
#
# Modes:
#   --check-chunk=<path>
#       Validate one chunk without consulting a task manifest.
#
#   Full merge mode
#       Validate the point-task manifest and all present expected chunks. With
#       --allow-incomplete=true, write inventories and stop successfully when
#       tasks are missing. Otherwise require complete coverage and write one
#       merged RDS object.

options(warn = 1)


additional_merge_parse_arguments <- function(arguments) {
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


additional_merge_validate_names <- function(arguments, allowed) {
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


additional_merge_parse_boolean <- function(
    arguments,
    name,
    default = FALSE) {
  value <- arguments[[name]]

  if (is.null(value)) {
    return(default)
  }

  additional_merge_boolean_value(value, paste0("--", name))
}


additional_merge_boolean_value <- function(value, label) {
  if (length(value) != 1L || is.na(value)) {
    stop(label, " must contain one Boolean value.", call. = FALSE)
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


additional_merge_resolve_path <- function(
    path,
    repo_root,
    must_work = FALSE) {
  candidate <- if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    path
  } else {
    file.path(repo_root, path)
  }

  normalizePath(candidate, winslash = "/", mustWork = must_work)
}


additional_merge_integer <- function(value, label, minimum = 0L) {
  numeric_value <- suppressWarnings(as.numeric(value))

  if (length(value) != 1L ||
      is.na(value) ||
      length(numeric_value) != 1L ||
      !is.finite(numeric_value) ||
      numeric_value != round(numeric_value) ||
      numeric_value < minimum ||
      numeric_value > .Machine$integer.max) {
    stop(label, " must be a valid integer.", call. = FALSE)
  }

  as.integer(numeric_value)
}


additional_merge_text <- function(value, label) {
  if (length(value) != 1L ||
      is.na(value) ||
      !nzchar(trimws(as.character(value)))) {
    stop(label, " must be one nonempty value.", call. = FALSE)
  }

  trimws(as.character(value))
}


additional_merge_assert_integer <- function(actual, expected, label) {
  parsed <- suppressWarnings(as.integer(actual))

  if (length(actual) != 1L ||
      is.na(actual) ||
      length(parsed) != 1L ||
      is.na(parsed) ||
      !identical(parsed, as.integer(expected))) {
    stop(
      label,
      " differs: expected ",
      expected,
      ", found ",
      paste(actual, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  invisible(actual)
}


additional_merge_assert_text <- function(actual, expected, label) {
  if (length(actual) != 1L ||
      is.na(actual) ||
      !identical(as.character(actual), as.character(expected))) {
    stop(
      label,
      " differs: expected `",
      expected,
      "`, found `",
      paste(actual, collapse = ", "),
      "`.",
      call. = FALSE
    )
  }

  invisible(actual)
}


additional_merge_assert_close <- function(
    actual,
    expected,
    label,
    tolerance = 1e-8) {
  if (!identical(dim(actual), dim(expected)) ||
      length(actual) != length(expected)) {
    stop(label, " has unexpected dimensions.", call. = FALSE)
  }

  difference <- suppressWarnings(max(abs(actual - expected)))

  if (!is.finite(difference) || difference > tolerance) {
    stop(
      label,
      " differs by ",
      format(difference, scientific = TRUE),
      "; tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(difference)
}


additional_merge_same_object <- function(actual, expected, label) {
  comparison <- all.equal(
    actual,
    expected,
    tolerance = 1e-12,
    check.attributes = FALSE
  )

  if (!isTRUE(comparison)) {
    stop(label, " differs across chunks.", call. = FALSE)
  }

  invisible(TRUE)
}


additional_merge_write_csv <- function(value, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  utils::write.csv(value, temporary_file, row.names = FALSE, na = "")
  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = dirname(path),
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop("Could not preserve existing CSV: ", path, call. = FALSE)
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop("Could not install completed CSV: ", path, call. = FALSE)
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


additional_merge_save_rds <- function(value, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop(
      "Merged point output already exists: ",
      path,
      "\nUse --overwrite=true only when deliberately replacing it.",
      call. = FALSE
    )
  }

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  saveRDS(value, temporary_file, compress = "xz")
  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = dirname(path),
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop("Could not preserve existing merged output: ", path, call. = FALSE)
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop("Could not install merged point output: ", path, call. = FALSE)
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


additional_merge_expected_filename <- function(
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end) {
  sprintf(
    "%s__%s__setting-%d__rep-%04d-%04d.rds",
    analysis,
    phase,
    as.integer(setting_index),
    as.integer(replicate_start),
    as.integer(replicate_end)
  )
}


additional_merge_expected_groups <- function() {
  rbind(
    data.frame(
      analysis = "rank_misspecification",
      phase = "standard",
      setting_index = seq_len(4L),
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis = "approximate_low_rank",
      phase = "standard",
      setting_index = seq_len(4L),
      stringsAsFactors = FALSE
    ),
    data.frame(
      analysis = "approximate_low_rank",
      phase = "mrdag",
      setting_index = seq_len(4L),
      stringsAsFactors = FALSE
    )
  )
}


additional_merge_manifest_boolean <- function(manifest, row_index, name) {
  additional_merge_boolean_value(
    manifest[[name]][[row_index]],
    paste0("Manifest row ", row_index, " field ", name)
  )
}


additional_merge_validate_manifest <- function(
    manifest,
    repo_root,
    core) {
  required_columns <- c(
    "manifest_schema_version",
    "task_id",
    "resource_class",
    "array_index",
    "analysis",
    "phase",
    "setting_index",
    "scenario",
    "measurement_error_weight",
    "genetic_effect_weight",
    "third_singular_value",
    "configuration_status",
    "chunk_index",
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

  if (nrow(manifest) == 0L) {
    stop("Point-task manifest contains no rows.", call. = FALSE)
  }

  task_ids <- suppressWarnings(as.integer(manifest$task_id))
  array_indices <- suppressWarnings(as.integer(manifest$array_index))

  if (anyNA(task_ids) ||
      anyNA(array_indices) ||
      !identical(task_ids, seq_len(nrow(manifest))) ||
      anyDuplicated(manifest$expected_output_file)) {
    stop("Point-task manifest task IDs or output files are invalid.", call. = FALSE)
  }

  for (resource_class in c("standard", "mrdag")) {
    selected_indices <- array_indices[
      manifest$resource_class == resource_class
    ]

    if (length(selected_indices) == 0L ||
        !identical(sort(selected_indices), seq_along(selected_indices))) {
      stop(
        "Point-task manifest has invalid ",
        resource_class,
        " array indices.",
        call. = FALSE
      )
    }
  }

  config_files <- c(
    rank_misspecification = file.path(
      repo_root,
      "paper",
      "config",
      "rank_misspecification_settings.csv"
    ),
    approximate_low_rank = file.path(
      repo_root,
      "paper",
      "config",
      "approximate_low_rank_settings.csv"
    )
  )
  configs <- list(
    rank_misspecification = core$additional_read_config(
      config_files[["rank_misspecification"]],
      "rank_misspecification"
    ),
    approximate_low_rank = core$additional_read_config(
      config_files[["approximate_low_rank"]],
      "approximate_low_rank"
    )
  )
  worker_path <- "paper/scripts/22_run_additional_simulation_chunk.R"
  allowed_analyses <- names(configs)

  for (row_index in seq_len(nrow(manifest))) {
    row <- manifest[row_index, , drop = FALSE]
    analysis <- additional_merge_text(row$analysis, "Manifest analysis")
    phase <- additional_merge_text(row$phase, "Manifest phase")
    resource_class <- additional_merge_text(
      row$resource_class,
      "Manifest resource class"
    )
    setting_index <- additional_merge_integer(
      row$setting_index,
      "Manifest setting index",
      1L
    )
    replicate_start <- additional_merge_integer(
      row$replicate_start,
      "Manifest replicate start",
      1L
    )
    replicate_end <- additional_merge_integer(
      row$replicate_end,
      "Manifest replicate end",
      1L
    )
    replicate_count <- additional_merge_integer(
      row$replicate_count,
      "Manifest replicate count",
      1L
    )
    allow_provisional <- additional_merge_manifest_boolean(
      manifest,
      row_index,
      "allow_provisional"
    )

    if (row$manifest_schema_version != "1.0.0" ||
        !analysis %in% allowed_analyses ||
        !phase %in% c("standard", "mrdag") ||
        !setting_index %in% seq_len(4L) ||
        replicate_end < replicate_start ||
        replicate_count != replicate_end - replicate_start + 1L) {
      stop(
        "Point-task manifest row ",
        row_index,
        " is internally invalid.",
        call. = FALSE
      )
    }

    expected_resource <- if (phase == "mrdag") "mrdag" else "standard"

    if (resource_class != expected_resource ||
        (analysis == "rank_misspecification" && phase != "standard") ||
        (phase == "mrdag" && analysis != "approximate_low_rank")) {
      stop(
        "Manifest analysis/phase/resource mismatch in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    config <- configs[[analysis]]
    selected_config <- config[
      config$setting_index == setting_index,
      ,
      drop = FALSE
    ]
    reference <- selected_config[1L, , drop = FALSE]

    if (row$scenario != reference$scenario ||
        abs(row$measurement_error_weight -
          reference$measurement_error_weight) > 1e-12 ||
        abs(row$genetic_effect_weight -
          reference$genetic_effect_weight) > 1e-12 ||
        abs(row$third_singular_value -
          reference$third_singular_value) > 1e-12 ||
        row$configuration_status != reference$status ||
        row$seed_base != reference$seed_base) {
      stop(
        "Manifest configuration metadata mismatch in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    expected_config_path <- substring(
      normalizePath(
        config_files[[analysis]],
        winslash = "/",
        mustWork = TRUE
      ),
      nchar(normalizePath(repo_root, winslash = "/", mustWork = TRUE)) + 2L
    )
    expected_file <- additional_merge_expected_filename(
      analysis,
      phase,
      setting_index,
      replicate_start,
      replicate_end
    )

    if (row$worker_script != worker_path ||
        row$configuration_file != expected_config_path ||
        row$expected_output_file != expected_file ||
        basename(row$expected_output_file) != row$expected_output_file ||
        !nzchar(trimws(row$output_directory))) {
      stop(
        "Manifest path or filename mismatch in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    hash_values <- unlist(
      row[, c("configuration_md5", "core_md5", "worker_md5")],
      use.names = FALSE
    )

    if (anyNA(hash_values) ||
        any(!grepl("^[0-9a-fA-F]{32}$", hash_values))) {
      stop("Manifest has an invalid MD5 in row ", row_index, ".", call. = FALSE)
    }

    if (reference$status != "locked") {
      if (analysis != "approximate_low_rank" ||
          !allow_provisional ||
          replicate_count > 5L) {
        stop(
          "Invalid provisional task in manifest row ",
          row_index,
          ".",
          call. = FALSE
        )
      }
    }

    if (analysis == "rank_misspecification" && allow_provisional) {
      stop("Rank task cannot be provisional in row ", row_index, ".", call. = FALSE)
    }

    if (phase == "mrdag") {
      niter <- additional_merge_integer(
        row$mrdag_niter,
        "Manifest MrDAG niter",
        1L
      )
      burnin <- additional_merge_integer(
        row$mrdag_burnin,
        "Manifest MrDAG burnin",
        0L
      )

      if (burnin >= niter) {
        stop("Manifest has invalid MrDAG iterations in row ", row_index, ".", call. = FALSE)
      }
    }
  }

  actual_groups <- unique(
    manifest[, c("analysis", "phase", "setting_index")]
  )
  expected_groups <- additional_merge_expected_groups()
  actual_keys <- paste(
    actual_groups$analysis,
    actual_groups$phase,
    actual_groups$setting_index,
    sep = "/"
  )
  expected_keys <- paste(
    expected_groups$analysis,
    expected_groups$phase,
    expected_groups$setting_index,
    sep = "/"
  )

  if (!setequal(actual_keys, expected_keys) ||
      length(actual_keys) != length(expected_keys)) {
    stop("Point-task manifest does not contain the expected 12 groups.", call. = FALSE)
  }

  group_totals <- integer(nrow(expected_groups))

  for (group_index in seq_len(nrow(expected_groups))) {
    key <- expected_groups[group_index, , drop = FALSE]
    selected <- manifest[
      manifest$analysis == key$analysis &
        manifest$phase == key$phase &
        manifest$setting_index == key$setting_index,
      ,
      drop = FALSE
    ]
    selected <- selected[order(selected$replicate_start), , drop = FALSE]

    if (selected$replicate_start[[1L]] != 1L ||
        (nrow(selected) > 1L && any(
          selected$replicate_start[-1L] !=
            selected$replicate_end[-nrow(selected)] + 1L
        ))) {
      stop(
        "Manifest has a gap or overlap for ",
        paste(key$analysis, key$phase, key$setting_index, sep = "/"),
        ".",
        call. = FALSE
      )
    }

    group_totals[[group_index]] <- max(selected$replicate_end)
  }

  total_replicates <- unique(group_totals)

  if (length(total_replicates) != 1L) {
    stop("Point-task groups have different total replicate counts.", call. = FALSE)
  }

  if (any(manifest$configuration_status != "locked") &&
      total_replicates > 5L) {
    stop("A provisional point manifest cannot exceed 5 replicates.", call. = FALSE)
  }

  list(
    total_replicates = as.integer(total_replicates),
    group_keys = expected_groups,
    configs = configs
  )
}


additional_merge_validate_matrix <- function(
    value,
    expected_rows,
    replicate_ids,
    label) {
  expected_dimension <- c(as.integer(expected_rows), length(replicate_ids))

  if (!is.matrix(value) || !identical(dim(value), expected_dimension)) {
    stop(label, " has unexpected dimensions.", call. = FALSE)
  }

  if (any(!is.finite(value))) {
    stop(label, " contains a non-finite value.", call. = FALSE)
  }

  if (!is.null(colnames(value)) &&
      !identical(colnames(value), as.character(replicate_ids))) {
    stop(label, " has incorrect replicate column names.", call. = FALSE)
  }

  invisible(value)
}


additional_merge_provenance_md5 <- function(input_manifest, path) {
  normalized_paths <- gsub("\\\\", "/", input_manifest$repository_path)
  selected <- input_manifest$md5[normalized_paths == path]

  if (length(selected) != 1L ||
      is.na(selected) ||
      !grepl("^[0-9a-fA-F]{32}$", selected)) {
    stop("Chunk provenance is missing or duplicates: ", path, call. = FALSE)
  }

  tolower(selected)
}


additional_merge_validate_chunk <- function(
    chunk,
    chunk_file,
    repo_root,
    core,
    worker,
    expected_row = NULL,
    tolerance = 1e-8) {
  required_components <- c(
    "metadata",
    "input_manifest",
    "configuration",
    "result_specification",
    "truth",
    "estimates",
    "biases",
    "successful_replicates",
    "errors",
    "session_info"
  )
  missing_components <- setdiff(required_components, names(chunk))

  if (length(missing_components) > 0L) {
    stop(
      basename(chunk_file),
      " is missing components: ",
      paste(missing_components, collapse = ", "),
      call. = FALSE
    )
  }

  metadata <- chunk$metadata
  required_metadata <- c(
    "schema_version",
    "analysis",
    "phase",
    "setting_index",
    "scenario",
    "replicate_start",
    "replicate_end",
    "replicate_ids",
    "replicate_count",
    "seed_base",
    "data_seeds",
    "estimator_seeds",
    "rng_kind",
    "mrdag_niter",
    "mrdag_burnin",
    "allow_provisional",
    "complete",
    "archived_simulation_result_rdata_read",
    "frozen_snapshot_written"
  )
  missing_metadata <- setdiff(required_metadata, names(metadata))

  if (!is.list(metadata) || length(missing_metadata) > 0L) {
    stop(
      basename(chunk_file),
      " has incomplete metadata: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }

  analysis <- additional_merge_text(metadata$analysis, "Chunk analysis")
  phase <- additional_merge_text(metadata$phase, "Chunk phase")
  setting_index <- additional_merge_integer(
    metadata$setting_index,
    "Chunk setting",
    1L
  )
  replicate_start <- additional_merge_integer(
    metadata$replicate_start,
    "Chunk replicate start",
    1L
  )
  replicate_end <- additional_merge_integer(
    metadata$replicate_end,
    "Chunk replicate end",
    1L
  )
  seed_base <- additional_merge_integer(
    metadata$seed_base,
    "Chunk seed base",
    1L
  )
  allow_provisional <- additional_merge_boolean_value(
    metadata$allow_provisional,
    "Chunk allow_provisional"
  )

  if (metadata$schema_version != "1.0.0" ||
      !analysis %in% c("rank_misspecification", "approximate_low_rank") ||
      !phase %in% c("standard", "mrdag") ||
      !setting_index %in% seq_len(4L) ||
      (analysis == "rank_misspecification" && phase != "standard") ||
      (phase == "mrdag" && analysis != "approximate_low_rank") ||
      replicate_end < replicate_start ||
      !isTRUE(metadata$complete) ||
      !identical(metadata$archived_simulation_result_rdata_read, FALSE) ||
      !identical(metadata$frozen_snapshot_written, FALSE)) {
    stop(basename(chunk_file), " has invalid task metadata.", call. = FALSE)
  }

  if (!is.data.frame(chunk$errors) || nrow(chunk$errors) != 0L) {
    stop(basename(chunk_file), " contains estimator errors.", call. = FALSE)
  }

  canonical_filename <- additional_merge_expected_filename(
    analysis,
    phase,
    setting_index,
    replicate_start,
    replicate_end
  )

  if (!identical(basename(chunk_file), canonical_filename)) {
    stop(basename(chunk_file), " has a noncanonical filename.", call. = FALSE)
  }

  replicate_ids <- as.integer(metadata$replicate_ids)
  expected_replicate_ids <- seq.int(replicate_start, replicate_end)

  if (!identical(replicate_ids, expected_replicate_ids) ||
      metadata$replicate_count != length(expected_replicate_ids)) {
    stop(basename(chunk_file), " has inconsistent replicate IDs.", call. = FALSE)
  }

  expected_data_seeds <- vapply(
    replicate_ids,
    function(replicate_id) {
      core$additional_make_replicate_seed(
        seed_base = seed_base,
        setting_index = setting_index,
        replicate_id = replicate_id
      )
    },
    integer(1)
  )

  if (!identical(as.integer(metadata$data_seeds), expected_data_seeds) ||
      anyDuplicated(metadata$data_seeds)) {
    stop(basename(chunk_file), " has inconsistent data seeds.", call. = FALSE)
  }

  if (!identical(
    as.character(metadata$rng_kind),
    c("L'Ecuyer-CMRG", "Inversion", "Rejection")
  )) {
    stop(basename(chunk_file), " has unexpected RNG metadata.", call. = FALSE)
  }

  config_file <- file.path(
    repo_root,
    "paper",
    "config",
    paste0(analysis, "_settings.csv")
  )
  full_config <- core$additional_read_config(config_file, analysis)
  expected_config <- full_config[
    full_config$setting_index == setting_index,
    ,
    drop = FALSE
  ]
  expected_specification <- worker$additional_chunk_result_specification(
    analysis,
    phase,
    expected_config
  )
  expected_config <- expected_specification$config
  expected_specification <- expected_specification$specification

  if (!identical(
    as.character(metadata$scenario),
    unique(as.character(expected_config$scenario))
  )) {
    stop(basename(chunk_file), " has incorrect scenario metadata.", call. = FALSE)
  }

  configuration_status <- unique(as.character(expected_config$status))

  if ((analysis == "rank_misspecification" && allow_provisional) ||
      (configuration_status != "locked" &&
        (!allow_provisional || length(replicate_ids) > 5L))) {
    stop(basename(chunk_file), " has an invalid provisional flag.", call. = FALSE)
  }

  additional_merge_same_object(
    chunk$configuration,
    expected_config,
    paste(basename(chunk_file), "configuration")
  )
  additional_merge_same_object(
    chunk$result_specification,
    expected_specification,
    paste(basename(chunk_file), "result specification")
  )

  delta <- unique(expected_config$third_singular_value)
  expected_truth <- core$additional_make_generic_truth(
    third_singular_value = delta,
    px = 9L,
    py = 3L,
    seed = 123L
  )

  if (!is.list(chunk$truth) ||
      !is.matrix(chunk$truth$C) ||
      !identical(dim(chunk$truth$C), c(3L, 9L)) ||
      any(!is.finite(chunk$truth$C))) {
    stop(basename(chunk_file), " has invalid truth metadata.", call. = FALSE)
  }

  additional_merge_assert_close(
    chunk$truth$C,
    expected_truth$C,
    paste(basename(chunk_file), "truth C"),
    tolerance
  )
  additional_merge_assert_close(
    as.numeric(chunk$truth$singular_values),
    as.numeric(expected_truth$singular_values),
    paste(basename(chunk_file), "truth singular values"),
    tolerance
  )
  additional_merge_assert_integer(
    chunk$truth$numerical_rank,
    expected_truth$numerical_rank,
    paste(basename(chunk_file), "truth rank")
  )

  result_keys <- expected_specification$result_key

  if (!identical(names(chunk$estimates), result_keys) ||
      !identical(names(chunk$biases), result_keys) ||
      !identical(names(chunk$successful_replicates), result_keys) ||
      !identical(names(metadata$estimator_seeds), result_keys)) {
    stop(basename(chunk_file), " has unexpected result keys.", call. = FALSE)
  }

  truth_vector <- as.vector(expected_truth$C)

  for (result_index in seq_len(nrow(expected_specification))) {
    result_key <- expected_specification$result_key[[result_index]]
    method <- expected_specification$method[[result_index]]
    working_rank <- expected_specification$working_rank[[result_index]]
    additional_merge_validate_matrix(
      chunk$estimates[[result_key]],
      27L,
      replicate_ids,
      paste(basename(chunk_file), result_key, "estimates")
    )
    additional_merge_validate_matrix(
      chunk$biases[[result_key]],
      27L,
      replicate_ids,
      paste(basename(chunk_file), result_key, "biases")
    )
    expected_bias <- sweep(
      chunk$estimates[[result_key]],
      1L,
      truth_vector,
      FUN = "-"
    )
    additional_merge_assert_close(
      chunk$biases[[result_key]],
      expected_bias,
      paste(basename(chunk_file), result_key, "bias reconstruction"),
      tolerance
    )
    expected_estimator_seeds <- vapply(
      expected_data_seeds,
      function(data_seed) {
        worker$additional_chunk_make_result_seed(
          data_seed,
          method,
          working_rank
        )
      },
      integer(1)
    )

    if (!identical(
      as.integer(metadata$estimator_seeds[[result_key]]),
      expected_estimator_seeds
    )) {
      stop(
        basename(chunk_file),
        " has inconsistent estimator seeds for ",
        result_key,
        ".",
        call. = FALSE
      )
    }
  }

  if (any(as.integer(chunk$successful_replicates) != length(replicate_ids))) {
    stop(basename(chunk_file), " has inconsistent success counts.", call. = FALSE)
  }

  required_input_columns <- c("repository_path", "md5")

  if (!is.data.frame(chunk$input_manifest) ||
      !all(required_input_columns %in% names(chunk$input_manifest)) ||
      nrow(chunk$input_manifest) == 0L ||
      anyNA(chunk$input_manifest[, required_input_columns]) ||
      any(!nzchar(chunk$input_manifest$repository_path)) ||
      any(!grepl("^[0-9a-fA-F]{32}$", chunk$input_manifest$md5))) {
    stop(basename(chunk_file), " has invalid input provenance.", call. = FALSE)
  }

  provenance_paths <- gsub("\\\\", "/", chunk$input_manifest$repository_path)

  if (anyDuplicated(provenance_paths) || any(grepl(
    "(^|/)results(/|$)|\\.[Rr][Dd]ata$|(^|/)output(/|$)",
    provenance_paths
  ))) {
    stop(
      basename(chunk_file),
      " lists an archived or generated result as an input.",
      call. = FALSE
    )
  }

  core_path <- "paper/scripts/20_additional_simulation_core.R"
  worker_path <- "paper/scripts/22_run_additional_simulation_chunk.R"
  config_path <- paste0("paper/config/", analysis, "_settings.csv")
  additional_merge_provenance_md5(chunk$input_manifest, core_path)
  additional_merge_provenance_md5(chunk$input_manifest, worker_path)
  additional_merge_provenance_md5(chunk$input_manifest, config_path)

  if (!is.null(expected_row)) {
    additional_merge_assert_text(
      analysis,
      expected_row$analysis,
      "Chunk analysis"
    )
    additional_merge_assert_text(phase, expected_row$phase, "Chunk phase")
    additional_merge_assert_integer(
      setting_index,
      expected_row$setting_index,
      "Chunk setting"
    )
    additional_merge_assert_text(
      metadata$scenario,
      expected_row$scenario,
      "Chunk scenario"
    )
    additional_merge_assert_integer(
      replicate_start,
      expected_row$replicate_start,
      "Chunk replicate start"
    )
    additional_merge_assert_integer(
      replicate_end,
      expected_row$replicate_end,
      "Chunk replicate end"
    )
    additional_merge_assert_integer(
      metadata$replicate_count,
      expected_row$replicate_count,
      "Chunk replicate count"
    )
    additional_merge_assert_integer(
      seed_base,
      expected_row$seed_base,
      "Chunk seed base"
    )
    expected_allow_provisional <- additional_merge_boolean_value(
      expected_row$allow_provisional,
      "Manifest allow_provisional"
    )

    if (!identical(allow_provisional, expected_allow_provisional)) {
      stop("Chunk provisional flag differs from its manifest row.", call. = FALSE)
    }

    if (basename(chunk_file) != expected_row$expected_output_file) {
      stop("Chunk filename differs from its manifest row.", call. = FALSE)
    }

    if (phase == "mrdag") {
      additional_merge_assert_integer(
        metadata$mrdag_niter,
        expected_row$mrdag_niter,
        "Chunk MrDAG niter"
      )
      additional_merge_assert_integer(
        metadata$mrdag_burnin,
        expected_row$mrdag_burnin,
        "Chunk MrDAG burnin"
      )
    }

    observed_hashes <- unname(c(
      configuration_md5 = additional_merge_provenance_md5(
        chunk$input_manifest,
        config_path
      ),
      core_md5 = additional_merge_provenance_md5(
        chunk$input_manifest,
        core_path
      ),
      worker_md5 = additional_merge_provenance_md5(
        chunk$input_manifest,
        worker_path
      )
    ))
    expected_hashes <- unname(tolower(c(
      configuration_md5 = expected_row$configuration_md5,
      core_md5 = expected_row$core_md5,
      worker_md5 = expected_row$worker_md5
    )))

    if (!identical(observed_hashes, expected_hashes)) {
      stop("Chunk input hashes differ from its manifest row.", call. = FALSE)
    }
  }

  invisible(chunk)
}


additional_merge_group <- function(
    manifest_group,
    chunks,
    total_replicates) {
  ordering <- order(manifest_group$replicate_start)
  manifest_group <- manifest_group[ordering, , drop = FALSE]
  chunks <- chunks[ordering]
  reference <- chunks[[1L]]

  if (length(chunks) > 1L) {
    for (chunk_index in seq.int(2L, length(chunks))) {
      chunk <- chunks[[chunk_index]]
      additional_merge_same_object(
        chunk$input_manifest,
        reference$input_manifest,
        "Input provenance"
      )
      additional_merge_same_object(
        chunk$configuration,
        reference$configuration,
        "Analysis configuration"
      )
      additional_merge_same_object(
        chunk$result_specification,
        reference$result_specification,
        "Result specification"
      )
      additional_merge_same_object(
        chunk$truth,
        reference$truth,
        "Simulation truth"
      )
    }
  }

  replicate_ids <- unlist(
    lapply(chunks, function(chunk) chunk$metadata$replicate_ids),
    use.names = FALSE
  )
  data_seeds <- unlist(
    lapply(chunks, function(chunk) chunk$metadata$data_seeds),
    use.names = FALSE
  )

  if (!identical(as.integer(replicate_ids), seq_len(total_replicates)) ||
      anyDuplicated(replicate_ids) ||
      anyDuplicated(data_seeds)) {
    stop(
      "Merged chunks have a replicate or seed gap/overlap for ",
      paste(
        reference$metadata$analysis,
        reference$metadata$phase,
        reference$metadata$setting_index,
        sep = "/"
      ),
      ".",
      call. = FALSE
    )
  }

  result_keys <- reference$result_specification$result_key
  merge_named_values <- function(component) {
    stats::setNames(
      lapply(result_keys, function(result_key) {
        value <- do.call(
          cbind,
          lapply(chunks, function(chunk) chunk[[component]][[result_key]])
        )
        colnames(value) <- as.character(replicate_ids)
        value
      }),
      result_keys
    )
  }
  estimator_seeds <- stats::setNames(
    lapply(result_keys, function(result_key) {
      as.integer(unlist(
        lapply(
          chunks,
          function(chunk) chunk$metadata$estimator_seeds[[result_key]]
        ),
        use.names = FALSE
      ))
    }),
    result_keys
  )

  list(
    analysis = reference$metadata$analysis,
    phase = reference$metadata$phase,
    setting_index = as.integer(reference$metadata$setting_index),
    scenario = reference$metadata$scenario,
    seed_base = as.integer(reference$metadata$seed_base),
    replicate_ids = as.integer(replicate_ids),
    data_seeds = as.integer(data_seeds),
    estimator_seeds = estimator_seeds,
    configuration = reference$configuration,
    result_specification = reference$result_specification,
    truth = reference$truth,
    estimates = merge_named_values("estimates"),
    biases = merge_named_values("biases"),
    input_manifest = reference$input_manifest,
    chunk_files = manifest_group$expected_output_file
  )
}


additional_merge_chunk_path <- function(
    row,
    repo_root,
    chunks_override = NULL) {
  if (!is.null(chunks_override)) {
    return(file.path(chunks_override, row$expected_output_file))
  }

  output_directory <- additional_merge_resolve_path(
    row$output_directory,
    repo_root,
    must_work = FALSE
  )
  file.path(output_directory, row$expected_output_file)
}


additional_merge_run_standalone <- function(
    chunk_argument,
    repo_root,
    core,
    worker) {
  chunk_file <- additional_merge_resolve_path(
    chunk_argument,
    repo_root,
    must_work = TRUE
  )
  chunk <- readRDS(chunk_file)
  additional_merge_validate_chunk(
    chunk,
    chunk_file,
    repo_root,
    core,
    worker
  )

  cat("MR-rr standalone additional point-chunk check\n")
  cat("Chunk:", chunk_file, "\n")
  cat("Analysis:", chunk$metadata$analysis, "\n")
  cat("Phase:", chunk$metadata$phase, "\n")
  cat("Setting:", chunk$metadata$setting_index, "\n")
  cat(
    "Replicates:",
    chunk$metadata$replicate_start,
    "through",
    chunk$metadata$replicate_end,
    "\n"
  )
  cat("Archived simulation-result RData read: no\n")
  cat("Frozen snapshot written: no\n")
  cat("Standalone additional point-chunk validation: PASS\n")
  invisible(chunk)
}


additional_merge_full <- function(
    arguments,
    repo_root,
    core,
    worker,
    allow_incomplete,
    overwrite) {
  manifest_argument <- arguments[["manifest"]]
  chunks_argument <- arguments[["chunks-dir"]]
  output_argument <- arguments[["output-dir"]]
  manifest_file <- if (is.null(manifest_argument) ||
      !nzchar(trimws(manifest_argument))) {
    file.path(repo_root, "paper", "config", "additional_point_tasks.csv")
  } else {
    additional_merge_resolve_path(manifest_argument, repo_root)
  }
  chunks_override <- if (is.null(chunks_argument) ||
      !nzchar(trimws(chunks_argument))) {
    NULL
  } else {
    additional_merge_resolve_path(chunks_argument, repo_root, must_work = TRUE)
  }
  output_directory <- if (is.null(output_argument) ||
      !nzchar(trimws(output_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "point_merged"
    )
  } else {
    additional_merge_resolve_path(output_argument, repo_root)
  }

  if (!file.exists(manifest_file)) {
    stop("Point-task manifest is missing: ", manifest_file, call. = FALSE)
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  manifest_file <- normalizePath(
    manifest_file,
    winslash = "/",
    mustWork = TRUE
  )
  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )
  manifest <- utils::read.csv(
    manifest_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  validation <- additional_merge_validate_manifest(manifest, repo_root, core)
  total_replicates <- validation$total_replicates
  chunk_paths <- vapply(
    seq_len(nrow(manifest)),
    function(row_index) {
      additional_merge_chunk_path(
        manifest[row_index, , drop = FALSE],
        repo_root,
        chunks_override
      )
    },
    character(1)
  )
  present <- file.exists(chunk_paths)
  missing_tasks <- manifest[!present, , drop = FALSE]
  missing_file <- file.path(output_directory, "missing_additional_point_tasks.csv")
  additional_merge_write_csv(missing_tasks, missing_file)
  inventory <- manifest
  inventory$chunk_path <- chunk_paths
  inventory$present <- present
  inventory$validated <- FALSE
  inventory$file_size_bytes <- NA_real_
  inventory$md5 <- NA_character_
  inventory$schema_version <- NA_character_
  inventory$complete <- NA
  inventory$archived_simulation_result_rdata_read <- NA
  inventory$frozen_snapshot_written <- NA
  chunks <- vector("list", nrow(manifest))

  cat("MR-rr additional point-chunk merger\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Expected chunks:", nrow(manifest), "\n")
  cat("Present chunks:", sum(present), "\n")
  cat("Missing chunks:", sum(!present), "\n")

  for (row_index in which(present)) {
    chunk <- readRDS(chunk_paths[[row_index]])
    additional_merge_validate_chunk(
      chunk,
      chunk_paths[[row_index]],
      repo_root,
      core,
      worker,
      expected_row = manifest[row_index, , drop = FALSE]
    )
    chunks[[row_index]] <- chunk
    file_information <- file.info(chunk_paths[[row_index]])
    inventory$validated[[row_index]] <- TRUE
    inventory$file_size_bytes[[row_index]] <- file_information$size
    inventory$md5[[row_index]] <- unname(tools::md5sum(chunk_paths[[row_index]]))
    inventory$schema_version[[row_index]] <- chunk$metadata$schema_version
    inventory$complete[[row_index]] <- chunk$metadata$complete
    inventory$archived_simulation_result_rdata_read[[row_index]] <-
      chunk$metadata$archived_simulation_result_rdata_read
    inventory$frozen_snapshot_written[[row_index]] <-
      chunk$metadata$frozen_snapshot_written
  }

  inventory_file <- file.path(
    output_directory,
    "additional_point_chunk_inventory.csv"
  )
  additional_merge_write_csv(inventory, inventory_file)

  if (nrow(missing_tasks) > 0L) {
    if (allow_incomplete) {
      cat("Existing expected chunks validated: PASS\n")
      cat("Missing-task inventory:", missing_file, "\n")
      cat("Full point merge deferred until every task is present.\n")
      cat("Incomplete additional point inventory: PASS\n")
      return(invisible(inventory))
    }

    stop(
      "Additional point output is incomplete. See: ",
      missing_file,
      "\nUse --allow-incomplete=true only for a progress check.",
      call. = FALSE
    )
  }

  if (!all(inventory$validated)) {
    stop("At least one expected point chunk was not validated.", call. = FALSE)
  }

  group_keys <- validation$group_keys
  groups <- vector("list", nrow(group_keys))

  for (group_index in seq_len(nrow(group_keys))) {
    key <- group_keys[group_index, , drop = FALSE]
    selected_indices <- which(
      manifest$analysis == key$analysis &
        manifest$phase == key$phase &
        manifest$setting_index == key$setting_index
    )
    cat(
      "Merging:",
      paste(key$analysis, key$phase, key$setting_index, sep = "/"),
      "\n"
    )
    groups[[group_index]] <- additional_merge_group(
      manifest[selected_indices, , drop = FALSE],
      chunks[selected_indices],
      total_replicates
    )
  }

  names(groups) <- vapply(
    groups,
    function(group) {
      paste(
        group$analysis,
        group$phase,
        group$setting_index,
        sep = "__"
      )
    },
    character(1)
  )

  for (setting_index in seq_len(4L)) {
    standard <- groups[[paste(
      "approximate_low_rank",
      "standard",
      setting_index,
      sep = "__"
    )]]
    mrdag <- groups[[paste(
      "approximate_low_rank",
      "mrdag",
      setting_index,
      sep = "__"
    )]]

    if (!identical(standard$replicate_ids, mrdag$replicate_ids) ||
        !identical(standard$data_seeds, mrdag$data_seeds)) {
      stop(
        "Approximate standard and MrDAG seed streams differ for setting ",
        setting_index,
        ".",
        call. = FALSE
      )
    }

    additional_merge_same_object(
      standard$truth,
      mrdag$truth,
      paste("Approximate-low-rank setting", setting_index, "truth")
    )
  }

  merged_file <- file.path(output_directory, "additional_point_results.rds")
  merged <- list(
    metadata = list(
      schema_version = "1.0.0",
      generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      total_replicates_per_setting = total_replicates,
      manifest_file = basename(manifest_file),
      manifest_md5 = unname(tools::md5sum(manifest_file)),
      inventory_file = basename(inventory_file),
      approximate_delta = unique(
        validation$configs$approximate_low_rank$third_singular_value
      ),
      approximate_configuration_status = unique(
        validation$configs$approximate_low_rank$status
      ),
      archived_simulation_result_rdata_read = FALSE,
      frozen_snapshot_written = FALSE,
      session_info = capture.output(utils::sessionInfo())
    ),
    manifest = manifest,
    inventory = inventory,
    groups = groups
  )
  additional_merge_save_rds(merged, merged_file, overwrite)

  cat("\nAll expected additional point chunks validated: PASS\n")
  cat("Replicate coverage 1 through", total_replicates, ": PASS\n")
  cat("Approximate standard/MrDAG shared seed streams: PASS\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Frozen snapshot written: no\n")
  cat("Merged output:", merged_file, "\n")
  cat("Additional point-chunk merge: PASS\n")
  invisible(merged)
}


additional_merge_run <- function() {
  arguments <- additional_merge_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "check-chunk",
    "manifest",
    "chunks-dir",
    "output-dir",
    "allow-incomplete",
    "overwrite"
  )
  additional_merge_validate_names(arguments, allowed_arguments)
  allow_incomplete <- additional_merge_parse_boolean(
    arguments,
    "allow-incomplete",
    FALSE
  )
  overwrite <- additional_merge_parse_boolean(arguments, "overwrite", FALSE)
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the merger script path.", call. = FALSE)
  }

  script_file <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  scripts_directory <- dirname(script_file)
  repo_root <- normalizePath(
    file.path(scripts_directory, "..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  core_file <- file.path(scripts_directory, "20_additional_simulation_core.R")
  worker_file <- file.path(
    scripts_directory,
    "22_run_additional_simulation_chunk.R"
  )

  if (!file.exists(file.path(repo_root, "DESCRIPTION")) ||
      !file.exists(core_file) ||
      !file.exists(worker_file)) {
    stop("Could not locate the repository, core, or point worker.", call. = FALSE)
  }

  core <- new.env(parent = globalenv())
  worker <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  sys.source(worker_file, envir = worker)
  chunk_argument <- arguments[["check-chunk"]]

  if (!is.null(chunk_argument) && nzchar(trimws(chunk_argument))) {
    incompatible <- intersect(
      names(arguments),
      c(
        "manifest",
        "chunks-dir",
        "output-dir",
        "allow-incomplete",
        "overwrite"
      )
    )

    if (length(incompatible) > 0L) {
      stop(
        "--check-chunk cannot be combined with full-merge arguments.",
        call. = FALSE
      )
    }

    return(additional_merge_run_standalone(
      chunk_argument,
      repo_root,
      core,
      worker
    ))
  }

  additional_merge_full(
    arguments,
    repo_root,
    core,
    worker,
    allow_incomplete,
    overwrite
  )
}


if (sys.nframe() == 0L) {
  additional_merge_run()
}
