#!/usr/bin/env Rscript

# Validate and merge bootstrap-inference chunks for the additional simulations.
#
# With --allow-incomplete=true, all present expected chunks are validated and
# missing-task inventories are written without producing a merged result.
# A complete run writes merged entrywise bootstrap SE/CI/coverage matrices and
# the median/IQR table summaries used by the manuscript tables.

options(warn = 1)


additional_bootstrap_merge_row_values <- function(row, point_dispatch) {
  integer_value <- point_dispatch$additional_dispatch_integer_value
  text_value <- point_dispatch$additional_dispatch_text_value
  boolean_value <- point_dispatch$additional_dispatch_boolean_value
  phase <- text_value(row, "phase")

  list(
    analysis = text_value(row, "analysis"),
    phase = phase,
    setting_index = integer_value(row, "setting_index", 1L),
    replicate_start = integer_value(row, "replicate_start", 1L),
    replicate_end = integer_value(row, "replicate_end", 1L),
    seed_base = integer_value(row, "seed_base", 1L),
    bootstrap_size = integer_value(row, "bootstrap_size", 2L),
    configured_bootstrap_size = integer_value(
      row,
      "configured_bootstrap_size",
      2L
    ),
    bootstrap_size_override = boolean_value(
      row$bootstrap_size_override,
      "Manifest field bootstrap_size_override"
    ),
    requested_cores = integer_value(row, "requested_cores", 1L),
    allow_provisional = boolean_value(
      row$allow_provisional,
      "Manifest field allow_provisional"
    ),
    mrdag_niter = integer_value(
      row,
      "mrdag_niter",
      minimum = 1L,
      allow_na = phase != "mrdag"
    ),
    mrdag_burnin = integer_value(
      row,
      "mrdag_burnin",
      minimum = 0L,
      allow_na = phase != "mrdag"
    )
  )
}


additional_bootstrap_merge_validate_manifest <- function(
    manifest,
    repo_root,
    core_file,
    point_worker_file,
    worker_file,
    core,
    manifest_generator,
    point_dispatch) {
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

  if (nrow(manifest) == 0L) {
    stop("Bootstrap-task manifest contains no rows.", call. = FALSE)
  }

  task_ids <- suppressWarnings(as.integer(manifest$task_id))
  array_indices <- suppressWarnings(as.integer(manifest$array_index))
  manifest_keys <- paste(
    manifest$resource_class,
    array_indices,
    sep = "/"
  )

  if (anyNA(task_ids) ||
      anyNA(array_indices) ||
      !identical(task_ids, seq_len(nrow(manifest))) ||
      any(!manifest$resource_class %in% c("standard", "mrdag")) ||
      anyDuplicated(manifest_keys) ||
      anyDuplicated(manifest$expected_output_file)) {
    stop("Bootstrap-task manifest keys are invalid.", call. = FALSE)
  }

  for (resource_class in c("standard", "mrdag")) {
    selected_indices <- array_indices[
      manifest$resource_class == resource_class
    ]

    if (length(selected_indices) == 0L ||
        !identical(sort(selected_indices), seq_along(selected_indices))) {
      stop(
        "Bootstrap-task manifest has invalid ",
        resource_class,
        " array indices.",
        call. = FALSE
      )
    }
  }

  expected_core_md5 <- unname(tools::md5sum(core_file))
  expected_point_worker_md5 <- unname(tools::md5sum(point_worker_file))
  expected_worker_md5 <- unname(tools::md5sum(worker_file))
  expected_worker_path <- normalizePath(
    worker_file,
    winslash = "/",
    mustWork = TRUE
  )
  expected_point_worker_path <- normalizePath(
    point_worker_file,
    winslash = "/",
    mustWork = TRUE
  )
  total_candidates <- integer(nrow(manifest))

  for (row_index in seq_len(nrow(manifest))) {
    row <- manifest[row_index, , drop = FALSE]
    values <- additional_bootstrap_merge_row_values(row, point_dispatch)
    resource_class <- point_dispatch$additional_dispatch_text_value(
      row,
      "resource_class"
    )
    configuration_status <- point_dispatch$additional_dispatch_text_value(
      row,
      "configuration_status"
    )
    replicate_count <- point_dispatch$additional_dispatch_integer_value(
      row,
      "replicate_count",
      1L
    )

    if (!identical(
          as.character(row$manifest_schema_version[[1L]]),
          "1.0.0"
        ) ||
        values$replicate_end < values$replicate_start ||
        replicate_count !=
          values$replicate_end - values$replicate_start + 1L ||
        !values$setting_index %in% seq_len(4L) ||
        !values$analysis %in% c(
          "rank_misspecification",
          "approximate_low_rank"
        ) ||
        !values$phase %in% c("standard", "mrdag")) {
      stop("Bootstrap manifest row ", row_index, " is invalid.", call. = FALSE)
    }

    if (!resource_class %in% c("standard", "mrdag") ||
        (resource_class == "standard" && values$phase != "standard") ||
        (resource_class == "mrdag" && values$phase != "mrdag") ||
        (values$analysis == "rank_misspecification" &&
          values$phase != "standard") ||
        (values$phase == "mrdag" &&
          values$analysis != "approximate_low_rank")) {
      stop(
        "Bootstrap resource mapping is invalid in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    if (values$bootstrap_size_override !=
        (values$bootstrap_size != values$configured_bootstrap_size) ||
        (values$bootstrap_size_override && values$bootstrap_size > 20L)) {
      stop(
        "Bootstrap-size metadata is invalid in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    if (values$phase == "mrdag" &&
        values$mrdag_burnin >= values$mrdag_niter) {
      stop("MrDAG settings are invalid in row ", row_index, ".", call. = FALSE)
    }

    if (values$analysis == "rank_misspecification" &&
        values$allow_provisional) {
      stop("Rank tasks cannot be provisional.", call. = FALSE)
    }

    if (configuration_status != "locked" &&
        (values$analysis != "approximate_low_rank" ||
          !values$allow_provisional ||
          replicate_count > 5L)) {
      stop(
        "A provisional bootstrap row is not a small validation task.",
        call. = FALSE
      )
    }

    configuration_file <- point_dispatch$additional_dispatch_resolve_path(
      point_dispatch$additional_dispatch_text_value(
        row,
        "configuration_file"
      ),
      repo_root
    )
    expected_configuration_file <- normalizePath(
      file.path(
        repo_root,
        "paper",
        "config",
        paste0(values$analysis, "_settings.csv")
      ),
      winslash = "/",
      mustWork = FALSE
    )

    if (!identical(configuration_file, expected_configuration_file)) {
      stop("Configuration path is invalid in row ", row_index, ".", call. = FALSE)
    }

    manifest_worker_path <- point_dispatch$additional_dispatch_resolve_path(
      point_dispatch$additional_dispatch_text_value(row, "worker_script"),
      repo_root
    )
    manifest_point_worker_path <-
      point_dispatch$additional_dispatch_resolve_path(
        point_dispatch$additional_dispatch_text_value(
          row,
          "point_worker_script"
        ),
        repo_root
      )

    if (!identical(manifest_worker_path, expected_worker_path) ||
        !identical(
          manifest_point_worker_path,
          expected_point_worker_path
        )) {
      stop("Worker path is invalid in row ", row_index, ".", call. = FALSE)
    }

    expected_output_name <- sprintf(
      "%s__%s__setting-%d__rep-%04d-%04d__bt-%04d.rds",
      values$analysis,
      values$phase,
      values$setting_index,
      values$replicate_start,
      values$replicate_end,
      values$bootstrap_size
    )
    manifest_output_name <- point_dispatch$additional_dispatch_text_value(
      row,
      "expected_output_file"
    )

    if (!identical(manifest_output_name, basename(manifest_output_name)) ||
        !identical(manifest_output_name, expected_output_name)) {
      stop("Output filename is invalid in row ", row_index, ".", call. = FALSE)
    }

    point_dispatch$additional_dispatch_text_value(row, "output_directory")

    point_dispatch$additional_dispatch_check_md5(
      configuration_file,
      row$configuration_md5,
      "Analysis configuration"
    )

    if (tolower(row$core_md5[[1L]]) != tolower(expected_core_md5) ||
        tolower(row$point_worker_md5[[1L]]) !=
          tolower(expected_point_worker_md5) ||
        tolower(row$worker_md5[[1L]]) != tolower(expected_worker_md5)) {
      stop("Script hashes differ in manifest row ", row_index, ".", call. = FALSE)
    }

    config <- core$additional_read_config(
      configuration_file,
      values$analysis
    )
    selected_config <- config[
      config$setting_index == values$setting_index,
      ,
      drop = FALSE
    ]

    if (nrow(selected_config) == 0L ||
        length(unique(selected_config$monte_carlo_replicates)) != 1L) {
      stop("Configuration setting is invalid in row ", row_index, ".", call. = FALSE)
    }

    reference <- selected_config[1L, , drop = FALSE]
    manifest_measurement_error <- suppressWarnings(as.numeric(
      row$measurement_error_weight[[1L]]
    ))
    manifest_genetic_effect <- suppressWarnings(as.numeric(
      row$genetic_effect_weight[[1L]]
    ))
    manifest_delta <- suppressWarnings(as.numeric(
      row$third_singular_value[[1L]]
    ))
    manifest_seed_base <- suppressWarnings(as.integer(
      row$seed_base[[1L]]
    ))

    if (anyNA(c(
          manifest_measurement_error,
          manifest_genetic_effect,
          manifest_delta,
          manifest_seed_base
        )) ||
        !identical(
          as.character(row$scenario[[1L]]),
          as.character(reference$scenario[[1L]])
        ) ||
        abs(
          manifest_measurement_error -
            reference$measurement_error_weight[[1L]]
        ) > 1e-12 ||
        abs(
          manifest_genetic_effect -
            reference$genetic_effect_weight[[1L]]
        ) > 1e-12 ||
        abs(
          manifest_delta - reference$third_singular_value[[1L]]
        ) > 1e-12 ||
        !identical(
          configuration_status,
          as.character(reference$status[[1L]])
        ) ||
        !identical(
          manifest_seed_base,
          as.integer(reference$seed_base[[1L]])
        ) ||
        values$configured_bootstrap_size !=
          as.integer(reference$bootstrap_size[[1L]])) {
      stop(
        "Configuration metadata differ in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    total_candidates[[row_index]] <- unique(
      selected_config$monte_carlo_replicates
    )
  }

  total_replicates <- max(as.integer(manifest$replicate_end))
  bootstrap_sizes <- unique(as.integer(manifest$bootstrap_size))
  configured_bootstrap_sizes <- unique(
    as.integer(manifest$configured_bootstrap_size)
  )

  if (length(bootstrap_sizes) != 1L ||
      anyNA(bootstrap_sizes) ||
      length(configured_bootstrap_sizes) != 1L ||
      anyNA(configured_bootstrap_sizes)) {
    stop(
      "Bootstrap-task groups must use one common bootstrap size.",
      call. = FALSE
    )
  }

  if (total_replicates > min(total_candidates)) {
    stop("Manifest exceeds the configured Monte Carlo total.", call. = FALSE)
  }

  manifest_generator$additional_bootstrap_manifest_validate_coverage(
    manifest,
    total_replicates
  )

  if ((any(manifest$configuration_status != "locked") ||
       bootstrap_sizes != configured_bootstrap_sizes) &&
      total_replicates > 5L) {
    stop(
      "A provisional or reduced-bootstrap manifest cannot exceed 5 replicates.",
      call. = FALSE
    )
  }

  list(
    total_replicates = total_replicates,
    bootstrap_size = bootstrap_sizes,
    group_keys = unique(
      manifest[, c("analysis", "phase", "setting_index")]
    )
  )
}


additional_bootstrap_merge_chunk_path <- function(
    row,
    repo_root,
    chunks_override,
    point_dispatch) {
  filename <- point_dispatch$additional_dispatch_text_value(
    row,
    "expected_output_file"
  )

  if (!is.null(chunks_override)) {
    return(file.path(chunks_override, filename))
  }

  directory <- point_dispatch$additional_dispatch_resolve_path(
    point_dispatch$additional_dispatch_text_value(
      row,
      "output_directory"
    ),
    repo_root
  )
  file.path(directory, filename)
}


additional_bootstrap_merge_validate_chunk <- function(
    chunk_path,
    row,
    repo_root,
    core,
    worker,
    dispatcher,
    point_dispatch,
    point_worker,
    point_merge) {
  values <- additional_bootstrap_merge_row_values(row, point_dispatch)
  chunk <- dispatcher$additional_bootstrap_dispatch_validate_completed_chunk(
    output_file = chunk_path,
    selected = row,
    repo_root = repo_root,
    analysis = values$analysis,
    phase = values$phase,
    setting_index = values$setting_index,
    replicate_start = values$replicate_start,
    replicate_end = values$replicate_end,
    seed_base = values$seed_base,
    bootstrap_size = values$bootstrap_size,
    configured_bootstrap_size = values$configured_bootstrap_size,
    bootstrap_size_override = values$bootstrap_size_override,
    requested_cores = values$requested_cores,
    allow_provisional = values$allow_provisional,
    mrdag_niter = values$mrdag_niter,
    mrdag_burnin = values$mrdag_burnin,
    core = core,
    worker = worker,
    point_dispatch = point_dispatch
  )

  configuration_file <- point_dispatch$additional_dispatch_resolve_path(
    as.character(row$configuration_file[[1L]]),
    repo_root
  )
  config <- core$additional_read_config(configuration_file, values$analysis)
  expected <- worker$additional_bootstrap_result_specification(
    analysis = values$analysis,
    phase = values$phase,
    selected_config = config[
      config$setting_index == values$setting_index, , drop = FALSE
    ],
    point = point_worker
  )
  point_merge$additional_merge_same_object(
    chunk$configuration, expected$config, "Canonical bootstrap configuration"
  )
  point_merge$additional_merge_same_object(
    chunk$result_specification, expected$specification,
    "Canonical bootstrap result specification"
  )
  truth <- core$additional_make_generic_truth(
    third_singular_value = unique(expected$config$third_singular_value),
    px = 9L, py = 3L, seed = 123L
  )
  point_merge$additional_merge_same_object(
    chunk$truth,
    truth[c("C", "singular_values", "numerical_rank")],
    "Canonical bootstrap truth"
  )
  expected_paths <- c(
    "freeze/current_analysis_20260823/project/data/dat_1e-4.csv",
    "freeze/current_analysis_20260823/project/data/rho_mat_1e-4.csv",
    "freeze/current_analysis_20260823/project/scripts/MR_rr_estimators.R",
    "paper/scripts/20_additional_simulation_core.R",
    "paper/scripts/22_run_additional_simulation_chunk.R",
    paste0("paper/config/", values$analysis, "_settings.csv"),
    "paper/scripts/26_run_additional_bootstrap_chunk.R"
  )
  input_paths <- gsub(
    "\\\\", "/", as.character(chunk$input_manifest$repository_path)
  )
  if (anyDuplicated(input_paths) ||
      !setequal(input_paths, expected_paths)) {
    stop("Bootstrap provenance does not contain the exact required inputs.",
         call. = FALSE)
  }
  minimum_successful <- max(
    2L, as.integer(ceiling(0.9 * values$bootstrap_size))
  )
  if (!identical(as.integer(chunk$metadata$minimum_successful_draws),
                 minimum_successful)) {
    stop("Bootstrap minimum-success threshold differs from the worker rule.",
         call. = FALSE)
  }
  invisible(chunk)
}


additional_bootstrap_merge_group <- function(
    manifest_group,
    chunks,
    total_replicates,
    point_merge) {
  ordering <- order(manifest_group$replicate_start)
  manifest_group <- manifest_group[ordering, , drop = FALSE]
  chunks <- chunks[ordering]
  reference <- chunks[[1L]]
  expected_keys <- names(reference$standard_errors)
  replicate_ids <- unlist(lapply(
    chunks,
    function(chunk) as.integer(chunk$metadata$replicate_ids)
  ))

  if (!identical(replicate_ids, seq_len(total_replicates))) {
    stop("Merged bootstrap replicate coverage is invalid.", call. = FALSE)
  }

  for (chunk_index in seq_along(chunks)) {
    chunk <- chunks[[chunk_index]]
    point_merge$additional_merge_same_object(
      chunk$configuration,
      reference$configuration,
      "Bootstrap configuration"
    )
    point_merge$additional_merge_same_object(
      chunk$result_specification,
      reference$result_specification,
      "Bootstrap result specification"
    )
    point_merge$additional_merge_same_object(
      chunk$truth,
      reference$truth,
      "Bootstrap truth"
    )

    if (!identical(names(chunk$standard_errors), expected_keys) ||
        chunk$metadata$bootstrap_size !=
          reference$metadata$bootstrap_size ||
        chunk$metadata$seed_base != reference$metadata$seed_base) {
      stop("Bootstrap chunk schema differs within a group.", call. = FALSE)
    }
  }

  bind_matrix_list <- function(name) {
    stats::setNames(
      lapply(expected_keys, function(key) {
        do.call(cbind, lapply(chunks, function(chunk) chunk[[name]][[key]]))
      }),
      expected_keys
    )
  }
  fit_seeds <- stats::setNames(
    lapply(expected_keys, function(key) {
      do.call(cbind, lapply(chunks, function(chunk) chunk$fit_seeds[[key]]))
    }),
    expected_keys
  )
  successful_draws <- do.call(
    cbind,
    lapply(chunks, function(chunk) chunk$successful_draws)
  )
  errors <- do.call(rbind, lapply(chunks, function(chunk) chunk$errors))
  rownames(errors) <- NULL

  list(
    analysis = reference$metadata$analysis,
    phase = reference$metadata$phase,
    setting_index = reference$metadata$setting_index,
    scenario = reference$metadata$scenario,
    replicate_ids = replicate_ids,
    data_seeds = unlist(lapply(
      chunks,
      function(chunk) as.integer(chunk$metadata$data_seeds)
    )),
    resample_seeds = unlist(lapply(
      chunks,
      function(chunk) as.integer(chunk$metadata$resample_seeds)
    )),
    bootstrap_size = reference$metadata$bootstrap_size,
    configured_bootstrap_size =
      reference$metadata$configured_bootstrap_size,
    configuration = reference$configuration,
    result_specification = reference$result_specification,
    truth = reference$truth,
    fit_seeds = fit_seeds,
    successful_draws = successful_draws,
    standard_errors = bind_matrix_list("standard_errors"),
    ci_lower = bind_matrix_list("ci_lower"),
    ci_upper = bind_matrix_list("ci_upper"),
    coverage = bind_matrix_list("coverage"),
    errors = errors,
    chunk_files = manifest_group$expected_output_file
  )
}


additional_bootstrap_merge_label <- function(result_key) {
  replacements <- c(
    ivw = "IVW",
    srivw = "SRIVW",
    naive_mr_rr = "Naive MR-rr",
    mr_rr = "MR-rr",
    regularized_mr_rr = "Reg. MR-rr",
    mrdag = "MrDAG"
  )
  base_key <- sub("^r[123]__", "", result_key)

  if (!base_key %in% names(replacements)) {
    stop("Unknown bootstrap result key: ", result_key, call. = FALSE)
  }

  unname(replacements[[base_key]])
}


additional_bootstrap_merge_summary <- function(groups) {
  rows <- list()
  row_index <- 0L

  summarize <- function(values, scale, digits) {
    values <- scale * values
    c(
      med = round(stats::median(values), digits),
      q1 = round(unname(stats::quantile(values, 0.25)), digits),
      q3 = round(unname(stats::quantile(values, 0.75)), digits)
    )
  }

  for (group in groups) {
    specification <- group$result_specification

    for (key in names(group$standard_errors)) {
      specification_row <- specification[
        specification$result_key == key,
        ,
        drop = FALSE
      ]
      average_se_by_entry <- rowMeans(group$standard_errors[[key]])
      cp_by_entry <- rowMeans(group$coverage[[key]])
      se <- summarize(average_se_by_entry, scale = 1, digits = 3L)
      cp <- summarize(cp_by_entry, scale = 100, digits = 1L)
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        analysis = group$analysis,
        setting_index = group$setting_index,
        scenario = group$scenario,
        phase = group$phase,
        result_key = key,
        method = specification_row$method[[1L]],
        estimator = additional_bootstrap_merge_label(key),
        working_rank = specification_row$working_rank[[1L]],
        monte_carlo_replicates = length(group$replicate_ids),
        bootstrap_size = group$bootstrap_size,
        configuration_status = unique(group$configuration$status),
        third_singular_value = unique(
          group$configuration$third_singular_value
        ),
        failed_bootstrap_fits = sum(
          group$bootstrap_size - group$successful_draws[key, ]
        ),
        SE_med = unname(se[["med"]]),
        SE_q1 = unname(se[["q1"]]),
        SE_q3 = unname(se[["q3"]]),
        CP_med = unname(cp[["med"]]),
        CP_q1 = unname(cp[["q1"]]),
        CP_q3 = unname(cp[["q3"]]),
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }

  output <- do.call(rbind, rows)
  ordering <- order(
    match(
      output$analysis,
      c("rank_misspecification", "approximate_low_rank")
    ),
    output$setting_index,
    match(output$phase, c("standard", "mrdag")),
    output$working_rank,
    match(
      output$method,
      c(
        "ivw",
        "srivw",
        "naive_mr_rr",
        "mr_rr",
        "regularized_mr_rr",
        "mrdag"
      )
    )
  )
  output <- output[ordering, , drop = FALSE]
  rownames(output) <- NULL
  output
}


additional_merge_bootstrap_chunks <- function() {
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the merger script path.", call. = FALSE)
  }

  merger_file <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  scripts_directory <- dirname(merger_file)
  dependency_names <- c(
    core = "20_additional_simulation_core.R",
    point_worker = "22_run_additional_simulation_chunk.R",
    point_dispatch = "24_dispatch_additional_point_task.R",
    point_merge = "25_merge_additional_point_chunks.R",
    worker = "26_run_additional_bootstrap_chunk.R",
    manifest_generator = "27_make_additional_bootstrap_manifest.R",
    dispatcher = "28_dispatch_additional_bootstrap_task.R"
  )
  dependencies <- stats::setNames(
    file.path(scripts_directory, dependency_names),
    names(dependency_names)
  )
  missing <- dependencies[!file.exists(dependencies)]

  if (length(missing) > 0L) {
    stop(
      "Bootstrap merger dependencies are missing: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  environments <- lapply(dependencies, function(path) {
    environment <- new.env(parent = globalenv())
    sys.source(path, envir = environment)
    environment
  })
  names(environments) <- names(dependencies)
  core <- environments$core
  point_worker <- environments$point_worker
  point_dispatch <- environments$point_dispatch
  point_merge <- environments$point_merge
  worker <- environments$worker
  manifest_generator <- environments$manifest_generator
  dispatcher <- environments$dispatcher
  repo_root <- core$additional_find_repo_root(dirname(merger_file))
  arguments <- point_merge$additional_merge_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "manifest",
    "chunks-dir",
    "output-dir",
    "allow-incomplete",
    "overwrite"
  )
  point_merge$additional_merge_validate_names(arguments, allowed_arguments)
  allow_incomplete <- point_merge$additional_merge_parse_boolean(
    arguments,
    "allow-incomplete",
    FALSE
  )
  overwrite <- point_merge$additional_merge_parse_boolean(
    arguments,
    "overwrite",
    FALSE
  )
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
    point_merge$additional_merge_resolve_path(
      manifest_argument,
      repo_root
    )
  }
  chunks_argument <- arguments[["chunks-dir"]]
  chunks_override <- if (is.null(chunks_argument) ||
      !nzchar(trimws(chunks_argument))) {
    NULL
  } else {
    point_merge$additional_merge_resolve_path(
      chunks_argument,
      repo_root,
      must_work = TRUE
    )
  }
  output_argument <- arguments[["output-dir"]]
  output_directory <- if (is.null(output_argument) ||
      !nzchar(trimws(output_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "bootstrap_merged"
    )
  } else {
    point_merge$additional_merge_resolve_path(
      output_argument,
      repo_root
    )
  }

  if (!file.exists(manifest_file)) {
    stop("Bootstrap manifest is missing: ", manifest_file, call. = FALSE)
  }

  frozen_directory <- normalizePath(
    file.path(repo_root, "freeze"), winslash = "/", mustWork = TRUE
  )
  output_directory <- normalizePath(
    output_directory, winslash = "/", mustWork = FALSE
  )
  if (identical(tolower(output_directory), tolower(frozen_directory)) ||
      startsWith(tolower(output_directory),
                 paste0(tolower(frozen_directory), "/"))) {
    stop("Bootstrap outputs must not be written under freeze/.", call. = FALSE)
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
  validation <- additional_bootstrap_merge_validate_manifest(
    manifest = manifest,
    repo_root = repo_root,
    core_file = dependencies[["core"]],
    point_worker_file = dependencies[["point_worker"]],
    worker_file = dependencies[["worker"]],
    core = core,
    manifest_generator = manifest_generator,
    point_dispatch = point_dispatch
  )
  chunk_paths <- vapply(
    seq_len(nrow(manifest)),
    function(row_index) {
      additional_bootstrap_merge_chunk_path(
        manifest[row_index, , drop = FALSE],
        repo_root,
        chunks_override,
        point_dispatch
      )
    },
    character(1)
  )
  present <- file.exists(chunk_paths)
  missing_tasks <- manifest[!present, , drop = FALSE]
  missing_file <- file.path(
    output_directory,
    "missing_additional_bootstrap_tasks.csv"
  )
  point_merge$additional_merge_write_csv(missing_tasks, missing_file)
  inventory <- manifest
  inventory$chunk_path <- chunk_paths
  inventory$present <- present
  inventory$validated <- FALSE
  inventory$file_size_bytes <- NA_real_
  inventory$md5 <- NA_character_
  inventory$complete <- NA
  inventory$all_bootstrap_draws_successful <- NA
  inventory$estimator_failure_count <- NA_integer_
  chunks <- vector("list", nrow(manifest))

  cat("MR-rr additional bootstrap-chunk merger\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Expected chunks:", nrow(manifest), "\n")
  cat("Present chunks:", sum(present), "\n")
  cat("Missing chunks:", sum(!present), "\n")

  for (row_index in which(present)) {
    chunk <- additional_bootstrap_merge_validate_chunk(
      chunk_path = chunk_paths[[row_index]],
      row = manifest[row_index, , drop = FALSE],
      repo_root = repo_root,
      core = core,
      worker = worker,
      dispatcher = dispatcher,
      point_dispatch = point_dispatch,
      point_worker = point_worker,
      point_merge = point_merge
    )
    chunks[[row_index]] <- chunk
    information <- file.info(chunk_paths[[row_index]])
    inventory$validated[[row_index]] <- TRUE
    inventory$file_size_bytes[[row_index]] <- information$size
    inventory$md5[[row_index]] <- unname(
      tools::md5sum(chunk_paths[[row_index]])
    )
    inventory$complete[[row_index]] <- chunk$metadata$complete
    inventory$all_bootstrap_draws_successful[[row_index]] <-
      chunk$metadata$all_bootstrap_draws_successful
    inventory$estimator_failure_count[[row_index]] <- nrow(chunk$errors)
  }

  inventory_file <- file.path(
    output_directory,
    "additional_bootstrap_chunk_inventory.csv"
  )
  point_merge$additional_merge_write_csv(inventory, inventory_file)

  if (nrow(missing_tasks) > 0L) {
    if (allow_incomplete) {
      cat("Existing expected bootstrap chunks validated: PASS\n")
      cat("Missing-task inventory:", missing_file, "\n")
      cat("Full bootstrap merge deferred until every task is present.\n")
      cat("Incomplete additional bootstrap inventory: PASS\n")
      return(invisible(inventory))
    }

    stop(
      "Additional bootstrap output is incomplete. See: ",
      missing_file,
      "\nUse --allow-incomplete=true only for a progress check.",
      call. = FALSE
    )
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
    groups[[group_index]] <- additional_bootstrap_merge_group(
      manifest[selected_indices, , drop = FALSE],
      chunks[selected_indices],
      validation$total_replicates,
      point_merge
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
        !identical(standard$data_seeds, mrdag$data_seeds) ||
        !identical(standard$resample_seeds, mrdag$resample_seeds)) {
      stop(
        "Approximate standard and MrDAG seed streams differ for setting ",
        setting_index,
        ".",
        call. = FALSE
      )
    }

    point_merge$additional_merge_same_object(
      standard$truth,
      mrdag$truth,
      paste("Approximate-low-rank setting", setting_index, "truth")
    )
  }

  table_summary <- additional_bootstrap_merge_summary(groups)
  summary_file <- file.path(
    output_directory,
    "additional_bootstrap_summary.csv"
  )
  merged_file <- file.path(
    output_directory,
    "additional_bootstrap_results.rds"
  )

  if ((file.exists(summary_file) || file.exists(merged_file)) && !overwrite) {
    stop(
      "Merged bootstrap output exists; use --overwrite=true to replace it.",
      call. = FALSE
    )
  }

  merged <- list(
    metadata = list(
      schema_version = "1.0.0",
      generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      total_replicates_per_setting = validation$total_replicates,
      bootstrap_size = validation$bootstrap_size,
      manifest_file = basename(manifest_file),
      manifest_md5 = unname(tools::md5sum(manifest_file)),
      inventory_file = basename(inventory_file),
      summary_file = basename(summary_file),
      approximate_delta = unique(
        manifest$third_singular_value[
          manifest$analysis == "approximate_low_rank"
        ]
      ),
      approximate_configuration_status = unique(
        manifest$configuration_status[
          manifest$analysis == "approximate_low_rank"
        ]
      ),
      point_result_files_read = FALSE,
      archived_simulation_result_rdata_read = FALSE,
      frozen_snapshot_written = FALSE,
      session_info = capture.output(utils::sessionInfo())
    ),
    manifest = manifest,
    inventory = inventory,
    groups = groups,
    table_summary = table_summary
  )
  point_merge$additional_merge_write_csv(table_summary, summary_file)
  point_merge$additional_merge_save_rds(merged, merged_file, overwrite)

  cat("\nAll expected additional bootstrap chunks validated: PASS\n")
  cat(
    "Replicate coverage 1 through",
    validation$total_replicates,
    ": PASS\n"
  )
  cat("Approximate standard/MrDAG common resamples: PASS\n")
  cat("Recorded bootstrap fit failures:",
      sum(table_summary$failed_bootstrap_fits), "\n")
  if (sum(table_summary$failed_bootstrap_fits) > 0L) {
    cat("Inspect failure counts before interpreting SE and coverage.\n")
  }
  cat("Point-result files read: no\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Frozen snapshot written: no\n")
  cat("Bootstrap summary:", summary_file, "\n")
  cat("Merged output:", merged_file, "\n")
  cat("Additional bootstrap-chunk merge: PASS\n")
  invisible(merged)
}


if (sys.nframe() == 0L) {
  additional_merge_bootstrap_chunks()
}
