#!/usr/bin/env Rscript

# Validate and merge portable MR-rr simulation chunks.
#
# This script has two modes:
#
# 1. --check-chunk=<path>
#    Validate one chunk produced by 16_run_simulation_chunk.R.
#
# 2. Full merge mode (the default)
#    Read paper/config/full_simulation_tasks.csv, require and validate every
#    expected chunk, verify complete replicate coverage and shared seed streams,
#    and write legacy-compatible result objects under
#    paper/output/full_run/merged.
#
# No archived simulation-result RData is read. The immutable freeze is used by
# the worker only for raw calibration inputs and the legacy estimator source.

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


resolve_repo_path <- function(path, repo_root, must_work = FALSE) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    candidate <- path
  } else {
    candidate <- file.path(repo_root, path)
  }

  normalizePath(
    candidate,
    winslash = "/",
    mustWork = must_work
  )
}


assert_scalar <- function(value, label) {
  if (length(value) != 1L || is.na(value)) {
    stop(label, " must be one non-missing value.", call. = FALSE)
  }

  invisible(value)
}


assert_identical_value <- function(actual, expected, label) {
  if (length(actual) != 1L || is.na(actual) ||
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


assert_integer_value <- function(actual, expected, label) {
  if (length(actual) != 1L || is.na(actual) ||
      as.integer(actual) != as.integer(expected)) {
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


assert_close <- function(
    actual,
    expected,
    label,
    tolerance = 1e-10) {
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


make_replicate_seed <- function(
    seed_base,
    design_index,
    setting_index,
    phase,
    replicate_id) {
  stream_index <- if (identical(phase, "sparse")) 2L else 1L

  seed <-
    as.double(seed_base) +
    (design_index - 1L) * 10000000 +
    (stream_index - 1L) * 5000000 +
    (setting_index - 1L) * 100000 +
    replicate_id

  if (!is.finite(seed) || seed > .Machine$integer.max) {
    stop(
      "A derived replicate seed exceeds the R integer range.",
      call. = FALSE
    )
  }

  as.integer(seed)
}


expected_methods <- function(phase) {
  switch(
    phase,
    main_no_mrdag = c(
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr",
      "ivw",
      "srivw"
    ),
    main_mrdag = "mrdag",
    sparse = "sparse_mr_rr",
    stop("Unknown simulation phase: ", phase, call. = FALSE)
  )
}


expected_chunk_filename <- function(
    design,
    phase,
    setting,
    replicate_start,
    replicate_end) {
  sprintf(
    "%s__%s__setting-%d__rep-%04d-%04d.rds",
    design,
    phase,
    as.integer(setting),
    as.integer(replicate_start),
    as.integer(replicate_end)
  )
}


validate_manifest <- function(manifest) {
  required_columns <- c(
    "task_id",
    "resource_class",
    "array_index",
    "design",
    "phase",
    "setting",
    "scenario",
    "me_weight",
    "effect_weight",
    "replicate_start",
    "replicate_end",
    "replicate_count",
    "seed_base",
    "mrdag_niter",
    "mrdag_burnin",
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

  if (nrow(manifest) == 0L) {
    stop("Task manifest contains no rows.", call. = FALSE)
  }

  if (anyDuplicated(manifest$task_id) ||
      anyDuplicated(manifest$expected_output_file)) {
    stop(
      "Task manifest contains duplicate task IDs or output files.",
      call. = FALSE
    )
  }

  allowed_designs <- c("generic_low_rank", "sparse_loading")
  allowed_phases <- c("main_no_mrdag", "sparse", "main_mrdag")
  expected_settings <- data.frame(
    setting = seq_len(4L),
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

  if (any(!manifest$design %in% allowed_designs) ||
      any(!manifest$phase %in% allowed_phases) ||
      any(!manifest$setting %in% expected_settings$setting)) {
    stop(
      "Task manifest contains an unknown design, phase, or setting.",
      call. = FALSE
    )
  }

  expected_resource <- ifelse(
    manifest$phase == "main_mrdag",
    "mrdag",
    "standard"
  )

  if (any(manifest$resource_class != expected_resource)) {
    stop(
      "A task manifest row has the wrong resource class for its phase.",
      call. = FALSE
    )
  }

  for (row_index in seq_len(nrow(manifest))) {
    row <- manifest[row_index, , drop = FALSE]
    expected_setting <- expected_settings[
      expected_settings$setting == row$setting,
      ,
      drop = FALSE
    ]

    if (row$scenario != expected_setting$scenario ||
        abs(row$me_weight - expected_setting$me_weight) > 1e-12 ||
        abs(row$effect_weight - expected_setting$effect_weight) > 1e-12) {
      stop(
        "Task manifest setting metadata is inconsistent in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    if (is.na(row$replicate_start) || is.na(row$replicate_end) ||
        row$replicate_start < 1L ||
        row$replicate_end < row$replicate_start ||
        row$replicate_count !=
          row$replicate_end - row$replicate_start + 1L) {
      stop(
        "Task manifest has an invalid replicate range in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    expected_file <- expected_chunk_filename(
      design = row$design,
      phase = row$phase,
      setting = row$setting,
      replicate_start = row$replicate_start,
      replicate_end = row$replicate_end
    )

    if (row$expected_output_file != expected_file) {
      stop(
        "Task manifest output filename is inconsistent in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }

    if (row$phase == "main_mrdag" &&
        (is.na(row$mrdag_niter) || is.na(row$mrdag_burnin) ||
        row$mrdag_burnin < 0L ||
        row$mrdag_burnin >= row$mrdag_niter)) {
      stop(
        "Task manifest has invalid MrDAG iterations in row ",
        row_index,
        ".",
        call. = FALSE
      )
    }
  }

  group_keys <- unique(
    manifest[, c("design", "phase", "setting")]
  )
  expected_group_count <-
    length(allowed_designs) * length(allowed_phases) * 4L

  if (nrow(group_keys) != expected_group_count) {
    stop(
      "Task manifest does not contain all 24 design/phase/setting groups.",
      call. = FALSE
    )
  }

  total_replicates <- unique(
    vapply(
      seq_len(nrow(group_keys)),
      function(group_index) {
        group <- group_keys[group_index, , drop = FALSE]
        selected <- manifest[
          manifest$design == group$design &
          manifest$phase == group$phase &
          manifest$setting == group$setting,
          ,
          drop = FALSE
        ]
        selected <- selected[
          order(selected$replicate_start),
          ,
          drop = FALSE
        ]

        if (selected$replicate_start[[1L]] != 1L ||
            any(
              selected$replicate_start[-1L] !=
              selected$replicate_end[-nrow(selected)] + 1L
            )) {
          stop(
            "Task manifest has a gap or overlap for ",
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

        max(selected$replicate_end)
      },
      integer(1)
    )
  )

  if (length(total_replicates) != 1L) {
    stop(
      "Task groups do not share one total replicate count.",
      call. = FALSE
    )
  }

  list(
    total_replicates = total_replicates,
    group_keys = group_keys
  )
}


validate_matrix_component <- function(
    matrix_value,
    expected_rows,
    replicate_ids,
    label) {
  expected_dimensions <- c(
    as.integer(expected_rows),
    length(replicate_ids)
  )

  if (!is.matrix(matrix_value) ||
      !identical(dim(matrix_value), expected_dimensions)) {
    stop(label, " has unexpected dimensions.", call. = FALSE)
  }

  if (any(!is.finite(matrix_value))) {
    stop(label, " contains a non-finite value.", call. = FALSE)
  }

  expected_column_names <- as.character(replicate_ids)

  if (!is.null(colnames(matrix_value)) &&
      !identical(colnames(matrix_value), expected_column_names)) {
    stop(label, " has incorrect replicate column names.", call. = FALSE)
  }

  invisible(matrix_value)
}


validate_chunk <- function(
    chunk,
    chunk_file,
    expected_row = NULL,
    tolerance = 1e-8) {
  required_components <- c(
    "metadata",
    "input_manifest",
    "setting",
    "regularization_rate",
    "sparse_lambda",
    "parameters",
    "truth",
    "methods",
    "estimates",
    "biases",
    "predictions",
    "sparse_B",
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
    "design",
    "design_index",
    "phase",
    "setting_index",
    "scenario",
    "replicate_start",
    "replicate_end",
    "replicate_ids",
    "replicate_count",
    "seed_base",
    "replicate_seeds",
    "rng_kind",
    "mrdag_niter",
    "mrdag_burnin",
    "complete",
    "archived_result_rdata_read"
  )
  missing_metadata <- setdiff(required_metadata, names(metadata))

  if (length(missing_metadata) > 0L) {
    stop(
      basename(chunk_file),
      " is missing metadata: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }

  assert_identical_value(
    metadata$schema_version,
    "1.0.0",
    "Chunk schema version"
  )

  allowed_designs <- c("generic_low_rank", "sparse_loading")
  allowed_phases <- c("main_no_mrdag", "main_mrdag", "sparse")

  if (!metadata$design %in% allowed_designs ||
      !metadata$phase %in% allowed_phases ||
      !metadata$setting_index %in% seq_len(4L)) {
    stop(
      basename(chunk_file),
      " has an unknown design, phase, or setting.",
      call. = FALSE
    )
  }

  if (!isTRUE(metadata$complete)) {
    stop(basename(chunk_file), " is not marked complete.", call. = FALSE)
  }

  if (!identical(metadata$archived_result_rdata_read, FALSE)) {
    stop(
      basename(chunk_file),
      " reports reading archived simulation-result RData.",
      call. = FALSE
    )
  }

  if (!is.data.frame(chunk$errors) || nrow(chunk$errors) != 0L) {
    stop(
      basename(chunk_file),
      " contains estimator errors.",
      call. = FALSE
    )
  }

  replicate_ids <- as.integer(metadata$replicate_ids)
  expected_ids <- seq.int(
    as.integer(metadata$replicate_start),
    as.integer(metadata$replicate_end)
  )

  if (!identical(replicate_ids, expected_ids) ||
      metadata$replicate_count != length(expected_ids)) {
    stop(
      basename(chunk_file),
      " has inconsistent replicate metadata.",
      call. = FALSE
    )
  }

  expected_design_index <- match(metadata$design, allowed_designs)
  assert_integer_value(
    metadata$design_index,
    expected_design_index,
    "Chunk design index"
  )

  expected_seeds <- vapply(
    replicate_ids,
    function(replicate_id) {
      make_replicate_seed(
        seed_base = as.integer(metadata$seed_base),
        design_index = expected_design_index,
        setting_index = as.integer(metadata$setting_index),
        phase = metadata$phase,
        replicate_id = replicate_id
      )
    },
    integer(1)
  )

  if (!identical(as.integer(metadata$replicate_seeds), expected_seeds) ||
      anyDuplicated(metadata$replicate_seeds)) {
    stop(
      basename(chunk_file),
      " has inconsistent deterministic replicate seeds.",
      call. = FALSE
    )
  }

  if (!identical(
    as.character(metadata$rng_kind),
    c("L'Ecuyer-CMRG", "Inversion", "Rejection")
  )) {
    stop(
      basename(chunk_file),
      " has unexpected RNG metadata.",
      call. = FALSE
    )
  }

  methods <- expected_methods(metadata$phase)

  if (!identical(as.character(chunk$methods), methods) ||
      !identical(names(chunk$estimates), methods) ||
      !identical(names(chunk$biases), methods) ||
      !identical(names(chunk$predictions), methods)) {
    stop(
      basename(chunk_file),
      " has unexpected method components.",
      call. = FALSE
    )
  }

  parameters <- chunk$parameters

  if (!is.list(parameters) ||
      !identical(as.integer(parameters$py), 3L) ||
      !identical(as.integer(parameters$px), 9L) ||
      !identical(as.integer(parameters$r_RR), 2L) ||
      !is.matrix(parameters$C) ||
      !identical(dim(parameters$C), c(3L, 9L)) ||
      any(!is.finite(parameters$C))) {
    stop(
      basename(chunk_file),
      " has invalid simulation parameters.",
      call. = FALSE
    )
  }

  if (length(parameters$iv_strength) != 1L ||
      !is.finite(parameters$iv_strength) ||
      parameters$iv_strength <= 0) {
    stop(
      basename(chunk_file),
      " has invalid population instrument strength.",
      call. = FALSE
    )
  }

  if (!is.list(chunk$truth) ||
      !is.matrix(chunk$truth$C) ||
      !identical(dim(chunk$truth$C), c(3L, 9L)) ||
      length(chunk$truth$fixed_exposure) != 9L ||
      length(chunk$truth$true_prediction) != 3L ||
      any(!is.finite(c(
        chunk$truth$C,
        chunk$truth$fixed_exposure,
        chunk$truth$true_prediction
      )))) {
    stop(
      basename(chunk_file),
      " has invalid truth metadata.",
      call. = FALSE
    )
  }

  assert_close(
    chunk$truth$C,
    parameters$C,
    "Chunk truth C",
    tolerance = tolerance
  )
  assert_close(
    as.numeric(chunk$truth$true_prediction),
    as.numeric(parameters$C %*% chunk$truth$fixed_exposure),
    "Chunk true prediction",
    tolerance = tolerance
  )

  for (method in methods) {
    validate_matrix_component(
      chunk$estimates[[method]],
      expected_rows = 27L,
      replicate_ids = replicate_ids,
      label = paste(basename(chunk_file), method, "estimates")
    )
    validate_matrix_component(
      chunk$biases[[method]],
      expected_rows = 27L,
      replicate_ids = replicate_ids,
      label = paste(basename(chunk_file), method, "biases")
    )
    validate_matrix_component(
      chunk$predictions[[method]],
      expected_rows = 3L,
      replicate_ids = replicate_ids,
      label = paste(basename(chunk_file), method, "predictions")
    )

    expected_bias <- sweep(
      chunk$estimates[[method]],
      1L,
      as.vector(parameters$C),
      FUN = "-"
    )
    assert_close(
      chunk$biases[[method]],
      expected_bias,
      paste(basename(chunk_file), method, "bias reconstruction"),
      tolerance = tolerance
    )

    reconstructed_prediction <- vapply(
      seq_len(ncol(chunk$estimates[[method]])),
      function(column_index) {
        effect <- matrix(
          chunk$estimates[[method]][, column_index],
          nrow = 3L,
          ncol = 9L
        )
        as.numeric(effect %*% chunk$truth$fixed_exposure)
      },
      numeric(3L)
    )
    reconstructed_prediction <- matrix(
      reconstructed_prediction,
      nrow = 3L,
      ncol = length(replicate_ids)
    )
    assert_close(
      chunk$predictions[[method]],
      reconstructed_prediction,
      paste(basename(chunk_file), method, "prediction reconstruction"),
      tolerance = tolerance
    )
  }

  if (metadata$phase == "sparse") {
    validate_matrix_component(
      chunk$sparse_B,
      expected_rows = 18L,
      replicate_ids = replicate_ids,
      label = paste(basename(chunk_file), "sparse B")
    )
  } else if (!is.null(chunk$sparse_B)) {
    stop(
      basename(chunk_file),
      " has sparse B output outside the sparse phase.",
      call. = FALSE
    )
  }

  if (!identical(names(chunk$successful_replicates), methods) ||
      any(as.integer(chunk$successful_replicates) != length(replicate_ids))) {
    stop(
      basename(chunk_file),
      " has inconsistent success counts.",
      call. = FALSE
    )
  }

  required_input_columns <- c("role", "repository_path", "md5")

  if (!is.data.frame(chunk$input_manifest) ||
      !all(required_input_columns %in% names(chunk$input_manifest)) ||
      nrow(chunk$input_manifest) == 0L ||
      anyNA(chunk$input_manifest[, required_input_columns]) ||
      any(!nzchar(chunk$input_manifest$md5))) {
    stop(
      basename(chunk_file),
      " has invalid input provenance.",
      call. = FALSE
    )
  }

  result_path_detected <- grepl(
    "(^|[/\\\\])results([/\\\\]|$)|\\.[Rr][Dd]ata$",
    chunk$input_manifest$repository_path
  )

  if (any(result_path_detected)) {
    stop(
      basename(chunk_file),
      " lists an archived result as an input.",
      call. = FALSE
    )
  }

  if (!is.null(expected_row)) {
    assert_identical_value(
      metadata$design,
      expected_row$design,
      "Chunk design"
    )
    assert_identical_value(
      metadata$phase,
      expected_row$phase,
      "Chunk phase"
    )
    assert_integer_value(
      metadata$setting_index,
      expected_row$setting,
      "Chunk setting"
    )
    assert_identical_value(
      metadata$scenario,
      expected_row$scenario,
      "Chunk scenario"
    )
    assert_integer_value(
      metadata$replicate_start,
      expected_row$replicate_start,
      "Chunk replicate start"
    )
    assert_integer_value(
      metadata$replicate_end,
      expected_row$replicate_end,
      "Chunk replicate end"
    )
    assert_integer_value(
      metadata$replicate_count,
      expected_row$replicate_count,
      "Chunk replicate count"
    )
    assert_integer_value(
      metadata$seed_base,
      expected_row$seed_base,
      "Chunk seed base"
    )

    if (basename(chunk_file) != expected_row$expected_output_file) {
      stop(
        "Chunk filename differs from its manifest row: ",
        basename(chunk_file),
        ".",
        call. = FALSE
      )
    }

    if (metadata$phase == "main_mrdag") {
      assert_integer_value(
        metadata$mrdag_niter,
        expected_row$mrdag_niter,
        "Chunk MrDAG niter"
      )
      assert_integer_value(
        metadata$mrdag_burnin,
        expected_row$mrdag_burnin,
        "Chunk MrDAG burnin"
      )
    }

    if (nrow(chunk$setting) != 1L ||
        chunk$setting$setting_index != expected_row$setting ||
        chunk$setting$scenario != expected_row$scenario ||
        abs(chunk$setting$me_weight - expected_row$me_weight) > 1e-12 ||
        abs(chunk$setting$effect_weight -
          expected_row$effect_weight) > 1e-12) {
      stop(
        basename(chunk_file),
        " has setting metadata inconsistent with the manifest.",
        call. = FALSE
      )
    }
  }

  invisible(chunk)
}


same_object <- function(actual, expected, label) {
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


merge_chunk_group <- function(
    manifest_group,
    chunk_directory,
    total_replicates) {
  manifest_group <- manifest_group[
    order(manifest_group$replicate_start),
    ,
    drop = FALSE
  ]
  chunks <- vector("list", nrow(manifest_group))

  for (row_index in seq_len(nrow(manifest_group))) {
    row <- manifest_group[row_index, , drop = FALSE]
    chunk_file <- file.path(
      chunk_directory,
      row$expected_output_file
    )
    chunk <- readRDS(chunk_file)
    validate_chunk(chunk, chunk_file, expected_row = row)
    chunks[[row_index]] <- chunk
  }

  reference <- chunks[[1L]]

  for (chunk_index in seq_along(chunks)[-1L]) {
    chunk <- chunks[[chunk_index]]
    same_object(
      chunk$input_manifest,
      reference$input_manifest,
      "Input provenance"
    )
    same_object(
      chunk$parameters,
      reference$parameters,
      "Simulation parameters"
    )
    same_object(chunk$truth, reference$truth, "Simulation truth")
    same_object(
      chunk$regularization_rate,
      reference$regularization_rate,
      "Regularization rate"
    )
    same_object(
      chunk$sparse_lambda,
      reference$sparse_lambda,
      "Sparse tuning parameter"
    )
  }

  replicate_ids <- unlist(
    lapply(chunks, function(chunk) chunk$metadata$replicate_ids),
    use.names = FALSE
  )
  replicate_seeds <- unlist(
    lapply(chunks, function(chunk) chunk$metadata$replicate_seeds),
    use.names = FALSE
  )

  if (!identical(as.integer(replicate_ids), seq_len(total_replicates)) ||
      anyDuplicated(replicate_ids) ||
      anyDuplicated(replicate_seeds)) {
    stop(
      "Merged chunks have a gap, overlap, or duplicate seed for ",
      paste(
        reference$metadata$design,
        reference$metadata$phase,
        reference$metadata$setting_index,
        sep = "/"
      ),
      ".",
      call. = FALSE
    )
  }

  methods <- reference$methods
  merge_method_matrices <- function(component) {
    stats::setNames(
      lapply(
        methods,
        function(method) {
          merged <- do.call(
            cbind,
            lapply(chunks, function(chunk) chunk[[component]][[method]])
          )
          colnames(merged) <- as.character(replicate_ids)
          merged
        }
      ),
      methods
    )
  }

  sparse_B <- if (reference$metadata$phase == "sparse") {
    merged <- do.call(
      cbind,
      lapply(chunks, function(chunk) chunk$sparse_B)
    )
    colnames(merged) <- as.character(replicate_ids)
    merged
  } else {
    NULL
  }

  list(
    design = reference$metadata$design,
    phase = reference$metadata$phase,
    setting_index = reference$metadata$setting_index,
    scenario = reference$metadata$scenario,
    setting = reference$setting,
    seed_base = reference$metadata$seed_base,
    replicate_ids = as.integer(replicate_ids),
    replicate_seeds = as.integer(replicate_seeds),
    input_manifest = reference$input_manifest,
    regularization_rate = reference$regularization_rate,
    sparse_lambda = reference$sparse_lambda,
    parameters = reference$parameters,
    truth = reference$truth,
    methods = methods,
    estimates = merge_method_matrices("estimates"),
    biases = merge_method_matrices("biases"),
    predictions = merge_method_matrices("predictions"),
    sparse_B = sparse_B,
    chunk_files = manifest_group$expected_output_file
  )
}


write_csv_atomically <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path),
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)

  utils::write.csv(
    data,
    temporary_file,
    row.names = FALSE,
    na = ""
  )

  if (file.exists(path) && !file.remove(path)) {
    stop("Could not replace CSV file: ", path, call. = FALSE)
  }

  if (!file.rename(temporary_file, path)) {
    stop("Could not move completed CSV to: ", path, call. = FALSE)
  }

  invisible(path)
}


save_objects_atomically <- function(objects, path, overwrite = FALSE) {
  if (file.exists(path) && !overwrite) {
    stop(
      "Merged output already exists: ",
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
  environment <- list2env(objects, parent = emptyenv())

  save(
    list = names(objects),
    file = temporary_file,
    envir = environment,
    compress = "xz"
  )

  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = dirname(path),
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop(
        "Could not preserve existing merged output: ",
        path,
        call. = FALSE
      )
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop("Could not move completed output to: ", path, call. = FALSE)
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


make_named_scenario_list <- function(groups, component, method) {
  scenarios <- c(
    "me_2.5_effect_0.25",
    "me_2.5_effect_1",
    "me_1_effect_0.25",
    "me_1_effect_1"
  )

  stats::setNames(
    lapply(
      scenarios,
      function(scenario) {
        selected <- groups[
          vapply(
            groups,
            function(group) identical(group$scenario, scenario),
            logical(1)
          )
        ]

        if (length(selected) != 1L) {
          stop(
            "Expected exactly one merged group for scenario ",
            scenario,
            ".",
            call. = FALSE
          )
        }

        selected[[1L]][[component]][[method]]
      }
    ),
    scenarios
  )
}


build_legacy_main_result <- function(main_groups, mrdag_groups) {
  scenarios <- c(
    "me_2.5_effect_0.25",
    "me_2.5_effect_1",
    "me_1_effect_0.25",
    "me_1_effect_1"
  )
  reference_groups <- stats::setNames(
    lapply(
      scenarios,
      function(scenario) {
        selected <- main_groups[
          vapply(
            main_groups,
            function(group) identical(group$scenario, scenario),
            logical(1)
          )
        ]

        if (length(selected) != 1L) {
          stop(
            "Expected one main group for scenario ",
            scenario,
            ".",
            call. = FALSE
          )
        }

        selected[[1L]]
      }
    ),
    scenarios
  )

  list(
    result_AB_list = make_named_scenario_list(
      main_groups, "biases", "naive_mr_rr"
    ),
    result_AB_d_list = make_named_scenario_list(
      main_groups, "biases", "mr_rr"
    ),
    result_AB_d_r_list = make_named_scenario_list(
      main_groups, "biases", "regularized_mr_rr"
    ),
    iv_strength_list = unname(vapply(
      reference_groups,
      function(group) group$parameters$iv_strength,
      numeric(1)
    )),
    parameters_list = unname(lapply(
      reference_groups,
      function(group) group$parameters
    )),
    result_C_ivw_list = make_named_scenario_list(
      main_groups, "biases", "ivw"
    ),
    result_C_adivw_list = make_named_scenario_list(
      main_groups, "biases", "srivw"
    ),
    result_MrDAG_list = make_named_scenario_list(
      mrdag_groups, "biases", "mrdag"
    ),
    AB_list = make_named_scenario_list(
      main_groups, "estimates", "naive_mr_rr"
    ),
    AB_d_list = make_named_scenario_list(
      main_groups, "estimates", "mr_rr"
    ),
    AB_d_r_list = make_named_scenario_list(
      main_groups, "estimates", "regularized_mr_rr"
    ),
    C_ivw_list = make_named_scenario_list(
      main_groups, "estimates", "ivw"
    ),
    C_adivw_list = make_named_scenario_list(
      main_groups, "estimates", "srivw"
    ),
    MrDAG_list = make_named_scenario_list(
      mrdag_groups, "estimates", "mrdag"
    ),
    Y_pred_AB_list = make_named_scenario_list(
      main_groups, "predictions", "naive_mr_rr"
    ),
    Y_pred_AB_d_list = make_named_scenario_list(
      main_groups, "predictions", "mr_rr"
    ),
    Y_pred_AB_d_r_list = make_named_scenario_list(
      main_groups, "predictions", "regularized_mr_rr"
    ),
    Y_pred_C_ivw_list = make_named_scenario_list(
      main_groups, "predictions", "ivw"
    ),
    Y_pred_C_adivw_list = make_named_scenario_list(
      main_groups, "predictions", "srivw"
    ),
    Y_pred_MrDAG_list = make_named_scenario_list(
      mrdag_groups, "predictions", "mrdag"
    )
  )
}


build_legacy_sparse_result <- function(sparse_groups) {
  scenarios <- c(
    "me_2.5_effect_0.25",
    "me_2.5_effect_1",
    "me_1_effect_0.25",
    "me_1_effect_1"
  )
  reference_groups <- stats::setNames(
    lapply(
      scenarios,
      function(scenario) {
        selected <- sparse_groups[
          vapply(
            sparse_groups,
            function(group) identical(group$scenario, scenario),
            logical(1)
          )
        ]

        if (length(selected) != 1L) {
          stop(
            "Expected one sparse group for scenario ",
            scenario,
            ".",
            call. = FALSE
          )
        }

        selected[[1L]]
      }
    ),
    scenarios
  )

  result_C_sparse_list <- make_named_scenario_list(
    sparse_groups,
    "biases",
    "sparse_mr_rr"
  )
  C_sparse_list <- make_named_scenario_list(
    sparse_groups,
    "estimates",
    "sparse_mr_rr"
  )
  Y_pred_C_sparse_list <- make_named_scenario_list(
    sparse_groups,
    "predictions",
    "sparse_mr_rr"
  )
  B_sparse_list <- stats::setNames(
    lapply(reference_groups, function(group) group$sparse_B),
    scenarios
  )
  sparse_success_rate_list <- stats::setNames(
    lapply(
      C_sparse_list,
      function(result) {
        mean(colSums(is.finite(result)) == nrow(result))
      }
    ),
    scenarios
  )

  list(
    result_C_sparse_list = result_C_sparse_list,
    C_sparse_list = C_sparse_list,
    Y_pred_C_sparse_list = Y_pred_C_sparse_list,
    B_sparse_list = B_sparse_list,
    sparse_success_rate_list = sparse_success_rate_list,
    iv_strength_list = unname(vapply(
      reference_groups,
      function(group) group$parameters$iv_strength,
      numeric(1)
    )),
    parameters_list = unname(lapply(
      reference_groups,
      function(group) group$parameters
    ))
  )
}


make_full_run_metadata <- function(
    design,
    manifest_file,
    inventory_file,
    total_replicates,
    groups) {
  list(
    schema_version = "1.0.0",
    generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    design = design,
    total_replicates_per_setting = total_replicates,
    manifest_file = basename(manifest_file),
    manifest_md5 = unname(tools::md5sum(manifest_file)),
    chunk_inventory_file = basename(inventory_file),
    archived_result_rdata_read = FALSE,
    rng_note = paste(
      "Deterministic per-replicate L'Ecuyer-CMRG streams.",
      "main_no_mrdag and main_mrdag share simulation seeds."
    ),
    input_manifest = groups[[1L]]$input_manifest,
    session_info = capture.output(utils::sessionInfo())
  )
}


run_standalone_chunk_check <- function(chunk_argument, repo_root) {
  chunk_file <- resolve_repo_path(
    chunk_argument,
    repo_root,
    must_work = TRUE
  )
  chunk <- readRDS(chunk_file)
  validate_chunk(chunk, chunk_file)

  cat("MR-rr standalone simulation chunk check\n")
  cat("Chunk:", chunk_file, "\n")
  cat("Design:", chunk$metadata$design, "\n")
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
  cat("Standalone chunk validation: PASS\n")
  invisible(chunk)
}


merge_full_simulation <- function(
    arguments,
    repo_root,
    allow_incomplete,
    overwrite) {
  manifest_argument <- arguments[["manifest"]]
  chunks_argument <- arguments[["chunks-dir"]]
  output_argument <- arguments[["output-dir"]]

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
  chunk_directory <- if (is.null(chunks_argument) ||
      !nzchar(trimws(chunks_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "full_run",
      "chunks"
    )
  } else {
    resolve_repo_path(chunks_argument, repo_root)
  }
  output_directory <- if (is.null(output_argument) ||
      !nzchar(trimws(output_argument))) {
    file.path(
      repo_root,
      "paper",
      "output",
      "full_run",
      "merged"
    )
  } else {
    resolve_repo_path(output_argument, repo_root)
  }

  if (!file.exists(manifest_file)) {
    stop("Task manifest is missing: ", manifest_file, call. = FALSE)
  }

  if (!dir.exists(chunk_directory)) {
    stop("Chunk directory is missing: ", chunk_directory, call. = FALSE)
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  manifest_file <- normalizePath(
    manifest_file,
    winslash = "/",
    mustWork = TRUE
  )
  chunk_directory <- normalizePath(
    chunk_directory,
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
  manifest_validation <- validate_manifest(manifest)
  total_replicates <- manifest_validation$total_replicates
  chunk_paths <- file.path(
    chunk_directory,
    manifest$expected_output_file
  )
  present <- file.exists(chunk_paths)
  missing_tasks <- manifest[!present, , drop = FALSE]
  missing_inventory_file <- file.path(
    output_directory,
    "missing_simulation_tasks.csv"
  )

  write_csv_atomically(missing_tasks, missing_inventory_file)

  cat("MR-rr full simulation chunk merger\n")
  cat("Manifest:", manifest_file, "\n")
  cat("Chunk directory:", chunk_directory, "\n")
  cat("Expected chunks:", nrow(manifest), "\n")
  cat("Present chunks:", sum(present), "\n")
  cat("Missing chunks:", sum(!present), "\n")

  inventory <- manifest
  inventory$chunk_path <- chunk_paths
  inventory$present <- present
  inventory$validated <- FALSE
  inventory$file_size_bytes <- NA_real_
  inventory$md5 <- NA_character_
  inventory$schema_version <- NA_character_
  inventory$complete <- NA
  inventory$archived_result_rdata_read <- NA

  for (row_index in which(present)) {
    chunk <- readRDS(chunk_paths[[row_index]])
    validate_chunk(
      chunk,
      chunk_paths[[row_index]],
      expected_row = manifest[row_index, , drop = FALSE]
    )
    file_info <- file.info(chunk_paths[[row_index]])
    inventory$validated[[row_index]] <- TRUE
    inventory$file_size_bytes[[row_index]] <- file_info$size
    inventory$md5[[row_index]] <- unname(
      tools::md5sum(chunk_paths[[row_index]])
    )
    inventory$schema_version[[row_index]] <-
      chunk$metadata$schema_version
    inventory$complete[[row_index]] <- chunk$metadata$complete
    inventory$archived_result_rdata_read[[row_index]] <-
      chunk$metadata$archived_result_rdata_read
  }

  inventory_file <- file.path(
    output_directory,
    "simulation_chunk_inventory.csv"
  )
  write_csv_atomically(inventory, inventory_file)

  if (nrow(missing_tasks) > 0L) {
    if (allow_incomplete) {
      cat("Existing expected chunks: PASS\n")
      cat("Missing-task inventory:", missing_inventory_file, "\n")
      cat("Full merge deferred until every manifest task is present.\n")
      cat("Incomplete full-run inventory: PASS\n")
      return(invisible(inventory))
    }

    stop(
      "Full simulation output is incomplete. See: ",
      missing_inventory_file,
      "\nUse --allow-incomplete=true only for a progress check.",
      call. = FALSE
    )
  }

  if (!all(inventory$validated)) {
    stop("At least one expected chunk was not validated.", call. = FALSE)
  }

  group_keys <- manifest_validation$group_keys
  merged_groups <- vector("list", nrow(group_keys))

  for (group_index in seq_len(nrow(group_keys))) {
    key <- group_keys[group_index, , drop = FALSE]
    manifest_group <- manifest[
      manifest$design == key$design &
      manifest$phase == key$phase &
      manifest$setting == key$setting,
      ,
      drop = FALSE
    ]

    cat(
      "Merging:",
      paste(key$design, key$phase, key$setting, sep = "/"),
      "\n"
    )
    merged_groups[[group_index]] <- merge_chunk_group(
      manifest_group = manifest_group,
      chunk_directory = chunk_directory,
      total_replicates = total_replicates
    )
  }

  names(merged_groups) <- vapply(
    merged_groups,
    function(group) {
      paste(
        group$design,
        group$phase,
        group$setting_index,
        sep = "__"
      )
    },
    character(1)
  )

  for (design in c("generic_low_rank", "sparse_loading")) {
    for (setting_index in seq_len(4L)) {
      main <- merged_groups[[paste(
        design,
        "main_no_mrdag",
        setting_index,
        sep = "__"
      )]]
      mrdag <- merged_groups[[paste(
        design,
        "main_mrdag",
        setting_index,
        sep = "__"
      )]]

      if (!identical(main$replicate_ids, mrdag$replicate_ids) ||
          !identical(main$replicate_seeds, mrdag$replicate_seeds)) {
        stop(
          "Main and MrDAG seed streams differ for ",
          design,
          ", setting ",
          setting_index,
          ".",
          call. = FALSE
        )
      }

      same_object(
        main$parameters,
        mrdag$parameters,
        paste(design, "setting", setting_index, "parameters")
      )
      same_object(
        main$truth,
        mrdag$truth,
        paste(design, "setting", setting_index, "truth")
      )
    }
  }

  output_specification <- list(
    generic_low_rank = list(
      main_file = "simulate_result_pred_260717_regularC.RData",
      sparse_file = "simulate_result_pred_sparse_260717_regularC.RData"
    ),
    sparse_loading = list(
      main_file = "simulate_result_pred_260718_sparseC.RData",
      sparse_file = "simulate_result_pred_sparse_260718_sparseC.RData"
    )
  )

  output_paths <- unlist(
    lapply(
      output_specification,
      function(specification) {
        file.path(
          output_directory,
          c(specification$main_file, specification$sparse_file)
        )
      }
    ),
    use.names = FALSE
  )

  existing_outputs <- output_paths[file.exists(output_paths)]

  if (length(existing_outputs) > 0L && !overwrite) {
    stop(
      paste(
        "Merged outputs already exist:",
        paste(existing_outputs, collapse = "\n  "),
        "Use --overwrite=true only when deliberately replacing them.",
        sep = "\n  "
      ),
      call. = FALSE
    )
  }

  for (design in names(output_specification)) {
    design_groups <- merged_groups[
      vapply(
        merged_groups,
        function(group) identical(group$design, design),
        logical(1)
      )
    ]
    main_groups <- design_groups[
      vapply(
        design_groups,
        function(group) identical(group$phase, "main_no_mrdag"),
        logical(1)
      )
    ]
    mrdag_groups <- design_groups[
      vapply(
        design_groups,
        function(group) identical(group$phase, "main_mrdag"),
        logical(1)
      )
    ]
    sparse_groups <- design_groups[
      vapply(
        design_groups,
        function(group) identical(group$phase, "sparse"),
        logical(1)
      )
    ]

    simulate_result_prediction <- build_legacy_main_result(
      main_groups,
      mrdag_groups
    )
    sparse_results <- build_legacy_sparse_result(sparse_groups)
    full_run_metadata <- make_full_run_metadata(
      design = design,
      manifest_file = manifest_file,
      inventory_file = inventory_file,
      total_replicates = total_replicates,
      groups = design_groups
    )

    specification <- output_specification[[design]]
    main_file <- file.path(
      output_directory,
      specification$main_file
    )
    sparse_file <- file.path(
      output_directory,
      specification$sparse_file
    )

    save_objects_atomically(
      list(
        simulate_result_prediction = simulate_result_prediction,
        full_run_metadata = full_run_metadata
      ),
      main_file,
      overwrite = overwrite
    )
    save_objects_atomically(
      list(
        sparse_results = sparse_results,
        full_run_metadata = full_run_metadata
      ),
      sparse_file,
      overwrite = overwrite
    )
  }

  cat("\nAll expected chunks validated: PASS\n")
  cat("Replicate coverage 1 through", total_replicates, ": PASS\n")
  cat("Main/MrDAG shared seed streams: PASS\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Merged output directory:", output_directory, "\n")
  cat("Full simulation chunk merge: PASS\n")

  invisible(
    list(
      manifest = manifest,
      inventory = inventory,
      groups = merged_groups,
      output_paths = output_paths
    )
  )
}


run_merge_script <- function() {
  arguments <- parse_named_arguments(
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
  validate_argument_names(arguments, allowed_arguments)

  allow_incomplete <- parse_boolean_argument(
    arguments,
    "allow-incomplete",
    default = FALSE
  )
  overwrite <- parse_boolean_argument(
    arguments,
    "overwrite",
    default = FALSE
  )
  repo_root <- locate_repo_root()
  chunk_argument <- arguments[["check-chunk"]]

  if (!is.null(chunk_argument) && nzchar(trimws(chunk_argument))) {
    incompatible_arguments <- intersect(
      names(arguments),
      c(
        "manifest",
        "chunks-dir",
        "output-dir",
        "allow-incomplete",
        "overwrite"
      )
    )

    if (length(incompatible_arguments) > 0L) {
      stop(
        "--check-chunk cannot be combined with full-merge arguments.",
        call. = FALSE
      )
    }

    return(run_standalone_chunk_check(chunk_argument, repo_root))
  }

  merge_full_simulation(
    arguments = arguments,
    repo_root = repo_root,
    allow_incomplete = allow_incomplete,
    overwrite = overwrite
  )
}


if (sys.nframe() == 0L) {
  run_merge_script()
}
