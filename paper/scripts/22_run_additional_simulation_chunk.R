#!/usr/bin/env Rscript

# Portable point-estimation chunk worker for the additional MR-rr simulations.
#
# Supported tasks:
#   rank_misspecification / standard
#   approximate_low_rank / standard
#   approximate_low_rank / mrdag
#
# The worker reads frozen raw inputs and the frozen estimator implementation,
# but never reads archived simulation-result RData and never writes to freeze/.

options(warn = 1)


additional_chunk_parse_arguments <- function(arguments) {
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


additional_chunk_validate_argument_names <- function(
    arguments,
    allowed_names) {
  unknown <- setdiff(names(arguments), allowed_names)

  if (length(unknown) > 0L) {
    stop(
      "Unknown arguments: ",
      paste(paste0("--", unknown), collapse = ", "),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


additional_chunk_require_argument <- function(arguments, name) {
  value <- arguments[[name]]

  if (is.null(value) || !nzchar(trimws(value))) {
    stop("Missing required argument --", name, ".", call. = FALSE)
  }

  trimws(value)
}


additional_chunk_parse_integer <- function(
    arguments,
    name,
    default = NULL,
    minimum = 1L) {
  value <- arguments[[name]]

  if (is.null(value)) {
    if (is.null(default)) {
      stop("Missing required argument --", name, ".", call. = FALSE)
    }

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


additional_chunk_parse_boolean <- function(
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


additional_chunk_resolve_path <- function(path, repo_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(
    file.path(repo_root, path),
    winslash = "/",
    mustWork = FALSE
  )
}


additional_chunk_safe_errors <- function() {
  data.frame(
    replicate_id = integer(),
    result_key = character(),
    message = character(),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


additional_chunk_append_error <- function(
    errors,
    replicate_id,
    result_key,
    message) {
  rbind(
    errors,
    data.frame(
      replicate_id = as.integer(replicate_id),
      result_key = as.character(result_key),
      message = as.character(message),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  )
}


additional_chunk_save_rds <- function(object, path) {
  output_directory <- dirname(path)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = output_directory,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)
  saveRDS(object, temporary_file, compress = "xz")

  backup_file <- NULL

  if (file.exists(path)) {
    backup_file <- tempfile(
      pattern = paste0(basename(path), "."),
      tmpdir = output_directory,
      fileext = ".bak"
    )

    if (!file.rename(path, backup_file)) {
      stop("Could not preserve existing output: ", path, call. = FALSE)
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop("Could not install completed output: ", path, call. = FALSE)
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


additional_chunk_make_result_seed <- function(
    data_seed,
    method,
    working_rank) {
  method_offsets <- c(
    ivw = 10000L,
    srivw = 20000L,
    naive_mr_rr = 30000L,
    mr_rr = 40000L,
    regularized_mr_rr = 50000L,
    sparse_mr_rr = 60000L,
    mrdag = 70000L
  )
  method_offset <- method_offsets[[method]]

  if (is.null(method_offset)) {
    stop("Unknown method: ", method, call. = FALSE)
  }

  seed <-
    as.double(data_seed) +
    method_offset +
    as.double(working_rank) * 1000

  if (!is.finite(seed) || seed > .Machine$integer.max) {
    stop("Derived estimator seed exceeds the R integer range.", call. = FALSE)
  }

  as.integer(seed)
}


additional_chunk_set_seed <- function(seed) {
  set.seed(
    seed,
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}


additional_chunk_fit_method <- function(
    method,
    simulated_data,
    parameters,
    regularization_rate,
    sparse_lambda,
    estimators,
    mrdag_niter,
    mrdag_burnin) {
  switch(
    method,
    naive_mr_rr = estimators$mr_rr_naive(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      r = parameters$r_RR,
      W = parameters$weight.matrix
    )$AB,
    mr_rr = estimators$mr_rr(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      W = parameters$weight.matrix
    )$AB,
    regularized_mr_rr = estimators$mr_rr_regularized(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      regularization_rate = regularization_rate,
      W = parameters$weight.matrix
    )$AB,
    ivw = estimators$ivw_multiple_outcomes(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    ),
    srivw = estimators$adivw_multiple_outcomes(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    ),
    sparse_mr_rr = estimators$mr_rr_sparse(
      GAMMA_hat = simulated_data$GAMMA_hat,
      gamma_hat = simulated_data$gamma_hat,
      W = parameters$weight.matrix,
      Sigma_X = parameters$Sigma_X,
      lambda = rep(sparse_lambda, parameters$px),
      r = parameters$r_RR,
      max_iter = 100L,
      tol = 1e-2
    ),
    mrdag = estimators$Mr_DAG(
      Y = simulated_data$GAMMA_hat,
      X = simulated_data$gamma_hat,
      niter = mrdag_niter,
      burnin = mrdag_burnin
    ),
    stop("Unknown method: ", method, call. = FALSE)
  )
}


additional_chunk_validate_fit <- function(
    fit,
    method,
    parameters,
    result_key,
    baseline) {
  if (method == "sparse_mr_rr") {
    effect_estimate <- fit$AB
    baseline$assert_dimensions(
      fit$A,
      c(parameters$py, parameters$r_RR),
      paste(result_key, "sparse A")
    )
    baseline$assert_dimensions(
      fit$B,
      c(parameters$r_RR, parameters$px),
      paste(result_key, "sparse B")
    )
    baseline$assert_finite(fit$A, paste(result_key, "sparse A"))
    baseline$assert_finite(fit$B, paste(result_key, "sparse B"))
    reconstruction_difference <- max(abs(
      effect_estimate - fit$A %*% fit$B
    ))

    if (!is.finite(reconstruction_difference) ||
        reconstruction_difference > 1e-8) {
      stop(
        result_key,
        " failed AB = A %*% B; maximum difference = ",
        format(reconstruction_difference, scientific = TRUE),
        ".",
        call. = FALSE
      )
    }
  } else {
    effect_estimate <- fit
  }

  baseline$validate_effect_estimate(
    estimate = effect_estimate,
    parameters = parameters,
    method = result_key,
    require_low_rank = method %in% c(
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr",
      "sparse_mr_rr"
    )
  )

  effect_estimate
}


additional_chunk_result_specification <- function(
    analysis,
    phase,
    selected_config) {
  if (analysis == "rank_misspecification") {
    rows <- selected_config[order(selected_config$working_rank), , drop = FALSE]
    specification <- do.call(
      rbind,
      lapply(seq_len(nrow(rows)), function(row_index) {
        data.frame(
          result_key = paste0(
            "r",
            rows$working_rank[[row_index]],
            "__",
            c("regularized_mr_rr", "sparse_mr_rr")
          ),
          config_row_index = row_index,
          working_rank = rows$working_rank[[row_index]],
          method = c("regularized_mr_rr", "sparse_mr_rr"),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      }))
    rownames(specification) <- NULL
    return(list(config = rows, specification = specification))
  }

  methods <- if (phase == "standard") {
    setdiff(
      strsplit(
        selected_config$method_set[[1L]],
        split = ";",
        fixed = TRUE
      )[[1L]],
      "mrdag"
    )
  } else {
    "mrdag"
  }

  specification <- data.frame(
    result_key = methods,
    config_row_index = rep(1L, length(methods)),
    working_rank = rep(selected_config$working_rank[[1L]], length(methods)),
    method = methods,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  list(config = selected_config, specification = specification)
}


additional_run_point_chunk <- function() {
  arguments <- additional_chunk_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "analysis",
    "phase",
    "setting",
    "start",
    "end",
    "seed-base",
    "mrdag-niter",
    "mrdag-burnin",
    "output-dir",
    "allow-provisional",
    "overwrite"
  )
  additional_chunk_validate_argument_names(arguments, allowed_arguments)
  analysis <- additional_chunk_require_argument(arguments, "analysis")
  phase <- additional_chunk_require_argument(arguments, "phase")
  setting_index <- additional_chunk_parse_integer(
    arguments,
    "setting",
    minimum = 1L
  )
  replicate_start <- additional_chunk_parse_integer(
    arguments,
    "start",
    minimum = 1L
  )
  replicate_end <- additional_chunk_parse_integer(
    arguments,
    "end",
    minimum = 1L
  )
  seed_base <- additional_chunk_parse_integer(
    arguments,
    "seed-base",
    default = 123L,
    minimum = 1L
  )
  mrdag_niter <- additional_chunk_parse_integer(
    arguments,
    "mrdag-niter",
    default = 1000L,
    minimum = 1L
  )
  mrdag_burnin <- additional_chunk_parse_integer(
    arguments,
    "mrdag-burnin",
    default = 200L,
    minimum = 0L
  )
  allow_provisional <- additional_chunk_parse_boolean(
    arguments,
    "allow-provisional",
    default = FALSE
  )
  overwrite <- additional_chunk_parse_boolean(
    arguments,
    "overwrite",
    default = FALSE
  )

  if (!analysis %in% c(
    "rank_misspecification",
    "approximate_low_rank"
  )) {
    stop(
      "--analysis must be rank_misspecification or approximate_low_rank.",
      call. = FALSE
    )
  }

  if (!phase %in% c("standard", "mrdag")) {
    stop("--phase must be standard or mrdag.", call. = FALSE)
  }

  if (analysis == "rank_misspecification" && phase != "standard") {
    stop(
      "Rank misspecification supports only --phase=standard.",
      call. = FALSE
    )
  }

  if (!setting_index %in% seq_len(4L)) {
    stop("--setting must be an integer from 1 through 4.", call. = FALSE)
  }

  if (replicate_end < replicate_start) {
    stop("--end must be greater than or equal to --start.", call. = FALSE)
  }

  if (mrdag_burnin >= mrdag_niter) {
    stop("--mrdag-burnin must be smaller than --mrdag-niter.", call. = FALSE)
  }

  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the worker script path.", call. = FALSE)
  }

  worker_file <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  core_file <- file.path(
    dirname(worker_file),
    "20_additional_simulation_core.R"
  )

  if (!file.exists(core_file)) {
    stop("Additional simulation core is missing: ", core_file, call. = FALSE)
  }

  core <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  repo_root <- core$additional_find_repo_root(dirname(worker_file))
  config_file <- file.path(
    repo_root,
    "paper",
    "config",
    paste0(analysis, "_settings.csv")
  )
  config <- core$additional_read_config(config_file, analysis)
  selected_config <- config[
    config$setting_index == setting_index,
    ,
    drop = FALSE
  ]

  if (nrow(selected_config) == 0L) {
    stop("Selected setting is absent from the configuration.", call. = FALSE)
  }

  if (replicate_end > unique(selected_config$monte_carlo_replicates)) {
    stop(
      "Requested replicate range exceeds the configured total.",
      call. = FALSE
    )
  }

  if (seed_base != unique(selected_config$seed_base)) {
    stop("--seed-base differs from the locked configuration.", call. = FALSE)
  }

  if (analysis == "approximate_low_rank" &&
      any(selected_config$status != "locked")) {
    if (!allow_provisional) {
      stop(
        paste(
          "Approximate-low-rank delta is still provisional.",
          "Use --allow-provisional=true only for a small validation task;",
          "lock the configuration before the full run."
        ),
        call. = FALSE
      )
    }

    if (replicate_end - replicate_start + 1L > 5L) {
      stop(
        "A provisional approximate-low-rank task is limited to 5 replicates.",
        call. = FALSE
      )
    }
  }

  specification_object <- additional_chunk_result_specification(
    analysis = analysis,
    phase = phase,
    selected_config = selected_config
  )
  selected_config <- specification_object$config
  result_specification <- specification_object$specification

  if (nrow(result_specification) == 0L ||
      anyDuplicated(result_specification$result_key)) {
    stop("Point-estimation result specification is invalid.", call. = FALSE)
  }

  baseline <- core$additional_load_baseline_core(repo_root)

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package `MASS` is required.", call. = FALSE)
  }

  estimators <- core$additional_load_estimators(repo_root)
  calibration <- core$additional_build_calibration(
    repo_root = repo_root,
    baseline = baseline,
    pz = 177L
  )
  delta <- unique(selected_config$third_singular_value)
  truth_object <- core$additional_make_generic_truth(
    third_singular_value = delta,
    px = calibration$px,
    py = calibration$py,
    seed = 123L
  )
  causal_effect <- truth_object$C
  replicate_ids <- seq.int(replicate_start, replicate_end)
  replicate_count <- length(replicate_ids)
  effect_dimension <- calibration$py * calibration$px
  result_keys <- result_specification$result_key
  estimates <- stats::setNames(
    lapply(result_keys, function(result_key) {
      matrix(
        NA_real_,
        nrow = effect_dimension,
        ncol = replicate_count,
        dimnames = list(NULL, as.character(replicate_ids))
      )
    }),
    result_keys
  )
  biases <- stats::setNames(
    lapply(result_keys, function(result_key) {
      matrix(
        NA_real_,
        nrow = effect_dimension,
        ncol = replicate_count,
        dimnames = list(NULL, as.character(replicate_ids))
      )
    }),
    result_keys
  )
  data_seeds <- integer(replicate_count)
  estimator_seeds <- stats::setNames(
    lapply(result_keys, function(result_key) integer(replicate_count)),
    result_keys
  )
  errors <- additional_chunk_safe_errors()
  output_argument <- arguments[["output-dir"]]

  if (is.null(output_argument) || !nzchar(trimws(output_argument))) {
    output_directory <- file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "point_chunks"
    )
  } else {
    output_directory <- additional_chunk_resolve_path(
      trimws(output_argument),
      repo_root
    )
  }

  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )
  output_file <- file.path(
    output_directory,
    sprintf(
      "%s__%s__setting-%d__rep-%04d-%04d.rds",
      analysis,
      phase,
      setting_index,
      replicate_start,
      replicate_end
    )
  )

  if (file.exists(output_file) && !overwrite) {
    stop(
      "Output already exists; use --overwrite=true to replace it: ",
      output_file,
      call. = FALSE
    )
  }

  cat("MR-rr additional-simulation point chunk worker\n")
  cat("Repository root:", repo_root, "\n")
  cat("Analysis:", analysis, "\n")
  cat("Phase:", phase, "\n")
  cat("Setting:", setting_index, "\n")
  cat("Replicates:", replicate_start, "through", replicate_end, "\n")
  cat("Third singular value:", delta, "\n")
  cat("Configuration status:", paste(unique(selected_config$status), collapse = ", "), "\n")
  cat("Result keys:", paste(result_keys, collapse = ", "), "\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Output:", output_file, "\n\n")
  start_time <- proc.time()[["elapsed"]]

  for (replicate_position in seq_along(replicate_ids)) {
    replicate_id <- replicate_ids[[replicate_position]]
    canonical_config <- selected_config[1L, , drop = FALSE]
    simulated <- core$additional_simulate_dataset(
      config_row = canonical_config,
      causal_effect = causal_effect,
      replicate_id = replicate_id,
      calibration = calibration,
      estimators = estimators,
      baseline = baseline
    )
    data_seeds[[replicate_position]] <- simulated$replicate_seed
    parameters_by_row <- lapply(
      seq_len(nrow(selected_config)),
      function(row_index) {
        core$additional_build_parameters(
          config_row = selected_config[row_index, , drop = FALSE],
          causal_effect = causal_effect,
          calibration = calibration,
          estimators = estimators,
          baseline = baseline
        )
      }
    )

    for (specification_index in seq_len(nrow(result_specification))) {
      specification <- result_specification[
        specification_index,
        ,
        drop = FALSE
      ]
      result_key <- specification$result_key[[1L]]
      method <- specification$method[[1L]]
      config_row_index <- specification$config_row_index[[1L]]
      config_row <- selected_config[config_row_index, , drop = FALSE]
      parameters <- parameters_by_row[[config_row_index]]
      estimator_seed <- additional_chunk_make_result_seed(
        data_seed = simulated$replicate_seed,
        method = method,
        working_rank = parameters$r_RR
      )
      estimator_seeds[[result_key]][[replicate_position]] <- estimator_seed
      additional_chunk_set_seed(estimator_seed)
      fit <- tryCatch(
        additional_chunk_fit_method(
          method = method,
          simulated_data = simulated$data,
          parameters = parameters,
          regularization_rate = config_row$regularization_rate,
          sparse_lambda = config_row$sparse_lambda,
          estimators = estimators,
          mrdag_niter = mrdag_niter,
          mrdag_burnin = mrdag_burnin
        ),
        error = function(error) error
      )

      if (inherits(fit, "error")) {
        errors <- additional_chunk_append_error(
          errors,
          replicate_id,
          result_key,
          conditionMessage(fit)
        )
        next
      }

      effect_estimate <- tryCatch(
        additional_chunk_validate_fit(
          fit = fit,
          method = method,
          parameters = parameters,
          result_key = paste0(
            analysis,
            "/",
            phase,
            "/setting=",
            setting_index,
            "/replicate=",
            replicate_id,
            "/",
            result_key
          ),
          baseline = baseline
        ),
        error = function(error) error
      )

      if (inherits(effect_estimate, "error")) {
        errors <- additional_chunk_append_error(
          errors,
          replicate_id,
          result_key,
          conditionMessage(effect_estimate)
        )
        next
      }

      estimates[[result_key]][, replicate_position] <-
        as.vector(effect_estimate)
      biases[[result_key]][, replicate_position] <-
        as.vector(effect_estimate - causal_effect)
    }

    cat(sprintf(
      "Completed replicate %d (%d/%d)\n",
      replicate_id,
      replicate_position,
      replicate_count
    ))
    flush.console()
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  successful_replicates <- vapply(
    estimates,
    function(value) sum(colSums(is.finite(value)) == nrow(value)),
    integer(1)
  )
  freeze_project <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project"
  )
  required_files <- c(
    file.path(freeze_project, "data", "dat_1e-4.csv"),
    file.path(freeze_project, "data", "rho_mat_1e-4.csv"),
    file.path(freeze_project, "scripts", "MR_rr_estimators.R"),
    core_file,
    config_file,
    worker_file
  )
  missing_inputs <- required_files[!file.exists(required_files)]

  if (length(missing_inputs) > 0L) {
    stop(
      "Required inputs disappeared during the task: ",
      paste(missing_inputs, collapse = ", "),
      call. = FALSE
    )
  }

  input_manifest <- data.frame(
    repository_path = substring(
      normalizePath(required_files, winslash = "/", mustWork = TRUE),
      nchar(normalizePath(repo_root, winslash = "/", mustWork = TRUE)) + 2L
    ),
    md5 = unname(tools::md5sum(required_files)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  point_chunk <- list(
    metadata = list(
      schema_version = "1.0.0",
      generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      analysis = analysis,
      phase = phase,
      setting_index = setting_index,
      scenario = unique(selected_config$scenario),
      replicate_start = replicate_start,
      replicate_end = replicate_end,
      replicate_ids = replicate_ids,
      replicate_count = replicate_count,
      seed_base = seed_base,
      data_seeds = data_seeds,
      estimator_seeds = estimator_seeds,
      rng_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      mrdag_niter = if (phase == "mrdag") mrdag_niter else NA_integer_,
      mrdag_burnin = if (phase == "mrdag") {
        mrdag_burnin
      } else {
        NA_integer_
      },
      allow_provisional = allow_provisional,
      elapsed_seconds = elapsed_seconds,
      complete = nrow(errors) == 0L,
      archived_simulation_result_rdata_read = FALSE,
      frozen_snapshot_written = FALSE
    ),
    input_manifest = input_manifest,
    configuration = selected_config,
    result_specification = result_specification,
    truth = list(
      C = causal_effect,
      singular_values = truth_object$singular_values,
      numerical_rank = truth_object$numerical_rank
    ),
    estimates = estimates,
    biases = biases,
    successful_replicates = successful_replicates,
    errors = errors,
    session_info = capture.output(utils::sessionInfo())
  )
  additional_chunk_save_rds(point_chunk, output_file)

  cat("\nChunk output written:", output_file, "\n")
  cat(
    "Successful replicates by result key:",
    paste(
      paste0(names(successful_replicates), "=", successful_replicates),
      collapse = ", "
    ),
    "\n"
  )
  cat("Elapsed time:", sprintf("%.1f seconds", elapsed_seconds), "\n")

  if (nrow(errors) > 0L) {
    print(errors, row.names = FALSE)
    stop(
      "Point-estimation chunk completed with ",
      nrow(errors),
      " failed fits.",
      call. = FALSE
    )
  }

  cat("Archived simulation-result RData independence: PASS\n")
  cat("Frozen snapshot written: no\n")
  cat("Additional-simulation point chunk: PASS\n")
  invisible(point_chunk)
}


if (sys.nframe() == 0L) {
  additional_run_point_chunk()
}
