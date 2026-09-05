#!/usr/bin/env Rscript

# Portable chunk worker for the full MR-rr simulation rerun.
#
# This worker reads only the frozen raw calibration inputs and the frozen
# estimator implementation. It does not read archived simulation-result RData.
# The same command can be called directly on a workstation or from a Slurm
# array task.

options(warn = 1)


locate_repository_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION")) &&
        dir.exists(file.path(current, "paper", "scripts")) &&
        dir.exists(file.path(
          current,
          "freeze",
          "current_analysis_20260823",
          "project"
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
      paste(
        "Every argument must use --name=value syntax.",
        "Example: --design=generic_low_rank --phase=main_no_mrdag"
      ),
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


require_argument <- function(arguments, name) {
  value <- arguments[[name]]

  if (is.null(value) || !nzchar(trimws(value))) {
    stop(
      "Missing required argument --",
      name,
      ".",
      call. = FALSE
    )
  }

  trimws(value)
}


parse_integer_argument <- function(
    arguments,
    name,
    default = NULL,
    minimum = 1L) {
  value <- arguments[[name]]

  if (is.null(value)) {
    if (is.null(default)) {
      stop(
        "Missing required argument --",
        name,
        ".",
        call. = FALSE
      )
    }

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
      "The derived replicate seed exceeds the R integer range.",
      call. = FALSE
    )
  }

  as.integer(seed)
}


safe_error_table <- function() {
  data.frame(
    replicate_id = integer(),
    method = character(),
    message = character(),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


append_error <- function(error_table, replicate_id, method, message) {
  rbind(
    error_table,
    data.frame(
      replicate_id = as.integer(replicate_id),
      method = as.character(method),
      message = as.character(message),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  )
}


fit_one_method <- function(
    method,
    simulated_data,
    parameters,
    regularization_rate,
    sparse_lambda,
    estimators,
    mrdag_niter,
    mrdag_burnin) {
  GAMMA_hat <- simulated_data$GAMMA_hat
  gamma_hat <- simulated_data$gamma_hat

  switch(
    method,
    naive_mr_rr = estimators$mr_rr_naive(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      W = parameters$weight.matrix
    )$AB,
    mr_rr = estimators$mr_rr(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      W = parameters$weight.matrix
    )$AB,
    regularized_mr_rr = estimators$mr_rr_regularized(
      Y = GAMMA_hat,
      X = gamma_hat,
      r = parameters$r_RR,
      Sigma_X = parameters$Sigma_X,
      regularization_rate = regularization_rate,
      W = parameters$weight.matrix
    )$AB,
    ivw = estimators$ivw_multiple_outcomes(
      Y = GAMMA_hat,
      X = gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    ),
    srivw = estimators$adivw_multiple_outcomes(
      Y = GAMMA_hat,
      X = gamma_hat,
      Sigma_X = parameters$Sigma_X,
      Sigma_Y = parameters$Sigma_Y
    ),
    mrdag = estimators$Mr_DAG(
      Y = GAMMA_hat,
      X = gamma_hat,
      niter = mrdag_niter,
      burnin = mrdag_burnin
    ),
    sparse_mr_rr = estimators$mr_rr_sparse(
      GAMMA_hat = GAMMA_hat,
      gamma_hat = gamma_hat,
      W = parameters$weight.matrix,
      Sigma_X = parameters$Sigma_X,
      lambda = rep(sparse_lambda, parameters$px),
      r = parameters$r_RR,
      max_iter = 100L,
      tol = 1e-2
    ),
    stop(
      "Unknown estimator method: ",
      method,
      call. = FALSE
    )
  )
}


save_rds_atomically <- function(object, path) {
  output_directory <- dirname(path)
  temporary_file <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = output_directory,
    fileext = ".tmp"
  )
  on.exit(unlink(temporary_file), add = TRUE)

  saveRDS(
    object,
    file = temporary_file,
    compress = "xz"
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
        "Could not preserve the existing chunk before replacement: ",
        path,
        call. = FALSE
      )
    }
  }

  if (!file.rename(temporary_file, path)) {
    if (!is.null(backup_file) && file.exists(backup_file)) {
      file.rename(backup_file, path)
    }

    stop(
      "Could not atomically move the completed chunk to: ",
      path,
      call. = FALSE
    )
  }

  if (!is.null(backup_file) && file.exists(backup_file)) {
    unlink(backup_file)
  }

  invisible(path)
}


run_simulation_chunk <- function() {
  arguments <- parse_named_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "design",
    "phase",
    "setting",
    "start",
    "end",
    "seed-base",
    "mrdag-niter",
    "mrdag-burnin",
    "output-dir",
    "overwrite"
  )
  validate_argument_names(arguments, allowed_arguments)

  design <- require_argument(arguments, "design")
  phase <- require_argument(arguments, "phase")
  setting_index <- parse_integer_argument(
    arguments,
    "setting",
    minimum = 1L
  )
  replicate_start <- parse_integer_argument(
    arguments,
    "start",
    minimum = 1L
  )
  replicate_end <- parse_integer_argument(
    arguments,
    "end",
    minimum = 1L
  )
  seed_base <- parse_integer_argument(
    arguments,
    "seed-base",
    default = 123L,
    minimum = 1L
  )
  mrdag_niter <- parse_integer_argument(
    arguments,
    "mrdag-niter",
    default = 1000L,
    minimum = 1L
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

  allowed_designs <- c(
    "generic_low_rank",
    "sparse_loading"
  )
  allowed_phases <- c(
    "main_no_mrdag",
    "main_mrdag",
    "sparse"
  )

  if (!design %in% allowed_designs) {
    stop(
      "--design must be generic_low_rank or sparse_loading.",
      call. = FALSE
    )
  }

  if (!phase %in% allowed_phases) {
    stop(
      paste(
        "--phase must be main_no_mrdag, main_mrdag,",
        "or sparse."
      ),
      call. = FALSE
    )
  }

  if (!setting_index %in% seq_len(4L)) {
    stop(
      "--setting must be an integer from 1 through 4.",
      call. = FALSE
    )
  }

  if (replicate_end < replicate_start) {
    stop(
      "--end must be greater than or equal to --start.",
      call. = FALSE
    )
  }

  if (phase == "main_mrdag" && mrdag_burnin >= mrdag_niter) {
    stop(
      "--mrdag-burnin must be smaller than --mrdag-niter.",
      call. = FALSE
    )
  }

  repo_root <- locate_repository_root()
  smoke_core_file <- file.path(
    repo_root,
    "paper",
    "scripts",
    "15_smoke_test_from_scratch_simulations.R"
  )

  if (!file.exists(smoke_core_file)) {
    stop(
      "Required simulation core is missing: ",
      smoke_core_file,
      call. = FALSE
    )
  }

  source(smoke_core_file, local = FALSE)

  freeze_root <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project"
  )
  estimator_file <- file.path(
    freeze_root,
    "scripts",
    "MR_rr_estimators.R"
  )
  data_file <- file.path(
    freeze_root,
    "data",
    "dat_1e-4.csv"
  )
  correlation_file <- file.path(
    freeze_root,
    "data",
    "rho_mat_1e-4.csv"
  )
  required_files <- c(
    estimator_file,
    data_file,
    correlation_file,
    smoke_core_file
  )
  invisible(lapply(required_files, require_file))

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop(
      "Package `MASS` is required.",
      call. = FALSE
    )
  }

  output_directory_argument <- arguments[["output-dir"]]

  if (is.null(output_directory_argument) ||
      !nzchar(trimws(output_directory_argument))) {
    output_directory <- file.path(
      repo_root,
      "paper",
      "output",
      "full_run",
      "chunks"
    )
  } else {
    output_directory <- output_directory_argument
  }

  dir.create(
    output_directory,
    recursive = TRUE,
    showWarnings = FALSE
  )
  output_directory <- normalizePath(
    output_directory,
    winslash = "/",
    mustWork = TRUE
  )

  output_filename <- sprintf(
    "%s__%s__setting-%d__rep-%04d-%04d.rds",
    design,
    phase,
    setting_index,
    replicate_start,
    replicate_end
  )
  output_file <- file.path(
    output_directory,
    output_filename
  )

  if (file.exists(output_file) && !overwrite) {
    stop(
      "Output already exists: ",
      output_file,
      "\nUse --overwrite=true only when deliberately rerunning this chunk.",
      call. = FALSE
    )
  }

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
  selected_setting <- settings[
    setting_index,
    ,
    drop = FALSE
  ]
  scenario <- selected_setting$scenario

  generic_regularization_rates <- c(
    me_2.5_effect_0.25 = 1.007845e-10,
    me_1_effect_0.25 = 1.578970e-12,
    me_2.5_effect_1 = 2.340371e-12,
    me_1_effect_1 = 3.103420e-15
  )
  sparse_regularization_rates <- c(
    me_2.5_effect_0.25 = 9.109067e-11,
    me_1_effect_0.25 = 1.578970e-12,
    me_2.5_effect_1 = 2.340371e-12,
    me_1_effect_1 = 3.105502e-15
  )
  sparse_lambdas <- c(
    me_2.5_effect_0.25 = 1e-4,
    me_1_effect_0.25 = 1e-3,
    me_2.5_effect_1 = 1e-3,
    me_1_effect_1 = 1e-3
  )

  estimators <- new.env(parent = globalenv())
  sys.source(estimator_file, envir = estimators)

  lip_data <- utils::read.csv(
    data_file,
    check.names = FALSE
  )
  lip_correlation <- utils::read.csv(
    correlation_file,
    check.names = FALSE
  )
  calibration <- build_calibration(
    lip_data = lip_data,
    lip_correlation = lip_correlation,
    pz = 177L
  )

  generic_effect <- generate_generic_effect(
    px = calibration$px,
    py = calibration$py,
    rank = 2L,
    seed = 123L
  )
  sparse_truth <- generate_sparse_loading_effect(
    calibration = calibration,
    estimators = estimators,
    rank = 2L,
    seed = 123L
  )
  causal_effect <- if (design == "generic_low_rank") {
    generic_effect
  } else {
    sparse_truth$C
  }
  regularization_rates <- if (design == "generic_low_rank") {
    generic_regularization_rates
  } else {
    sparse_regularization_rates
  }
  regularization_rate <- regularization_rates[[scenario]]
  sparse_lambda <- sparse_lambdas[[scenario]]

  parameters <- build_parameters(
    causal_effect = causal_effect,
    exposure_error_weight = selected_setting$me_weight,
    genetic_effect_weight = selected_setting$effect_weight,
    calibration = calibration,
    estimators = estimators,
    rank = 2L
  )

  set.seed(123L)
  fixed_exposure <- stats::rnorm(
    calibration$px,
    mean = 0,
    sd = 1
  )
  true_prediction <- as.vector(
    causal_effect %*% fixed_exposure
  )

  methods <- switch(
    phase,
    main_no_mrdag = c(
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr",
      "ivw",
      "srivw"
    ),
    main_mrdag = "mrdag",
    sparse = "sparse_mr_rr"
  )
  replicate_ids <- seq.int(
    replicate_start,
    replicate_end
  )
  replicate_count <- length(replicate_ids)
  effect_dimension <- parameters$py * parameters$px

  estimates <- setNames(
    lapply(
      methods,
      function(method) {
        matrix(
          NA_real_,
          nrow = effect_dimension,
          ncol = replicate_count,
          dimnames = list(
            NULL,
            as.character(replicate_ids)
          )
        )
      }
    ),
    methods
  )
  biases <- setNames(
    lapply(
      methods,
      function(method) {
        matrix(
          NA_real_,
          nrow = effect_dimension,
          ncol = replicate_count,
          dimnames = list(
            NULL,
            as.character(replicate_ids)
          )
        )
      }
    ),
    methods
  )
  predictions <- setNames(
    lapply(
      methods,
      function(method) {
        matrix(
          NA_real_,
          nrow = parameters$py,
          ncol = replicate_count,
          dimnames = list(
            NULL,
            as.character(replicate_ids)
          )
        )
      }
    ),
    methods
  )
  sparse_B <- if (phase == "sparse") {
    matrix(
      NA_real_,
      nrow = parameters$r_RR * parameters$px,
      ncol = replicate_count,
      dimnames = list(
        NULL,
        as.character(replicate_ids)
      )
    )
  } else {
    NULL
  }
  replicate_seeds <- integer(replicate_count)
  errors <- safe_error_table()

  design_index <- match(design, allowed_designs)

  cat("MR-rr full simulation chunk worker\n")
  cat("Repository root:", repo_root, "\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Design:", design, "\n")
  cat("Phase:", phase, "\n")
  cat(
    "Setting:",
    setting_index,
    sprintf(
      "(me=%s, effect=%s)",
      selected_setting$me_weight,
      selected_setting$effect_weight
    ),
    "\n"
  )
  cat(
    "Replicates:",
    replicate_start,
    "through",
    replicate_end,
    "\n"
  )
  cat("Methods:", paste(methods, collapse = ", "), "\n")
  cat("Output:", output_file, "\n\n")

  chunk_start_time <- proc.time()[["elapsed"]]

  for (replicate_position in seq_along(replicate_ids)) {
    replicate_id <- replicate_ids[[replicate_position]]
    replicate_seed <- make_replicate_seed(
      seed_base = seed_base,
      design_index = design_index,
      setting_index = setting_index,
      phase = phase,
      replicate_id = replicate_id
    )
    replicate_seeds[[replicate_position]] <- replicate_seed

    set.seed(
      replicate_seed,
      kind = "L'Ecuyer-CMRG",
      normal.kind = "Inversion",
      sample.kind = "Rejection"
    )
    simulated_data <- simulate_dataset(parameters)

    for (method in methods) {
      fit <- tryCatch(
        fit_one_method(
          method = method,
          simulated_data = simulated_data,
          parameters = parameters,
          regularization_rate = regularization_rate,
          sparse_lambda = sparse_lambda,
          estimators = estimators,
          mrdag_niter = mrdag_niter,
          mrdag_burnin = mrdag_burnin
        ),
        error = function(error) error
      )

      if (inherits(fit, "error")) {
        errors <- append_error(
          errors,
          replicate_id = replicate_id,
          method = method,
          message = conditionMessage(fit)
        )
        next
      }

      if (method == "sparse_mr_rr") {
        effect_estimate <- fit$AB
        loading_estimate <- fit$B
      } else {
        effect_estimate <- fit
        loading_estimate <- NULL
      }

      validation_error <- tryCatch(
        {
          validate_effect_estimate(
            estimate = effect_estimate,
            parameters = parameters,
            method = paste(
              design,
              scenario,
              method,
              "replicate",
              replicate_id
            ),
            require_low_rank = method %in% c(
              "naive_mr_rr",
              "mr_rr",
              "regularized_mr_rr",
              "sparse_mr_rr"
            )
          )
          NULL
        },
        error = function(error) error
      )

      if (inherits(validation_error, "error")) {
        errors <- append_error(
          errors,
          replicate_id = replicate_id,
          method = method,
          message = conditionMessage(validation_error)
        )
        next
      }

      estimates[[method]][, replicate_position] <-
        as.vector(effect_estimate)
      biases[[method]][, replicate_position] <-
        as.vector(effect_estimate - causal_effect)
      predictions[[method]][, replicate_position] <-
        as.vector(effect_estimate %*% fixed_exposure)

      if (method == "sparse_mr_rr") {
        if (!identical(
          dim(loading_estimate),
          c(parameters$r_RR, parameters$px)
        ) || any(!is.finite(loading_estimate))) {
          errors <- append_error(
            errors,
            replicate_id = replicate_id,
            method = method,
            message = "Sparse loading matrix has invalid dimensions or values."
          )
          estimates[[method]][, replicate_position] <- NA_real_
          biases[[method]][, replicate_position] <- NA_real_
          predictions[[method]][, replicate_position] <- NA_real_
          next
        }

        sparse_B[, replicate_position] <-
          as.vector(loading_estimate)
      }
    }

    cat(
      sprintf(
        "Completed replicate %d (%d/%d)\n",
        replicate_id,
        replicate_position,
        replicate_count
      )
    )
    flush.console()
  }

  elapsed_seconds <-
    proc.time()[["elapsed"]] - chunk_start_time
  successful_columns <- vapply(
    estimates,
    function(value) sum(colSums(is.finite(value)) == nrow(value)),
    integer(1)
  )

  input_manifest <- data.frame(
    role = c(
      "raw simulation data",
      "raw correlation matrix",
      "frozen legacy estimator",
      "portable simulation core"
    ),
    repository_path = c(
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "data",
        "dat_1e-4.csv"
      ),
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "data",
        "rho_mat_1e-4.csv"
      ),
      file.path(
        "freeze",
        "current_analysis_20260823",
        "project",
        "scripts",
        "MR_rr_estimators.R"
      ),
      file.path(
        "paper",
        "scripts",
        "15_smoke_test_from_scratch_simulations.R"
      )
    ),
    md5 = unname(tools::md5sum(required_files)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  simulation_chunk <- list(
    metadata = list(
      schema_version = "1.0.0",
      generated_at_utc = format(
        Sys.time(),
        tz = "UTC",
        usetz = TRUE
      ),
      design = design,
      design_index = design_index,
      phase = phase,
      setting_index = setting_index,
      scenario = scenario,
      replicate_start = replicate_start,
      replicate_end = replicate_end,
      replicate_ids = replicate_ids,
      replicate_count = replicate_count,
      seed_base = seed_base,
      replicate_seeds = replicate_seeds,
      rng_kind = c(
        "L'Ecuyer-CMRG",
        "Inversion",
        "Rejection"
      ),
      rng_note = paste(
        "Each replicate has an independent deterministic seed.",
        "Results are invariant to chunk size and Slurm scheduling.",
        "Aggregate summaries, rather than archived replicate columns,",
        "are the comparison target."
      ),
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
      elapsed_seconds = elapsed_seconds,
      complete = nrow(errors) == 0L,
      archived_result_rdata_read = FALSE
    ),
    input_manifest = input_manifest,
    setting = selected_setting,
    regularization_rate = regularization_rate,
    sparse_lambda = sparse_lambda,
    parameters = parameters,
    truth = list(
      C = causal_effect,
      B = if (design == "sparse_loading") {
        sparse_truth$B
      } else {
        NULL
      },
      fixed_exposure = fixed_exposure,
      true_prediction = true_prediction
    ),
    methods = methods,
    estimates = estimates,
    biases = biases,
    predictions = predictions,
    sparse_B = sparse_B,
    successful_replicates = successful_columns,
    errors = errors,
    session_info = capture.output(utils::sessionInfo())
  )

  save_rds_atomically(
    simulation_chunk,
    output_file
  )

  cat("\nChunk output written:", output_file, "\n")
  cat(
    "Successful replicates by method:",
    paste(
      paste0(names(successful_columns), "=", successful_columns),
      collapse = ", "
    ),
    "\n"
  )
  cat(
    "Elapsed time:",
    sprintf("%.1f seconds", elapsed_seconds),
    "\n"
  )

  if (nrow(errors) > 0L) {
    print(errors, row.names = FALSE)
    stop(
      "The chunk completed with ",
      nrow(errors),
      " estimator error(s). Partial output was saved for diagnosis.",
      call. = FALSE
    )
  }

  cat("Portable simulation chunk: PASS\n")
  invisible(simulation_chunk)
}


if (sys.nframe() == 0L) {
  run_simulation_chunk()
}
