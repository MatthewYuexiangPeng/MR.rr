#!/usr/bin/env Rscript

# Run one bootstrap-inference chunk for the additional MR-rr simulations.
#
# Each Monte Carlo replicate is regenerated from its locked replicate seed.
# Within a replicate, every method and working rank uses the same SNP-level
# nonparametric bootstrap samples. The worker stores only the entrywise
# bootstrap SE, percentile interval, and coverage—not the full draw matrices.
# It reads the frozen calibration data and estimator implementation but never
# reads archived simulation-result RData and never writes under freeze/.

options(warn = 1)


additional_bootstrap_parse_arguments <- function(arguments) {
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


additional_bootstrap_validate_argument_names <- function(
    arguments,
    allowed) {
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


additional_bootstrap_require_argument <- function(arguments, name) {
  value <- arguments[[name]]

  if (is.null(value) || !nzchar(trimws(value))) {
    stop("Missing required argument --", name, ".", call. = FALSE)
  }

  trimws(value)
}


additional_bootstrap_parse_integer <- function(
    arguments,
    name,
    default = NULL,
    minimum = 0L) {
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


additional_bootstrap_parse_boolean <- function(
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


additional_bootstrap_resolve_path <- function(path, repo_root) {
  if (grepl("^(/|[A-Za-z]:[/\\\\])", path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }

  normalizePath(
    file.path(repo_root, path),
    winslash = "/",
    mustWork = FALSE
  )
}


additional_bootstrap_repo_path <- function(path, repo_root) {
  normalized_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  normalized_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  prefix <- paste0(normalized_root, "/")

  if (!startsWith(normalized_path, prefix)) {
    stop("Required input is outside the repository: ", path, call. = FALSE)
  }

  substring(normalized_path, nchar(prefix) + 1L)
}


additional_bootstrap_set_seed <- function(seed) {
  set.seed(
    seed,
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}


additional_bootstrap_make_resample_seed <- function(data_seed) {
  seed <- as.double(data_seed) + 900000

  if (!is.finite(seed) || seed > .Machine$integer.max) {
    stop("Derived bootstrap-resample seed is invalid.", call. = FALSE)
  }

  as.integer(seed)
}


additional_bootstrap_make_fit_seed <- function(
    data_seed,
    method,
    working_rank,
    bootstrap_id) {
  method_index <- match(
    method,
    c(
      "ivw",
      "srivw",
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr",
      "sparse_mr_rr",
      "mrdag"
    )
  )

  if (is.na(method_index)) {
    stop("Unknown bootstrap method: ", method, call. = FALSE)
  }

  seed <-
    as.double(data_seed) +
    2000000 +
    as.double(method_index) * 100000 +
    as.double(working_rank) * 10000 +
    as.double(bootstrap_id)

  if (!is.finite(seed) || seed > .Machine$integer.max) {
    stop("Derived bootstrap-fit seed is invalid.", call. = FALSE)
  }

  as.integer(seed)
}


additional_bootstrap_result_specification <- function(
    analysis,
    phase,
    selected_config,
    point) {
  point_specification <- point$additional_chunk_result_specification(
    analysis = analysis,
    phase = phase,
    selected_config = selected_config
  )
  specification <- point_specification$specification

  if (analysis == "rank_misspecification") {
    specification <- specification[
      specification$method == "regularized_mr_rr",
      ,
      drop = FALSE
    ]
  } else if (phase == "standard") {
    specification <- specification[
      specification$method != "sparse_mr_rr",
      ,
      drop = FALSE
    ]
  }

  rownames(specification) <- NULL

  if (nrow(specification) == 0L ||
      anyDuplicated(specification$result_key)) {
    stop("Bootstrap result specification is invalid.", call. = FALSE)
  }

  if (analysis == "rank_misspecification" &&
      !identical(
        specification$result_key,
        paste0("r", 1:3, "__regularized_mr_rr")
      )) {
    stop(
      "Rank-misspecification bootstrap methods are incorrect.",
      call. = FALSE
    )
  }

  expected_approximate <- if (phase == "standard") {
    c(
      "ivw",
      "srivw",
      "naive_mr_rr",
      "mr_rr",
      "regularized_mr_rr"
    )
  } else {
    "mrdag"
  }

  if (analysis == "approximate_low_rank" &&
      !identical(specification$result_key, expected_approximate)) {
    stop(
      "Approximate-low-rank bootstrap methods are incorrect.",
      call. = FALSE
    )
  }

  list(
    config = point_specification$config,
    specification = specification
  )
}


additional_bootstrap_empty_errors <- function() {
  data.frame(
    replicate_id = integer(0),
    bootstrap_id = integer(0),
    result_key = character(0),
    message = character(0),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


additional_bootstrap_draw <- function(
    bootstrap_id,
    indices,
    simulated_data,
    parameters_by_row,
    selected_config,
    result_specification,
    data_seed,
    estimators,
    baseline,
    point,
    mrdag_niter,
    mrdag_burnin) {
  sampled_data <- simulated_data
  sampled_data$gamma_hat <- simulated_data$gamma_hat[indices, , drop = FALSE]
  sampled_data$GAMMA_hat <- simulated_data$GAMMA_hat[indices, , drop = FALSE]
  result_keys <- result_specification$result_key
  estimates <- stats::setNames(
    lapply(result_keys, function(key) NULL),
    result_keys
  )
  errors <- additional_bootstrap_empty_errors()

  for (specification_index in seq_len(nrow(result_specification))) {
    specification <- result_specification[
      specification_index,
      ,
      drop = FALSE
    ]
    result_key <- specification$result_key[[1L]]
    method <- specification$method[[1L]]
    config_row_index <- specification$config_row_index[[1L]]
    parameters <- parameters_by_row[[config_row_index]]
    config_row <- selected_config[config_row_index, , drop = FALSE]
    fit_seed <- additional_bootstrap_make_fit_seed(
      data_seed = data_seed,
      method = method,
      working_rank = parameters$r_RR,
      bootstrap_id = bootstrap_id
    )
    additional_bootstrap_set_seed(fit_seed)
    fit <- tryCatch(
      point$additional_chunk_fit_method(
        method = method,
        simulated_data = sampled_data,
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
      errors <- rbind(
        errors,
        data.frame(
          replicate_id = NA_integer_,
          bootstrap_id = bootstrap_id,
          result_key = result_key,
          message = conditionMessage(fit),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      )
      next
    }

    effect_estimate <- tryCatch(
      point$additional_chunk_validate_fit(
        fit = fit,
        method = method,
        parameters = parameters,
        result_key = paste0(
          "bootstrap=",
          bootstrap_id,
          "/",
          result_key
        ),
        baseline = baseline
      ),
      error = function(error) error
    )

    if (inherits(effect_estimate, "error")) {
      errors <- rbind(
        errors,
        data.frame(
          replicate_id = NA_integer_,
          bootstrap_id = bootstrap_id,
          result_key = result_key,
          message = conditionMessage(effect_estimate),
          stringsAsFactors = FALSE,
          row.names = NULL
        )
      )
      next
    }

    estimates[[result_key]] <- as.vector(effect_estimate)
  }

  list(estimates = estimates, errors = errors)
}


additional_bootstrap_map_draws <- function(
    bootstrap_ids,
    run_draw,
    n_cores) {
  effective_cores <- min(length(bootstrap_ids), n_cores)

  if (.Platform$OS.type != "windows" && effective_cores > 1L) {
    return(parallel::mclapply(
      bootstrap_ids,
      run_draw,
      mc.cores = effective_cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    ))
  }

  lapply(bootstrap_ids, run_draw)
}


additional_bootstrap_summarize_draws <- function(
    draw_results,
    result_key,
    truth_vector,
    effect_dimension) {
  draw_matrix <- matrix(
    NA_real_,
    nrow = length(draw_results),
    ncol = effect_dimension
  )

  for (bootstrap_id in seq_along(draw_results)) {
    estimate <- draw_results[[bootstrap_id]]$estimates[[result_key]]

    if (!is.null(estimate)) {
      if (length(estimate) != effect_dimension ||
          any(!is.finite(estimate))) {
        stop(
          "A bootstrap draw has invalid dimensions or values for ",
          result_key,
          ".",
          call. = FALSE
        )
      }

      draw_matrix[bootstrap_id, ] <- estimate
    }
  }

  successful <- rowSums(is.finite(draw_matrix)) == effect_dimension
  successful_count <- sum(successful)

  if (successful_count < 2L) {
    return(list(
      standard_error = rep(NA_real_, effect_dimension),
      ci_lower = rep(NA_real_, effect_dimension),
      ci_upper = rep(NA_real_, effect_dimension),
      coverage = rep(NA_real_, effect_dimension),
      successful_draws = successful_count
    ))
  }

  successful_matrix <- draw_matrix[successful, , drop = FALSE]
  standard_error <- apply(
    successful_matrix,
    2L,
    stats::sd,
    na.rm = TRUE
  )
  interval <- apply(
    successful_matrix,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE,
    type = 7L
  )

  if (effect_dimension == 1L) {
    interval <- matrix(interval, nrow = 2L)
  }

  ci_lower <- interval[1L, ]
  ci_upper <- interval[2L, ]

  list(
    standard_error = standard_error,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    coverage = as.numeric(
      truth_vector >= ci_lower & truth_vector <= ci_upper
    ),
    successful_draws = successful_count
  )
}


additional_bootstrap_matrix_list <- function(
    result_keys,
    effect_dimension,
    replicate_ids,
    mode = "numeric") {
  stats::setNames(
    lapply(result_keys, function(result_key) {
      initial <- if (mode == "integer") NA_integer_ else NA_real_
      matrix(
        initial,
        nrow = effect_dimension,
        ncol = length(replicate_ids),
        dimnames = list(NULL, as.character(replicate_ids))
      )
    }),
    result_keys
  )
}


additional_bootstrap_seed_list <- function(
    result_specification,
    data_seeds,
    bootstrap_size,
    replicate_ids) {
  result_keys <- result_specification$result_key
  seeds <- stats::setNames(
    lapply(result_keys, function(result_key) {
      matrix(
        NA_integer_,
        nrow = bootstrap_size,
        ncol = length(replicate_ids),
        dimnames = list(
          as.character(seq_len(bootstrap_size)),
          as.character(replicate_ids)
        )
      )
    }),
    result_keys
  )

  for (specification_index in seq_len(nrow(result_specification))) {
    specification <- result_specification[
      specification_index,
      ,
      drop = FALSE
    ]
    key <- specification$result_key[[1L]]

    for (replicate_position in seq_along(replicate_ids)) {
      seeds[[key]][, replicate_position] <- vapply(
        seq_len(bootstrap_size),
        function(bootstrap_id) {
          additional_bootstrap_make_fit_seed(
            data_seed = data_seeds[[replicate_position]],
            method = specification$method[[1L]],
            working_rank = specification$working_rank[[1L]],
            bootstrap_id = bootstrap_id
          )
        },
        integer(1)
      )
    }
  }

  seeds
}


additional_bootstrap_run <- function() {
  arguments <- additional_bootstrap_parse_arguments(
    commandArgs(trailingOnly = TRUE)
  )
  allowed_arguments <- c(
    "analysis",
    "phase",
    "setting",
    "start",
    "end",
    "seed-base",
    "bootstrap-size",
    "allow-bootstrap-override",
    "n-cores",
    "mrdag-niter",
    "mrdag-burnin",
    "output-dir",
    "allow-provisional",
    "overwrite"
  )
  additional_bootstrap_validate_argument_names(
    arguments,
    allowed_arguments
  )
  analysis <- additional_bootstrap_require_argument(arguments, "analysis")
  phase <- additional_bootstrap_require_argument(arguments, "phase")
  setting_index <- additional_bootstrap_parse_integer(
    arguments,
    "setting",
    minimum = 1L
  )
  replicate_start <- additional_bootstrap_parse_integer(
    arguments,
    "start",
    minimum = 1L
  )
  replicate_end <- additional_bootstrap_parse_integer(
    arguments,
    "end",
    minimum = 1L
  )
  seed_base <- additional_bootstrap_parse_integer(
    arguments,
    "seed-base",
    default = 123L,
    minimum = 1L
  )
  requested_bootstrap_size <- if (is.null(arguments[["bootstrap-size"]])) {
    NULL
  } else {
    additional_bootstrap_parse_integer(
      arguments,
      "bootstrap-size",
      minimum = 2L
    )
  }
  allow_bootstrap_override <- additional_bootstrap_parse_boolean(
    arguments,
    "allow-bootstrap-override",
    default = FALSE
  )
  default_cores <- suppressWarnings(as.integer(Sys.getenv(
    "SLURM_CPUS_PER_TASK",
    unset = "1"
  )))

  if (length(default_cores) != 1L ||
      is.na(default_cores) ||
      default_cores < 1L) {
    default_cores <- 1L
  }

  n_cores <- additional_bootstrap_parse_integer(
    arguments,
    "n-cores",
    default = default_cores,
    minimum = 1L
  )
  mrdag_niter <- additional_bootstrap_parse_integer(
    arguments,
    "mrdag-niter",
    default = 1000L,
    minimum = 1L
  )
  mrdag_burnin <- additional_bootstrap_parse_integer(
    arguments,
    "mrdag-burnin",
    default = 200L,
    minimum = 0L
  )
  allow_provisional <- additional_bootstrap_parse_boolean(
    arguments,
    "allow-provisional",
    default = FALSE
  )
  overwrite <- additional_bootstrap_parse_boolean(
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
  scripts_directory <- dirname(worker_file)
  core_file <- file.path(
    scripts_directory,
    "20_additional_simulation_core.R"
  )
  point_file <- file.path(
    scripts_directory,
    "22_run_additional_simulation_chunk.R"
  )

  if (!file.exists(core_file) || !file.exists(point_file)) {
    stop("Additional simulation dependencies are missing.", call. = FALSE)
  }

  core <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  point <- new.env(parent = globalenv())
  sys.source(point_file, envir = point)
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

  configured_total <- unique(selected_config$monte_carlo_replicates)

  if (replicate_end > configured_total) {
    stop(
      "Requested replicate range exceeds the configured total.",
      call. = FALSE
    )
  }

  if (seed_base != unique(selected_config$seed_base)) {
    stop("--seed-base differs from the configuration.", call. = FALSE)
  }

  configured_bootstrap_size <- unique(selected_config$bootstrap_size)
  bootstrap_size <- if (is.null(requested_bootstrap_size)) {
    configured_bootstrap_size
  } else {
    requested_bootstrap_size
  }

  if (bootstrap_size != configured_bootstrap_size &&
      !allow_bootstrap_override) {
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

  specification_object <- additional_bootstrap_result_specification(
    analysis = analysis,
    phase = phase,
    selected_config = selected_config,
    point = point
  )
  selected_config <- specification_object$config
  result_specification <- specification_object$specification
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
  truth_vector <- as.vector(causal_effect)
  replicate_ids <- seq.int(replicate_start, replicate_end)
  replicate_count <- length(replicate_ids)
  effect_dimension <- calibration$py * calibration$px
  result_keys <- result_specification$result_key
  standard_errors <- additional_bootstrap_matrix_list(
    result_keys,
    effect_dimension,
    replicate_ids
  )
  ci_lower <- additional_bootstrap_matrix_list(
    result_keys,
    effect_dimension,
    replicate_ids
  )
  ci_upper <- additional_bootstrap_matrix_list(
    result_keys,
    effect_dimension,
    replicate_ids
  )
  coverage <- additional_bootstrap_matrix_list(
    result_keys,
    effect_dimension,
    replicate_ids,
    mode = "integer"
  )
  data_seeds <- vapply(
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
  resample_seeds <- vapply(
    data_seeds,
    additional_bootstrap_make_resample_seed,
    integer(1)
  )
  fit_seeds <- additional_bootstrap_seed_list(
    result_specification = result_specification,
    data_seeds = data_seeds,
    bootstrap_size = bootstrap_size,
    replicate_ids = replicate_ids
  )
  successful_draws <- matrix(
    NA_integer_,
    nrow = length(result_keys),
    ncol = replicate_count,
    dimnames = list(result_keys, as.character(replicate_ids))
  )
  errors <- additional_bootstrap_empty_errors()
  minimum_successful_draws <- max(
    2L,
    as.integer(ceiling(0.9 * bootstrap_size))
  )
  output_argument <- arguments[["output-dir"]]

  if (is.null(output_argument) || !nzchar(trimws(output_argument))) {
    output_directory <- file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "bootstrap_chunks"
    )
  } else {
    output_directory <- additional_bootstrap_resolve_path(
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
      "%s__%s__setting-%d__rep-%04d-%04d__bt-%04d.rds",
      analysis,
      phase,
      setting_index,
      replicate_start,
      replicate_end,
      bootstrap_size
    )
  )

  if (file.exists(output_file) && !overwrite) {
    stop(
      "Output already exists; use --overwrite=true to replace it: ",
      output_file,
      call. = FALSE
    )
  }

  effective_cores <- if (.Platform$OS.type == "windows") {
    1L
  } else {
    min(n_cores, bootstrap_size)
  }

  cat("MR-rr additional-simulation bootstrap chunk worker\n")
  cat("Repository root:", repo_root, "\n")
  cat("Analysis:", analysis, "\n")
  cat("Phase:", phase, "\n")
  cat("Setting:", setting_index, "\n")
  cat("Replicates:", replicate_start, "through", replicate_end, "\n")
  cat("Bootstrap draws per replicate:", bootstrap_size, "\n")
  cat("Configured bootstrap draws:", configured_bootstrap_size, "\n")
  cat("Parallel bootstrap workers:", effective_cores, "\n")
  cat("Third singular value:", delta, "\n")
  cat(
    "Configuration status:",
    paste(unique(selected_config$status), collapse = ", "),
    "\n"
  )
  cat("Result keys:", paste(result_keys, collapse = ", "), "\n")
  cat("Point-result files read: no\n")
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

    if (simulated$replicate_seed != data_seeds[[replicate_position]]) {
      stop("Regenerated data seed is inconsistent.", call. = FALSE)
    }

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
    additional_bootstrap_set_seed(
      resample_seeds[[replicate_position]]
    )
    bootstrap_indices <- lapply(
      seq_len(bootstrap_size),
      function(bootstrap_id) {
        sample.int(
          n = nrow(simulated$data$gamma_hat),
          size = nrow(simulated$data$gamma_hat),
          replace = TRUE
        )
      }
    )
    run_draw <- function(bootstrap_id) {
      tryCatch(
        additional_bootstrap_draw(
          bootstrap_id = bootstrap_id,
          indices = bootstrap_indices[[bootstrap_id]],
          simulated_data = simulated$data,
          parameters_by_row = parameters_by_row,
          selected_config = selected_config,
          result_specification = result_specification,
          data_seed = simulated$replicate_seed,
          estimators = estimators,
          baseline = baseline,
          point = point,
          mrdag_niter = mrdag_niter,
          mrdag_burnin = mrdag_burnin
        ),
        error = function(error) {
          list(
            fatal_error = conditionMessage(error),
            bootstrap_id = bootstrap_id
          )
        }
      )
    }
    draw_results <- additional_bootstrap_map_draws(
      bootstrap_ids = seq_len(bootstrap_size),
      run_draw = run_draw,
      n_cores = effective_cores
    )
    fatal <- vapply(
      draw_results,
      function(result) !is.null(result$fatal_error),
      logical(1)
    )

    if (any(fatal)) {
      first <- draw_results[[which(fatal)[[1L]]]]
      stop(
        "Bootstrap draw ",
        first$bootstrap_id,
        " failed outside estimator error handling: ",
        first$fatal_error,
        call. = FALSE
      )
    }

    replicate_errors <- do.call(
      rbind,
      lapply(draw_results, function(result) result$errors)
    )

    if (nrow(replicate_errors) > 0L) {
      replicate_errors$replicate_id <- replicate_id
      errors <- rbind(errors, replicate_errors)
    }

    for (result_key in result_keys) {
      summary <- additional_bootstrap_summarize_draws(
        draw_results = draw_results,
        result_key = result_key,
        truth_vector = truth_vector,
        effect_dimension = effect_dimension
      )
      standard_errors[[result_key]][, replicate_position] <-
        summary$standard_error
      ci_lower[[result_key]][, replicate_position] <- summary$ci_lower
      ci_upper[[result_key]][, replicate_position] <- summary$ci_upper
      coverage[[result_key]][, replicate_position] <- summary$coverage
      successful_draws[result_key, replicate_position] <-
        summary$successful_draws
    }

    cat(
      sprintf(
        "Completed replicate %d (%d/%d): %s\n",
        replicate_id,
        replicate_position,
        replicate_count,
        paste(
          paste0(
            result_keys,
            "=",
            successful_draws[, replicate_position],
            "/",
            bootstrap_size
          ),
          collapse = ", "
        )
      )
    )
    flush.console()
  }

  elapsed_seconds <- proc.time()[["elapsed"]] - start_time
  enough_draws <- successful_draws >= minimum_successful_draws
  output_values_finite <- vapply(
    c(standard_errors, ci_lower, ci_upper, coverage),
    function(value) all(is.finite(value)),
    logical(1)
  )
  complete <- all(enough_draws) && all(output_values_finite)
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
    point_file,
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
    repository_path = vapply(
      required_files,
      additional_bootstrap_repo_path,
      character(1),
      repo_root = repo_root
    ),
    md5 = unname(tools::md5sum(required_files)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  bootstrap_chunk <- list(
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
      resample_seeds = resample_seeds,
      rng_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
      bootstrap_size = bootstrap_size,
      configured_bootstrap_size = configured_bootstrap_size,
      bootstrap_size_override = bootstrap_size != configured_bootstrap_size,
      minimum_successful_draws = minimum_successful_draws,
      requested_cores = n_cores,
      effective_cores = effective_cores,
      parallel_strategy = if (effective_cores > 1L) {
        "forked_bootstrap_draws"
      } else {
        "sequential_bootstrap_draws"
      },
      mrdag_niter = if (phase == "mrdag") mrdag_niter else NA_integer_,
      mrdag_burnin = if (phase == "mrdag") {
        mrdag_burnin
      } else {
        NA_integer_
      },
      allow_provisional = allow_provisional,
      elapsed_seconds = elapsed_seconds,
      complete = complete,
      all_bootstrap_draws_successful = all(
        successful_draws == bootstrap_size
      ),
      point_result_files_read = FALSE,
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
    fit_seeds = fit_seeds,
    successful_draws = successful_draws,
    standard_errors = standard_errors,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    coverage = coverage,
    errors = errors,
    session_info = capture.output(utils::sessionInfo())
  )
  point$additional_chunk_save_rds(bootstrap_chunk, output_file)

  cat("\nBootstrap chunk output written:", output_file, "\n")
  cat(
    "Minimum successful draws required:",
    minimum_successful_draws,
    "of",
    bootstrap_size,
    "\n"
  )
  cat("Recorded estimator failures:", nrow(errors), "\n")
  cat("Elapsed time:", sprintf("%.1f seconds", elapsed_seconds), "\n")

  if (!complete) {
    stop(
      "Bootstrap chunk failed completeness or finite-output checks.",
      call. = FALSE
    )
  }

  cat("Common resamples across methods and working ranks: PASS\n")
  cat("Point-result files read: no\n")
  cat("Archived simulation-result RData independence: PASS\n")
  cat("Frozen snapshot written: no\n")
  cat("Additional-simulation bootstrap chunk: PASS\n")
  invisible(bootstrap_chunk)
}


if (sys.nframe() == 0L) {
  additional_bootstrap_run()
}
