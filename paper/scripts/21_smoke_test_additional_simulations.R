#!/usr/bin/env Rscript

# From-scratch smoke test for the two additional MR-rr simulation designs.
#
# The script reads only the locked configuration files, frozen raw calibration
# inputs, and frozen estimator implementation. It runs a small number of point
# estimator fits, writes diagnostics outside the frozen snapshot, and never
# reads archived simulation-result RData.

options(warn = 1)


additional_smoke_require_file <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }

  invisible(path)
}


additional_smoke_snapshot_files <- function(root) {
  relative_paths <- sort(list.files(
    root,
    all.files = TRUE,
    full.names = FALSE,
    recursive = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  ))
  paths <- file.path(root, relative_paths)
  information <- file.info(paths)

  data.frame(
    path = relative_paths,
    size = as.double(information$size),
    modified = as.double(information$mtime),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


additional_smoke_parse_replicates <- function(value) {
  parsed <- suppressWarnings(as.integer(value))

  if (length(parsed) != 1L ||
      is.na(parsed) ||
      parsed < 1L ||
      parsed > 5L) {
    stop(
      "MRRR_ADDITIONAL_SMOKE_REPS must be an integer from 1 through 5.",
      call. = FALSE
    )
  }

  parsed
}


additional_smoke_make_method_seed <- function(
    data_seed,
    method) {
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
    stop("Unknown smoke-test method: ", method, call. = FALSE)
  }

  seed <- as.double(data_seed) + method_offset

  if (!is.finite(seed) || seed > .Machine$integer.max) {
    stop("Derived method seed exceeds the R integer range.", call. = FALSE)
  }

  as.integer(seed)
}


additional_smoke_set_seed <- function(seed) {
  set.seed(
    seed,
    kind = "L'Ecuyer-CMRG",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}


additional_smoke_fit_regularized <- function(
    simulated_data,
    parameters,
    regularization_rate,
    estimators) {
  estimators$mr_rr_regularized(
    Y = simulated_data$GAMMA_hat,
    X = simulated_data$gamma_hat,
    r = parameters$r_RR,
    Sigma_X = parameters$Sigma_X,
    regularization_rate = regularization_rate,
    W = parameters$weight.matrix
  )$AB
}


additional_smoke_fit_sparse <- function(
    simulated_data,
    parameters,
    sparse_lambda,
    estimators) {
  estimators$mr_rr_sparse(
    GAMMA_hat = simulated_data$GAMMA_hat,
    gamma_hat = simulated_data$gamma_hat,
    W = parameters$weight.matrix,
    Sigma_X = parameters$Sigma_X,
    lambda = rep(sparse_lambda, parameters$px),
    r = parameters$r_RR,
    max_iter = 100L,
    tol = 1e-2
  )
}


additional_smoke_fit_main_method <- function(
    method,
    simulated_data,
    parameters,
    regularization_rate,
    estimators) {
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
    regularized_mr_rr = additional_smoke_fit_regularized(
      simulated_data = simulated_data,
      parameters = parameters,
      regularization_rate = regularization_rate,
      estimators = estimators
    ),
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
    stop("Unknown main estimator method: ", method, call. = FALSE)
  )
}


additional_smoke_validate_estimate <- function(
    estimate,
    parameters,
    method,
    baseline,
    require_low_rank) {
  baseline$validate_effect_estimate(
    estimate = estimate,
    parameters = parameters,
    method = method,
    require_low_rank = require_low_rank
  )
}


additional_smoke_validate_sparse_fit <- function(
    sparse_fit,
    parameters,
    label,
    baseline) {
  additional_smoke_validate_estimate(
    estimate = sparse_fit$AB,
    parameters = parameters,
    method = label,
    baseline = baseline,
    require_low_rank = TRUE
  )
  baseline$assert_dimensions(
    sparse_fit$A,
    c(parameters$py, parameters$r_RR),
    paste(label, "A")
  )
  baseline$assert_dimensions(
    sparse_fit$B,
    c(parameters$r_RR, parameters$px),
    paste(label, "B")
  )

  reconstruction_difference <- max(abs(
    sparse_fit$AB - sparse_fit$A %*% sparse_fit$B
  ))

  if (!is.finite(reconstruction_difference) ||
      reconstruction_difference > 1e-8) {
    stop(
      label,
      " failed AB = A %*% B; maximum difference = ",
      format(reconstruction_difference, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  reconstruction_difference
}


additional_smoke_summary_row <- function(
    analysis_id,
    config_row,
    replicate_id,
    method,
    data_seed,
    method_seed,
    estimate,
    truth,
    sparse_reconstruction_difference = NA_real_) {
  data.frame(
    analysis_id = analysis_id,
    setting_index = config_row$setting_index,
    scenario = config_row$scenario,
    replicate_id = as.integer(replicate_id),
    truth_rank = config_row$truth_rank,
    working_rank = config_row$working_rank,
    third_singular_value = config_row$third_singular_value,
    method = method,
    data_seed = as.integer(data_seed),
    method_seed = as.integer(method_seed),
    estimated_rank = qr(estimate, tol = 1e-8)$rank,
    maximum_absolute_error = max(abs(estimate - truth)),
    causal_rmse = sqrt(mean((estimate - truth)^2)),
    sparse_reconstruction_difference =
      sparse_reconstruction_difference,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}


additional_smoke_run_rank_design <- function(
    config,
    replicates,
    truth,
    calibration,
    estimators,
    baseline,
    core) {
  summary_rows <- list()
  data_checks <- list()
  summary_index <- 0L
  check_index <- 0L

  cat("Running working-rank misspecification smoke test.\n")

  for (setting_index in seq_len(4L)) {
    selected <- config[
      config$setting_index == setting_index,
      ,
      drop = FALSE
    ]
    selected <- selected[order(selected$working_rank), , drop = FALSE]

    for (replicate_id in seq_len(replicates)) {
      generated <- lapply(
        seq_len(nrow(selected)),
        function(row_index) {
          core$additional_simulate_dataset(
            config_row = selected[row_index, , drop = FALSE],
            causal_effect = truth,
            replicate_id = replicate_id,
            calibration = calibration,
            estimators = estimators,
            baseline = baseline
          )
        }
      )
      reference_data <- generated[[1L]]$data
      gamma_difference <- max(vapply(
        generated,
        function(value) max(abs(
          value$data$gamma_hat - reference_data$gamma_hat
        )),
        numeric(1)
      ))
      outcome_difference <- max(vapply(
        generated,
        function(value) max(abs(
          value$data$GAMMA_hat - reference_data$GAMMA_hat
        )),
        numeric(1)
      ))

      if (!is.finite(gamma_difference) ||
          !is.finite(outcome_difference) ||
          gamma_difference != 0 ||
          outcome_difference != 0) {
        stop(
          "Working ranks did not receive identical data in setting ",
          setting_index,
          ", replicate ",
          replicate_id,
          ".",
          call. = FALSE
        )
      }

      check_index <- check_index + 1L
      data_checks[[check_index]] <- data.frame(
        setting_index = setting_index,
        replicate_id = replicate_id,
        data_seed = generated[[1L]]$replicate_seed,
        maximum_gamma_difference = gamma_difference,
        maximum_outcome_difference = outcome_difference,
        stringsAsFactors = FALSE,
        row.names = NULL
      )

      for (row_index in seq_len(nrow(selected))) {
        config_row <- selected[row_index, , drop = FALSE]
        simulated <- generated[[row_index]]
        parameters <- simulated$parameters
        rank_label <- paste0(
          config_row$scenario,
          "/r=",
          config_row$working_rank,
          "/rep=",
          replicate_id
        )

        method <- "regularized_mr_rr"
        method_seed <- additional_smoke_make_method_seed(
          simulated$replicate_seed,
          method
        )
        additional_smoke_set_seed(method_seed)
        regularized_estimate <- additional_smoke_fit_regularized(
          simulated_data = simulated$data,
          parameters = parameters,
          regularization_rate = config_row$regularization_rate,
          estimators = estimators
        )
        additional_smoke_validate_estimate(
          estimate = regularized_estimate,
          parameters = parameters,
          method = paste(rank_label, method),
          baseline = baseline,
          require_low_rank = TRUE
        )
        summary_index <- summary_index + 1L
        summary_rows[[summary_index]] <- additional_smoke_summary_row(
          analysis_id = "rank_misspecification",
          config_row = config_row,
          replicate_id = replicate_id,
          method = method,
          data_seed = simulated$replicate_seed,
          method_seed = method_seed,
          estimate = regularized_estimate,
          truth = truth
        )

        method <- "sparse_mr_rr"
        method_seed <- additional_smoke_make_method_seed(
          simulated$replicate_seed,
          method
        )
        additional_smoke_set_seed(method_seed)
        sparse_fit <- additional_smoke_fit_sparse(
          simulated_data = simulated$data,
          parameters = parameters,
          sparse_lambda = config_row$sparse_lambda,
          estimators = estimators
        )
        reconstruction_difference <-
          additional_smoke_validate_sparse_fit(
            sparse_fit = sparse_fit,
            parameters = parameters,
            label = paste(rank_label, method),
            baseline = baseline
          )
        summary_index <- summary_index + 1L
        summary_rows[[summary_index]] <- additional_smoke_summary_row(
          analysis_id = "rank_misspecification",
          config_row = config_row,
          replicate_id = replicate_id,
          method = method,
          data_seed = simulated$replicate_seed,
          method_seed = method_seed,
          estimate = sparse_fit$AB,
          truth = truth,
          sparse_reconstruction_difference = reconstruction_difference
        )
      }
    }

    cat("  Setting", setting_index, ": PASS\n")
  }

  list(
    summary = do.call(rbind, summary_rows),
    common_data_checks = do.call(rbind, data_checks)
  )
}


additional_smoke_run_approximate_design <- function(
    config,
    replicates,
    truth,
    calibration,
    estimators,
    baseline,
    core,
    include_mrdag,
    mrdag_niter,
    mrdag_burnin) {
  summary_rows <- list()
  summary_index <- 0L

  cat("Running approximate-low-rank smoke test.\n")

  for (setting_index in seq_len(4L)) {
    config_row <- config[
      config$setting_index == setting_index,
      ,
      drop = FALSE
    ]

    for (replicate_id in seq_len(replicates)) {
      simulated <- core$additional_simulate_dataset(
        config_row = config_row,
        causal_effect = truth,
        replicate_id = replicate_id,
        calibration = calibration,
        estimators = estimators,
        baseline = baseline
      )
      parameters <- simulated$parameters
      label <- paste0(
        config_row$scenario,
        "/rep=",
        replicate_id
      )

      main_methods <- setdiff(
        core$additional_split_methods(config_row$method_set[[1L]]),
        c("sparse_mr_rr", "mrdag")
      )
      fits <- stats::setNames(vector("list", length(main_methods)), main_methods)

      for (method in main_methods) {
        method_seed <- additional_smoke_make_method_seed(
          simulated$replicate_seed,
          method
        )
        additional_smoke_set_seed(method_seed)
        fits[[method]] <- additional_smoke_fit_main_method(
          method = method,
          simulated_data = simulated$data,
          parameters = parameters,
          regularization_rate = config_row$regularization_rate,
          estimators = estimators
        )
      }

      if (include_mrdag) {
        method_seed <- additional_smoke_make_method_seed(
          simulated$replicate_seed,
          "mrdag"
        )
        additional_smoke_set_seed(method_seed)
        fits$mrdag <- estimators$Mr_DAG(
          Y = simulated$data$GAMMA_hat,
          X = simulated$data$gamma_hat,
          niter = mrdag_niter,
          burnin = mrdag_burnin
        )
      }

      for (method in names(fits)) {
        method_seed <- additional_smoke_make_method_seed(
          simulated$replicate_seed,
          method
        )
        require_low_rank <- method %in% c(
          "naive_mr_rr",
          "mr_rr",
          "regularized_mr_rr"
        )
        additional_smoke_validate_estimate(
          estimate = fits[[method]],
          parameters = parameters,
          method = paste(label, method),
          baseline = baseline,
          require_low_rank = require_low_rank
        )
        summary_index <- summary_index + 1L
        summary_rows[[summary_index]] <- additional_smoke_summary_row(
          analysis_id = "approximate_low_rank",
          config_row = config_row,
          replicate_id = replicate_id,
          method = method,
          data_seed = simulated$replicate_seed,
          method_seed = method_seed,
          estimate = fits[[method]],
          truth = truth
        )
      }

      method <- "sparse_mr_rr"
      method_seed <- additional_smoke_make_method_seed(
        simulated$replicate_seed,
        method
      )
      additional_smoke_set_seed(method_seed)
      sparse_fit <- additional_smoke_fit_sparse(
        simulated_data = simulated$data,
        parameters = parameters,
        sparse_lambda = config_row$sparse_lambda,
        estimators = estimators
      )
      reconstruction_difference <-
        additional_smoke_validate_sparse_fit(
          sparse_fit = sparse_fit,
          parameters = parameters,
          label = paste(label, method),
          baseline = baseline
        )
      summary_index <- summary_index + 1L
      summary_rows[[summary_index]] <- additional_smoke_summary_row(
        analysis_id = "approximate_low_rank",
        config_row = config_row,
        replicate_id = replicate_id,
        method = method,
        data_seed = simulated$replicate_seed,
        method_seed = method_seed,
        estimate = sparse_fit$AB,
        truth = truth,
        sparse_reconstruction_difference = reconstruction_difference
      )
    }

    cat("  Setting", setting_index, ": PASS\n")
  }

  do.call(rbind, summary_rows)
}


additional_run_simulation_smoke <- function() {
  script_argument <- grep(
    "^--file=",
    commandArgs(trailingOnly = FALSE),
    value = TRUE
  )

  if (length(script_argument) == 0L) {
    stop("Could not determine the smoke-test script path.", call. = FALSE)
  }

  script_path <- normalizePath(
    sub("^--file=", "", script_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  core_file <- file.path(
    dirname(script_path),
    "20_additional_simulation_core.R"
  )
  additional_smoke_require_file(core_file, "Additional simulation core")
  core <- new.env(parent = globalenv())
  sys.source(core_file, envir = core)
  repo_root <- core$additional_find_repo_root(dirname(script_path))
  freeze_root <- file.path(
    repo_root,
    "freeze",
    "current_analysis_20260823",
    "project"
  )
  freeze_before <- additional_smoke_snapshot_files(freeze_root)
  replicates <- additional_smoke_parse_replicates(Sys.getenv(
    "MRRR_ADDITIONAL_SMOKE_REPS",
    unset = "2"
  ))
  baseline <- core$additional_load_baseline_core(repo_root)
  include_mrdag <- baseline$parse_boolean(
    Sys.getenv(
      "MRRR_ADDITIONAL_SMOKE_INCLUDE_MRDAG",
      unset = "true"
    ),
    "MRRR_ADDITIONAL_SMOKE_INCLUDE_MRDAG"
  )
  mrdag_niter <- baseline$parse_integer(
    Sys.getenv(
      "MRRR_ADDITIONAL_SMOKE_MRDAG_NITER",
      unset = "1000"
    ),
    "MRRR_ADDITIONAL_SMOKE_MRDAG_NITER",
    minimum = 1L
  )
  mrdag_burnin <- baseline$parse_integer(
    Sys.getenv(
      "MRRR_ADDITIONAL_SMOKE_MRDAG_BURNIN",
      unset = "200"
    ),
    "MRRR_ADDITIONAL_SMOKE_MRDAG_BURNIN",
    minimum = 0L
  )

  if (include_mrdag && mrdag_burnin >= mrdag_niter) {
    stop(
      "MRRR_ADDITIONAL_SMOKE_MRDAG_BURNIN must be smaller than ",
      "MRRR_ADDITIONAL_SMOKE_MRDAG_NITER.",
      call. = FALSE
    )
  }

  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package `MASS` is required.", call. = FALSE)
  }

  rank_config <- core$additional_read_config(
    file.path(
      repo_root,
      "paper",
      "config",
      "rank_misspecification_settings.csv"
    ),
    "rank_misspecification"
  )
  approximate_config <- core$additional_read_config(
    file.path(
      repo_root,
      "paper",
      "config",
      "approximate_low_rank_settings.csv"
    ),
    "approximate_low_rank"
  )
  estimators <- core$additional_load_estimators(repo_root)
  calibration <- core$additional_build_calibration(
    repo_root = repo_root,
    baseline = baseline,
    pz = 177L
  )
  exact_truth <- core$additional_make_generic_truth(
    third_singular_value = 0,
    px = calibration$px,
    py = calibration$py,
    seed = 123L
  )$C
  delta <- unique(approximate_config$third_singular_value)
  approximate_truth <- core$additional_make_generic_truth(
    third_singular_value = delta,
    px = calibration$px,
    py = calibration$py,
    seed = 123L
  )$C
  output_setting <- trimws(Sys.getenv(
    "MRRR_ADDITIONAL_OUTPUT_DIR",
    unset = ""
  ))

  if (nzchar(output_setting)) {
    output_root <- normalizePath(
      output_setting,
      winslash = "/",
      mustWork = FALSE
    )
  } else {
    output_root <- file.path(
      repo_root,
      "paper",
      "output",
      "additional_simulations",
      "smoke"
    )
  }

  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  output_root <- normalizePath(
    output_root,
    winslash = "/",
    mustWork = TRUE
  )

  cat("MR-rr additional-simulation from-scratch smoke test\n")
  cat("Repository root:", repo_root, "\n")
  cat("Replicates per setting:", replicates, "\n")
  cat("Configured approximate-low-rank delta:", delta, "\n")
  cat("MrDAG included:", include_mrdag, "\n")
  cat("Archived simulation-result RData read: no\n")
  cat("Output directory:", output_root, "\n")

  rank_result <- additional_smoke_run_rank_design(
    config = rank_config,
    replicates = replicates,
    truth = exact_truth,
    calibration = calibration,
    estimators = estimators,
    baseline = baseline,
    core = core
  )
  approximate_summary <- additional_smoke_run_approximate_design(
    config = approximate_config,
    replicates = replicates,
    truth = approximate_truth,
    calibration = calibration,
    estimators = estimators,
    baseline = baseline,
    core = core,
    include_mrdag = include_mrdag,
    mrdag_niter = mrdag_niter,
    mrdag_burnin = mrdag_burnin
  )
  summary_table <- rbind(
    rank_result$summary,
    approximate_summary
  )
  rownames(summary_table) <- NULL
  expected_rank_rows <- 4L * replicates * 3L * 2L
  expected_approximate_methods <- if (include_mrdag) 7L else 6L
  expected_approximate_rows <-
    4L * replicates * expected_approximate_methods
  result_key_columns <- c(
    "analysis_id",
    "setting_index",
    "replicate_id",
    "working_rank",
    "method"
  )
  result_keys <- do.call(
    paste,
    c(summary_table[, result_key_columns, drop = FALSE], sep = "/")
  )

  if (nrow(rank_result$summary) != expected_rank_rows ||
      nrow(approximate_summary) != expected_approximate_rows ||
      any(!is.finite(summary_table$maximum_absolute_error)) ||
      any(!is.finite(summary_table$causal_rmse)) ||
      anyDuplicated(result_keys)) {
    stop("Smoke-test result inventory is invalid.", call. = FALSE)
  }

  expected_approximate_methods <- core$additional_split_methods(
    approximate_config$method_set[[1L]]
  )

  if (!include_mrdag) {
    expected_approximate_methods <- setdiff(
      expected_approximate_methods,
      "mrdag"
    )
  }

  if (!setequal(
    unique(approximate_summary$method),
    expected_approximate_methods
  )) {
    stop(
      "Approximate-low-rank smoke test did not run the expected methods.",
      call. = FALSE
    )
  }

  summary_file <- file.path(
    output_root,
    "additional_simulation_smoke_summary.csv"
  )
  data_check_file <- file.path(
    output_root,
    "rank_common_data_checks.csv"
  )
  session_file <- file.path(output_root, "sessionInfo.txt")
  utils::write.csv(
    summary_table,
    summary_file,
    row.names = FALSE,
    na = ""
  )
  utils::write.csv(
    rank_result$common_data_checks,
    data_check_file,
    row.names = FALSE
  )
  writeLines(
    capture.output(utils::sessionInfo()),
    session_file,
    useBytes = TRUE
  )
  freeze_after <- additional_smoke_snapshot_files(freeze_root)

  if (!identical(freeze_before, freeze_after)) {
    stop(
      "The frozen snapshot changed during the smoke test.",
      call. = FALSE
    )
  }

  cat("Rank-misspecification estimator rows:", expected_rank_rows, "\n")
  cat(
    "Approximate-low-rank estimator rows:",
    expected_approximate_rows,
    "\n"
  )
  cat("Identical simulated data across working ranks: PASS\n")
  cat("Finite estimator outputs and rank constraints: PASS\n")
  cat("Frozen snapshot unchanged: PASS\n")
  cat("Additional-simulation from-scratch smoke test: PASS\n")

  invisible(list(
    summary = summary_table,
    common_data_checks = rank_result$common_data_checks,
    output_files = c(summary_file, data_check_file, session_file)
  ))
}


if (sys.nframe() == 0L) {
  additional_run_simulation_smoke()
}
