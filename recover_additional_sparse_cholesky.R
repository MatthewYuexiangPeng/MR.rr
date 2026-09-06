#!/usr/bin/env Rscript

# Recover the three failed sparse fits for setting 1, replicate 840.
# Original chunks, manifests, frozen inputs, and bootstrap results remain intact.
# The original input_manifest describes initial production; numerical_recovery
# records the additional code and computations used for the recovered estimates.
# Version 2 fixes truth reconstruction after the worker switches to L'Ecuyer-CMRG.

options(warn = 1)

recovery_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

recovery_hash <- function(paths) {
  recovery_assert(all(file.exists(paths)), "A required provenance file is missing.")
  stats::setNames(unname(tools::md5sum(paths)), paths)
}

recovery_with_truth_rng <- function(expression) {
  # Production constructs truth in a fresh --vanilla session, before the worker
  # switches RNG kind for simulation and estimator streams. Reproduce that
  # truth explicitly, without changing the caller's kind or random state.
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  old_kind <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  force(expression)
}

recovery_truth_function <- function(original_truth) {
  force(original_truth)
  function(third_singular_value = 0, px = 9L, py = 3L, seed = 123L) {
    recovery_with_truth_rng(original_truth(third_singular_value, px, py, seed))
  }
}

recovery_test_truth_rng <- function(original_truth, safe_truth) {
  recovery_with_truth_rng({
    for (delta in c(0, 0.1)) {
      RNGkind("Mersenne-Twister", "Inversion", "Rejection")
      expected <- original_truth(third_singular_value = delta, seed = 123L)
      for (kind in c("Mersenne-Twister", "L'Ecuyer-CMRG")) {
        RNGkind(kind, "Inversion", "Rejection")
        set.seed(840L)
        before_kind <- RNGkind()
        before_seed <- get(".Random.seed", envir = .GlobalEnv)
        actual <- safe_truth(third_singular_value = delta, seed = 123L)
        recovery_assert(identical(actual, expected),
                        paste("Truth depends on caller RNG kind:", kind, "delta", delta))
        recovery_assert(identical(RNGkind(), before_kind) &&
                          identical(get(".Random.seed", envir = .GlobalEnv), before_seed),
                        "Truth reconstruction changed the caller's RNG stream.")
        error <- tryCatch(safe_truth(third_singular_value = -1), error = identity)
        recovery_assert(inherits(error, "error") && identical(RNGkind(), before_kind) &&
                          identical(get(".Random.seed", envir = .GlobalEnv), before_seed),
                        "Failed truth reconstruction changed the caller's RNG stream.")
      }
    }
    rm(".Random.seed", envir = .GlobalEnv)
    invisible(safe_truth(seed = 123L))
    recovery_assert(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE),
                    "Truth reconstruction did not restore the absent-seed state.")
  })
  cat("Truth RNG isolation (both analyses, both RNG kinds, error restoration): PASS\n")
  invisible(TRUE)
}

recovery_chol <- function(x) {
  tryCatch(base::chol(x), error = function(e) e)
}

recovery_matrix_diagnostic <- function(x, label, old_is_psd) {
  ev <- eigen((x + t(x)) / 2, symmetric = TRUE, only.values = TRUE)$values
  factor <- recovery_chol(x)
  data.frame(
    matrix = label,
    min_eigenvalue = min(ev),
    max_eigenvalue = max(ev),
    relative_min_eigenvalue = min(ev) / max(abs(ev)),
    old_is_psd = old_is_psd(x),
    chol_success = !inherits(factor, "error"),
    chol_error = if (inherits(factor, "error")) conditionMessage(factor) else "",
    stringsAsFactors = FALSE
  )
}

recovery_make_sparse <- function(estimators) {
  original <- estimators$mr_rr_sparse
  statements <- as.list(body(original))[-1L]
  expected_guard <- quote(if (!is_psd(matrix_part2)) {
    matrix_part2 <- .nearest_psd(matrix_part2, epsilon = 1e-6)
  })
  recovery_assert(sum(vapply(statements, identical, logical(1),
                            y = expected_guard)) == 1L,
                  "Frozen sparse implementation differs from the reviewed projection rule.")
  for (statement in list(quote(R <- chol(matrix_part1)),
                         quote(Q <- chol(matrix_part2)))) {
    recovery_assert(sum(vapply(statements, identical, logical(1),
                              y = statement)) == 1L,
                    "Frozen sparse Cholesky expressions differ from the reviewed code.")
  }
  recovery_assert(identical(formals(estimators$is_psd)$tol, 1e-8),
                  "The frozen PSD tolerance differs from 1e-8.")

  original_environment <- environment(original)
  old_is_psd <- estimators$is_psd
  state <- new.env(parent = emptyenv())
  state$fallbacks <- 0L
  scope <- new.env(parent = original_environment)
  scope$is_psd <- function(A, tol = 1e-8) {
    if (!old_is_psd(A, tol = tol)) return(FALSE)
    factor <- tryCatch(base::chol(A), error = function(e) e)
    if (inherits(factor, "error")) {
      state$fallbacks <- state$fallbacks + 1L
      return(FALSE)
    }
    TRUE
  }
  patched <- original
  environment(patched) <- scope
  recovery_assert(identical(environment(original), original_environment),
                  "Original estimator environment was modified.")
  recovery_assert(identical(body(patched), body(original)),
                  "The estimator body was unexpectedly changed.")
  list(fit = patched, guard = scope$is_psd, state = state)
}

recover_additional_sparse <- function() {
  started <- proc.time()[["elapsed"]]
  root <- normalizePath(Sys.getenv("MRRR_REPO_ROOT", unset = getwd()),
                        winslash = "/", mustWork = TRUE)
  setwd(root)
  full <- file.path(root, "paper/output/additional_simulations/cluster_full_b300_delta01")
  manifest_file <- file.path(full, "point_tasks.csv")
  audit_dir <- file.path(full, "sparse_cholesky_recovery")
  recovered_dir <- file.path(full, "point_chunks_recovered")
  merged_dir <- file.path(full, "point_merged_recovered")
  final_file <- file.path(merged_dir, "additional_point_results.rds")
  recovery_assert(!dir.exists(recovered_dir) && !file.exists(final_file),
                  "Recovery outputs already exist; this script will not overwrite them.")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

  script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  recovery_assert(length(script_arg) == 1L, "Cannot identify this recovery script.")
  recovery_script <- normalizePath(sub("^--file=", "", script_arg),
                                  winslash = "/", mustWork = TRUE)
  load_script <- function(name) {
    environment <- new.env(parent = globalenv())
    sys.source(file.path(root, "paper/scripts", name), envir = environment)
    environment
  }
  core <- load_script("20_additional_simulation_core.R")
  worker <- load_script("22_run_additional_simulation_chunk.R")
  merger <- load_script("25_merge_additional_point_chunks.R")
  # Scope the RNG correction to this recovery process. Every existing merger
  # truth check uses this same core environment; the on-disk inputs stay intact.
  original_truth <- core$additional_make_generic_truth
  core$additional_make_generic_truth <- recovery_truth_function(original_truth)
  recovery_test_truth_rng(original_truth, core$additional_make_generic_truth)
  manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  validation <- merger$additional_merge_validate_manifest(manifest, root, core)
  recovery_assert(validation$total_replicates == 1000L && nrow(manifest) == 240L,
                  "This recovery expects the completed 1,000-replicate point manifest.")
  selected <- which(manifest$resource_class == "standard" & manifest$array_index == 9L)
  recovery_assert(length(selected) == 1L, "Point task standard/9 is not unique.")
  row <- manifest[selected, , drop = FALSE]
  recovery_assert(row$analysis == "rank_misspecification" && row$setting_index == 1L &&
                    row$replicate_start == 801L && row$replicate_end == 900L,
                  "Point task standard/9 differs from the reported failed chunk.")
  original_paths <- vapply(seq_len(nrow(manifest)), function(i) {
    merger$additional_merge_chunk_path(manifest[i, , drop = FALSE], root)
  }, character(1))
  original_hashes <- recovery_hash(original_paths)
  original_file <- original_paths[selected]
  original <- readRDS(original_file)
  keys <- paste0("r", 1:3, "__sparse_mr_rr")
  error_rows <- original$errors
  recovery_assert(is.data.frame(error_rows) && nrow(error_rows) == 3L &&
                    all(error_rows$replicate_id == 840L) &&
                    setequal(error_rows$result_key, keys) &&
                    !anyDuplicated(error_rows$result_key) &&
                    all(grepl("leading minor of order 9 is not positive", error_rows$message)),
                  "Observed errors differ from the three reported sparse Cholesky failures.")
  position <- match(840L, original$metadata$replicate_ids)
  recovery_assert(!is.na(position) && !isTRUE(original$metadata$complete),
                  "Unexpected original completion status or replicate indices.")
  for (key in names(original$estimates)) {
    keep <- if (key %in% keys) setdiff(seq_len(ncol(original$estimates[[key]])), position) else
      seq_len(ncol(original$estimates[[key]]))
    recovery_assert(all(is.finite(original$estimates[[key]][, keep, drop = FALSE])) &&
                      all(is.finite(original$biases[[key]][, keep, drop = FALSE])),
                    paste("Additional missing/nonfinite results found for", key))
    if (key %in% keys) {
      recovery_assert(all(is.na(original$estimates[[key]][, position])) &&
                        all(is.na(original$biases[[key]][, position])),
                      paste("The failed result is not an entirely missing column:", key))
    }
  }

  initial_inputs <- file.path(root, original$input_manifest$repository_path)
  recovery_assert(identical(unname(recovery_hash(initial_inputs)),
                            tolower(as.character(original$input_manifest$md5))),
                  "Original chunk inputs changed after production.")
  bootstrap_files <- file.path(full, "bootstrap_merged",
                               c("additional_bootstrap_results.rds", "additional_bootstrap_summary.csv"))
  protected_paths <- unique(c(initial_inputs, original_paths, manifest_file,
    bootstrap_files, recovery_script,
    file.path(root, "paper/scripts/15_smoke_test_from_scratch_simulations.R"),
    file.path(root, "paper/scripts/25_merge_additional_point_chunks.R")))
  protected_hashes <- recovery_hash(protected_paths)

  baseline <- core$additional_load_baseline_core(root)
  estimators <- core$additional_load_estimators(root)
  calibration <- core$additional_build_calibration(root, baseline, pz = 177L)
  config <- core$additional_read_config(file.path(root,
    "paper/config/rank_misspecification_settings.csv"), "rank_misspecification")
  config <- config[config$setting_index == 1L, , drop = FALSE]
  specification <- worker$additional_chunk_result_specification(
    "rank_misspecification", "standard", config)
  config <- specification$config
  recovery_assert(identical(config, original$configuration) &&
                    identical(specification$specification, original$result_specification),
                  "Recovered configuration or result specification differs from the original.")
  truth <- core$additional_make_generic_truth(0, calibration$px, calibration$py, seed = 123L)
  merger$additional_merge_assert_close(truth$C, original$truth$C, "Truth reconstruction", 1e-12)
  simulated <- core$additional_simulate_dataset(config[1L, , drop = FALSE], truth$C,
    840L, calibration, estimators, baseline)
  recovery_assert(identical(as.integer(simulated$replicate_seed),
                            as.integer(original$metadata$data_seeds[position])),
                  "Regenerated data seed differs from the original.")

  # Reproduce the original matrix expressions, including multiplication order.
  GAMMA_hat <- simulated$data$GAMMA_hat
  gamma_hat <- simulated$data$gamma_hat
  n <- nrow(gamma_hat)
  P_GAMMA <- GAMMA_hat %*% solve(t(GAMMA_hat) %*% GAMMA_hat) %*% t(GAMMA_hat)
  Sigma_gammahat <- t(gamma_hat) %*% gamma_hat / n
  part1 <- Sigma_gammahat - t(P_GAMMA %*% gamma_hat) %*% (P_GAMMA %*% gamma_hat) / n
  part2 <- part1 - simulated$parameters$Sigma_X
  diagnostics <- rbind(
    recovery_matrix_diagnostic(part1, "matrix_part1", estimators$is_psd),
    recovery_matrix_diagnostic(part2, "matrix_part2_before_projection", estimators$is_psd))
  write.csv(diagnostics, file.path(audit_dir, "matrix_diagnostics.csv"), row.names = FALSE)
  print(diagnostics, row.names = FALSE, digits = 16)
  recovery_assert(diagnostics$chol_success[1L] && diagnostics$old_is_psd[2L] &&
                    !diagnostics$chol_success[2L] && diagnostics$min_eigenvalue[2L] >= -1e-8,
                  "Failure is outside the PSD-tolerance/Cholesky mismatch; recovery stopped.")
  projected <- estimators$.nearest_psd(part2, epsilon = 1e-6)
  recovery_assert(!inherits(recovery_chol(projected), "error"),
                  "The existing projection does not yield a valid Cholesky factor.")
  diagnostics <- rbind(diagnostics,
    recovery_matrix_diagnostic(projected, "matrix_part2_after_existing_projection", estimators$is_psd))
  write.csv(diagnostics, file.path(audit_dir, "matrix_diagnostics.csv"), row.names = FALSE)

  patch <- recovery_make_sparse(estimators)
  patched_estimators <- new.env(parent = emptyenv())
  patched_estimators$mr_rr_sparse <- patch$fit
  recovered <- original
  reproduction_differences <- numeric(3L)
  names(reproduction_differences) <- paste0("r", 1:3)
  fit_elapsed <- numeric(3L)
  names(fit_elapsed) <- keys
  for (rank in 1:3) {
    cfg <- config[config$working_rank == rank, , drop = FALSE]
    recovery_assert(nrow(cfg) == 1L, "Working rank configuration is not unique.")
    parameters <- core$additional_build_parameters(cfg, truth$C, calibration, estimators, baseline)
    regular_key <- paste0("r", rank, "__regularized_mr_rr")
    worker$additional_chunk_set_seed(original$metadata$estimator_seeds[[regular_key]][position])
    regular_fit <- worker$additional_chunk_fit_method("regularized_mr_rr", simulated$data,
      parameters, cfg$regularization_rate, cfg$sparse_lambda, estimators, 1000L, 200L)
    reference <- original$estimates[[regular_key]][, position]
    reproduction_differences[rank] <- max(abs(as.vector(regular_fit) - reference))
    recovery_assert(reproduction_differences[rank] <= 1e-10 * max(1, max(abs(reference))),
                    paste("Existing regularized fit could not be reproduced for rank", rank))

    key <- keys[rank]
    worker$additional_chunk_set_seed(original$metadata$estimator_seeds[[key]][position])
    original_error <- tryCatch(worker$additional_chunk_fit_method("sparse_mr_rr", simulated$data,
      parameters, cfg$regularization_rate, cfg$sparse_lambda, estimators, 1000L, 200L),
      error = function(e) e)
    recovery_assert(inherits(original_error, "error") &&
                      grepl("leading minor of order 9 is not positive", conditionMessage(original_error)),
                    paste("The original sparse failure did not reproduce for rank", rank))
    worker$additional_chunk_set_seed(original$metadata$estimator_seeds[[key]][position])
    before_fit <- proc.time()[["elapsed"]]
    fit <- worker$additional_chunk_fit_method("sparse_mr_rr", simulated$data,
      parameters, cfg$regularization_rate, cfg$sparse_lambda, patched_estimators, 1000L, 200L)
    estimate <- worker$additional_chunk_validate_fit(fit, "sparse_mr_rr", parameters, key, baseline)
    fit_elapsed[rank] <- proc.time()[["elapsed"]] - before_fit
    recovered$estimates[[key]][, position] <- as.vector(estimate)
    recovered$biases[[key]][, position] <- as.vector(estimate - truth$C)
    cat("Recovered replicate 840:", key, "PASS\n")
  }
  recovery_assert(patch$state$fallbacks == 3L, "Expected exactly three added projection fallbacks.")
  for (key in names(original$estimates)) {
    keep <- if (key %in% keys) setdiff(seq_len(ncol(original$estimates[[key]])), position) else
      seq_len(ncol(original$estimates[[key]]))
    for (component in c("estimates", "biases")) {
      recovery_assert(identical(original[[component]][[key]][, keep, drop = FALSE],
                                recovered[[component]][[key]][, keep, drop = FALSE]),
                      paste("A previously successful result changed:", component, key))
    }
  }

  audit <- list(
    schema_version = "1.1.0",
    recovery_version = "2",
    recovered_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    truth_rng_kind = c("Mersenne-Twister", "Inversion", "Rejection"),
    replicate_rng_kind = c("L'Ecuyer-CMRG", "Inversion", "Rejection"),
    truth_rng_policy = "Match fresh production-session truth and restore caller RNG kind and state.",
    policy = "Apply the existing eigenvalue floor 1e-6 also when is_psd passes but chol fails.",
    scope = "Same rule for sparse MR-rr; original successful fits follow the identical computation path.",
    initial_provenance_note = "input_manifest describes initial production; this record describes the recovery.",
    original_chunk = original_file,
    original_chunk_md5 = unname(original_hashes[selected]),
    recovery_script = recovery_script,
    recovery_script_md5 = unname(recovery_hash(recovery_script)),
    additional_source_hashes = protected_hashes[!names(protected_hashes) %in% original_paths],
    analysis = "rank_misspecification", setting_index = 1L,
    replicate_id = 840L, data_seed = simulated$replicate_seed,
    result_keys = keys,
    original_errors = error_rows,
    original_estimator_seeds = lapply(original$metadata$estimator_seeds[keys], `[`, position),
    matrix_diagnostics = diagnostics,
    matrices = list(matrix_part1 = part1, matrix_part2 = part2, projected_part2 = projected),
    projection_frobenius_change = norm(projected - part2, "F"),
    projection_relative_frobenius_change = norm(projected - part2, "F") / norm(part2, "F"),
    regularized_reproduction_max_differences = reproduction_differences,
    changed_fit_count = 3L, previously_successful_entries_identical = TRUE,
    guard_definition = deparse(patch$guard), fit_elapsed_seconds = fit_elapsed,
    session_info = capture.output(sessionInfo())
  )
  recovered$errors <- original$errors[0L, , drop = FALSE]
  recovered$successful_replicates <- vapply(recovered$estimates, function(x) {
    as.integer(sum(colSums(is.finite(x)) == nrow(x)))
  }, integer(1))
  recovered$metadata$complete <- TRUE
  recovered$metadata$numerical_recovery_applied <- TRUE
  recovered$numerical_recovery <- audit
  new_file <- file.path(recovered_dir, basename(original_file))
  merger$additional_merge_validate_chunk(recovered, new_file, root, core, worker, expected_row = row)
  recovery_assert(identical(recovery_hash(protected_paths), protected_hashes),
                  "An original input/result changed during recovery.")

  dir.create(recovered_dir, recursive = TRUE)
  others <- setdiff(seq_along(original_paths), selected)
  recovery_assert(all(file.copy(original_paths[others],
    file.path(recovered_dir, basename(original_paths[others])), overwrite = FALSE)),
    "Could not copy the previously successful chunks.")
  copied_hashes <- recovery_hash(file.path(recovered_dir, basename(original_paths[others])))
  recovery_assert(identical(unname(copied_hashes), unname(original_hashes[others])),
                  "A copied successful chunk differs from its original.")
  worker$additional_chunk_save_rds(recovered, new_file)
  saveRDS(audit, file.path(audit_dir, "recovery_audit.rds"))
  write.csv(error_rows, file.path(audit_dir, "original_failed_fits.csv"), row.names = FALSE)

  # Extend output provenance while retaining EVERY existing merger validation.
  original_save <- merger$additional_merge_save_rds
  merger$additional_merge_save_rds <- function(object, path, overwrite) {
    recovery_assert(basename(path) == "additional_point_results.rds",
                    "Unexpected merger output path.")
    object$metadata$numerical_recovery_applied <- TRUE
    object$metadata$recovered_fit_count <- 3L
    object$numerical_recovery <- audit
    for (i in seq_along(object$groups)) {
      group <- object$groups[[i]]
      if (group$analysis == "rank_misspecification" && group$setting_index == 1L) {
        object$groups[[i]]$numerical_recovery <- audit
      }
    }
    original_save(object, path, overwrite)
  }
  merger$additional_merge_full(
    arguments = list(manifest = manifest_file, `chunks-dir` = recovered_dir,
                     `output-dir` = merged_dir),
    repo_root = root, core = core, worker = worker,
    allow_incomplete = FALSE, overwrite = FALSE)
  final <- readRDS(final_file)
  recovery_assert(isTRUE(final$metadata$numerical_recovery_applied) &&
                    final$metadata$recovered_fit_count == 3L &&
                    identical(final$numerical_recovery, audit),
                  "Merged output lost the numerical recovery provenance.")
  recovery_assert(identical(recovery_hash(protected_paths), protected_hashes),
                  "An original input/result changed during the merge.")
  writeLines(c("RECOVERY AND POINT MERGE: PASS", paste("Merged output:", final_file)),
             file.path(audit_dir, "SUCCESS.txt"))
  cat("\nOriginal chunks, frozen inputs, and bootstrap results unchanged: PASS\n")
  cat("Previously successful estimates unchanged: PASS\n")
  cat("Recovered fits: 3 (replicate 840, sparse ranks 1/2/3)\n")
  cat("Recovery provenance retained in merged RDS: PASS\n")
  cat("Total elapsed seconds:", proc.time()[["elapsed"]] - started, "\n")
  cat("RECOVERY AND POINT MERGE: PASS\n")
}

if (sys.nframe() == 0L) recover_additional_sparse()
