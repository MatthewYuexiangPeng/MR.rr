#!/usr/bin/env Rscript
# Targeted core validation. Default mode includes actual CVXR sparse fits.
# --base-only is a partial check for environments without a native CVXR solver.

paper_validate_spectral_engine <- function(repo_root = getwd(), output_dir = NULL,
                                           base_only = FALSE) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  if (is.null(output_dir)) output_dir <- file.path(root, "paper/output/spectral_rebuild/core_validation")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  api <- new.env(parent = baseenv())
  sys.source(file.path(root, "paper/lib/paper_engine.R"), envir = api)
  engine <- api$paper_engine_load(root)
  e <- engine$estimators
  rows <- list()
  check <- function(label, expression) {
    started <- proc.time()[[3L]]
    error <- tryCatch({ force(expression); NULL }, error = function(err) conditionMessage(err))
    passed <- is.null(error)
    rows[[length(rows) + 1L]] <<- data.frame(
      check = label, status = if (passed) "PASS" else "FAIL",
      elapsed_seconds = proc.time()[[3L]] - started,
      detail = if (passed) "" else error, stringsAsFactors = FALSE
    )
    cat(if (passed) "PASS: " else "FAIL: ", label,
        if (!passed) paste0(" - ", error) else "", "\n", sep = "")
    invisible(passed)
  }
  close_to <- function(actual, expected, tolerance = 1e-10) {
    stopifnot(identical(dim(actual), dim(expected)), all(is.finite(actual)),
              max(abs(actual - expected)) <= tolerance * max(1, max(abs(expected))))
  }
  throws <- function(expression) {
    stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
  }
  original_surrogate <- function(Y, X, Sigma_X) {
    n <- nrow(X)
    P <- Y %*% solve(t(Y) %*% Y) %*% t(Y)
    first <- t(X) %*% X / n - t(P %*% X) %*% (P %*% X) / n
    second <- first - Sigma_X
    if (!e$.is_psd(second)) second <- e$.nearest_psd(second, 1e-6)
    P %*% X + (diag(1, n) - P) %*% X %*% solve(chol(first)) %*% chol(second)
  }
  cat("MR.rr spectral paper core validation\nR:", R.version.string,
      "\nMode:", if (base_only) "base-only (partial)" else "full core (with CVXR)", "\n")
  api$paper_with_seed(9062026L, {
    X <- matrix(stats::rnorm(280), 70, 4)
    Y <- X %*% matrix(c(.5, .1, .2, 0, 0, -.3, .2, .4), 4, 2) +
      matrix(stats::rnorm(140), 70, 2) * .3
    Sigma_X <- diag(rep(1e-3, 4))
    W <- matrix(c(2, .2, .2, 1.5), 2, 2)
    check("Current source entry and source fingerprints", {
      stopifnot(identical(parent.env(e), baseenv()), environmentIsLocked(e))
      api$paper_engine_check_sources(engine)
    })
    check("Spectral/legacy equivalence, weighted rank 1 and 2", {
      for (r in 1:2) {
        fit <- api$paper_engine_fit(engine, "regularized", "simulation", Y, X, r,
                                    Sigma_X, W, regularization_rate = .01)
        old <- e$mr_rr_regularized(Y, X, r, Sigma_X, .01, W, implementation = "legacy")
        close_to(fit$AB, old$AB)
        close_to(crossprod(fit$A, W %*% fit$A), diag(r))
        stopifnot(identical(fit$paper$regularized_implementation, "spectral"))
      }
    })
    check("Singular/indefinite spectral inverse retains signed eigenvalues", {
      S <- diag(c(1, 0, -1))
      throws(solve(S))
      close_to(e$.regularized_inverse(S, .1), diag(c(1 / 1.1, 0, -1 / 1.1)))
    })
    check("Zero phi matches corrected MR-rr for nonsingular input", {
      fit <- api$paper_engine_fit(engine, "regularized", "simulation", Y, X, 2,
                                  Sigma_X, W, regularization_rate = 0)
      direct <- e$mr_rr(Y, X, 2, Sigma_X, W)
      close_to(fit$AB, direct$AB)
    })
    check("Unspecified tuning and incompatible support are rejected", {
      throws(api$paper_engine_fit(engine, "regularized", "simulation", Y, X, 1, Sigma_X))
      throws(api$paper_engine_fit(engine, "sparse", "simulation", Y, X, 1, Sigma_X))
      throws(api$paper_engine_fit(engine, "sparse", "real_bootstrap", Y, X, 1, Sigma_X))
      throws(api$paper_engine_fit(engine, "naive", "simulation", Y, X, 1,
                                  sparse_support = matrix(TRUE, 1, 4)))
    })
    check("Successful sparse Cholesky and surrogate are unchanged", {
      M <- matrix(c(2, .1, .1, 1), 2, 2)
      corrected <- e$.sparse_corrected_chol(M)
      stopifnot(identical(corrected$factor, chol(M)), identical(corrected$covariance, M),
                !corrected$diagnostics$corrected_covariance_projected)
      close_to(e$.construct_gamma_tilde(Y, X, Sigma_X), original_surrogate(Y, X, Sigma_X))
    })
    check("PSD-pass/Cholesky-fail uses the recorded 1e-6 projection", {
      M <- diag(c(2e-4, 1e-4, -3.990090310468490e-9))
      stopifnot(e$.is_psd(M))
      throws(chol(M))
      corrected <- e$.sparse_corrected_chol(M)
      stopifnot(corrected$diagnostics$projection_reason == "chol_failed",
                corrected$diagnostics$eigenvalue_floor == 1e-6)
      close_to(crossprod(corrected$factor), diag(c(2e-4, 1e-4, 1e-6)), 1e-12)
      singular <- e$.sparse_corrected_chol(diag(c(1, 0)))
      stopifnot(singular$diagnostics$projection_reason == "chol_failed")
    })
    check("Clearly indefinite covariance retains the historical projection", {
      M <- diag(c(2e-4, -2e-5))
      corrected <- e$.sparse_corrected_chol(M)
      stopifnot(corrected$diagnostics$projection_reason == "not_psd")
      close_to(crossprod(corrected$factor), diag(c(2e-4, 1e-6)), 1e-12)
    })
    # An exact construction: residual X'X/n = 0.5 I; the corrected fourth
    # eigenvalue is slightly negative and lies inside the historical PSD tolerance.
    bad_Y <- diag(8)[, 1:2, drop = FALSE]
    bad_X <- matrix(0, 8, 4)
    bad_X[3:6, ] <- 2 * diag(4)
    bad_X[1:2, ] <- matrix(c(.2, .1, .1, .3, .2, -.1, .4, .2), 2, 4)
    bad_Sigma <- diag(c(.49, .48, .47, .5 + 4e-9))
    check("Recovered surrogate satisfies its second-moment identity", {
      throws(original_surrogate(bad_Y, bad_X, bad_Sigma))
      recovered <- e$.construct_gamma_tilde(bad_Y, bad_X, bad_Sigma, diagnostics = TRUE)
      P <- bad_Y %*% t(bad_Y)
      target <- crossprod(P %*% bad_X) / 8 + diag(c(.01, .02, .03, 1e-6))
      close_to(crossprod(recovered$value) / 8, target, 1e-12)
      stopifnot(recovered$diagnostics$projection_reason == "chol_failed")
    })
    check("Uncorrected singular surrogate is not silently repaired", {
      throws(e$.construct_gamma_tilde(bad_Y, matrix(0, 8, 4), diag(4)))
    })
    check("Fixed-support real bootstrap preserves rank and exact zeros", {
      init <- e$mr_rr(Y, X, 1, Sigma_X, W)
      support <- matrix(c(TRUE, FALSE, TRUE, FALSE), 1, 4)
      fit <- api$paper_engine_fit(engine, "sparse", "real_bootstrap", Y, X, 1,
        Sigma_X, W, sparse_support = support, sparse_A_init = init$A, sparse_B_init = init$B)
      direct <- e$mr_rr_sparse_refit(Y, X, Sigma_X, support, W, init$A, init$B)
      close_to(fit$AB, direct$AB)
      stopifnot(fit$paper$sparse_refitted, identical(fit$support, support),
                all(fit$B[!support] == 0))
      throws(api$paper_engine_fit(engine, "sparse", "real_bootstrap", Y, X, 2,
        Sigma_X, W, sparse_support = support, sparse_A_init = init$A, sparse_B_init = init$B))
    })
    check("Scoped RNG restores kind/state on success and failure", {
      before <- .Random.seed
      kind <- RNGkind()
      a <- api$paper_with_seed(123L, stats::rnorm(12), kind = "Mersenne-Twister")
      stopifnot(identical(before, .Random.seed), identical(kind, RNGkind()))
      b <- api$paper_with_seed(123L, stats::rnorm(12), kind = "Mersenne-Twister")
      stopifnot(identical(a, b))
      throws(api$paper_with_seed(99L, stop("intentional")))
      stopifnot(identical(before, .Random.seed), identical(kind, RNGkind()))
      rm(".Random.seed", envir = .GlobalEnv)
      api$paper_with_seed(55L, stats::runif(3))
      stopifnot(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      assign(".Random.seed", before, envir = .GlobalEnv)
    })

    if (!base_only) {
      available <- check("Native CVXR and OSQP available", {
        stopifnot(requireNamespace("CVXR", quietly = TRUE))
        stopifnot("OSQP" %in% CVXR::installed_solvers())
        cat("CVXR version:", as.character(utils::packageVersion("CVXR")), "\n")
      })
      if (available) {
        check("Simulation/tuning sparse fits do not perform post-selection refit", {
          selected <- e$mr_rr_sparse(Y, X, 1, Sigma_X, lambda = .05, W = W)
          for (stage in c("simulation", "tuning")) {
            fit <- api$paper_engine_fit(engine, "sparse", stage, Y, X, 1, Sigma_X, W,
                                        sparse_lambda = .05)
            close_to(fit$AB, selected$AB, 1e-8)
            stopifnot(!fit$paper$sparse_refitted, is.null(fit$bootstrap_state))
          }
        })
        check("Real point selection/refit supplies a reusable bootstrap state", {
          fit <- api$paper_engine_fit(engine, "sparse", "real_point", Y, X, 1,
                                      Sigma_X, W, sparse_lambda = .05)
          state <- fit$bootstrap_state
          stopifnot(fit$paper$sparse_refitted,
                    identical(state$support, fit$sparse_selection$B != 0),
                    identical(state$A_init, fit$sparse_selection$A))
          indices <- api$paper_with_seed(123L, sample.int(nrow(X), replace = TRUE))
          boot <- api$paper_engine_fit(engine, "sparse", "real_bootstrap",
            Y[indices, , drop = FALSE], X[indices, , drop = FALSE], 1, Sigma_X, W,
            sparse_support = state$support, sparse_A_init = state$A_init,
            sparse_B_init = state$B_init)
          stopifnot(identical(boot$support, state$support), all(boot$B[!state$support] == 0))
        })
        check("Sparse solver handles the PSD/Cholesky regression fixture", {
          fit <- api$paper_engine_fit(engine, "sparse", "simulation", bad_Y, bad_X, 1,
                                      bad_Sigma, sparse_lambda = .001)
          stopifnot(all(is.finite(fit$AB)),
                    fit$numerical_diagnostics$projection_reason == "chol_failed")
        })
      }
    }
  }, kind = "Mersenne-Twister")

  table <- do.call(rbind, rows)
  utils::write.csv(table, file.path(output_dir, "checks.csv"), row.names = FALSE)
  utils::write.csv(engine$source_manifest, file.path(output_dir, "engine_sources.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "session_info.txt"))
  success <- all(table$status == "PASS")
  status <- if (!success) "FAIL" else if (base_only) "PARTIAL_PASS" else "PASS"
  writeLines(c(paste("Status:", status), paste("Checks:", nrow(table)),
    "No Monte Carlo production jobs were run by this check.",
    if (base_only) "CVXR sparse selection was not tested in base-only mode." else
      "Full core mode: see checks.csv for the actual executed CVXR/OSQP checks."),
    file.path(output_dir, "STATUS.txt"))
  cat("SPECTRAL CORE VALIDATION:", status, "\nReport:", normalizePath(output_dir, winslash = "/"), "\n")
  if (!success) stop("Spectral engine validation failed; see checks.csv.", call. = FALSE)
  invisible(table)
}

if (sys.nframe() == 0L) {
  arguments <- commandArgs(trailingOnly = TRUE)
  allowed <- grepl("^--(repo|output-dir)=.+$", arguments) | arguments == "--base-only"
  if (any(!allowed)) stop("Unknown argument: ", paste(arguments[!allowed], collapse = ", "))
  value <- function(key, default = NULL) {
    selected <- arguments[startsWith(arguments, paste0("--", key, "="))]
    if (length(selected) > 1L) stop("Duplicate argument: ", key)
    if (length(selected)) sub(paste0("^--", key, "="), "", selected) else default
  }
  paper_validate_spectral_engine(value("repo", getwd()), value("output-dir"),
                                 base_only = "--base-only" %in% arguments)
}
