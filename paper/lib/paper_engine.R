# Shared MR-rr computation entry for the spectral paper rebuild.
# Sourcing this file does not load frozen code, fit models, or change RNG state.

paper_engine_load <- function(repo_root = getwd()) {
  root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  paths <- c(
    "R/input-validation.R", "R/matrix-utils.R", "R/cvxr-utils.R",
    "R/mr-rr.R", "R/mr-rr-naive.R", "R/mr-rr-regularized.R",
    "R/mr-rr-sparse.R", "R/mr-rr-sparse-refit.R"
  )
  tracked <- c("DESCRIPTION", "paper/lib/paper_engine.R", paths)
  absolute <- file.path(root, tracked)
  if (any(!file.exists(absolute))) {
    stop("Missing engine source: ", paste(tracked[!file.exists(absolute)], collapse = ", "),
         call. = FALSE)
  }
  if (read.dcf(file.path(root, "DESCRIPTION"), fields = "Package")[1L] != "MR.rr") {
    stop("repo_root must contain the MR.rr package.", call. = FALSE)
  }
  # A base-environment parent avoids resolving estimators from .GlobalEnv or a
  # previously installed copy of MR.rr. The submitted R/ sources are the core.
  estimators <- new.env(parent = baseenv())
  for (path in paths) sys.source(file.path(root, path), envir = estimators)
  lockEnvironment(estimators, bindings = TRUE)
  structure(list(
    api_version = "spectral-paper-core-1",
    root = root,
    estimators = estimators,
    source_manifest = data.frame(
      path = tracked, md5 = unname(tools::md5sum(absolute)),
      stringsAsFactors = FALSE
    )
  ), class = "mrrr_paper_engine")
}


paper_engine_check_sources <- function(engine) {
  stopifnot(inherits(engine, "mrrr_paper_engine"))
  actual <- unname(tools::md5sum(file.path(engine$root, engine$source_manifest$path)))
  if (!identical(actual, engine$source_manifest$md5)) {
    stop("Engine source files changed after loading. Start a new versioned run.",
         call. = FALSE)
  }
  invisible(TRUE)
}


paper_with_seed <- function(seed, expression, kind = "L'Ecuyer-CMRG") {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != floor(seed)) {
    stop("seed must be one non-negative integer.", call. = FALSE)
  }
  if (!is.character(kind) || length(kind) != 1L || is.na(kind) ||
      !kind %in% c("L'Ecuyer-CMRG", "Mersenne-Twister")) {
    stop("Unsupported RNG kind for the paper pipeline.", call. = FALSE)
  }
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  RNGkind(kind, "Inversion", "Rejection")
  set.seed(as.integer(seed))
  force(expression)
}


paper_engine_fit <- function(
    engine, method, stage, Y, X, rank, Sigma_X = NULL, W = NULL,
    regularization_rate = NULL, sparse_lambda = NULL,
    sparse_threshold = 1e-2, sparse_max_iter = 100L, sparse_tol = 1e-2,
    sparse_solver = "OSQP", sparse_support = NULL,
    sparse_A_init = NULL, sparse_B_init = NULL) {
  stopifnot(inherits(engine, "mrrr_paper_engine"))
  method <- match.arg(method, c("naive", "corrected", "regularized", "sparse"))
  stage <- match.arg(stage, c("simulation", "real_point", "real_bootstrap", "tuning"))
  e <- engine$estimators
  inputs <- e$.validate_estimator_inputs(
    Y, X, rank, Sigma_X, W, require_sigma_x = method != "naive"
  )
  Y <- inputs$Y
  X <- inputs$X
  rank <- inputs$r
  Sigma_X <- inputs$Sigma_X
  W <- inputs$W
  if (method == "regularized" && is.null(regularization_rate)) {
    stop("Supply regularization_rate from the analysis configuration or tuning step.",
         call. = FALSE)
  }
  has_support_inputs <- !is.null(sparse_support) || !is.null(sparse_A_init) ||
    !is.null(sparse_B_init)
  if (has_support_inputs && !(method == "sparse" && stage == "real_bootstrap")) {
    stop("Fixed support/initial values are only accepted for sparse real_bootstrap.",
         call. = FALSE)
  }
  if (method == "sparse") {
    if (stage == "real_bootstrap") {
      if (is.null(sparse_support) || is.null(sparse_A_init) || is.null(sparse_B_init)) {
        stop("real_bootstrap requires original-data sparse_support, sparse_A_init and sparse_B_init.",
             call. = FALSE)
      }
      if (!is.matrix(sparse_support) || nrow(sparse_support) != rank) {
        stop("Fixed sparse support must have exactly rank rows.", call. = FALSE)
      }
    } else if (is.null(sparse_lambda)) {
      stop("Supply sparse_lambda from the analysis configuration.", call. = FALSE)
    }
  }

  warnings <- character()
  selection <- NULL
  fit <- withCallingHandlers({
    switch(method,
      naive = e$mr_rr_naive(Y, X, r = rank, W = W),
      corrected = e$mr_rr(Y, X, r = rank, Sigma_X = Sigma_X, W = W),
      regularized = e$mr_rr_regularized(
        Y, X, r = rank, Sigma_X = Sigma_X, W = W,
        regularization_rate = regularization_rate, implementation = "spectral"
      ),
      sparse = {
        if (stage != "real_bootstrap") {
          selection <- e$mr_rr_sparse(
            Y, X, r = rank, Sigma_X = Sigma_X, W = W,
            lambda = sparse_lambda, threshold = sparse_threshold,
            max_iter = sparse_max_iter, tol = sparse_tol, solver = sparse_solver
          )
        }
        if (stage %in% c("real_point", "real_bootstrap")) {
          if (stage == "real_point") {
            sparse_support <- selection$B != 0
            sparse_A_init <- selection$A
            sparse_B_init <- selection$B
          }
          e$mr_rr_sparse_refit(
            Y, X, Sigma_X = Sigma_X, W = W, support = sparse_support,
            A_init = sparse_A_init, B_init = sparse_B_init,
            max_iter = sparse_max_iter, tol = sparse_tol
          )
        } else {
          selection
        }
      }
    )
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    # Warnings remain visible in the job log and are also returned to the worker.
  })

  expected_dimensions <- list(A = c(ncol(Y), rank), B = c(rank, ncol(X)),
                              AB = c(ncol(Y), ncol(X)))
  for (name in names(expected_dimensions)) {
    value <- fit[[name]]
    if (!is.matrix(value) || !is.numeric(value) || is.complex(value) ||
        !identical(dim(value), expected_dimensions[[name]]) || any(!is.finite(value))) {
      stop("Invalid estimator output: ", name, call. = FALSE)
    }
  }
  scale <- max(1, max(abs(fit$AB)))
  if (max(abs(fit$AB - fit$A %*% fit$B)) > 1e-10 * scale) {
    stop("Estimator returned inconsistent A, B and AB.", call. = FALSE)
  }
  fit$paper <- list(
    engine_version = engine$api_version, method = method, stage = stage,
    working_rank = rank, warnings = warnings,
    regularization_rate = if (method == "regularized") regularization_rate else NULL,
    regularized_implementation = if (method == "regularized") "spectral" else NULL,
    sparse_refitted = method == "sparse" && stage %in% c("real_point", "real_bootstrap"),
    sparse_lambda = if (method == "sparse" && stage != "real_bootstrap") sparse_lambda else NULL,
    sparse_threshold = if (method == "sparse" && stage != "real_bootstrap") sparse_threshold else NULL,
    sparse_max_iter = if (method == "sparse") sparse_max_iter else NULL,
    sparse_tol = if (method == "sparse") sparse_tol else NULL,
    sparse_solver = if (!is.null(selection)) sparse_solver else NULL
  )
  if (method == "sparse" && stage == "real_point") {
    fit$sparse_selection <- selection
    # Preserve selection initial values as used in the manuscript's bootstrap.
    fit$bootstrap_state <- list(
      support = sparse_support, A_init = sparse_A_init, B_init = sparse_B_init
    )
  }
  fit
}
