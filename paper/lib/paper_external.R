# External comparator adapters for the unified spectral paper run.

paper_external_require <- function(methods) {
  required <- character()
  if (any(methods %in% c("ivw", "srivw"))) required <- c(required, "mr.divw")
  if ("mrdag" %in% methods) required <- c(required, "MrDAG")
  if ("sparse_mr_rr" %in% methods) required <- c(required, "CVXR", "osqp")
  for (p in required) if (!requireNamespace(p, quietly = TRUE))
    stop("Required package is unavailable: ", p, call. = FALSE)
  expected <- list("mr.divw" = c("mvmr.ivw", "mvmr.divw"),
                   MrDAG = c("MrDAG", "get_causaleffects"))
  for (p in intersect(required, names(expected))) for (f in expected[[p]])
    if (!f %in% getNamespaceExports(p))
      stop("This run requires the audited ", p, "::", f,
           " interface. A different package API must be reviewed explicitly.", call. = FALSE)
  if ("CVXR" %in% required && !"OSQP" %in% CVXR::installed_solvers())
    stop("CVXR cannot find OSQP.", call. = FALSE)
  invisible(required)
}

paper_external_se <- function(Sigma, n) {
  stopifnot(is.matrix(Sigma), all(is.finite(Sigma)), all(diag(Sigma) > 0), n >= 1L)
  # Each trait has its own column; its calibrated SE is constant across SNPs.
  matrix(rep(sqrt(diag(Sigma)), each = n), nrow = n, ncol = nrow(Sigma))
}

paper_external_fit <- function(method, Y, X, Sigma_X, Sigma_Y, configuration) {
  method <- match.arg(method, c("ivw", "srivw", "mrdag"))
  if (method %in% c("ivw", "srivw")) {
    f <- getExportedValue("mr.divw", if (method == "ivw") "mvmr.ivw" else "mvmr.divw")
    sx <- paper_external_se(Sigma_X, nrow(X))
    sy <- paper_external_se(Sigma_Y, nrow(Y))
    result <- matrix(NA_real_, ncol(Y), ncol(X))
    for (j in seq_len(ncol(Y))) {
      args <- list(beta.exposure = X, se.exposure = sx,
        beta.outcome = as.vector(Y[, j]), se.outcome = sy[, j])
      if (method == "srivw") args["phi_cand"] <- list(NULL)
      fit <- do.call(f, args)
      if (length(fit$beta.hat) != ncol(X)) stop("Unexpected mr.divw coefficient shape.")
      result[j, ] <- fit$beta.hat
    }
    return(list(AB = result))
  }
  sampler <- getExportedValue("MrDAG", "MrDAG")
  effects <- getExportedValue("MrDAG", "get_causaleffects")
  py <- ncol(Y); px <- ncol(X)
  fit <- NULL
  # Keep package progress output out of the per-bootstrap log; warnings are
  # handled by the caller and retained in the result diagnostics.
  invisible(utils::capture.output({
    fit <- sampler(data = data.frame(Y, X),
      niter = as.integer(configuration$mrdag_niter),
      burnin = as.integer(configuration$mrdag_burnin), thin = 5L,
      MrDAGcheck = list(Y_idx = seq_len(py), X_idx = py + seq_len(px)), fileName = NULL)
    causal <- effects(fit, ord = c(py + seq_len(px), seq_len(py)))$causalEffects
  }))
  list(AB = t(as.matrix(causal)))
}
