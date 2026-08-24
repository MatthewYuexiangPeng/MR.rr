# ============================================================ #
# bootstrap_core.R
# Function definitions only — no data loading, no execution.
# Loaded by run_bootstrap_one.R
# ============================================================ #

suppressPackageStartupMessages({
  library(parallel)
  library(mr.divw)
})

# ============================================================ #
# Universal setting lists (must match main simulation)
# ============================================================ #
me_weight_list     <- c("2.5", "1")
effect_weight_list <- c("0.25", "1")

# from main simulation (matches setting 1: simulate_result_pred_260712_177.RData or setting 2: simulate_result_pred_260714_sparseC.RData)
# sim setting 1: regular C
# regularization_rate_list <- c(1.007845e-10,   # idx 1: me=2.5, effect=0.25
#                               2.340371e-12,   # idx 2: me=2.5, effect=1
#                               1.578970e-12,   # idx 3: me=1,   effect=0.25
#                               3.103420e-15)   # idx 4: me=1,   effect=1

# sim setting 2: C with sparse B to be updated
regularization_rate_list <- c(9.109067e-11,   # idx 1: me=2.5, effect=0.25
                              2.340371e-12,   # idx 2: me=2.5, effect=1
                              1.578970e-12,   # idx 3: me=1,   effect=0.25
                              3.105502e-15)   # idx 4: me=1,   effect=1


# ============================================================ #
# Helpers
# ============================================================ #
.get_sim_index <- function(me_weight_index, effect_weight_index, len_me = 2) {
  (me_weight_index - 1) * len_me + effect_weight_index
}

.simulation_data_only <- function(parameters) {
  n <- 177
  py <- parameters$py;  px <- parameters$px
  VX_tilde <- parameters$VX_tilde
  Sigma_X  <- parameters$Sigma_X
  Sigma_Y  <- parameters$Sigma_Y
  C        <- parameters$C
  
  gamma_j_star <- MASS::mvrnorm(n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star <- gamma_j_star %*% t(C)
  
  x_j_hat <- matrix(0, n, px)
  y_j_hat <- matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j, ] <- MASS::mvrnorm(1, gamma_j_star[j, ], Sigma_X, tol = 100)
    y_j_hat[j, ] <- MASS::mvrnorm(1, Gamma_j_star[j, ], Sigma_Y, tol = 100)
  }
  list(x_j_hat = x_j_hat, y_j_hat = y_j_hat)
}

# ============================================================ #
# Main: no-MrDAG parallel bootstrap
# ============================================================ #
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = NULL,
                                       bootstrap_size = 500,
                                       iteration = 1000, r_rank = 2, n_cores = NULL,
                                       estimator_set = "mr_only",
                                       save_estimates = TRUE) {
  
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  
  me_index     <- match(as.character(me_weight), me_weight_list)
  effect_index <- match(as.character(effect_weight), effect_weight_list)
  param_index  <- .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters   <- parameters_list[[param_index]]
  
  if (is.null(regularization_rate)) {
    regularization_rate <- regularization_rate_list[param_index]
  }
  
  C           <- parameters$C
  py          <- parameters$py;   px <- parameters$px
  Sigma_X_hat <- parameters$Sigma_X
  Sigma_Y_hat <- parameters$Sigma_Y
  C_vec       <- as.vector(C)
  d           <- px * py
  pz          <- 177
  W_precomputed <- solve(Sigma_Y_hat)
  
  estimator_order <- switch(estimator_set,
                            "all"      = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG"),
                            "no_mrdag" = c("IVW", "adIVW", "Naive", "MR", "MR_r"),
                            "mr_only"  = c("Naive", "MR", "MR_r"),
                            stop("estimator_set must be one of: 'all', 'no_mrdag', 'mr_only'"))
  
  start_time <- Sys.time()
  cat(sprintf("\n[no_mrdag] me=%s effect=%s (idx %d) | regu=%.3e | iter=%d bt=%d cores=%d\n",
              me_weight, effect_weight, param_index, regularization_rate,
              iteration, bootstrap_size, n_cores))
  flush.console()
  
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  
  clusterExport(cl,
                varlist = c("parameters", "regularization_rate", "r_rank",
                            "estimator_order", "C_vec", "Sigma_X_hat", "Sigma_Y_hat",
                            "pz", "px", "py", "d", "bootstrap_size",
                            "estimator_set", "save_estimates", "W_precomputed",
                            ".simulation_data_only"),
                envir = environment())
  
  clusterEvalQ(cl, {
    .libPaths(c("/home/peng.1276/R/library", .libPaths()))
    library(MASS)
    library(mr.divw)
    source("scripts/MR_rr_estimators.R")
  })
  
  run_one_iter <- function(loop_id) {
    sim <- .simulation_data_only(parameters)
    x_j_hat <- sim$x_j_hat;  y_j_hat <- sim$y_j_hat
    
    result_mats <- lapply(estimator_order, function(i) matrix(NA_real_, bootstrap_size, d))
    names(result_mats) <- estimator_order
    
    for (bt in 1:bootstrap_size) {
      idx  <- sample.int(pz, pz, replace = TRUE)
      x_bt <- x_j_hat[idx, ];  y_bt <- y_j_hat[idx, ]
      
      if (estimator_set %in% c("all", "no_mrdag", "mr_only")) {
        result_mats$Naive[bt, ] <- as.vector(mr_rr_naive(y_bt, x_bt, r = r_rank, W = W_precomputed)$AB)
        result_mats$MR[bt, ]    <- as.vector(mr_rr(y_bt, x_bt, r = r_rank, W = W_precomputed, Sigma_X = Sigma_X_hat)$AB)
        result_mats$MR_r[bt, ]  <- as.vector(mr_rr_regularized(y_bt, x_bt, r = r_rank, W = W_precomputed,
                                                               Sigma_X = Sigma_X_hat,
                                                               regularization_rate = regularization_rate)$AB)
      }
      if (estimator_set %in% c("all", "no_mrdag")) {
        result_mats$IVW[bt, ]   <- as.vector(ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat))
        result_mats$adIVW[bt, ] <- as.vector(adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat))
      }
      if (estimator_set == "all") {
        result_mats$MrDAG[bt, ] <- as.vector(Mr_DAG(y_bt, x_bt))
      }
    }
    
    coverage_vec <- list();  ci_length_vec <- list()
    for (name in estimator_order) {
      ci <- apply(result_mats[[name]], 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
      coverage_vec[[name]]  <- as.numeric(C_vec >= ci[1, ] & C_vec <= ci[2, ])
      ci_length_vec[[name]] <- ci[2, ] - ci[1, ]
    }
    list(coverage = coverage_vec, ci_length = ci_length_vec,
         estimates = if (save_estimates) result_mats else NULL)
  }
  
  clusterExport(cl, varlist = "run_one_iter", envir = environment())
  
  have_pb <- requireNamespace("pbapply", quietly = TRUE)
  coverage_results <- if (have_pb) pbapply::pblapply(1:iteration, run_one_iter, cl = cl)
  else         parLapply(cl, 1:iteration, run_one_iter)
  
  stopCluster(cl);  on.exit()
  
  # ===== Aggregate =====
  coverage_list <- lapply(estimator_order, function(name) {
    mat <- sapply(coverage_results, function(r) as.numeric(r$coverage[[name]]))
    if (is.vector(mat)) mat <- matrix(mat, nrow = 1) else if (nrow(mat) != d) mat <- t(mat)
    mat
  }); names(coverage_list) <- estimator_order
  
  ci_length_list <- lapply(estimator_order, function(name) {
    mat <- sapply(coverage_results, function(r) as.numeric(r$ci_length[[name]]))
    if (is.vector(mat)) mat <- matrix(mat, nrow = 1) else if (nrow(mat) != d) mat <- t(mat)
    mat
  }); names(ci_length_list) <- estimator_order
  
  estimates_list <- NULL
  if (save_estimates) {
    estimates_list <- lapply(estimator_order, function(name)
      lapply(coverage_results, function(r) r$estimates[[name]]))
    names(estimates_list) <- estimator_order
  }
  
  med_coverage <- lapply(coverage_list, function(m) median(rowMeans(m, na.rm = TRUE), na.rm = TRUE))
  ci_length_summary <- lapply(ci_length_list, function(m) {
    ml <- rowMeans(m, na.rm = TRUE)
    list(median = median(ml, na.rm = TRUE), mean = mean(ml, na.rm = TRUE),
         sd = sd(ml, na.rm = TRUE), min = min(ml, na.rm = TRUE), max = max(ml, na.rm = TRUE))
  })
  
  total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("[no_mrdag] done in %.1f min\n", total / 60))
  for (n in estimator_order)
    cat(sprintf("  %s: coverage=%.4f, CI_median=%.4f\n",
                n, med_coverage[[n]], ci_length_summary[[n]]$median))
  flush.console()
  
  list(
    setting = list(me_weight = me_weight, effect_weight = effect_weight,
                   param_index = param_index, regularization_rate = regularization_rate,
                   r_rank = r_rank),
    med_coverage = med_coverage, ci_length_summary = ci_length_summary,
    coverage_matrix = coverage_list, ci_length_matrix = ci_length_list,
    estimates = estimates_list
  )
}

# ============================================================ #
# Main: MrDAG-only parallel bootstrap
# ============================================================ #
nonpara_bootstrap_parallel_mrdag <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                             iteration = 300, bootstrap_size = 100,
                                             r_rank = 2, n_cores = NULL,
                                             mrdag_niter = 500, mrdag_burnin = 100,
                                             save_estimates = TRUE) {
  
  if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
  
  me_index     <- match(as.character(me_weight), me_weight_list)
  effect_index <- match(as.character(effect_weight), effect_weight_list)
  param_index  <- .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters   <- parameters_list[[param_index]]
  
  C_vec <- as.vector(parameters$C)
  px    <- parameters$px;  py <- parameters$py
  d     <- px * py;  pz <- 177
  estimator_order <- c("MrDAG")
  
  start_time <- Sys.time()
  cat(sprintf("\n[mrdag_only] me=%s effect=%s (idx %d) | niter=%d burnin=%d | iter=%d bt=%d cores=%d\n",
              me_weight, effect_weight, param_index,
              mrdag_niter, mrdag_burnin, iteration, bootstrap_size, n_cores))
  flush.console()
  
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  
  clusterExport(cl,
                varlist = c("parameters", "C_vec", "px", "py", "d", "pz",
                            "bootstrap_size", "mrdag_niter", "mrdag_burnin",
                            "save_estimates", "estimator_order",
                            ".simulation_data_only", "Mr_DAG"),
                envir = environment())
  
  clusterEvalQ(cl, {
    .libPaths(c("/home/peng.1276/R/library", .libPaths()))
    library(MASS)
    library(MrDAG)
  })
  
  run_one_iter <- function(loop_id) {
    sim <- .simulation_data_only(parameters)
    x_j_hat <- sim$x_j_hat;  y_j_hat <- sim$y_j_hat
    
    draws_mat <- matrix(NA_real_, bootstrap_size, d)
    for (bt in 1:bootstrap_size) {
      idx  <- sample.int(pz, pz, replace = TRUE)
      x_bt <- x_j_hat[idx, , drop = FALSE]
      y_bt <- y_j_hat[idx, , drop = FALSE]
      est_mat <- tryCatch(
        Mr_DAG(y_bt, x_bt, niter = mrdag_niter, burnin = mrdag_burnin),
        error = function(e) matrix(NA_real_, py, px))
      draws_mat[bt, ] <- as.vector(est_mat)
    }
    
    ci <- apply(draws_mat, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    list(coverage  = list(MrDAG = as.numeric(C_vec >= ci[1, ] & C_vec <= ci[2, ])),
         ci_length = list(MrDAG = ci[2, ] - ci[1, ]),
         estimates = if (save_estimates) list(MrDAG = draws_mat) else NULL)
  }
  
  clusterExport(cl, "run_one_iter", envir = environment())
  
  have_pb <- requireNamespace("pbapply", quietly = TRUE)
  coverage_results <- if (have_pb) pbapply::pblapply(1:iteration, run_one_iter, cl = cl)
  else         parLapply(cl, 1:iteration, run_one_iter)
  
  stopCluster(cl);  on.exit()
  
  # Aggregate (same shape as no_mrdag version)
  coverage_list  <- list(MrDAG = { m <- sapply(coverage_results, function(r) r$coverage$MrDAG)
  if (is.vector(m)) matrix(m, nrow = 1) else if (nrow(m) != d) t(m) else m })
  ci_length_list <- list(MrDAG = { m <- sapply(coverage_results, function(r) r$ci_length$MrDAG)
  if (is.vector(m)) matrix(m, nrow = 1) else if (nrow(m) != d) t(m) else m })
  estimates_list <- if (save_estimates) list(MrDAG = lapply(coverage_results, function(r) r$estimates$MrDAG)) else NULL
  
  med_coverage <- list(MrDAG = median(rowMeans(coverage_list$MrDAG, na.rm = TRUE), na.rm = TRUE))
  ml <- rowMeans(ci_length_list$MrDAG, na.rm = TRUE)
  ci_length_summary <- list(MrDAG = list(
    median = median(ml, na.rm = TRUE), mean = mean(ml, na.rm = TRUE),
    sd = sd(ml, na.rm = TRUE), min = min(ml, na.rm = TRUE), max = max(ml, na.rm = TRUE)))
  
  total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("[mrdag_only] done in %.1f min | coverage=%.4f CI_median=%.4f\n",
              total / 60, med_coverage$MrDAG, ci_length_summary$MrDAG$median))
  flush.console()
  
  list(
    setting = list(me_weight = me_weight, effect_weight = effect_weight,
                   param_index = param_index, r_rank = r_rank,
                   mrdag_niter = mrdag_niter, mrdag_burnin = mrdag_burnin),
    med_coverage = med_coverage, ci_length_summary = ci_length_summary,
    coverage_matrix = coverage_list, ci_length_matrix = ci_length_list,
    estimates = estimates_list
  )
}