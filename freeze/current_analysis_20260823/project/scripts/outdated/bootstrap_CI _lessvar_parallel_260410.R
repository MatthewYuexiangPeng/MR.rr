# ============================================================ #
# 1. loading environment and data ####
# ============================================================ #
rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
library(GRAPPLE)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap::pheatmap)
# library(readxl)
# library(glmnet)
library(foreach)
library(doParallel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(latex2exp)


set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
# setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
setwd("~/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

exp_name = read.csv('data/traits_1e-4.csv')$x
out_name = c("LAS", "CES", "SVS")
names(out_name) <- 0:2
names(exp_name) <- 0:8

lip_data_real = read.csv('data/dat_1e-4.csv')
lip_corr_real = read.csv('data/rho_mat_1e-4.csv')
n_Y = c(1241207,1245612,1241619)

source("scripts/MR_rr_estimators.R")

# # regress lip_data$gamma_out on lip_data$gamma_exp1 to lip_data$gamma_exp24
# mean(abs(lm(lip_data$gamma_out ~ ., data = lip_data[,paste0('gamma_exp',1:9)])$coefficients))
# # the scale is similar to the generated C in .simulation


# ============================================================ #
# 2. hidden functions ####
# ============================================================ #
.sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
}


.nuclear_norm <- function(A) {
  sum(svd(A, nu = 0, nv = 0)$d)
}


.get_index <- function(i, py) {
  return(c(i - ((i-1) %/% py) * py, (i-1) %/% py + 1))
}


.plot_entrybias = function (data) {
  ggplot2::ggplot()+
    ggplot2::geom_boxplot(ggplot2::aes(y = data))+
    ggplot2::labs(title = "Naive MR-rr", y = "entrywise bias")
}


.plot_entrybias_d = function (data) {
  ggplot2::ggplot()+
    ggplot2::geom_boxplot(ggplot2::aes(y = data))+
    ggplot2::labs(title = "MR-rr", y = "entrywise bias")
}


.plot_entrybias_d_r = function (data) {
  ggplot2::ggplot()+
    ggplot2::geom_boxplot(ggplot2::aes(y = data))+
    ggplot2::labs(title = "MR-rr with regularization", y = "entrywise bias")
}


# new ver
rank_test_M1 = function(W, y_j_hat, x_j_hat, Sigma_X, print = TRUE, min_rank = 1, bt_loop = 5000){
  # # test 250826
  # min_rank = 1
  # bt_loop = 100
  
  # use bootstrap to estimate var(\hat C_MRrr)
  # store C hats
  px = ncol(x_j_hat)
  py = ncol(y_j_hat)
  pz = nrow(y_j_hat)
  
  p_value = c()
  for (r in min_rank:(min(px,py)-1)) {
    C_hat_MRrr_vec = as.vector(mr_rr(y_j_hat, x_j_hat, r=r, W = W, Sigma_X = Sigma_X)$AB)
    
    C_hat_matrix = matrix(0, nrow = bt_loop, ncol = px * py)
    for (b in 1:bt_loop) {
      sample_idx = sample(1:nrow(y_j_hat), replace = TRUE)
      y_j_hat_bt = y_j_hat[sample_idx, ]
      x_j_hat_bt = x_j_hat[sample_idx, ]
      C_hat_MRrr_bt_vec = as.vector(mr_rr(y_j_hat_bt, x_j_hat_bt, r=r, W = W, Sigma_X = Sigma_X)$AB)
      C_hat_matrix[b, ] = C_hat_MRrr_bt_vec
    }
    
    # variance of C_hat
    var_C_hat = cov(C_hat_matrix)
    # add regularized on inverse of var_C_hat
    # C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% C_hat_MRrr_vec
    C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% C_hat_MRrr_vec
    
    # # consider centering c hat
    # mu_hat = colMeans(C_hat_matrix)
    # C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% (C_hat_MRrr_vec - mu_hat)
    
    C_tilde = matrix(C_tilde_vec, nrow = py, ncol = px)
    
    # singular value of C_tilde
    svd_C_tilde = svd(C_tilde)
    lambda_tilde = svd_C_tilde$d
    
    # M1
    k = min(px, py)
    M1 = sum(lambda_tilde[(r+1):k]^2)
    # compare to chi-square distribution (py-r,px-r)
    p_value = c(p_value, 1 - pchisq(M1, df = (py - r) * (px - r)))
  }
  
  # 第一个为False的位置
  res_idx <- which(p_value >= 0.05)[1]
  
  if (is.na(res_idx)) {
    if (print) {
      cat("No rank passes the test; all p-values are < 0.05.\n")
      cat("Defaulting to full rank:", min(px, py), "\n")
    }
    return(min(px, py))
  } else {
    if (print) {
      cat("The lowest rank satisfying the test is:", res_idx+min_rank-1, 
          "; p-value =", round(p_value[res_idx], 3), "\n")
    }
    return(res_idx+min_rank-1)
  }
}


# me_weight: measurement error (Sigma_X) weight; effect_weight: random effect (Sigma_gammagamma) weight
.get_parameters = function(C, me_weight, effect_weight, r_RR = 2){
  # ## test
  # me_weight =1
  # r_RR = 2
  # px=9
  # r_approx = 2
  
  # r_RR is the rank chose by user when performing RRR
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 1000 # pz can be changed to any number
  px = 9
  py = 3 # py can be changed to any number from 1 to 9
  gamma_j = as.matrix(lip_data_real[,paste0('gamma_exp',1:9)])
  var_Z = 2 * lip_data_real$ImpMAF * (1 - lip_data_real$ImpMAF) # each Z is sum of two alleles, so var(Z) = var(Z^1+Z^2) = 2*var(Z^1), where Z1, Z2 ~ Binomial(eaf.outcome)
  var_Z_raw = var_Z
  
  
  z_index = sample(1:length(var_Z), pz, replace = TRUE)
  var_Z = var_Z[z_index]
  sqrt_var_Z = sqrt(var_Z)
  
  # TODO: consider use the same index here?
  # lip_data_index_z <- sample(1:114, pz, replace = TRUE) # Question: Is this p_Z too large?
  gamma_j_sample <- gamma_j[z_index,]
  gamma_j_star_temp <- gamma_j_sample * sqrt_var_Z # Q 2025: suppose to be devided instead of multiply? no should be multi
  
  
  # Sigma_Xj & Sigma_X & VX_tilde ----------------------------------------------------------------
  sigma_gamma_j = as.matrix(lip_data_real[,paste0('se_exp',1:9)])
  sigma_gamma_j = sigma_gamma_j * sqrt(var_Z_raw)
  Corr_X = lip_corr_real[1:9,1:9]
  Corr_X = as.matrix(Corr_X)
  pz_lip_data = nrow(gamma_j)
  Sigma_Xj = lapply(1:pz_lip_data, function(j)
    diag(sigma_gamma_j[j,]) %*% Corr_X %*% diag(sigma_gamma_j[j,]))
  Sigma_Xj_sample = Sigma_Xj[z_index] # Q 2025: same index
  Sigma_X_temp = lapply(1:pz, function(j) Sigma_Xj_sample[[j]]) # Q 2025: changed from lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d <- array(unlist(Sigma_X_temp), dim = c(9, 9, length(Sigma_X_temp)))
  Sigma_X <- apply(array_3d, c(1, 2), mean)
  
  VX_tilde <- cov(gamma_j_star_temp) - Sigma_X
  
  # weight to adjust IV strength
  Sigma_X = me_weight * Sigma_X
  VX_tilde = effect_weight * VX_tilde
  
  
  # Sigma_Y & weight.matrix -----------------------------------------------------------
  # var_Y = sample((lip_data_real$se_out1^2 * var_Z_raw) * n_Y,py) # Q 2025: devided by var z?
  # Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))
  
  sigma_Gamma_j = as.matrix(lip_data_real[,paste0('se_out',2:(py+1))])
  sigma_Gamma_j = sigma_Gamma_j * sqrt(var_Z_raw)
  Corr_Y = lip_corr_real[11:13,11:13]
  Corr_Y = as.matrix(Corr_Y)
  
  Sigma_Yj = lapply(1:pz_lip_data, function(j)
    diag(sigma_Gamma_j[j,]) %*% Corr_Y %*% diag(sigma_Gamma_j[j,]))
  Sigma_Yj_sample = Sigma_Yj[z_index] # Q 2025: same index
  Sigma_Y_temp = lapply(1:pz, function(j) Sigma_Yj_sample[[j]]) # Q 2025: changed from lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d_Y <- array(unlist(Sigma_Y_temp), dim = c(3, 3, length(Sigma_Y_temp)))
  Sigma_Y <- apply(array_3d_Y, c(1, 2), mean)
  
  # Sigma_Y = me_weight * Sigma_Y # TODO: should not weight Y?
  weight.matrix = solve(Sigma_Y)
  
  
  # SigmaXX, SigmaYX, SigmaXY ----------------------------------------------------------------
  r =  min(px,py) # maximum rank of true C
  # U = matrix(rnorm(py*r), py, r)
  # V = matrix(rnorm(px*r), r, px)
  # # assume true C to be rank r=3 for now. C=C^(r)
  # # use modified THM 2.1 to define true C^(r), (it is equivalent to define it as A_d*B_d)
  # eignvalue_matrix = diag(c(sample(c(sqrt(0.3), sqrt(0.2), sqrt(0.2)), r_approx, replace = TRUE),
  #                           rep(0, r-r_approx))) # set at 0.001 if want C to be rank > r
  # C <- U %*% eignvalue_matrix %*% V
  
  # SigmaXX, SigmaYX and SigmaXY
  SigmaXX <- Sigma_X + VX_tilde
  SigmaYX <- C %*% VX_tilde
  # SigmaXX <- Sigma_X_upscale + VX_tilde_upscale
  # SigmaYX <- C %*% VX_tilde_upscale
  SigmaXY <- t(SigmaYX)
  sqrt_Gamma <- .sqrt_matrix(weight.matrix)
  sqrt_Gamma_inv <- solve(sqrt_Gamma)
  
  # mr_rr_naive population level estimator
  M = sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX) %*% SigmaXY %*% sqrt_Gamma
  V = eigen(M)$vec[, 1:r_RR, drop = FALSE]
  A = sqrt_Gamma_inv %*% V
  B = t(V) %*% sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX)
  
  # debiased RRR population level estimator
  M_d = sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde) %*% SigmaXY %*% sqrt_Gamma
  # M_d = sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde_upscale) %*% SigmaXY %*% sqrt_Gamma
  V_d = eigen(M_d)$vec[, 1:r_RR, drop = FALSE]
  V_d = matrix(as.numeric(V_d),py,r_RR)
  A_d = sqrt_Gamma_inv %*% V_d
  B_d = t(V_d) %*% sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde)
  # B_d = t(V_d) %*% sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde_upscale)
  
  
  # VY_tilde ----------------------------------------------------------------
  VY_tilde = C %*% VX_tilde %*% t(C)
  # VY_tilde_upscale = C %*% VX_tilde_upscale %*% t(C)
  
  parameters = list(py = py, px = px, var_Z = var_Z, VX_tilde = VX_tilde,
                    Sigma_X = Sigma_X, Sigma_Y = Sigma_Y,
                    weight.matrix = weight.matrix, SigmaXX = SigmaXX,
                    SigmaYX = SigmaYX, SigmaXY = SigmaXY, C = C, A = A, B = B,
                    A_d = A_d, B_d = B_d, C_r = A_d %*% B_d, r_RR = r_RR,
                    VY_tilde = VY_tilde
                    # ,true.A_sparse = true.A_sparse, true.B_sparse = true.B_sparse
  )
  return(parameters)
}


.simulation = function(parameters, regularization_rate = 1e-13) {
  n = 1000
  # n: number of samples (pz)
  # get parameters
  py = parameters$py
  px = parameters$px
  var_Z = parameters$var_Z
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  SigmaXX = parameters$SigmaXX
  SigmaYX = parameters$SigmaYX
  SigmaXY = parameters$SigmaXY
  C = parameters$C
  r_RR = parameters$r_RR
  VY_tilde = parameters$VY_tilde
  W = parameters$weight.matrix
  W_sqrt = .sqrt_matrix(parameters$weight.matrix)
  
  # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
  gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  # sample x_j_hat, y_j_hat
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }
  
  # if we detect correct rank
  # r_tested <- rank_test(
  #   W = W,
  #   y_j_hat = y_j_hat,
  #   x_j_hat = x_j_hat,
  #   Sigma_X = Sigma_X,
  #   print = FALSE,
  #   bt_loop = 1000,
  #   alpha = 0.05,
  #   seed = 123,
  #   return_details = FALSE
  # )
  
  # M0 or M1 
  r_tested <- rank_test_M1(W, y_j_hat, x_j_hat, Sigma_X, print = FALSE, min_rank = 1)
  
  if (r_RR < r_tested) {
    rank_correct <- "small"
  } else if (r_RR > r_tested) {
    rank_correct <- "large"
  } else {
    rank_correct <- "eq"
  }
  
  # compute A_hat, B_hat
  result <- mr_rr_naive(y_j_hat, x_j_hat, r=r_RR, W=W) # TODO: changed here need check the result
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB # sample level estimator A_hat * B_hat
  
  result_d <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  result_d_r <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)
  
  result_ivw <- ivw_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)
  result_adivw <- adivw_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)
  
  # result_GRAPPLE <- GRAPPLE_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)
  result_MrDAG <- Mr_DAG(y_j_hat, x_j_hat)
  
  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB
  
  A_d_hat_r = result_d_r$A
  B_d_hat_r = result_d_r$B
  AB_d_hat_r = result_d_r$AB
  
  return(list(A_hat, B_hat, AB_hat,
              A_d_hat, B_d_hat, AB_d_hat,
              A_d_hat_r, B_d_hat_r, AB_d_hat_r,
              x_j_hat, y_j_hat,
              result_ivw, 
              result_adivw, 
              result_MrDAG,
              # result_GRAPPLE,
              rank_correct,
              r_tested
  ))
}


.get_sim_index = function(me_weight_index, effect_weight_index, len_me = 2){
  sim_index = (me_weight_index-1)*len_me + effect_weight_index
  return(sim_index)
}


.simulation_mrdag <- function(parameters, n = 1000) {
  # n = pz
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  C = parameters$C
  
  # gamma_j_star and Gamma_j_star
  gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  # sample x_j_hat, y_j_hat
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j, ] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j, ], Sigma = Sigma_X, tol = 100)
    y_j_hat[j, ] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j, ], Sigma = Sigma_Y, tol = 100)
  }
  
  list(x_j_hat = x_j_hat, y_j_hat = y_j_hat)
}


# ============================================================ #
# 3. generate low rank true C ####
# ============================================================ #
set.seed(123)
px = 9
py = 3
r_approx = 2
r =  min(px,py) # maximum rank of true C

C = matrix(rnorm(py*px), py, px)
U = svd(C)$u
V = svd(C)$v
C = U %*% diag(c(rep(1, r_approx),
                 rep(0, r-r_approx))) %*% t(V)


# ============================================================ #
# 4. set universal parameters ####
# ============================================================ #
me_weight_list = c(2.5, 1)
effect_weight_list = c(0.25, 1)
# library(mr.divw)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap)
# library(readxl)
# library(glmnet)


# source("R/MR_rr_simulation.R")
# source("R/MR_rr_estimators.R")
# data("lip_data")
# data("lip_corr")
# data("lip_samplesize")

set.seed(123)
# setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")

source("scripts/MR_rr_estimators.R")

####
.simulation_data_only <- function(parameters) {
  n = 1000
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  C = parameters$C
  
  # true latent gamma and Gamma
  gamma_j_star = MASS::mvrnorm(n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j, ] = MASS::mvrnorm(1, gamma_j_star[j, ], Sigma_X, tol = 100)
    y_j_hat[j, ] = MASS::mvrnorm(1, Gamma_j_star[j, ], Sigma_Y, tol = 100)
  }
  
  return(list(x_j_hat = x_j_hat, y_j_hat = y_j_hat))
}


me_weight_list = c("2.5", "1")
effect_weight_list = c("0.25", "1")

# from main simulation
regularization_rate_list = c(1.514199e-10, 5.241443e-13, 1.214013e-12, 8.647016e-17)


#### 260409 parallel bootstrap - except MrDAG ####
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = 1e-13, bootstrap_size = 500,
                                       iteration = 1000, r_rank = 2, n_cores = NULL,
                                       estimator_set = "mr_only",
                                       save_estimates = TRUE) {
  
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }
  
  me_index     = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index  = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters   = parameters_list[[param_index]]
  
  C            = parameters$C
  py           = parameters$py
  px           = parameters$px
  Sigma_X_hat  = parameters$Sigma_X
  Sigma_Y_hat  = parameters$Sigma_Y
  C_vec        = as.vector(C)
  d            = px * py
  pz           = 1000
  
  # ===== Pre-compute =====
  W_precomputed = solve(Sigma_Y_hat)
  
  estimator_order = switch(estimator_set,
                           "all"      = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG"),
                           "no_mrdag" = c("IVW", "adIVW", "Naive", "MR", "MR_r"),
                           "mr_only"  = c("Naive", "MR", "MR_r"),
                           stop("estimator_set must be one of: 'all', 'no_mrdag', 'mr_only'")
  )
  
  set_name = switch(estimator_set,
                    "all"      = "All estimators (including MrDAG)",
                    "no_mrdag" = "All estimators EXCEPT MrDAG",
                    "mr_only"  = "Only MR family (Naive, MR, MR_r)"
  )
  
  start_time = Sys.time()
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Starting Parallel Bootstrap\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("Estimators: %s\n", set_name))
  cat(sprintf("ME weight: %s, Effect weight: %s\n", me_weight, effect_weight))
  cat(sprintf("Iterations: %d, Bootstrap size: %d\n", iteration, bootstrap_size))
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Started at: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("========================================\n\n"))
  flush.console()
  
  # ===== Cluster =====
  cl <- makeCluster(n_cores)
  on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
  
  clusterExport(cl,
                varlist = c("parameters", "regularization_rate", "r_rank",
                            "estimator_order", "C_vec", "Sigma_X_hat", "Sigma_Y_hat",
                            "pz", "px", "py", "d", "bootstrap_size",
                            "estimator_set", "save_estimates", "W_precomputed"),
                envir = environment())
  
  clusterEvalQ(cl, {
    library(MASS)
    library(mr.divw)
    source("scripts/MR_rr_estimators.R")
  })
  
  # ===== Worker function: runs ONE iteration =====
  run_one_iter <- function(loop_id) {
    sim     = .simulation(parameters, regularization_rate)
    x_j_hat = sim[[10]]
    y_j_hat = sim[[11]]
    
    # matrix storage: rows = bootstrap draws, cols = entries of C
    result_mats <- lapply(estimator_order, function(i) matrix(NA_real_, bootstrap_size, d))
    names(result_mats) <- estimator_order
    
    for (bt in 1:bootstrap_size) {
      idx  = sample.int(pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]
      
      if (estimator_set %in% c("all", "no_mrdag", "mr_only")) {
        r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = W_precomputed)
        r_mr    = mr_rr(y_bt, x_bt, r = r_rank, W = W_precomputed, Sigma_X = Sigma_X_hat)
        r_mr_r  = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = W_precomputed,
                                    Sigma_X = Sigma_X_hat,
                                    regularization_rate = regularization_rate)
        result_mats$Naive[bt, ] = as.vector(r_naive$AB)
        result_mats$MR[bt, ]    = as.vector(r_mr$AB)
        result_mats$MR_r[bt, ]  = as.vector(r_mr_r$AB)
      }
      
      if (estimator_set %in% c("all", "no_mrdag")) {
        r_ivw   = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        result_mats$IVW[bt, ]   = as.vector(r_ivw)
        result_mats$adIVW[bt, ] = as.vector(r_adivw)
      }
      
      if (estimator_set == "all") {
        r_mrdag = Mr_DAG(y_bt, x_bt)
        result_mats$MrDAG[bt, ] = as.vector(r_mrdag)
      }
    }
    
    # ===== vectorized CI / coverage =====
    coverage_vec  = list()
    ci_length_vec = list()
    for (name in estimator_order) {
      ci <- apply(result_mats[[name]], 2, quantile,
                  probs = c(0.025, 0.975), na.rm = TRUE)
      coverage_vec[[name]]  = as.numeric(C_vec >= ci[1, ] & C_vec <= ci[2, ])
      ci_length_vec[[name]] = ci[2, ] - ci[1, ]
    }
    
    list(
      coverage  = coverage_vec,
      ci_length = ci_length_vec,
      estimates = if (save_estimates) result_mats else NULL
    )
  }
  
  clusterExport(cl, varlist = "run_one_iter", envir = environment())
  
  # ===== Run all iterations in parallel =====
  have_pb <- requireNamespace("pbapply", quietly = TRUE)
  if (have_pb) {
    coverage_results <- pbapply::pblapply(1:iteration, run_one_iter, cl = cl)
  } else {
    cat("(install 'pbapply' for a progress bar)\n"); flush.console()
    coverage_results <- parLapply(cl, 1:iteration, run_one_iter)
  }
  
  stopCluster(cl)
  on.exit()
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  fmt <- function(s) if (s < 60) sprintf("%.1f secs", s)
  else if (s < 3600) sprintf("%.1f mins", s/60)
  else sprintf("%dh %dm", floor(s/3600), floor((s %% 3600)/60))
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Bootstrap completed!\n"))
  cat(sprintf("Total time: %s\n", fmt(total_time)))
  cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("========================================\n\n"))
  flush.console()
  
  # ===== Aggregate =====
  coverage_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(r) as.numeric(r$coverage[[name]]))
    if (is.vector(mat)) mat = matrix(mat, nrow = 1)
    else if (nrow(mat) != d) mat = t(mat)
    mat
  })
  names(coverage_list) = estimator_order
  
  ci_length_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(r) as.numeric(r$ci_length[[name]]))
    if (is.vector(mat)) mat = matrix(mat, nrow = 1)
    else if (nrow(mat) != d) mat = t(mat)
    mat
  })
  names(ci_length_list) = estimator_order
  
  estimates_list = NULL
  if (save_estimates) {
    estimates_list = lapply(estimator_order, function(name) {
      lapply(coverage_results, function(r) r$estimates[[name]])
    })
    names(estimates_list) = estimator_order
  }
  
  med_coverage = lapply(coverage_list, function(mat) {
    median(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
  })
  
  ci_length_summary = lapply(ci_length_list, function(mat) {
    mean_lengths = rowMeans(mat, na.rm = TRUE)
    list(
      median = median(mean_lengths, na.rm = TRUE),
      mean   = mean(mean_lengths,   na.rm = TRUE),
      sd     = sd(mean_lengths,     na.rm = TRUE),
      min    = min(mean_lengths,    na.rm = TRUE),
      max    = max(mean_lengths,    na.rm = TRUE)
    )
  })
  
  cat("\nCoverage (median over entries):\n")
  for (name in estimator_order) {
    cat(sprintf("  %s: %.4f\n", name, med_coverage[[name]]))
  }
  cat("\nCI Length Summary:\n")
  for (name in estimator_order) {
    cat(sprintf("  %s: Median = %.4f, Mean = %.4f, SD = %.4f\n",
                name,
                ci_length_summary[[name]]$median,
                ci_length_summary[[name]]$mean,
                ci_length_summary[[name]]$sd))
  }
  flush.console()
  
  list(
    med_coverage      = med_coverage,
    ci_length_summary = ci_length_summary,
    coverage_matrix   = coverage_list,
    ci_length_matrix  = ci_length_list,
    estimates         = estimates_list
  )
}


set.seed(123)

# 方式1：运行单个任务
test_sim_pred_filename = "results/simulate_result_pred_260212_1000_newC_M1.RData"
load(test_sim_pred_filename)


result <- nonpara_bootstrap_parallel(
  parameters_list = simulate_result_prediction$parameters_list,
  me_weight = "1",
  effect_weight = "1",
  regularization_rate = 1e-13,
  bootstrap_size = 100, #150
  iteration = 1000, #300
  n_cores = 10,
  estimator_set = "no_mrdag",
  save_estimates = TRUE
)


save(result, file = "results/CI_coverage_parallel_no_mrdag_260409_1_1.RData")


#single setting result
load("results/CI_coverage_parallel_no_mrdag_260409_1_1.RData")
CI_cov_result = result

summary_table = data.frame(
  Estimator = character(),
  Coverage = numeric(),
  CI_Median_Length = numeric(),
  CI_Mean_Length = numeric(),
  stringsAsFactors = FALSE
)
for (est_name in names(CI_cov_result$med_coverage)) {
  summary_table = rbind(summary_table, data.frame(
    Estimator = est_name,
    Coverage = round(CI_cov_result$med_coverage[[est_name]],4),
    CI_Median_Length = round(CI_cov_result$ci_length_summary[[est_name]]$median,3),
    CI_Mean_Length = round(CI_cov_result$ci_length_summary[[est_name]]$mean,3)
  ))
}
print(summary_table)

# 60224 result
# 2.5 0.25
# Estimator Coverage CI_Median_Length CI_Mean_Length
# 1       IVW   0.3902            0.149          0.147
# 2     adIVW   0.9315            0.985          2.822
# 3     Naive   0.2952            0.110          0.117
# 4        MR   0.9580            1.060          1.381
# 5      MR_r   0.9577            0.916          1.092

# 1 0.25
# Estimator Coverage CI_Median_Length CI_Mean_Length
# 1       IVW   0.6702            0.176          0.185
# 2     adIVW   0.9206            0.330          0.445
# 3     Naive   0.6123            0.134          0.147
# 4        MR   0.9253            0.279          0.312
# 5      MR_r   0.9253            0.279          0.312

# 2.5 1
# Estimator Coverage CI_Median_Length CI_Mean_Length
# 1       IVW   0.5730            0.109          0.111
# 2     adIVW   0.8932            0.219          0.212
# 3     Naive   0.4740            0.079          0.088
# 4        MR   0.9240            0.124          0.150
# 5      MR_r   0.9240            0.124          0.150

#11
# Estimator Coverage CI_Median_Length CI_Mean_Length
# 1       IVW   0.7765            0.114          0.116
# 2     adIVW   0.9184            0.160          0.152
# 3     Naive   0.7159            0.085          0.091
# 4        MR   0.9178            0.111          0.115
# 5      MR_r   0.9178            0.111          0.115



#### 260409 MrDAG-only parallel bootstrap ####
nonpara_bootstrap_parallel_mrdag <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                             iteration = 1000, bootstrap_size = 100,
                                             n_cores = NULL,
                                             mrdag_niter = 500, mrdag_burnin = 200,
                                             save_estimates = TRUE) {
  
  if (is.null(n_cores)) n_cores = parallel::detectCores() - 1
  
  me_index     = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index  = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters   = parameters_list[[param_index]]
  
  C_vec = as.vector(parameters$C)
  px    = parameters$px
  py    = parameters$py
  d     = px * py
  pz    = 1000
  
  start_time = Sys.time()
  cat(sprintf("\n[MrDAG only] iter=%d, bt=%d, cores=%d | me=%s effect=%s\n",
              iteration, bootstrap_size, n_cores, me_weight, effect_weight))
  flush.console()
  
  cl <- parallel::makeCluster(n_cores)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
  
  wd0 <- getwd()
  parallel::clusterExport(cl, "wd0", envir = environment())
  parallel::clusterEvalQ(cl, setwd(wd0))
  
  parallel::clusterEvalQ(cl, {
    library(MASS)
    library(MrDAG)
    NULL
  })
  
  parallel::clusterExport(
    cl,
    varlist = c("parameters", "C_vec", "px", "py", "d", "pz",
                "bootstrap_size", "mrdag_niter", "mrdag_burnin", "save_estimates"),
    envir = environment()
  )
  
  for (fn in c(".get_sim_index", ".simulation_mrdag", "Mr_DAG")) {
    if (exists(fn, envir = environment(), inherits = TRUE)) {
      parallel::clusterExport(cl, fn, envir = environment())
    } else if (exists(fn, envir = .GlobalEnv, inherits = TRUE)) {
      parallel::clusterExport(cl, fn, envir = .GlobalEnv)
    } else {
      stop(paste("Missing function in master session:", fn))
    }
  }
  
  # ===== Worker: one iteration =====
  run_one_iter <- function(loop_id) {
    sim     = .simulation_mrdag(parameters, n = pz)
    x_j_hat = sim$x_j_hat
    y_j_hat = sim$y_j_hat
    
    # bootstrap_size x d matrix
    draws_mat <- matrix(NA_real_, bootstrap_size, d)
    
    for (bt in 1:bootstrap_size) {
      idx  = sample.int(pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, , drop = FALSE]
      y_bt = y_j_hat[idx, , drop = FALSE]
      
      est_mat = Mr_DAG(y_bt, x_bt, niter = mrdag_niter, burnin = mrdag_burnin)
      draws_mat[bt, ] = as.vector(est_mat)
    }
    
    ci <- apply(draws_mat, 2, quantile,
                probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE)
    
    list(
      coverage  = as.numeric(C_vec >= ci[1, ] & C_vec <= ci[2, ]),
      ci_length = ci[2, ] - ci[1, ],
      estimates = if (save_estimates) draws_mat else NULL
    )
  }
  
  parallel::clusterExport(cl, "run_one_iter", envir = environment())
  
  # ===== Run =====
  have_pb <- requireNamespace("pbapply", quietly = TRUE)
  coverage_results <- if (have_pb) {
    pbapply::pblapply(1:iteration, run_one_iter, cl = cl)
  } else {
    parallel::parLapply(cl, 1:iteration, run_one_iter)
  }
  
  parallel::stopCluster(cl)
  on.exit()
  
  # ===== Aggregate =====
  coverage_mat = sapply(coverage_results, function(z) z$coverage)
  if (is.vector(coverage_mat)) coverage_mat = matrix(coverage_mat, nrow = d)
  if (nrow(coverage_mat) != d)  coverage_mat = t(coverage_mat)
  
  ci_length_mat = sapply(coverage_results, function(z) z$ci_length)
  if (is.vector(ci_length_mat)) ci_length_mat = matrix(ci_length_mat, nrow = d)
  if (nrow(ci_length_mat) != d)  ci_length_mat = t(ci_length_mat)
  
  med_coverage = median(rowMeans(coverage_mat, na.rm = TRUE), na.rm = TRUE)
  
  mean_lengths = rowMeans(ci_length_mat, na.rm = TRUE)
  ci_length_summary = list(
    median = median(mean_lengths, na.rm = TRUE),
    mean   = mean(mean_lengths,   na.rm = TRUE),
    sd     = sd(mean_lengths,     na.rm = TRUE),
    min    = min(mean_lengths,    na.rm = TRUE),
    max    = max(mean_lengths,    na.rm = TRUE)
  )
  
  estimates_list = NULL
  if (save_estimates) {
    estimates_list = lapply(coverage_results, function(z) z$estimates)
  }
  
  total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("[MrDAG only] done | total=%0.1fs | coverage=%0.4f | CI_len_median=%0.4f\n",
              total_time, med_coverage, ci_length_summary$median))
  flush.console()
  
  list(
    med_coverage      = list(MrDAG = med_coverage),
    ci_length_summary = list(MrDAG = ci_length_summary),
    coverage_matrix   = list(MrDAG = coverage_mat),
    ci_length_matrix  = list(MrDAG = ci_length_mat),
    estimates         = list(MrDAG = estimates_list)
  )
}


## 用法（对齐你之前的格式）
set.seed(123)
load("results/simulate_result_pred_260409_1000_newC_M1.RData")

result_mrdag <- nonpara_bootstrap_parallel_mrdag(
  parameters_list = simulate_result_prediction$parameters_list,
  me_weight = "2.5",
  effect_weight = "0.25",
  iteration = 1000,
  bootstrap_size = 100,
  # niter_mrdag = 10000,   # 真实跑建议更大
  # burnin_mrdag = 2000,
  n_cores = 10,
  save_estimates = TRUE
)

save(result_mrdag, file = "results/CI_coverage_parallel_only_mrdag_260409_2.5_0.25.RData")

## summary_table 对齐
CI_cov_result <- result_mrdag
summary_table <- data.frame(
  Estimator = "MrDAG",
  Coverage = round(CI_cov_result$med_coverage$MrDAG, 4),
  CI_Median_Length = round(CI_cov_result$ci_length_summary$MrDAG$median, 3),
  CI_Mean_Length = round(CI_cov_result$ci_length_summary$MrDAG$mean, 3)
)
print(summary_table)
