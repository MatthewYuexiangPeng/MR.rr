#### loading environment and data ####
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


#### hidden functions ####
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


#### generate low rank true C ####
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

#### set universal parameters ####
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


#### 260223 parallel_onebyone_no MrDAG ####
# pre-calculate W
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                                 regularization_rate = 1e-13, bootstrap_size = 500,
                                                 iteration = 1000, r_rank = 2, n_cores = NULL,
                                                 estimator_set = "mr_only",
                                                 save_estimates = TRUE) {
  
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }
  
  me_index = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters = parameters_list[[param_index]]
  
  C = parameters$C
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X
  Sigma_Y_hat = parameters$Sigma_Y
  C_vec = as.vector(C)
  pz = 1000
  
  # ===== 关键优化：预计算 =====
  W_precomputed = solve(Sigma_Y_hat)  # 只计算一次！
  Sigma_X_inv_precomputed = solve(Sigma_X_hat)  # 如果需要的话
  
  # 预计算SE（如果用IVW/adIVW）
  se_exposure_row = NULL
  se_outcome_row = NULL
  if (estimator_set %in% c("all", "no_mrdag")) {
    se_exposure_row = t(as.matrix(sqrt(diag(Sigma_X_hat))))
    se_outcome_row = sqrt(diag(Sigma_Y_hat))
  }
  # ===============================
  
  estimator_order = switch(estimator_set,
                           "all" = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG"),
                           "no_mrdag" = c("IVW", "adIVW", "Naive", "MR", "MR_r"),
                           "mr_only" = c("Naive", "MR", "MR_r"),
                           stop("estimator_set must be one of: 'all', 'no_mrdag', 'mr_only'")
  )
  
  set_name = switch(estimator_set,
                    "all" = "All estimators (including MrDAG)",
                    "no_mrdag" = "All estimators EXCEPT MrDAG",
                    "mr_only" = "Only MR family (Naive, MR, MR_r)"
  )
  
  start_time = Sys.time()
  
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Starting OPTIMIZED Parallel Bootstrap\n"))
  cat(sprintf("========================================\n"))
  cat(sprintf("Estimators: %s\n", set_name))
  cat(sprintf("ME weight: %s, Effect weight: %s\n", me_weight, effect_weight))
  cat(sprintf("Iterations: %d, Bootstrap size: %d\n", iteration, bootstrap_size))
  cat(sprintf("Cores: %d\n", n_cores))
  cat(sprintf("Started at: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("========================================\n\n"))
  flush.console()
  
  # 创建集群
  cl <- makeCluster(n_cores)
  
  # 导出必要的变量和函数
  clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  
  helper_functions <- c(".simulation", ".simulation_data_only", ".sqrt_matrix", ".get_sim_index",
                        "mr_rr_naive", "mr_rr", "mr_rr_regularized",
                        "ivw_multiple_outcomes", "adivw_multiple_outcomes", "Mr_DAG")
  
  # ===== 关键：导出预计算的矩阵 =====
  clusterExport(cl, varlist = c("parameters", "regularization_rate", "r_rank",
                                "estimator_order", "C_vec", "Sigma_X_hat",
                                "Sigma_Y_hat", "pz", "px", "py", "bootstrap_size",
                                "estimator_set", "save_estimates",
                                "W_precomputed",  # 新增
                                "se_exposure_row", "se_outcome_row"),  # 新增
                envir = environment())
  
  for (func in helper_functions) {
    tryCatch({
      if (exists(func, envir = .GlobalEnv)) {
        clusterExport(cl, varlist = func, envir = .GlobalEnv)
      }
    }, error = function(e) {
      warning(paste("Could not export function:", func))
    })
  }
  
  clusterEvalQ(cl, {
    library(MASS)
    library(mr.divw)
  })
  
  clusterEvalQ(cl, {
    source("scripts/MR_rr_estimators.R")
  })
  
  # 进度跟踪
  progress_env <- new.env()
  progress_env$completed <- 0
  progress_env$start_time <- Sys.time()
  progress_env$iter_times <- numeric(0)
  
  format_time <- function(seconds) {
    if (is.na(seconds) || seconds < 0) return("N/A")
    if (seconds < 60) {
      return(sprintf("%.1f secs", seconds))
    } else if (seconds < 3600) {
      return(sprintf("%.1f mins", seconds / 60))
    } else {
      hours = floor(seconds / 3600)
      mins = floor((seconds %% 3600) / 60)
      return(sprintf("%dh %dm", hours, mins))
    }
  }
  
  coverage_results <- vector("list", iteration)
  
  for (loop in 1:iteration) {
    iter_start <- Sys.time()
    
    coverage_results[[loop]] <- parLapply(cl, list(loop), function(loop_id) {
      
      sim = .simulation(parameters, regularization_rate)
      x_j_hat = sim[[10]]
      y_j_hat = sim[[11]]
      
      result_lists = lapply(estimator_order, function(i) {
        replicate(px*py, numeric(bootstrap_size), simplify = FALSE)
      })
      names(result_lists) = estimator_order
      
      # ===== 优化后的Bootstrap循环 =====
      for (bt in 1:bootstrap_size) {
        idx = sample(1:pz, pz, replace = TRUE)
        x_bt = x_j_hat[idx, ]
        y_bt = y_j_hat[idx, ]
        
        # 使用预计算的W！
        if (estimator_set %in% c("all", "no_mrdag", "mr_only")) {
          r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = W_precomputed)  # 不再重复计算solve()
          r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = W_precomputed, Sigma_X = Sigma_X_hat)
          r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = W_precomputed,
                                     Sigma_X = Sigma_X_hat,
                                     regularization_rate = regularization_rate)
        }
        
        if (estimator_set %in% c("all", "no_mrdag")) {
          # 预计算SE，避免重复计算
          se_exposure = matrix(rep(se_exposure_row, each = pz), nrow = pz, byrow = TRUE)
          se_outcome = matrix(rep(se_outcome_row, each = pz), nrow = pz, byrow = TRUE)
          
          r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
          r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        }
        
        if (estimator_set == "all") {
          r_mrdag = Mr_DAG(y_bt, x_bt)
        }
        
        res_list = list()
        if (estimator_set %in% c("all", "no_mrdag")) {
          res_list$IVW = as.vector(r_ivw)
          res_list$adIVW = as.vector(r_adivw)
        }
        res_list$Naive = as.vector(r_naive$AB)
        res_list$MR = as.vector(r_mr$AB)
        res_list$MR_r = as.vector(r_mr_r$AB)
        if (estimator_set == "all") {
          res_list$MrDAG = as.vector(r_mrdag)
        }
        
        for (name in estimator_order) {
          for (j in 1:(px*py)) {
            result_lists[[name]][[j]][bt] = res_list[[name]][j]
          }
        }
      }
      
      coverage_vec = list()
      ci_length_vec = list()
      estimates_stored = list()
      
      for (name in estimator_order) {
        cov_vec = numeric(px*py)
        length_vec = numeric(px*py)
        
        if (save_estimates) {
          estimates_stored[[name]] = result_lists[[name]]
        }
        
        for (j in 1:(px*py)) {
          ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
          cov_vec[j] = as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
          length_vec[j] = ci[2] - ci[1]
        }
        
        coverage_vec[[name]] = cov_vec
        ci_length_vec[[name]] = length_vec
      }
      
      return(list(
        coverage = coverage_vec,
        ci_length = ci_length_vec,
        estimates = if(save_estimates) estimates_stored else NULL
      ))
    })[[1]]
    
    iter_end <- Sys.time()
    iter_duration <- as.numeric(difftime(iter_end, iter_start, units = "secs"))
    progress_env$iter_times <- c(progress_env$iter_times, iter_duration)
    progress_env$completed <- loop
    
    avg_time_per_iter <- mean(progress_env$iter_times)
    remaining_iters <- iteration - loop
    estimated_remaining_secs <- avg_time_per_iter * remaining_iters
    estimated_completion <- Sys.time() + estimated_remaining_secs
    
    progress_pct <- (loop / iteration) * 100
    cat(sprintf("[%3d/%d] (%.1f%%) | Time: %s | Avg: %s/iter | Remaining: %s | ETA: %s\n",
                loop, iteration, progress_pct,
                format_time(iter_duration),
                format_time(avg_time_per_iter),
                format_time(estimated_remaining_secs),
                format(estimated_completion, "%H:%M:%S")))
    flush.console()
    
    if (loop %% 10 == 0 || loop == iteration) {
      elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      cat(sprintf("  >> Progress: %.1f%% | Elapsed: %s | ETA: %s (at %s)\n",
                  progress_pct,
                  format_time(elapsed_total),
                  format_time(estimated_remaining_secs),
                  format(estimated_completion, "%Y-%m-%d %H:%M:%S")))
      flush.console()
    }
  }
  
  stopCluster(cl)
  
  total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Bootstrap completed!\n"))
  cat(sprintf("Total time: %s\n", format_time(total_time)))
  cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("========================================\n\n"))
  flush.console()
  
  # 整理结果（同之前）
  coverage_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result$coverage[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(coverage_list) = estimator_order
  
  ci_length_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result$ci_length[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(ci_length_list) = estimator_order
  
  estimates_list = NULL
  if (save_estimates) {
    estimates_list = lapply(estimator_order, function(name) {
      lapply(coverage_results, function(iter_result) {
        iter_result$estimates[[name]]
      })
    })
    names(estimates_list) = estimator_order
  }
  
  avg_coverage = lapply(coverage_list, function(mat) {
    mean(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
  })
  
  ci_length_summary = lapply(ci_length_list, function(mat) {
    mean_lengths = rowMeans(mat, na.rm = TRUE)
    
    list(
      median = median(mean_lengths, na.rm = TRUE),
      mean = mean(mean_lengths, na.rm = TRUE),
      sd = sd(mean_lengths, na.rm = TRUE),
      min = min(mean_lengths, na.rm = TRUE),
      max = max(mean_lengths, na.rm = TRUE)
    )
  })
  
  cat("\nCI Length Summary:\n")
  for (name in estimator_order) {
    cat(sprintf("  %s: Median = %.4f, Mean = %.4f, SD = %.4f\n", 
                name, 
                ci_length_summary[[name]]$median,
                ci_length_summary[[name]]$mean,
                ci_length_summary[[name]]$sd))
  }
  flush.console()
  
  return(list(
    avg_coverage = avg_coverage,
    ci_length_summary = ci_length_summary,
    coverage_matrix = coverage_list,
    ci_length_matrix = ci_length_list,
    estimates = estimates_list
  ))
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
  bootstrap_size = 300, #150
  iteration = 500, #300
  n_cores = 6,
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
for (est_name in names(CI_cov_result$avg_coverage)) {
  summary_table = rbind(summary_table, data.frame(
    Estimator = est_name,
    Coverage = round(CI_cov_result$avg_coverage[[est_name]],4),
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



#### 260226 parallel_only MrDAG ####
# 只跑 MrDAG：产出 avg_coverage / ci_length_summary / coverage_matrix / ci_length_matrix / estimates（可选）
nonpara_bootstrap_parallel_mrdag <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                             iteration = 300, bootstrap_size = 200,
                                             n_cores = NULL,
                                             mrdag_niter = 1000, mrdag_burnin = 200,
                                             save_estimates = TRUE) {
  
  if (is.null(n_cores)) n_cores = parallel::detectCores() - 1
  
  me_index = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters = parameters_list[[param_index]]
  
  C_vec = as.vector(parameters$C)
  px = parameters$px
  py = parameters$py
  pz = 1000
  
  start_time = Sys.time()
  cat(sprintf("\n[MrDAG only] iter=%d, bt=%d, cores=%d | me=%s effect=%s\n",
              iteration, bootstrap_size, n_cores, me_weight, effect_weight))
  
  cl <- parallel::makeCluster(n_cores)
  
  # 保证路径一致（防止 source 相对路径挂）
  wd0 <- getwd()
  parallel::clusterExport(cl, "wd0", envir = environment())
  parallel::clusterEvalQ(cl, setwd(wd0))
  
  # worker load
  parallel::clusterEvalQ(cl, {
    library(MASS)
    # MrDAG 包在你的 wrapper 里会用到
    # 如果 wrapper 直接调用 MrDAG::MrDAG / get_causaleffects，需要确保库可用
    library(MrDAG)
    NULL
  })
  
  # export objects
  parallel::clusterExport(
    cl,
    varlist = c("parameters", "C_vec", "px", "py", "pz", "bootstrap_size",
                "mrdag_niter", "mrdag_burnin", "save_estimates"),
    envir = environment()
  )
  
  # export functions（别只在 .GlobalEnv 找）
  for (fn in c(".get_sim_index", ".simulation_mrdag", "Mr_DAG")) {
    if (exists(fn, envir = environment(), inherits = TRUE)) {
      parallel::clusterExport(cl, fn, envir = environment())
    } else if (exists(fn, envir = .GlobalEnv, inherits = TRUE)) {
      parallel::clusterExport(cl, fn, envir = .GlobalEnv)
    } else {
      stop(paste("Missing function in master session:", fn))
    }
  }
  
  # main
  coverage_results <- parallel::parLapply(cl, 1:iteration, function(loop_id) {
    
    sim = .simulation_mrdag(parameters, n = pz)
    x_j_hat = sim$x_j_hat
    y_j_hat = sim$y_j_hat
    
    # store bootstrap draws for each entry
    # list of length (px*py); each element numeric(bootstrap_size)
    draws = replicate(px * py, numeric(bootstrap_size), simplify = FALSE)
    
    for (bt in 1:bootstrap_size) {
      idx = sample.int(pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, , drop = FALSE]
      y_bt = y_j_hat[idx, , drop = FALSE]
      
      # MrDAG estimate (py x px matrix)
      est_mat = Mr_DAG(y_bt, x_bt, niter = mrdag_niter, burnin = mrdag_burnin)
      
      est_vec = as.vector(est_mat)
      for (j in 1:(px * py)) {
        draws[[j]][bt] = est_vec[j]
      }
    }
    
    # coverage and ci length
    cov_vec = numeric(px * py)
    len_vec = numeric(px * py)
    for (j in 1:(px * py)) {
      ci = quantile(draws[[j]], c(0.025, 0.975), na.rm = TRUE, names = FALSE)
      cov_vec[j] = as.numeric(C_vec[j] >= ci[1] && C_vec[j] <= ci[2])
      len_vec[j] = ci[2] - ci[1]
    }
    
    list(
      coverage = cov_vec,
      ci_length = len_vec,
      estimates = if (save_estimates) draws else NULL
    )
  })
  
  parallel::stopCluster(cl)
  
  # coverage matrix: (px*py) x iteration
  coverage_mat = sapply(coverage_results, function(z) z$coverage)
  if (is.vector(coverage_mat)) coverage_mat = matrix(coverage_mat, nrow = px * py)
  if (nrow(coverage_mat) != px * py) coverage_mat = t(coverage_mat)
  
  ci_length_mat = sapply(coverage_results, function(z) z$ci_length)
  if (is.vector(ci_length_mat)) ci_length_mat = matrix(ci_length_mat, nrow = px * py)
  if (nrow(ci_length_mat) != px * py) ci_length_mat = t(ci_length_mat)
  
  avg_coverage = mean(rowMeans(coverage_mat, na.rm = TRUE), na.rm = TRUE)
  
  mean_lengths = rowMeans(ci_length_mat, na.rm = TRUE)
  ci_length_summary = list(
    median = median(mean_lengths, na.rm = TRUE),
    mean = mean(mean_lengths, na.rm = TRUE),
    sd = sd(mean_lengths, na.rm = TRUE),
    min = min(mean_lengths, na.rm = TRUE),
    max = max(mean_lengths, na.rm = TRUE)
  )
  
  estimates_list = NULL
  if (save_estimates) {
    estimates_list = lapply(coverage_results, function(z) z$estimates)
  }
  
  total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("[MrDAG only] done | total=%0.1fs | coverage=%0.4f | CI_len_median=%0.4f\n",
              total_time, avg_coverage, ci_length_summary$median))
  
  list(
    avg_coverage = list(MrDAG = avg_coverage),
    ci_length_summary = list(MrDAG = ci_length_summary),
    coverage_matrix = list(MrDAG = coverage_mat),
    ci_length_matrix = list(MrDAG = ci_length_mat),
    estimates = list(MrDAG = estimates_list)
  )
}

## 用法（对齐你之前的格式）
set.seed(123)
load("results/simulate_result_pred_260212_1000_newC_M1.RData")

result_mrdag <- nonpara_bootstrap_parallel_mrdag(
  parameters_list = simulate_result_prediction$parameters_list,
  me_weight = "2.5",
  effect_weight = "0.25",
  iteration = 300,
  bootstrap_size = 100,
  # niter_mrdag = 10000,   # 真实跑建议更大
  # burnin_mrdag = 2000,
  n_cores = 6,
  save_estimates = TRUE
)

save(result_mrdag, file = "results/CI_coverage_parallel_only_mrdag_260301_2.5_0.25.RData")

## summary_table 对齐
CI_cov_result <- result_mrdag
summary_table <- data.frame(
  Estimator = "MrDAG",
  Coverage = round(CI_cov_result$avg_coverage$MrDAG, 4),
  CI_Median_Length = round(CI_cov_result$ci_length_summary$MrDAG$median, 3),
  CI_Mean_Length = round(CI_cov_result$ci_length_summary$MrDAG$mean, 3)
)
print(summary_table)

#### 方式2：运行所有任务 ####
# 主程序：带总体进度的批量运行
run_all_parallel_bootstrap <- function(test_sim_pred_filename = "results/simulate_result_pred_260212_1000_newC_M1.RData",
                                       bootstrap_size = 300,
                                       iteration = 500,
                                       n_cores = 3,
                                       estimator_set = "no_mrdag",
                                       save_estimates = TRUE,
                                       temp_save_file = "results/CI_coverage_parallel_temp.RData") {
  load(test_sim_pred_filename)
  
  # 计算总任务数
  total_tasks <- length(me_weight_list) * length(effect_weight_list)
  current_task <- 0
  overall_start <- Sys.time()
  task_times <- numeric(0)
  
  CI_cov_result <- list()
  
  cat(sprintf("\n###########################################\n"))
  cat(sprintf("Running ALL Parallel Bootstrap Tasks\n"))
  cat(sprintf("###########################################\n"))
  cat(sprintf("Estimator set: %s\n", estimator_set))
  cat(sprintf("Total tasks: %d\n", total_tasks))
  cat(sprintf("Bootstrap size: %d, Iterations: %d\n", bootstrap_size, iteration))
  cat(sprintf("Cores per task: %d\n", n_cores))
  cat(sprintf("Started at: %s\n", format(overall_start, "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("###########################################\n\n"))
  flush.console()
  
  # 时间格式化函数
  format_time <- function(seconds) {
    if (is.na(seconds) || seconds < 0) return("N/A")
    if (seconds < 60) {
      return(sprintf("%.1f secs", seconds))
    } else if (seconds < 3600) {
      return(sprintf("%.1f mins", seconds / 60))
    } else {
      hours = floor(seconds / 3600)
      mins = floor((seconds %% 3600) / 60)
      return(sprintf("%dh %dm", hours, mins))
    }
  }
  
  for (i in seq_along(me_weight_list)) {
    for (j in seq_along(effect_weight_list)) {
      current_task <- current_task + 1
      task_start <- Sys.time()
      
      me <- me_weight_list[i]
      eff <- effect_weight_list[j]
      reg <- regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
      
      cat(sprintf("\n╔══════════════════════════════════════╗\n"))
      cat(sprintf("║ Task %d/%d: me=%s, effect=%s\n", current_task, total_tasks, me, eff))
      cat(sprintf("╚══════════════════════════════════════╝\n"))
      flush.console()
      
      result <- tryCatch({
        nonpara_bootstrap_parallel(
          parameters_list = simulate_result_prediction$parameters_list,
          me_weight = me,
          effect_weight = eff,
          regularization_rate = reg,
          bootstrap_size = bootstrap_size,
          iteration = iteration,
          n_cores = n_cores,
          estimator_set = estimator_set,
          save_estimates = save_estimates
        )
      }, error = function(e) {
        cat(sprintf("ERROR in task %d: %s\n", current_task, e$message))
        flush.console()
        return(NULL)
      })
      
      CI_cov_result[[paste0("me_", me, "_eff_", eff)]] <- result
      
      # 计算任务耗时和总体进度
      task_end <- Sys.time()
      task_duration <- as.numeric(difftime(task_end, task_start, units = "secs"))
      task_times <- c(task_times, task_duration)
      
      # 预估剩余时间
      avg_task_time <- mean(task_times)
      remaining_tasks <- total_tasks - current_task
      estimated_remaining <- avg_task_time * remaining_tasks
      eta <- Sys.time() + estimated_remaining
      
      # 打印任务完成信息
      cat(sprintf("\n>>> Task %d completed in %s\n", current_task, format_time(task_duration)))
      cat(sprintf(">>> Overall Progress: %d/%d (%.1f%%)\n", 
                  current_task, total_tasks, 100 * current_task / total_tasks))
      cat(sprintf(">>> Avg task time: %s | Remaining: %s | ETA: %s\n\n",
                  format_time(avg_task_time),
                  format_time(estimated_remaining),
                  format(eta, "%Y-%m-%d %H:%M:%S")))
      flush.console()
      
      # 每完成一个任务就保存
      save(CI_cov_result, file = temp_save_file)
      
      # 强制垃圾回收
      gc()
    }
  }
  
  # 最终保存
  timestamp <- format(Sys.time(), "%m%d")
  final_file <- sprintf("results/CI_coverage_%s_%d_%d_%s.RData", 
                        estimator_set, bootstrap_size, iteration, timestamp)
  save(CI_cov_result, file = final_file)
  
  # 总结信息
  total_elapsed <- as.numeric(difftime(Sys.time(), overall_start, units = "hours"))
  cat(sprintf("\n###########################################\n"))
  cat(sprintf("ALL TASKS COMPLETED!\n"))
  cat(sprintf("###########################################\n"))
  cat(sprintf("Total time: %.2f hours\n", total_elapsed))
  cat(sprintf("Average per task: %s\n", format_time(mean(task_times))))
  cat(sprintf("Final results saved to:\n  %s\n", final_file))
  cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("###########################################\n\n"))
  flush.console()
  
  return(CI_cov_result)
}


# all_results <- run_all_parallel_bootstrap(
#   test_sim_pred_filename = "results/simulate_result_pred_260212_1000_newC_M1.RData",
#   bootstrap_size = 100,
#   iteration = 300,
#   n_cores = 6,
#   estimator_set = "no_mrdag",
#   save_estimates = TRUE,
#   temp_save_file = "results/CI_coverage_parallel_temp_no_mrdag.RData"
# )


# load data
# load("results/results1104/1104.RData")
# load("results/results1104/CI_coverage_no_mrdag_150_300_1104.RData")
# load("results/results1104/CI_coverage_parallel_temp_no_mrdag.RData")

# or
# CI_cov_result = all_results
# CI_cov_result = result
# 
# cat("\n\n=== Overall Results Summary ===\n\n")
# 
# # 创建汇总表
# summary_table = data.frame(
#   Setting = character(),
#   Estimator = character(),
#   Coverage = numeric(),
#   CI_Median_Length = numeric(),
#   CI_Mean_Length = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# for (setting_name in names(CI_cov_result)) {
#   result = CI_cov_result[[setting_name]]
#   
#   if (!is.null(result)) {
#     for (est_name in names(result$avg_coverage)) {
#       summary_table = rbind(summary_table, data.frame(
#         Setting = setting_name,
#         Estimator = est_name,
#         Coverage = round(result$avg_coverage[[est_name]],4),
#         CI_Median_Length = round(result$ci_length_summary[[est_name]]$median,3),
#         CI_Mean_Length = round(result$ci_length_summary[[est_name]]$mean,3)
#       ))
#     }
#   }
# }
# 
# print(summary_table)
# 
# # 保存汇总表
# write.csv(summary_table, file = "results/summary_table_150_300_251104.csv", row.names = FALSE)
# 
# cat("\n\nResults saved to:\n")
# cat("  - results/CI_coverage_parallel_final_500_1000_mr_only.RData\n")
# cat("  - results/summary_table_mr_only.csv\n")


# #### non parallel_1102_outdated ####
# # 修改后的函数：添加进度显示与预估完成时间
# nonpara_bootstrap <- function(parameters_list, me_weight = "1", effect_weight = "1", 
#                               regularization_rate = 1e-13, bootstrap_size = 100, 
#                               iteration = 100, r_rank = 2) {
#   me_index = match(as.character(me_weight), me_weight_list)
#   effect_index = match(as.character(effect_weight), effect_weight_list)
#   param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
#   parameters = parameters_list[[param_index]]
#   
#   C = parameters$C
#   py = parameters$py
#   px = parameters$px
#   Sigma_X_hat = parameters$Sigma_X
#   Sigma_Y_hat = parameters$Sigma_Y
#   C_vec = as.vector(C)
#   pz = 1000
#   
#   estimator_order = c("IVW", "adIVW", "Naive", "MR", "MR_r")
#   coverage_list = lapply(estimator_order, function(i) matrix(NA, px*py, iteration))
#   names(coverage_list) = estimator_order
#   
#   # 初始化时间记录
#   start_time = Sys.time()
#   iteration_times = numeric(0)  # 存储每个iteration的耗时
#   
#   cat(sprintf("\n========================================\n"))
#   cat(sprintf("Starting bootstrap: me = %s, effect = %s\n", me_weight, effect_weight))
#   cat(sprintf("Total iterations: %d\n", iteration))
#   cat(sprintf("Bootstrap size: %d\n", bootstrap_size))
#   cat(sprintf("Started at: %s\n", format(start_time, "%Y-%m-%d %H:%M:%S")))
#   cat(sprintf("========================================\n\n"))
#   flush.console()  # 强制刷新输出
#   
#   for (loop in 1:iteration) {
#     iter_start = Sys.time()
#     
#     sim = .simulation_data_only(parameters)
#     x_j_hat = sim$x_j_hat
#     y_j_hat = sim$y_j_hat
#     
#     result_lists = lapply(estimator_order, function(i) replicate(px*py, numeric(bootstrap_size), simplify = FALSE))
#     names(result_lists) = estimator_order
#     
#     for (bt in 1:bootstrap_size) {
#       idx = sample(1:pz, pz, replace = TRUE)
#       x_bt = x_j_hat[idx, ]
#       y_bt = y_j_hat[idx, ]
#       
#       r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
#       r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
#       r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat,
#                                  regularization_rate = regularization_rate)
#       r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
#       r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
#       
#       res_list = list(
#         IVW = as.vector(r_ivw),
#         adIVW = as.vector(r_adivw),
#         Naive = as.vector(r_naive$AB),
#         MR = as.vector(r_mr$AB),
#         MR_r = as.vector(r_mr_r$AB)
#       )
#       
#       for (name in estimator_order) {
#         for (j in 1:(px*py)) {
#           result_lists[[name]][[j]][bt] = res_list[[name]][j]
#         }
#       }
#     }
#     
#     for (name in estimator_order) {
#       for (j in 1:(px*py)) {
#         ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
#         coverage_list[[name]][j, loop] = as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
#       }
#     }
#     
#     # 计算本次iteration耗时
#     iter_end = Sys.time()
#     iter_duration = as.numeric(difftime(iter_end, iter_start, units = "secs"))
#     iteration_times = c(iteration_times, iter_duration)
#     
#     # 计算平均耗时和预估剩余时间
#     avg_time_per_iter = mean(iteration_times)
#     remaining_iters = iteration - loop
#     estimated_remaining_secs = avg_time_per_iter * remaining_iters
#     estimated_completion = Sys.time() + estimated_remaining_secs
#     
#     # 格式化时间显示
#     format_time = function(seconds) {
#       if (seconds < 60) {
#         return(sprintf("%.1f secs", seconds))
#       } else if (seconds < 3600) {
#         return(sprintf("%.1f mins", seconds / 60))
#       } else {
#         hours = floor(seconds / 3600)
#         mins = floor((seconds %% 3600) / 60)
#         return(sprintf("%dh %dm", hours, mins))
#       }
#     }
#     
#     # 打印进度信息
#     progress_pct = (loop / iteration) * 100
#     cat(sprintf("[%3d/%d] (%.1f%%) | Time: %s | Avg: %s/iter | Remaining: %s | ETA: %s\n",
#                 loop, iteration, progress_pct,
#                 format_time(iter_duration),
#                 format_time(avg_time_per_iter),
#                 format_time(estimated_remaining_secs),
#                 format(estimated_completion, "%H:%M:%S")))
#     flush.console()  # 强制刷新输出
#     
#     # 每10次iteration显示一次详细信息
#     if (loop %% 10 == 0) {
#       elapsed_total = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
#       cat(sprintf("  >> Progress: %.1f%% | Elapsed: %s | ETA: %s (at %s)\n",
#                   progress_pct,
#                   format_time(elapsed_total),
#                   format_time(estimated_remaining_secs),
#                   format(estimated_completion, "%Y-%m-%d %H:%M:%S")))
#       flush.console()  # 强制刷新输出
#     }
#   }
#   
#   # 完成信息
#   total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
#   cat(sprintf("\n========================================\n"))
#   cat(sprintf("Bootstrap completed!\n"))
#   cat(sprintf("Total time: %s\n", format_time(total_time)))
#   cat(sprintf("Completed at: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
#   cat(sprintf("========================================\n\n"))
#   flush.console()  # 强制刷新输出
#   
#   avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
#   return(avg)
# }
# 
# 
# # 修改后的运行脚本：单独选择weight组合
# run_single_bootstrap <- function(me_weight, effect_weight, 
#                                  test_sim_pred_filename = "results/simulate_result_pred_260212_1000_newC_M1.RData",
#                                  bootstrap_size = 5, #50
#                                  iteration = 5,#1000
#                                  save_result = TRUE,
#                                  output_filename = NULL) {
#   
#   # 加载数据
#   load(test_sim_pred_filename)
#   
#   # 找到对应的索引
#   me_index = match(as.character(me_weight), me_weight_list)
#   effect_index = match(as.character(effect_weight), effect_weight_list)
#   
#   if (is.na(me_index) || is.na(effect_index)) {
#     stop(sprintf("Invalid weight combination: me = %s, effect = %s", me_weight, effect_weight))
#   }
#   
#   param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
#   reg = regularization_rate_list[param_index]
#   
#   cat(sprintf("\n###########################################\n"))
#   cat(sprintf("Running Bootstrap Analysis\n"))
#   cat(sprintf("###########################################\n"))
#   cat(sprintf("ME weight: %s\n", me_weight))
#   cat(sprintf("Effect weight: %s\n", effect_weight))
#   cat(sprintf("Regularization rate: %e\n", reg))
#   cat(sprintf("Bootstrap size: %d\n", bootstrap_size))
#   cat(sprintf("Iterations: %d\n", iteration))
#   cat(sprintf("###########################################\n"))
#   flush.console()  # 强制刷新输出
#   
#   # 运行bootstrap
#   result = nonpara_bootstrap(
#     parameters_list = simulate_result_prediction$parameters_list,
#     me_weight = me_weight,
#     effect_weight = effect_weight,
#     regularization_rate = reg,
#     bootstrap_size = bootstrap_size, 
#     iteration = iteration
#   )
#   
#   # 打印结果
#   cat("\n>>> Coverage Results:\n")
#   for (estimator in names(result)) {
#     cat(sprintf("  %s: %.4f\n", estimator, result[[estimator]]))
#   }
#   flush.console()  # 强制刷新输出
#   
#   # 保存结果
#   if (save_result) {
#     if (is.null(output_filename)) {
#       timestamp = format(Sys.time(), "%y%m%d_%H%M")
#       output_filename = sprintf("results/CI_cov_me%s_eff%s_%s.RData", 
#                                 me_weight, effect_weight, timestamp)
#     }
#     
#     CI_cov_result = list()
#     CI_cov_result[[paste0("me_", me_weight, "_eff_", effect_weight)]] = result
#     
#     save(CI_cov_result, file = output_filename)
#     cat(sprintf("\n>>> Results saved to: %s\n\n", output_filename))
#     flush.console()  # 强制刷新输出
#   }
#   
#   return(result)
# }
# 
# 
# # 使用示例
# # 
# 
# # set.seed(123)
# # # 单独运行一个weight组合
# # result = run_single_bootstrap(
# #   me_weight = "1",
# #   effect_weight = "0.25",
# #   bootstrap_size = 300,
# #   iteration = 500
# # )
# 
# set.seed(123)
# # 或者批量运行所有组合（带进度）
# run_all_bootstraps <- function(test_sim_pred_filename = "results/simulate_result_pred_260212_1000_newC_M1.RData",
#                                bootstrap_size = 50,
#                                iteration = 1000) {
#   
#   load(test_sim_pred_filename)
#   
#   total_combinations = length(me_weight_list) * length(effect_weight_list)
#   current_combo = 0
#   overall_start = Sys.time()
#   
#   CI_cov_result = list()
#   
#   cat(sprintf("\n###########################################\n"))
#   cat(sprintf("Running ALL Bootstrap Combinations\n"))
#   cat(sprintf("###########################################\n"))
#   cat(sprintf("Total combinations: %d\n", total_combinations))
#   cat(sprintf("Started at: %s\n", format(overall_start, "%Y-%m-%d %H:%M:%S")))
#   cat(sprintf("###########################################\n\n"))
#   flush.console()  # 强制刷新输出
#   
#   for (i in seq_along(me_weight_list)) {
#     for (j in seq_along(effect_weight_list)) {
#       current_combo = current_combo + 1
#       me = me_weight_list[i]
#       eff = effect_weight_list[j]
#       
#       cat(sprintf("\n>>> Combination %d/%d: me = %s, effect = %s\n",
#                   current_combo, total_combinations, me, eff))
#       flush.console()  # 强制刷新输出
#       
#       result = run_single_bootstrap(
#         me_weight = me,
#         effect_weight = eff,
#         test_sim_pred_filename = test_sim_pred_filename,
#         bootstrap_size = bootstrap_size,
#         iteration = iteration,
#         save_result = FALSE
#       )
#       
#       CI_cov_result[[paste0("me_", me, "_eff_", eff)]] = result
#       
#       # 估算总体剩余时间
#       elapsed = as.numeric(difftime(Sys.time(), overall_start, units = "secs"))
#       avg_per_combo = elapsed / current_combo
#       remaining_combos = total_combinations - current_combo
#       remaining_time = avg_per_combo * remaining_combos
#       eta = Sys.time() + remaining_time
#       
#       cat(sprintf("\n>>> Overall Progress: %d/%d (%.1f%%) | ETA: %s\n\n",
#                   current_combo, total_combinations,
#                   100 * current_combo / total_combinations,
#                   format(eta, "%Y-%m-%d %H:%M:%S")))
#       flush.console()  # 强制刷新输出
#     }
#   }
#   
#   # 保存所有结果
#   timestamp = format(Sys.time(), "%y%m%d_%H%M")
#   output_file = sprintf("results/CI_cov_all_combinations_%s.RData", timestamp)
#   save(CI_cov_result, file = output_file)
#   
#   total_time = as.numeric(difftime(Sys.time(), overall_start, units = "hours"))
#   cat(sprintf("\n###########################################\n"))
#   cat(sprintf("ALL COMBINATIONS COMPLETED!\n"))
#   cat(sprintf("Total time: %.2f hours\n", total_time))
#   cat(sprintf("Results saved to: %s\n", output_file))
#   cat(sprintf("###########################################\n\n"))
#   flush.console()  # 强制刷新输出
#   
#   return(CI_cov_result)
# }
# 
# # 运行所有组合
# all_results = run_all_bootstraps(bootstrap_size = 150, iteration = 300)
# 
