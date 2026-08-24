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
setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

exp_name = read.csv('data/traits_1e-4.csv')$x

# delete first out 
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
  r_tested <- rank_test(
    W = W,
    y_j_hat = y_j_hat,
    x_j_hat = x_j_hat,
    Sigma_X = Sigma_X,
    print = FALSE,
    bt_loop = 1000,
    alpha = 0.05,
    seed = 123,
    return_details = FALSE
  )
  # M0 or M1 
  # r_tested_M1 <- rank_test(W, y_j_hat, x_j_hat, Sigma_X, print = FALSE, min_rank = 1)
  
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


#### functions for rank selection ####
.make_psd <- function(S, eps = 1e-10){
  S <- (S + t(S)) / 2
  e <- eigen(S, symmetric = TRUE)
  vals <- pmax(e$values, eps)
  e$vectors %*% diag(vals, length(vals)) %*% t(e$vectors)
}

# estimate Sigma_{gg} and Sigma_Y from observed hats (optional, but useful defaults)
.estimate_Sigma_gg <- function(x_j_hat, Sigma_X, eps = 1e-10){
  Sx <- cov(x_j_hat)
  Sgg <- Sx - Sigma_X
  .make_psd(Sgg, eps = eps)
}

.estimate_Sigma_Y <- function(y_j_hat, C_hat0, Sigma_gg_hat, eps = 1e-10){
  Sy <- cov(y_j_hat)
  SY <- Sy - C_hat0 %*% Sigma_gg_hat %*% t(C_hat0)
  .make_psd(SY, eps = eps)
}

# unconstrained (full-rank) and constrained rank-r estimators via your existing mr_rr()
.fit_C0 <- function(y_j_hat, x_j_hat, W, Sigma_X){
  k <- min(ncol(x_j_hat), ncol(y_j_hat))
  mr_rr(y_j_hat, x_j_hat, r = k, W = W, Sigma_X = Sigma_X)$AB
}

.fit_Cr <- function(y_j_hat, x_j_hat, r, W, Sigma_X){
  mr_rr(y_j_hat, x_j_hat, r = r, W = W, Sigma_X = Sigma_X)$AB
}


rank_test <- function(W, y_j_hat, x_j_hat, Sigma_X,
                         Sigma_Y = NULL, Sigma_gg = NULL,
                         print = TRUE, min_rank = 1,
                         bt_loop = 5000, alpha = 0.05,
                         eps_psd = 1e-10, seed = NULL,
                         return_details = TRUE){
  
  if (!is.null(seed)) set.seed(seed)
  
  px <- ncol(x_j_hat)
  py <- ncol(y_j_hat)
  pz <- nrow(y_j_hat)
  k  <- min(px, py)
  
  # Step 1: fit unconstrained C^0
  C0_hat <- .fit_C0(y_j_hat, x_j_hat, W = W, Sigma_X = Sigma_X)
  
  # estimate Sigma_gg if not given
  if (is.null(Sigma_gg)) {
    Sigma_gg_hat <- .estimate_Sigma_gg(x_j_hat, Sigma_X = Sigma_X, eps = eps_psd)
  } else {
    Sigma_gg_hat <- .make_psd(Sigma_gg, eps = eps_psd)
  }
  
  # estimate Sigma_Y if not given
  if (is.null(Sigma_Y)) {
    Sigma_Y_hat <- .estimate_Sigma_Y(y_j_hat, C_hat0 = C0_hat,
                                     Sigma_gg_hat = Sigma_gg_hat, eps = eps_psd)
  } else {
    Sigma_Y_hat <- .make_psd(Sigma_Y, eps = eps_psd)
  }
  
  # precompute chol for fast simulation
  chol_gg <- chol(Sigma_gg_hat + diag(eps_psd, px))
  chol_X  <- chol(Sigma_X      + diag(eps_psd, px))
  chol_Y  <- chol(Sigma_Y_hat  + diag(eps_psd, py))
  
  out <- data.frame(r = min_rank:(k-1),
                    T_obs = NA_real_,
                    q_crit = NA_real_,
                    reject = NA,
                    stringsAsFactors = FALSE)
  
  for (idx in seq_len(nrow(out))) {
    r <- out$r[idx]
    
    # Step 2: fit rank-r constrained C^(r)
    Cr_hat <- .fit_Cr(y_j_hat, x_j_hat, r = r, W = W, Sigma_X = Sigma_X)
    
    # residual and observed statistic
    E_r_hat <- C0_hat - Cr_hat
    T_obs <- sum(E_r_hat^2)  # ||E_r||_F^2
    
    # Step 4: constrained bootstrap under H0(r)
    T_star <- numeric(bt_loop)
    
    for (b in seq_len(bt_loop)) {
      
      # gamma* (pz x px)
      Zg <- matrix(rnorm(pz * px), pz, px)
      gamma_star <- Zg %*% t(chol_gg)
      
      # Gamma* = C^(r) gamma*  (pz x py)
      Gamma_star <- gamma_star %*% t(Cr_hat)
      
      # epsX, epsY
      ZX <- matrix(rnorm(pz * px), pz, px)
      ZY <- matrix(rnorm(pz * py), pz, py)
      epsX <- ZX %*% t(chol_X)
      epsY <- ZY %*% t(chol_Y)
      
      # hat stats under null
      x_star <- gamma_star + epsX
      y_star <- Gamma_star + epsY
      
      # refit C0* and Cr* on bootstrap sample
      C0_star <- .fit_C0(y_star, x_star, W = W, Sigma_X = Sigma_X)
      Cr_star <- .fit_Cr(y_star, x_star, r = r, W = W, Sigma_X = Sigma_X)
      
      # bootstrap statistic
      D_star <- C0_star - Cr_star
      T_star[b] <- sum(D_star^2)
    }
    
    # Step 5: critical value and decision
    q_crit <- unname(quantile(T_star, probs = 1 - alpha, type = 8, names = FALSE))
    reject <- (T_obs > q_crit)
    
    out$T_obs[idx] <- T_obs
    out$q_crit[idx] <- q_crit
    out$reject[idx] <- reject
  }
  
  # Rank selection: smallest r not rejected
  idx_keep <- which(out$reject == FALSE)[1]
  if (is.na(idx_keep)) {
    r_sel <- k
    if (print) {
      cat("All ranks rejected at alpha =", alpha, "; default to full rank:", k, "\n")
    }
  } else {
    r_sel <- out$r[idx_keep]
    if (print) {
      cat("Selected rank =", r_sel,
          " (first not rejected), with T_obs =", signif(out$T_obs[idx_keep], 4),
          " and q_crit =", signif(out$q_crit[idx_keep], 4), "\n")
    }
  }
  
  if (return_details) {
    return(list(r_sel = r_sel,
                table = out,
                C0_hat = C0_hat,
                Sigma_gg_hat = Sigma_gg_hat,
                Sigma_Y_hat = Sigma_Y_hat))
  } else {
    return(r_sel)
  }
}


# old versions
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


rank_test_M0 = function(W, y_j_hat, x_j_hat, Sigma_X, print = TRUE, min_rank = 1){
  px = ncol(x_j_hat)
  py = ncol(y_j_hat)
  pz = nrow(y_j_hat)
  sigmaxy = crossprod(x_j_hat, y_j_hat) / pz
  debiased_Sigma_xx = crossprod(x_j_hat) / pz - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  W_sqrt = .sqrt_matrix(W)
  V = W_sqrt %*% t(sigmaxy) %*% debiased_Sigma_xx_inv %*% sigmaxy %*% W_sqrt
  lambda = eigen(V)$values
  
  p_value = c()
  for (r in min_rank:(min(px,py)-1)) {
    log_sum_tail <- sum(log(1 + lambda[(r + 1):length(lambda)]))
    M = (pz-(px+py+1)/2) * log_sum_tail
    # 在卡方分布的百分比
    p_value = c(p_value, 1 - pchisq(M, df = (py - r) * (px - r)))
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
    return(res_idx)
  }
}


#### generate low rank true C ####
set.seed(123)
px = 9
py = 3
r_approx = 2
r =  min(px,py) # maximum rank of true C

# U = matrix(rnorm(py*r), py, r)
# V = matrix(rnorm(px*r), r, px)
# # assume true C to be rank r=3 for now. C=C^(r)
# # use modified THM 2.1 to define true C^(r), (it is equivalent to define it as A_d*B_d)
# eignvalue_matrix = diag(c(rep(0.25, r_approx),
#                           rep(0, r-r_approx))) # set at 0.001 if want C to be rank > r
# C <- U %*% eignvalue_matrix %*% V
# svd(C)$d
# #mean(abs(C))
# norm(C, "F")


C = matrix(rnorm(py*px), py, px)
U = svd(C)$u
V = svd(C)$v
C = U %*% diag(c(rep(1, r_approx),
                 rep(0, r-r_approx))) %*% t(V)

svd(C)$d
norm(C, "F")


#### set universal parameters ####
me_weight_list = c(2.5, 1)
effect_weight_list = c(0.25, 1)
# save_sim_filename = "results/simulate_result_pred_250826.RData"
# regularization_rate_list = c(1.422259e-11, 3.600086e-12, 1.335627e-12, 3.706175e-14)

#### choose regularization rate \lambda ####
# 0.3 15 [1] 1.502516e-10 good 3.716852e-13 1.523810e-12 3.777015e-17

# 0.75 10 [1] 2.900664e-10 1.032470e-17 4.833967e-12 1.643551e-32
# 0.5 12 [1] 1.977613e-10 2.100946e-15 4.502632e-12 1.265314e-23


# 0.3 26 [1] 3.239274e-09 3.249018e-12 1.050097e-10 9.435030e-18

# 0.25 [1] 1.476369e-09 4.629265e-12 5.508758e-11 6.620537e-17
# 0.3 [1] 4.372563e-09 4.385715e-12 1.417482e-10 1.273596e-17

# 0-35，0.19，0.1
# [1] 1.420287e-09 2.559184e-11 3.420459e-11 1.740193e-14
# c=0.3, 27, [1] 1.728388e-09 3.852705e-12 1.635144e-11 3.624514e-16
# 1.635144e-11 too small
# c=1/3, 26, [1] 2.424390e-09 2.962031e-12 1.779574e-11 1.555238e-16
# c=0.5 [1] 2.000784e-09 1.442122e-13 5.150635e-12 3.806087e-19
# 2.000784e-09 good can not be bigger; 5.150635e-12 too small
# [1] 9.585637e-09 4.048477e-14 7.567364e-11 6.649455e-24
# comment: 9.585637e-09 too big. 7.567364e-11 good, can be bigger.
set.seed(123)

eloop = 100
result_matrix = min_eigen_matrix = matrix(NA, nrow = eloop, ncol = length(me_weight_list)*length(effect_weight_list))
for (i in 1:eloop){
  result_list = list()
  
  # diagnose purpose: record min eigen values of Sigma_gammagamma_hat each simulation
  min_eigen_list = list() 
  
  for (me_weight in me_weight_list){
    for (effect_weight in effect_weight_list){
      # #test
      # me_weight = 2.5
      # effect_weight = 0.25
  
      parameters = .get_parameters(C, me_weight=me_weight, effect_weight=effect_weight, r_RR = 2)
      n = 1000
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
      
      Sigma_hgammahgamma_hat = cov(x_j_hat)
      Sigma_gammagamma_hat = Sigma_hgammahgamma_hat - parameters$Sigma_X
      mu_min = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% Sigma_gammagamma_hat %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
      
      min_eigen_list = c(min_eigen_list,min(eigen(Sigma_gammagamma_hat)$values))
      # min_eigen_list = c(min_eigen_list,max(eigen(Sigma_gammagamma_hat)$values))
      
      D_list = seq(0, 15, 1)
      regu_rate_list = c()
      obj_value = c()
      for (D in D_list) {
        sigma_y2 = mean(eigen(parameters$Sigma_X)$values)
        
        # regu_rate = sigma_y2^2 * exp(c * (q/sqrt(n)-(mu_min+1)))/n
        regu_rate = sigma_y2^2 * exp(0.3*(D - sqrt(n) * mu_min))/n
      
        regu_rate_list = c(regu_rate_list, regu_rate)
        result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate)
        C_hat = result_d$AB
        # objective function value
        obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
        # Q 2025: add debias term
        debias_term = sum(diag(W_sqrt %*% C_hat %*% Sigma_X %*% t(C_hat) %*% W_sqrt)) # tarce
        obj = obj - debias_term
        obj_value = c(obj_value, obj)
      }
      
      opt_rate = regu_rate_list[order(obj_value)][1]
      
      result_list= c(result_list, opt_rate)
    }
  }
  result_matrix[i,] = unlist(result_list)
  min_eigen_matrix[i,] = unlist(min_eigen_list)
}
# average
result_matrix_mean = apply(result_matrix, 2, mean)

min_eigen_matrix_mean = apply(min_eigen_matrix,2,mean)

regularization_rate_list=result_matrix_mean

print('regularization rate list is: ')
print(result_matrix_mean)

# # 26.1.19 compute relative rate (devided by sigma 2)
# set.seed(123)
# 
# eloop = 100
# n_combo = length(me_weight_list) * length(effect_weight_list)
# 
# result_matrix = matrix(NA, nrow = eloop, ncol = n_combo)
# min_eigen_matrix = matrix(NA, nrow = eloop, ncol = n_combo)
# 
# # NEW: store sigma_y2 for each weight-combo in each replicate
# sigma_y2_matrix = matrix(NA, nrow = eloop, ncol = n_combo)
# 
# for (i in 1:eloop){
#   result_list = c()
#   min_eigen_list = c()
#   sigma_y2_list = c()   # NEW: collect sigma_y2 for this replicate across combos
#   
#   for (me_weight in me_weight_list){
#     for (effect_weight in effect_weight_list){
#       
#       parameters = .get_parameters(C, me_weight=me_weight, effect_weight=effect_weight, r_RR = 2)
#       n = 1000
#       py = parameters$py
#       px = parameters$px
#       Sigma_X = parameters$Sigma_X
#       Sigma_Y = parameters$Sigma_Y
#       C = parameters$C
#       r_RR = parameters$r_RR
#       W = parameters$weight.matrix
#       W_sqrt = .sqrt_matrix(parameters$weight.matrix)
#       
#       # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
#       gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = parameters$VX_tilde, tol = 100)
#       Gamma_j_star = gamma_j_star %*% t(C)
#       
#       # sample x_j_hat, y_j_hat
#       x_j_hat = matrix(0, n, px)
#       y_j_hat = matrix(0, n, py)
#       for (j in 1:n) {
#         x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#         y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
#       }
#       
#       Sigma_hgammahgamma_hat = cov(x_j_hat)
#       Sigma_gammagamma_hat = Sigma_hgammahgamma_hat - parameters$Sigma_X
#       mu_min = min(
#         eigen(
#           solve(.sqrt_matrix(parameters$Sigma_X)) %*% Sigma_gammagamma_hat %*% solve(.sqrt_matrix(parameters$Sigma_X))
#         )$values
#       )
#       
#       min_eigen_list = c(min_eigen_list, min(eigen(Sigma_gammagamma_hat)$values))
#       
#       # NEW: compute sigma_y2 ONCE per (me_weight, effect_weight) combo
#       sigma_y2 = mean(eigen(parameters$Sigma_X)$values)
#       sigma_y2_list = c(sigma_y2_list, sigma_y2)
#       
#       D_list = seq(0, 15, 1)
#       regu_rate_list = c()
#       obj_value = c()
#       
#       for (D in D_list) {
#         # keep original formula unchanged
#         regu_rate = sigma_y2^2 * exp(0.3*(D - sqrt(n) * mu_min))/n
#         
#         regu_rate_list = c(regu_rate_list, regu_rate)
#         result_d <- mr_rr_regularized(
#           y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate
#         )
#         C_hat = result_d$AB
#         
#         # objective function value
#         obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
#         # Q 2025: add debias term
#         debias_term = sum(diag(W_sqrt %*% C_hat %*% Sigma_X %*% t(C_hat) %*% W_sqrt)) # trace
#         obj = obj - debias_term
#         obj_value = c(obj_value, obj)
#       }
#       
#       opt_rate = regu_rate_list[order(obj_value)][1]
#       result_list = c(result_list, opt_rate)
#     }
#   }
#   
#   result_matrix[i,] = result_list
#   min_eigen_matrix[i,] = min_eigen_list
#   sigma_y2_matrix[i,] = sigma_y2_list   # NEW
# }
# 
# # average opt rate per combo
# result_matrix_mean = apply(result_matrix, 2, mean)
# min_eigen_matrix_mean = apply(min_eigen_matrix, 2, mean)
# 
# # NEW: average sigma_y2 per combo
# sigma_y2_mean = apply(sigma_y2_matrix, 2, mean)
# 
# regularization_rate_list = result_matrix_mean
# 
# # NEW: relative regularization rate
# relative_rate_list = regularization_rate_list / sigma_y2_mean
# # If you instead want to normalize by sigma_y2^2, use:
# # relative_rate_list = regularization_rate_list / (sigma_y2_mean^2)
# 
# print('regularization rate list (mean opt_rate) is: ')
# print(regularization_rate_list)
# 
# print('sigma_y2 mean (per weight-combo) is: ')
# print(sigma_y2_mean)
# 
# print('relative regularization rate list is: ')
# print(relative_rate_list)



# order: me 2.5, effect 0.25; me 2.5, effect 1; me 1, effect 0.25; me 1, effect 1
# print("mean min eigen for sigma gammagamma hat:")
# print(min_eigen_matrix_mean)

regularization_rate_list = result_matrix_mean

#### prediction problem (can only run this instead of the standard sim. it included the standard sim) ####
set.seed(123)
# generate an X for some individual
X = rnorm(px, mean = 0, sd = 1)
Y = C %*% X
round(X,3)

# new: -0.560 -0.230  1.559  0.071  0.129  1.715  0.461 -1.265 -0.687

# # iv strength on x direction
# iv_strength_x = t(X) %*% solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)) %*% X / (t(X) %*% X)
# round(iv_strength_x,2)


#### main ####

set.seed(123)
save_sim_pred_filename = "results/simulate_result_pred_260214_1000_newC_testbt.RData"

run_simulation_prediction <- function(regularization_rate_list, eloop = 500){
  comb_names <- sprintf("me_%s_effect_%s", rep(me_weight_list, each = 2), effect_weight_list)
  template_list <- setNames(vector("list", length(comb_names)), comb_names)
  
  # bias
  result_AB_list = result_AB_d_list = result_AB_d_r_list = 
    result_C_ivw_list = result_C_adivw_list = result_MrDAG_list = # result_C_GRAPPLE_list =
    # C hat
    AB_list = AB_d_list = AB_d_r_list = 
    C_ivw_list = C_adivw_list = MrDAG_list =
    # C hat %*% X = Y_pred
    Y_pred_AB_list = Y_pred_AB_d_list = Y_pred_AB_d_r_list =
    Y_pred_C_ivw_list = Y_pred_C_adivw_list = Y_pred_MrDAG_list = # Y_pred_C_GRAPPLE_list =
    # rank correct
    rank_correct_list = selected_rank_list = template_list
  
  iv_strength_list = c()
  parameters_list = c()
  for (me_weight_index in seq_along(me_weight_list)){
    me_weight = me_weight_list[me_weight_index]
    for (effect_weight_index in seq_along(effect_weight_list)){
      effect_weight = effect_weight_list[effect_weight_index]
      
      # # test 250826
      # me_weight = me_weight_list[1]
      # effect_weight = effect_weight_list[1]
      # me_weight_index = 1
      # effect_weight_index = 1
      
      
      regularization_rate = regularization_rate_list[.get_sim_index(me_weight_index, effect_weight_index)]
      sim_name = paste0("me_", me_weight, "_effect_", effect_weight)
      
      parameters = .get_parameters(C, me_weight, effect_weight, r_RR = 2)
      C = parameters$C
      A = parameters$A
      B = parameters$B
      A_d = parameters$A_d
      B_d = parameters$B_d
      C_r = parameters$C_r
      py = parameters$py
      px = parameters$px
      
      iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values) * sqrt(1000) # Q 2025: times sqrt(pz) while presenting IV strength 
      iv_strength_list = c(iv_strength_list, iv_strength)
      
      parameters_list = c(parameters_list, list(c(parameters, iv_strength = iv_strength)))
      
      bias_AB_matrix = bias_AB_d_matrix = bias_AB_d_r_matrix = 
        bias_C_ivw_matrix = bias_C_adivw_matrix = bias_C_MrDAG_matrix =
        # bias_C_GRAPPLE_matrix =
        AB_matrix = AB_d_matrix = AB_d_r_matrix = 
        C_ivw_matrix = C_adivw_matrix = C_MrDAG_matrix =
        matrix(NA, nrow = py*px, ncol = eloop)
      
      Y_pred_AB_matrix = Y_pred_AB_d_matrix = Y_pred_AB_d_r_matrix =
        Y_pred_C_ivw_matrix = Y_pred_C_adivw_matrix = Y_pred_MrDAG_matrix = # Y_pred_C_GRAPPLE_matrix =
        matrix(NA, nrow = py, ncol = eloop)
      
      rank_correct = rep(NA, eloop)
      selected_rank = rep(NA, eloop)
      
      B_star = t(B) # in order to compute the norm
      B_d_star = t(B_d)
      for (i in 1:eloop) {
        simulation_result = .simulation(parameters, 
                                        regularization_rate=regularization_rate)
        simulation_result[[length(simulation_result)]]
        # if rank correct
        rank_correct[i] = simulation_result[[length(simulation_result)-1]]
        selected_rank[i] = simulation_result[[length(simulation_result)]]
        
        A_hat = simulation_result[[1]]
        B_hat = simulation_result[[2]]
        B_hat_star = t(simulation_result[[2]])
        AB_hat = simulation_result[[3]]
        
        bias_AB_vectorized = as.vector(AB_hat - C)
        # bias_AB_vectorized = as.vector(AB_hat - C_r)
        bias_AB_matrix[,i] = bias_AB_vectorized
        AB_matrix[,i] = as.vector(AB_hat)
        Y_pred_AB_matrix[,i] = as.vector(AB_hat %*% X)
        
        A_d_hat = simulation_result[[4]]
        B_d_hat = simulation_result[[5]]
        B_d_hat_star = t(simulation_result[[5]])
        AB_d_hat = simulation_result[[6]]
        
        # bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
        bias_AB_d_vectorized = as.vector(AB_d_hat - C)
        bias_AB_d_matrix[,i] = bias_AB_d_vectorized
        AB_d_matrix[,i] = as.vector(AB_d_hat)
        Y_pred_AB_d_matrix[,i] = as.vector(AB_d_hat %*% X)
        
        A_d_r_hat = simulation_result[[7]]
        B_d_r_hat = simulation_result[[8]]
        B_d_r_hat_star = t(simulation_result[[8]])
        AB_d_r_hat = simulation_result[[9]]
        
        # bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
        bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C)
        bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
        AB_d_r_matrix[,i] = as.vector(AB_d_r_hat)
        Y_pred_AB_d_r_matrix[,i] = as.vector(AB_d_r_hat %*% X)
        
        C_ivw = simulation_result[[12]]
        bias_C_ivw_vectorized = as.vector(C_ivw - C)
        bias_C_ivw_matrix[,i] = bias_C_ivw_vectorized
        C_ivw_matrix[,i] = as.vector(C_ivw)
        Y_pred_C_ivw_matrix[,i] = as.vector(C_ivw %*% X)
        
        C_adivw = simulation_result[[13]]
        bias_C_adivw_vectorized = as.vector(C_adivw - C)
        bias_C_adivw_matrix[,i] = bias_C_adivw_vectorized
        C_adivw_matrix[,i] = as.vector(C_adivw)
        Y_pred_C_adivw_matrix[,i] = as.vector(C_adivw %*% X)
        
        C_MrDAG = simulation_result[[14]]
        bias_C_MrDAG_vectorized = as.vector(C_MrDAG - C)
        bias_C_MrDAG_matrix[,i] = bias_C_MrDAG_vectorized
        C_MrDAG_matrix[,i] = as.vector(C_MrDAG)
        Y_pred_MrDAG_matrix[,i] = as.vector(C_MrDAG %*% X)
        
        # C_GRAPPLE = simulation_result[[15]]
        # bias_C_GRAPPLE_vectorized = as.vector(C_GRAPPLE - C)
        # bias_C_GRAPPLE_matrix[,i] = bias_C_GRAPPLE_vectorized
        
        # print "eloop" and eloop number in a line
        print(paste("eloop:", i, "me_weight:", me_weight, "effect_weight:", effect_weight))
      }
      
      # store the proportion of correctly specified rank
      rank_correct_list[[sim_name]] = prop.table(table(rank_correct))
      
      # store the proprtion of selected rank
      selected_rank_list[[sim_name]] = table(selected_rank)
      
      # store estimator bias
      result_AB_list[[sim_name]] = bias_AB_matrix
      result_AB_d_list[[sim_name]] = bias_AB_d_matrix
      result_AB_d_r_list[[sim_name]] = bias_AB_d_r_matrix
      result_C_ivw_list[[sim_name]] = bias_C_ivw_matrix
      result_C_adivw_list[[sim_name]] = bias_C_adivw_matrix
      result_MrDAG_list[[sim_name]] = bias_C_MrDAG_matrix
      # result_C_GRAPPLE_list[[sim_name]] = bias_C_GRAPPLE_matrix
      # store estimator value
      AB_list[[sim_name]] = AB_matrix
      AB_d_list[[sim_name]] = AB_d_matrix
      AB_d_r_list[[sim_name]] = AB_d_r_matrix
      C_ivw_list[[sim_name]] = C_ivw_matrix
      C_adivw_list[[sim_name]] = C_adivw_matrix
      MrDAG_list[[sim_name]] = C_MrDAG_matrix
      # store Y hat prediction value
      Y_pred_AB_list[[sim_name]] = Y_pred_AB_matrix
      Y_pred_AB_d_list[[sim_name]] = Y_pred_AB_d_matrix
      Y_pred_AB_d_r_list[[sim_name]] = Y_pred_AB_d_r_matrix
      Y_pred_C_ivw_list[[sim_name]] = Y_pred_C_ivw_matrix
      Y_pred_C_adivw_list[[sim_name]] = Y_pred_C_adivw_matrix
      Y_pred_MrDAG_list[[sim_name]] = Y_pred_MrDAG_matrix
    }
  }
  
  return(list(result_AB_list=result_AB_list, 
              result_AB_d_list=result_AB_d_list, 
              result_AB_d_r_list=result_AB_d_r_list,
              iv_strength_list=iv_strength_list, parameters_list=parameters_list,
              result_C_ivw_list=result_C_ivw_list, 
              result_C_adivw_list=result_C_adivw_list,
              result_MrDAG_list=result_MrDAG_list,
              # result_C_GRAPPLE_list=result_C_GRAPPLE_list
              AB_list=AB_list,
              AB_d_list=AB_d_list,
              AB_d_r_list=AB_d_r_list,
              C_ivw_list=C_ivw_list,
              C_adivw_list=C_adivw_list,
              MrDAG_list=MrDAG_list,
              # prediction
              Y_pred_AB_list=Y_pred_AB_list,
              Y_pred_AB_d_list=Y_pred_AB_d_list,
              Y_pred_AB_d_r_list=Y_pred_AB_d_r_list,
              Y_pred_C_ivw_list=Y_pred_C_ivw_list,
              Y_pred_C_adivw_list=Y_pred_C_adivw_list,
              Y_pred_MrDAG_list=Y_pred_MrDAG_list,
              # correct rank
              rank_correct_list = rank_correct_list,
              selected_rank_list = selected_rank_list
  ))
}

simulate_result_prediction = run_simulation_prediction(regularization_rate_list = regularization_rate_list, 
                                                       eloop = 1000)
save(simulate_result_prediction, file = save_sim_pred_filename)
round(simulate_result_prediction$iv_strength_list, 3)

load(save_sim_pred_filename)

# # use this result:
# test_sim_pred_filename = "results/simulate_result_pred_260210.RData"
# load(test_sim_pred_filename)
# use this one, perform better. 0.3 15 100loop


#### 1106 with sparse ####
set.seed(123)
.simulation_sparse <- function(parameters, sparse_lambda = rep(2e-3, 9)) {
  n = 1000
  
  # 获取参数
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  C = parameters$C
  r_RR = parameters$r_RR
  W = parameters$weight.matrix
  
  # 生成数据
  gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }
  
  # 计算Sparse MR-RR
  result_sparse <- tryCatch({
    mr_rr_sparse(
      GAMMA_hat = y_j_hat, 
      gamma_hat = x_j_hat, 
      W = W, 
      Sigma_X = Sigma_X, 
      lambda = sparse_lambda,
      r = r_RR,
      max_iter = 100,
      tol = 1e-2
    )
  }, error = function(e) {
    list(
      A = matrix(NA, py, r_RR), 
      B = matrix(NA, r_RR, px), 
      AB = matrix(NA, py, px),
      error = e$message
    )
  })
  
  # 检查是否成功
  success = !any(is.na(result_sparse$AB))
  
  return(list(
    C_sparse = result_sparse$AB,
    B_sparse = result_sparse$B,
    A_sparse = result_sparse$A,
    success = success,
    error = if(!success) result_sparse$error else NULL
  ))
}


run_simulation_sparse_prediction <- function(regularization_rate_list, 
                                             sparse_lambda_list = NULL,
                                             eloop = 500) {
  comb_names <- sprintf("me_%s_effect_%s", rep(me_weight_list, each = 2), effect_weight_list)
  template_list <- setNames(vector("list", length(comb_names)), comb_names)
  
  # 只需要sparse相关的结果
  result_C_sparse_list = C_sparse_list = Y_pred_C_sparse_list = 
    B_sparse_list = sparse_success_rate_list = template_list
  
  iv_strength_list = c()
  parameters_list = c()
  
  for (me_weight_index in seq_along(me_weight_list)) {
    me_weight = me_weight_list[me_weight_index]
    
    for (effect_weight_index in seq_along(effect_weight_list)) {
      effect_weight = effect_weight_list[effect_weight_index]
      
      regularization_rate = regularization_rate_list[.get_sim_index(me_weight_index, effect_weight_index)]
      sim_name = paste0("me_", me_weight, "_effect_", effect_weight)
      
      # 获取参数
      parameters = .get_parameters(C, me_weight, effect_weight, r_RR = 2)
      C = parameters$C
      A = parameters$A
      B = parameters$B
      py = parameters$py
      px = parameters$px
      
      # 计算IV strength
      iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% 
                                parameters$VX_tilde %*% 
                                solve(.sqrt_matrix(parameters$Sigma_X)))$values) * sqrt(1000)
      iv_strength_list = c(iv_strength_list, iv_strength)
      parameters_list = c(parameters_list, list(c(parameters, iv_strength = iv_strength)))
      
      # 设置sparse lambda
      if (is.null(sparse_lambda_list)) {
        sparse_lambda = rep(2e-3, px)
      } else {
        sparse_lambda = sparse_lambda_list[[sim_name]]
      }
      
      # 初始化结果矩阵
      bias_C_sparse_matrix = C_sparse_matrix = matrix(NA, nrow = py * px, ncol = eloop)
      Y_pred_C_sparse_matrix = matrix(NA, nrow = py, ncol = eloop)
      B_sparse_matrix = matrix(NA, nrow = 2 * px, ncol = eloop)  # r=2的情况
      sparse_success = rep(FALSE, eloop)
      
      for (i in 1:eloop) {
        # 生成数据
        simulation_result = .simulation_sparse(
          parameters = parameters,
          sparse_lambda = sparse_lambda
        )
        
        # 提取结果
        C_sparse = simulation_result$C_sparse
        B_sparse = simulation_result$B_sparse
        success = simulation_result$success
        
        # 记录是否成功
        sparse_success[i] = success
        
        if (success) {
          # 计算bias
          bias_C_sparse_matrix[, i] = as.vector(C_sparse - C)
          
          # 存储C估计值
          C_sparse_matrix[, i] = as.vector(C_sparse)
          
          # 存储B估计值（用于检查sparsity）
          B_sparse_matrix[, i] = as.vector(B_sparse)
          
          # 计算prediction
          Y_pred_C_sparse_matrix[, i] = as.vector(C_sparse %*% X)
        }
        
        print(paste("Sparse simulation - eloop:", i, 
                    "me_weight:", me_weight, 
                    "effect_weight:", effect_weight,
                    "success:", success))
      }
      
      # 存储结果
      result_C_sparse_list[[sim_name]] = bias_C_sparse_matrix
      C_sparse_list[[sim_name]] = C_sparse_matrix
      Y_pred_C_sparse_list[[sim_name]] = Y_pred_C_sparse_matrix
      B_sparse_list[[sim_name]] = B_sparse_matrix
      sparse_success_rate_list[[sim_name]] = mean(sparse_success)
      
      cat(sprintf("\n%s: Success rate = %.2f%%\n", 
                  sim_name, mean(sparse_success) * 100))
    }
  }
  
  return(list(
    result_C_sparse_list = result_C_sparse_list,
    C_sparse_list = C_sparse_list,
    Y_pred_C_sparse_list = Y_pred_C_sparse_list,
    B_sparse_list = B_sparse_list,
    sparse_success_rate_list = sparse_success_rate_list,
    iv_strength_list = iv_strength_list,
    parameters_list = parameters_list
  ))
}

sparse_results = run_simulation_sparse_prediction(
  regularization_rate_list = regularization_rate_list,
  sparse_lambda_list = NULL,  # 使用默认值 rep(2e-3, 9)
  eloop = 1000
)

save(sparse_results, file = "results/simulate_result_pred_sparse_1000_260210.RData")


#### plot ####
# plot main boxplot
plot_simulation_result <- function(sim_result, 
                                   me_weight, effect_weight,
                                   true_C,
                                   type = c("box", "ci"),
                                   width = 16, height = 8) {
  type <- match.arg(type)
  
  sim_name = paste0("me_", me_weight, "_effect_", effect_weight)
  today_str <- format(Sys.Date(), "%y%m%d")
  filename = sprintf("results/lessvar_plot_%s_type_%s_%s.png", sim_name, type, today_str)
  
  estimator_names <- c("C_ivw_list", "C_adivw_list", "AB_list", "AB_d_list", "AB_d_r_list", "MrDAG_list")
  pos_ids <- expand.grid(Row = 0:2, Col = 0:8)
  
  # ===== 整理数据 =====
  if (type == "box") {
    df_long <- purrr::map_dfr(estimator_names, function(est) {
      matrix_data <- sim_result[[est]][[sim_name]]
      n_sim <- ncol(matrix_data)
      
      purrr::map_dfr(1:(px*py), function(i) {
        tibble::tibble(
          Estimator = est,
          Row = pos_ids$Row[i],
          Col = pos_ids$Col[i],
          Value = matrix_data[i, ]
        )
      })
    })
    
    df_long$Estimator <- factor(df_long$Estimator,
                                levels = estimator_names,
                                labels = c("IVW", "adIVW", "Naive MR-rr", "MR-rr", "regularized MR-rr", "MrDAG"))
  } else {
    df_long <- purrr::map_dfr(estimator_names, function(est) {
      matrix_data <- sim_result[[est]][[sim_name]]
      
      purrr::map_dfr(1:(px*py), function(i) {
        vals <- matrix_data[i, ]
        tibble::tibble(
          Estimator = est,
          Row = pos_ids$Row[i],
          Col = pos_ids$Col[i],
          Median = median(vals, na.rm = TRUE),
          Lower = quantile(vals, 0.025, na.rm = TRUE),
          Upper = quantile(vals, 0.975, na.rm = TRUE)
        )
      })
    })
    
    df_long$Estimator <- factor(df_long$Estimator,
                                levels = estimator_names,
                                labels = c("IVW", "adIVW", "Naive MR-rr", "MR-rr", "regularized MR-rr", "MrDAG"))
  }
  
  # ===== true C 横线数据 =====
  df_true <- expand.grid(Row = 0:2, Col = 0:8) %>%
    dplyr::mutate(TrueValue = as.vector(true_C))
  
  # ===== 图层构建 =====
  plot_layers <- if (type == "box") {
    list(
      geom_boxplot(aes(y = Value, fill = Estimator), outlier.size = 0.3),
      geom_hline(data = df_true, aes(yintercept = TrueValue),
                 color = "red", linetype = "dashed", linewidth = 0.4,
                 inherit.aes = FALSE, show.legend = FALSE)
    )
  } else {
    list(
      geom_point(aes(y = Median, color = Estimator), position = position_dodge(width = 0.5)),
      geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Estimator),
                    width = 0.2, position = position_dodge(width = 0.5)),
      geom_hline(data = df_true, aes(yintercept = TrueValue),
                 color = "red", linetype = "dashed", linewidth = 0.4,
                 inherit.aes = FALSE, show.legend = FALSE)
    )
  }
  
  # ===== 绘图 =====
  p <- ggplot(df_long, aes(x = Estimator)) +
    plot_layers +
    facet_grid(
      rows = vars(Row), 
      cols = vars(Col),
      labeller = labeller(Row = out_name, Col = exp_name)
    ) +
    coord_cartesian(ylim = c(-0.5, 0.5)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 8),
          legend.position = "right") +
    ggtitle(paste0(
      if (type == "box") {
        "Simulation Boxplots of C Estimates (Red Dashed Line = True C)"
      } else {
        "Median and 95% CI of C Estimates (Red Dashed Line = True C)"
      },
      "\nmeasurement error weighted by ", me_weight, ", true effect weighted by ", effect_weight
    )) +
    ylab("Estimated Value") +
    xlab("")
  
  # 保存图
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}


plot_simulation_result(sim_result = simulate_result_prediction,
                       me_weight = "1",
                       effect_weight = "1",
                       true_C = C,
                       type = "box",
                       width = 10, height = 5)


# pred plot with sparse
plot_simulation_result_pred <- function(sim_result, 
                                        me_weight, effect_weight,
                                        Y,  # true Y = C %*% X
                                        sparse_result = NULL,  # 可选：sparse simulation 结果
                                        type = c("box", "ci"),
                                        width = 12, height = 6) {
  type <- match.arg(type)
  
  sim_name = paste0("me_", me_weight, "_effect_", effect_weight)
  today_str <- format(Sys.Date(), "%y%m%d")
  filename <- sprintf("results/pred_plot_with_sparse_%s_type_%s_%s.png", sim_name, type, today_str)
  
  estimator_names <- c("Y_pred_C_ivw_list", "Y_pred_C_adivw_list", 
                       "Y_pred_AB_list", "Y_pred_AB_d_list", 
                       "Y_pred_AB_d_r_list", "Y_pred_MrDAG_list")
  estimator_labels <- c("IVW", "adIVW", "Naive MR-rr", "MR-rr", 
                        "regularized MR-rr", "MrDAG")
  
  py <- nrow(Y)
  pos_ids <- data.frame(Row = 0:(py - 1))  # 每行是一个 Y_k
  
  # ===== 整理数据 =====
  if (type == "box") {
    df_long <- purrr::map2_dfr(estimator_names, estimator_labels, function(est, label) {
      matrix_data <- sim_result[[est]][[sim_name]]
      n_sim <- ncol(matrix_data)
      purrr::map_dfr(1:py, function(i) {
        tibble::tibble(
          Estimator = label,
          Row = pos_ids$Row[i],
          Value = matrix_data[i, ]
        )
      })
    })
    
    # 加入 sparse estimator 的预测结果
    if (!is.null(sparse_result)) {
      matrix_data <- sparse_result$Y_pred_C_sparse_list[[sim_name]]
      n_sim <- ncol(matrix_data)
      df_sparse <- purrr::map_dfr(1:py, function(i) {
        tibble::tibble(
          Estimator = "Sparse MR-rr",
          Row = pos_ids$Row[i],
          Value = matrix_data[i, ]
        )
      })
      df_long <- bind_rows(df_long, df_sparse)
    }
    
    # ✅ 设置显示顺序
    df_long$Estimator <- factor(df_long$Estimator, levels = c(
      "IVW", "adIVW", "Naive MR-rr", "MR-rr", 
      "regularized MR-rr", "MrDAG", "Sparse MR-rr"
    ))
    
  } else {
    df_long <- purrr::map2_dfr(estimator_names, estimator_labels, function(est, label) {
      matrix_data <- sim_result[[est]][[sim_name]]
      purrr::map_dfr(1:py, function(i) {
        vals <- matrix_data[i, ]
        tibble::tibble(
          Estimator = label,
          Row = pos_ids$Row[i],
          Median = median(vals, na.rm = TRUE),
          Lower = quantile(vals, 0.025, na.rm = TRUE),
          Upper = quantile(vals, 0.975, na.rm = TRUE)
        )
      })
    })
    
    # 加入 sparse estimator 的 CI
    if (!is.null(sparse_result)) {
      matrix_data <- sparse_result$Y_pred_C_sparse_list[[sim_name]]
      df_sparse <- purrr::map_dfr(1:py, function(i) {
        vals <- matrix_data[i, ]
        tibble::tibble(
          Estimator = "Sparse MR-rr",
          Row = pos_ids$Row[i],
          Median = median(vals, na.rm = TRUE),
          Lower = quantile(vals, 0.025, na.rm = TRUE),
          Upper = quantile(vals, 0.975, na.rm = TRUE)
        )
      })
      df_long <- bind_rows(df_long, df_sparse)
    }
    
    # ✅ 设置显示顺序
    df_long$Estimator <- factor(df_long$Estimator, levels = c(
      "IVW", "adIVW", "Naive MR-rr", "MR-rr", 
      "regularized MR-rr", "MrDAG", "Sparse MR-rr"
    ))
  }
  
  # ===== true Y 横线数据 =====
  df_true <- tibble::tibble(Row = 0:(py - 1), TrueValue = as.vector(Y))
  
  # ===== 图层构建 =====
  plot_layers <- if (type == "box") {
    list(
      geom_boxplot(aes(y = Value, fill = Estimator), outlier.size = 0.3),
      geom_hline(data = df_true, aes(yintercept = TrueValue),
                 color = "red", linetype = "dashed", linewidth = 0.8,
                 inherit.aes = FALSE, show.legend = FALSE)
    )
  } else {
    list(
      geom_point(aes(y = Median, color = Estimator), position = position_dodge(width = 0.5)),
      geom_errorbar(aes(ymin = Lower, ymax = Upper, color = Estimator),
                    width = 0.2, position = position_dodge(width = 0.5)),
      geom_hline(data = df_true, aes(yintercept = TrueValue),
                 color = "red", linetype = "dashed", linewidth = 0.8,
                 inherit.aes = FALSE, show.legend = FALSE)
    )
  }
  
  # ===== 绘图 =====
  p <- ggplot(df_long, aes(x = Estimator)) +
    plot_layers +
    facet_wrap(~ Row, ncol = 3, labeller = labeller(Row = out_name)) +
    coord_cartesian(ylim = c(-1, 1)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 10),
          legend.position = "right") +
    ggtitle(paste0(
      if (type == "box") {
        "Simulation Boxplots of Predicted Y (Red Dashed Line = True Y)"
      } else {
        "Median and 95% CI of Predicted Y (Red Dashed Line = True Y)"
      },
      "\nmeasurement error weighted by ", me_weight, 
      ", true effect weighted by ", effect_weight
    )) +
    ylab("Estimated Y value") +
    xlab("")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}



plot_simulation_result_pred(
  sim_result = simulate_result_prediction,
  sparse_result = sparse_results,  # 加入 sparse 结果
  me_weight = 2.5,
  effect_weight = 0.25,
  Y = Y,
  type = "box"
)

# generate table
gen_table <- function(simulate_result_prediction,
                      me_weight_to_plot, 
                      effect_weight_to_plot)
{
  parameters_list = simulate_result_prediction$parameters_list
  result_naive_MRrr_list = simulate_result_prediction$result_AB_list
  result_MRrr_list = simulate_result_prediction$result_AB_d_list
  result_r_MRrr_list = simulate_result_prediction$result_AB_d_r_list
  result_ivw_list = simulate_result_prediction$result_C_ivw_list
  result_adivw_list= simulate_result_prediction$result_C_adivw_list
  result_MrDAG_list = simulate_result_prediction$result_MrDAG_list
  
  px <- 9
  py <- 3
  sim_name <- paste0("me_", me_weight_to_plot, "_effect_", effect_weight_to_plot)
  
  ## --- 逐 entry 的 bias 与 sd（六种方法） ---
  # naive MR-rr
  abs_mean_entry_bias     <- mean_entry_bias     <- sd_entry_bias     <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(mean(result_naive_MRrr_list[[sim_name]][i,]))
    mean_entry_bias[i]     <- mean(result_naive_MRrr_list[[sim_name]][i,])
    sd_entry_bias[i]       <- sd(result_naive_MRrr_list[[sim_name]][i,])
  }
  
  # MR-rr
  abs_mean_entry_bias_d   <- mean_entry_bias_d   <- sd_entry_bias_d   <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(mean(result_MRrr_list[[sim_name]][i,]))
    mean_entry_bias_d[i]     <- mean(result_MRrr_list[[sim_name]][i,])
    sd_entry_bias_d[i]       <- sd(result_MRrr_list[[sim_name]][i,])
  }
  
  # MR-rr (regularized)
  abs_mean_entry_bias_d_r <- mean_entry_bias_d_r <- sd_entry_bias_d_r <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d_r[i] <- abs(mean(result_r_MRrr_list[[sim_name]][i,]))
    mean_entry_bias_d_r[i]     <- mean(result_r_MRrr_list[[sim_name]][i,])
    sd_entry_bias_d_r[i]       <- sd(result_r_MRrr_list[[sim_name]][i,])
  }
  
  # IVW
  abs_mean_entry_bias_ivw <- mean_entry_bias_ivw <- sd_entry_bias_ivw <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_ivw[i] <- abs(mean(result_ivw_list[[sim_name]][i,]))
    mean_entry_bias_ivw[i]     <- mean(result_ivw_list[[sim_name]][i,])
    sd_entry_bias_ivw[i]       <- sd(result_ivw_list[[sim_name]][i,])
  }
  
  # adIVW
  abs_mean_entry_bias_adivw <- mean_entry_bias_adivw <- sd_entry_bias_adivw <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_adivw[i] <- abs(mean(result_adivw_list[[sim_name]][i,]))
    mean_entry_bias_adivw[i]     <- mean(result_adivw_list[[sim_name]][i,])
    sd_entry_bias_adivw[i]       <- sd(result_adivw_list[[sim_name]][i,])
  }
  
  # MrDAG
  abs_mean_entry_bias_MrDAG <- mean_entry_bias_MrDAG <- sd_entry_bias_MrDAG <- rep(NA, px*py)
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_MrDAG[i] <- abs(mean(result_MrDAG_list[[sim_name]][i,]))
    mean_entry_bias_MrDAG[i]     <- mean(result_MrDAG_list[[sim_name]][i,])
    sd_entry_bias_MrDAG[i]       <- sd(result_MrDAG_list[[sim_name]][i,])
  }
  
  ## --- 组织成矩阵（py x px），便于后续整体统计 ---
  avg_bias      <- matrix(abs_mean_entry_bias,      nrow = px, ncol = py) |> t()
  avg_bias_d    <- matrix(abs_mean_entry_bias_d,    nrow = px, ncol = py) |> t()
  avg_bias_d_r  <- matrix(abs_mean_entry_bias_d_r,  nrow = px, ncol = py) |> t()
  avg_bias_ivw  <- matrix(abs_mean_entry_bias_ivw,  nrow = px, ncol = py) |> t()
  avg_bias_adivw<- matrix(abs_mean_entry_bias_adivw,nrow = px, ncol = py) |> t()
  avg_bias_MrDAG<- matrix(abs_mean_entry_bias_MrDAG,nrow = px, ncol = py) |> t()
  
  avg_sd        <- matrix(sd_entry_bias,      nrow = px, ncol = py) |> t()
  avg_sd_d      <- matrix(sd_entry_bias_d,    nrow = px, ncol = py) |> t()
  avg_sd_d_r    <- matrix(sd_entry_bias_d_r,  nrow = px, ncol = py) |> t()
  avg_sd_ivw    <- matrix(sd_entry_bias_ivw,  nrow = px, ncol = py) |> t()
  avg_sd_adivw  <- matrix(sd_entry_bias_adivw,nrow = px, ncol = py) |> t()
  avg_sd_MrDAG  <- matrix(sd_entry_bias_MrDAG,nrow = px, ncol = py) |> t()
  
  ## --- 与真 C_r 的规模对比 ---
  # 这里沿用你原来的取法：从 parameters_list 中取出对应情形的 C_r
  # 注意：这依赖于你已有的 .get_sim_index() / me_weight_list / effect_weight_list
  me_weight_index     <- match(me_weight_to_plot, me_weight_list)
  effect_weight_index <- match(effect_weight_to_plot, effect_weight_list)
  parameters <- parameters_list[[.get_sim_index(me_weight_index, effect_weight_index)]]
  C_r <- parameters$C_r
  
  ## --- 汇总数值（mean / median / IQR） ---
  mean_abs_C_entry <- round(mean(abs(C_r)), 3)
  
  mean_abs_bias_naive <- round(mean(abs(avg_bias)), 3)
  mean_abs_bias       <- round(mean(abs(avg_bias_d)), 3)
  mean_abs_bias_r     <- round(mean(abs(avg_bias_d_r)), 3)
  mean_abs_bias_ivw   <- round(mean(abs(avg_bias_ivw)), 3)
  mean_abs_bias_adivw <- round(mean(abs(avg_bias_adivw)), 3)
  mean_abs_bias_MrDAG <- round(mean(abs(avg_bias_MrDAG)), 3)
  
  mean_sd_naive <- round(mean(avg_sd), 3)
  mean_sd       <- round(mean(avg_sd_d), 3)
  mean_sd_r     <- round(mean(avg_sd_d_r), 3)
  mean_sd_ivw   <- round(mean(avg_sd_ivw), 3)
  mean_sd_adivw <- round(mean(avg_sd_adivw), 3)
  mean_sd_MrDAG <- round(mean(avg_sd_MrDAG), 3)
  
  median_abs_bias_naive <- round(median(abs(avg_bias)), 3)
  median_abs_bias       <- round(median(abs(avg_bias_d)), 3)
  median_abs_bias_r     <- round(median(abs(avg_bias_d_r)), 3)
  median_abs_bias_ivw   <- round(median(abs(avg_bias_ivw)), 3)
  median_abs_bias_adivw <- round(median(abs(avg_bias_adivw)), 3)
  median_abs_bias_MrDAG <- round(median(abs(avg_bias_MrDAG)), 3)
  
  median_sd_naive <- round(median(avg_sd), 3)
  median_sd       <- round(median(avg_sd_d), 3)
  median_sd_r     <- round(median(avg_sd_d_r), 3)
  median_sd_ivw   <- round(median(avg_sd_ivw), 3)
  median_sd_adivw <- round(median(avg_sd_adivw), 3)
  median_sd_MrDAG <- round(median(avg_sd_MrDAG), 3)
  
  quantile_abs_bias_naive <- round(quantile(abs(avg_bias),      c(0.25, 0.75)), 3)
  quantile_abs_bias       <- round(quantile(abs(avg_bias_d),    c(0.25, 0.75)), 3)
  quantile_abs_bias_r     <- round(quantile(abs(avg_bias_d_r),  c(0.25, 0.75)), 3)
  quantile_abs_bias_ivw   <- round(quantile(abs(avg_bias_ivw),  c(0.25, 0.75)), 3)
  quantile_abs_bias_adivw <- round(quantile(abs(avg_bias_adivw),c(0.25, 0.75)), 3)
  quantile_abs_bias_MrDAG <- round(quantile(abs(avg_bias_MrDAG),c(0.25, 0.75)), 3)
  
  quantile_sd_naive <- round(quantile(avg_sd,      c(0.25, 0.75)), 3)
  quantile_sd       <- round(quantile(avg_sd_d,    c(0.25, 0.75)), 3)
  quantile_sd_r     <- round(quantile(avg_sd_d_r,  c(0.25, 0.75)), 3)
  quantile_sd_ivw   <- round(quantile(avg_sd_ivw,  c(0.25, 0.75)), 3)
  quantile_sd_adivw <- round(quantile(avg_sd_adivw,c(0.25, 0.75)), 3)
  quantile_sd_MrDAG <- round(quantile(avg_sd_MrDAG,c(0.25, 0.75)), 3)
  
  ## --- 只返回数字 ---
  return(list(
    # mean_abs_C_entry = mean_abs_C_entry,
    # mean
    mean_abs_bias_ivw = mean_abs_bias_ivw, mean_sd_ivw = mean_sd_ivw,
    mean_abs_bias_adivw = mean_abs_bias_adivw, mean_sd_adivw = mean_sd_adivw,
    mean_abs_bias_naive = mean_abs_bias_naive, mean_sd_naive = mean_sd_naive,
    mean_abs_bias = mean_abs_bias, mean_sd = mean_sd,
    mean_abs_bias_r = mean_abs_bias_r, mean_sd_r = mean_sd_r,
    mean_abs_bias_MrDAG = mean_abs_bias_MrDAG, mean_sd_MrDAG = mean_sd_MrDAG,
    # median
    median_abs_bias_ivw = median_abs_bias_ivw, median_sd_ivw = median_sd_ivw,
    median_abs_bias_adivw = median_abs_bias_adivw, median_sd_adivw = median_sd_adivw,
    median_abs_bias_naive = median_abs_bias_naive, median_sd_naive = median_sd_naive,
    median_abs_bias = median_abs_bias, median_sd = median_sd,
    median_abs_bias_r = median_abs_bias_r, median_sd_r = median_sd_r,
    median_abs_bias_MrDAG = median_abs_bias_MrDAG, median_sd_MrDAG = median_sd_MrDAG,
    # IQR (0.25, 0.75)
    quantile_abs_bias_ivw = quantile_abs_bias_ivw, quantile_sd_ivw = quantile_sd_ivw,
    quantile_abs_bias_adivw = quantile_abs_bias_adivw, quantile_sd_adivw = quantile_sd_adivw,
    quantile_abs_bias_naive = quantile_abs_bias_naive, quantile_sd_naive = quantile_sd_naive,
    quantile_abs_bias = quantile_abs_bias, quantile_sd = quantile_sd,
    quantile_abs_bias_r = quantile_abs_bias_r, quantile_sd_r = quantile_sd_r,
    quantile_abs_bias_MrDAG = quantile_abs_bias_MrDAG, quantile_sd_MrDAG = quantile_sd_MrDAG
  ))
}


for (re in effect_weight_list){
  for (me in me_weight_list){
    print(me)
    print(re)
    print(gen_table(simulate_result_prediction,
              me_weight_to_plot = me, 
              effect_weight_to_plot = re)
    )
  }
}


