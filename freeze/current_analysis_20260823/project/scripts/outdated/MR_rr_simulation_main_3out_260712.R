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
# install.packages("Cairo")
# library(Cairo)


set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("D:/Users/YuexiangPeng/Documents/UW/Research/Ye Ting/sim_ArBr_bias")
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


# fix rank =2
.simulation <- function(parameters, regularization_rate = 1e-13) {
  n = 1000
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  C = parameters$C
  r_RR = parameters$r_RR
  W = parameters$weight.matrix
  
  # sample true effect
  gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  # sample x_j_hat, y_j_hat
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }
  
  # estimators (all with fixed r_RR)
  result     <- mr_rr_naive(y_j_hat, x_j_hat, r = r_RR, W = W)
  result_d   <- mr_rr(y_j_hat, x_j_hat, r = r_RR, W = W, Sigma_X = Sigma_X)
  result_d_r <- mr_rr_regularized(y_j_hat, x_j_hat, r = r_RR, W = W,
                                  Sigma_X = Sigma_X,
                                  regularization_rate = regularization_rate)
  result_ivw   <- ivw_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)
  result_adivw <- adivw_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)
  result_MrDAG <- Mr_DAG(y_j_hat, x_j_hat)
  
  return(list(
    result$A, result$B, result$AB,                          # 1-3: naive
    result_d$A, result_d$B, result_d$AB,                    # 4-6: MR-rr
    result_d_r$A, result_d_r$B, result_d_r$AB,              # 7-9: regularized
    x_j_hat, y_j_hat,                                       # 10-11
    result_ivw,                                              # 12
    result_adivw,                                            # 13
    result_MrDAG                                             # 14
  ))
}


.get_sim_index = function(me_weight_index, effect_weight_index, len_me = 2){
  sim_index = (me_weight_index-1)*len_me + effect_weight_index
  return(sim_index)
}


# ============================================================ #
# 3. functions for rank selection ####
# ============================================================ #
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


# Use this one eventrually in paper - 260603
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


# bootstrap old ver
rank_test_old <- function(W, y_j_hat, x_j_hat, Sigma_X,
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


# bootstrap new ver 260416
rank_test <- function(W, y_j_hat, x_j_hat, Sigma_X,
                         print = TRUE, min_rank = 1,
                         bt_loop = 500, alpha = 0.05) {
  px <- ncol(x_j_hat)
  py <- ncol(y_j_hat)
  pz <- nrow(y_j_hat)
  k  <- min(px, py)
  
  # 拟合 unconstrained (full-rank) 估计量，只需算一次
  C_hat_full <- mr_rr(y_j_hat, x_j_hat, r = k, W = W, Sigma_X = Sigma_X)$AB
  sv_full    <- svd(C_hat_full, nu = py, nv = px)
  
  out <- data.frame(r      = min_rank:(k - 1),
                    M1_obs = NA_real_,
                    q_crit = NA_real_,
                    reject = NA,
                    stringsAsFactors = FALSE)
  
  for (idx in seq_len(nrow(out))) {
    r <- out$r[idx]
    
    # ---- 1. 观测统计量：full-rank 估计量的尾部奇异值 ----
    M1_obs <- sum(sv_full$d[(r + 1):k]^2)
    
    # U2, V2 从 full-rank SVD 取正交补
    U2 <- sv_full$u[, (r + 1):py, drop = FALSE]
    V2 <- sv_full$v[, (r + 1):px, drop = FALSE]
    
    # ---- 2. Bootstrap ----
    M1_star <- numeric(bt_loop)
    for (b in 1:bt_loop) {
      idx_b <- sample.int(pz, pz, replace = TRUE)
      C_hat_star <- mr_rr(y_j_hat[idx_b, ], x_j_hat[idx_b, ],
                          r = k, W = W, Sigma_X = Sigma_X)$AB
      
      Delta <- C_hat_star - C_hat_full
      Proj  <- t(U2) %*% Delta %*% V2
      M1_star[b] <- sum(Proj^2)
    }
    
    # ---- 3. Critical value ----
    q_crit <- quantile(M1_star, probs = 1 - alpha, type = 8)
    reject <- (M1_obs > q_crit)
    
    out$M1_obs[idx] <- M1_obs
    out$q_crit[idx] <- q_crit
    out$reject[idx] <- reject
    
    if (print) {
      cat(sprintf("r=%d: M1=%.6f, q_crit=%.6f, reject=%s\n",
                  r, M1_obs, q_crit, reject))
    }
  }
  
  # ---- 4. 选最小不被 reject 的 r ----
  idx_keep <- which(out$reject == FALSE)[1]
  if (is.na(idx_keep)) {
    r_sel <- k
    if (print) cat("All rejected; default to full rank:", k, "\n")
  } else {
    r_sel <- out$r[idx_keep]
    if (print) cat("Selected rank =", r_sel, "\n")
  }
  
  return(r_sel)
}


# use this test M1 for best performance
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


# ============================================================ #
# 4. generate low rank true C ####
# ============================================================ #
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


# ============================================================ #
# 5. set universal parameters ####
# ============================================================ #
me_weight_list = c(2.5, 1)
effect_weight_list = c(0.25, 1)
# save_sim_filename = "results/simulate_result_pred_250826.RData"
# regularization_rate_list = c(1.422259e-11, 3.600086e-12, 1.335627e-12, 3.706175e-14)

# X
set.seed(123)
# generate an X for some individual
X = rnorm(px, mean = 0, sd = 1)
Y = C %*% X
round(X,3)

# new: -0.560 -0.230  1.559  0.071  0.129  1.715  0.461 -1.265 -0.687

# # iv strength on x direction
# iv_strength_x = t(X) %*% solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)) %*% X / (t(X) %*% X)
# round(iv_strength_x,2)

# ============================================================ #
# 6. choose regularization rate ####
# ============================================================ #
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

regularization_rate_list = result_matrix_mean


# ============================================================ #
# 7. simulation - main + prediction - other estimators 260416 ####
# ============================================================ #
# shared simulation function for both main + pred results
set.seed(123)
save_sim_pred_filename = "results/simulate_result_pred_260416_1000.RData"

run_simulation_prediction <- function(regularization_rate_list, eloop = 500){
  comb_names <- sprintf("me_%s_effect_%s", rep(me_weight_list, each = 2), effect_weight_list)
  template_list <- setNames(vector("list", length(comb_names)), comb_names)
  
  # bias
  result_AB_list = result_AB_d_list = result_AB_d_r_list = 
    result_C_ivw_list = result_C_adivw_list = result_MrDAG_list =
    # C hat
    AB_list = AB_d_list = AB_d_r_list = 
    C_ivw_list = C_adivw_list = MrDAG_list =
    # C hat %*% X = Y_pred
    Y_pred_AB_list = Y_pred_AB_d_list = Y_pred_AB_d_r_list =
    Y_pred_C_ivw_list = Y_pred_C_adivw_list = Y_pred_MrDAG_list = template_list
  
  iv_strength_list = c()
  parameters_list = c()
  for (me_weight_index in seq_along(me_weight_list)){
    me_weight = me_weight_list[me_weight_index]
    for (effect_weight_index in seq_along(effect_weight_list)){
      effect_weight = effect_weight_list[effect_weight_index]
      
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
      
      iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values) * sqrt(1000)
      iv_strength_list = c(iv_strength_list, iv_strength)
      parameters_list = c(parameters_list, list(c(parameters, iv_strength = iv_strength)))
      
      bias_AB_matrix = bias_AB_d_matrix = bias_AB_d_r_matrix = 
        bias_C_ivw_matrix = bias_C_adivw_matrix = bias_C_MrDAG_matrix =
        AB_matrix = AB_d_matrix = AB_d_r_matrix = 
        C_ivw_matrix = C_adivw_matrix = C_MrDAG_matrix =
        matrix(NA, nrow = py*px, ncol = eloop)
      
      Y_pred_AB_matrix = Y_pred_AB_d_matrix = Y_pred_AB_d_r_matrix =
        Y_pred_C_ivw_matrix = Y_pred_C_adivw_matrix = Y_pred_MrDAG_matrix =
        matrix(NA, nrow = py, ncol = eloop)
      
      B_star = t(B)
      B_d_star = t(B_d)
      for (i in 1:eloop) {
        simulation_result = .simulation(parameters, 
                                        regularization_rate=regularization_rate)
        
        A_hat = simulation_result[[1]]
        B_hat = simulation_result[[2]]
        B_hat_star = t(simulation_result[[2]])
        AB_hat = simulation_result[[3]]
        
        bias_AB_vectorized = as.vector(AB_hat - C)
        bias_AB_matrix[,i] = bias_AB_vectorized
        AB_matrix[,i] = as.vector(AB_hat)
        Y_pred_AB_matrix[,i] = as.vector(AB_hat %*% X)
        
        A_d_hat = simulation_result[[4]]
        B_d_hat = simulation_result[[5]]
        B_d_hat_star = t(simulation_result[[5]])
        AB_d_hat = simulation_result[[6]]
        
        bias_AB_d_vectorized = as.vector(AB_d_hat - C)
        bias_AB_d_matrix[,i] = bias_AB_d_vectorized
        AB_d_matrix[,i] = as.vector(AB_d_hat)
        Y_pred_AB_d_matrix[,i] = as.vector(AB_d_hat %*% X)
        
        A_d_r_hat = simulation_result[[7]]
        B_d_r_hat = simulation_result[[8]]
        B_d_r_hat_star = t(simulation_result[[8]])
        AB_d_r_hat = simulation_result[[9]]
        
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
        
        print(paste("eloop:", i, "me_weight:", me_weight, "effect_weight:", effect_weight))
      }
      
      # store estimator bias
      result_AB_list[[sim_name]] = bias_AB_matrix
      result_AB_d_list[[sim_name]] = bias_AB_d_matrix
      result_AB_d_r_list[[sim_name]] = bias_AB_d_r_matrix
      result_C_ivw_list[[sim_name]] = bias_C_ivw_matrix
      result_C_adivw_list[[sim_name]] = bias_C_adivw_matrix
      result_MrDAG_list[[sim_name]] = bias_C_MrDAG_matrix
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
              AB_list=AB_list,
              AB_d_list=AB_d_list,
              AB_d_r_list=AB_d_r_list,
              C_ivw_list=C_ivw_list,
              C_adivw_list=C_adivw_list,
              MrDAG_list=MrDAG_list,
              Y_pred_AB_list=Y_pred_AB_list,
              Y_pred_AB_d_list=Y_pred_AB_d_list,
              Y_pred_AB_d_r_list=Y_pred_AB_d_r_list,
              Y_pred_C_ivw_list=Y_pred_C_ivw_list,
              Y_pred_C_adivw_list=Y_pred_C_adivw_list,
              Y_pred_MrDAG_list=Y_pred_MrDAG_list
  ))
}

simulate_result_prediction = run_simulation_prediction(regularization_rate_list = regularization_rate_list, 
                                                       eloop = 1000)
save(simulate_result_prediction, file = save_sim_pred_filename)
round(simulate_result_prediction$iv_strength_list, 2)
# 7.05 28.51 19.45 82.25

load(save_sim_pred_filename)


# ============================================================ #
# 8. simulation - for rank test only ####
# ============================================================ #
# 260506 use M0 test
run_rank_selection_sim <- function(parameters, eloop = 300) {
  n  = 1000
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X  = parameters$Sigma_X
  Sigma_Y  = parameters$Sigma_Y
  C  = parameters$C
  r_RR = parameters$r_RR
  W  = parameters$weight.matrix
  
  selected_ranks <- rep(NA, eloop)
  
  for (i in 1:eloop) {
    set.seed(i)
    
    # 生成数据（和 main sim 相同的 DGP）
    gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
    Gamma_j_star = gamma_j_star %*% t(C)
    
    x_j_hat = matrix(0, n, px)
    y_j_hat = matrix(0, n, py)
    for (j in 1:n) {
      x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
      y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
    }
    
    # M0 rank test（chi-square based, no bootstrap）
    selected_ranks[i] <- rank_test_M0(W, y_j_hat, x_j_hat, Sigma_X,
                                      print = FALSE, min_rank = 1)
    
    # 判断正确性
    label <- ifelse(selected_ranks[i] == r_RR, "eq",
                    ifelse(selected_ranks[i] < r_RR, "small", "large"))
    cat(sprintf("Rank sim %d/%d: selected=%d, true=%d (%s)\n",
                i, eloop, selected_ranks[i], r_RR, label))
  }
  
  # 汇总
  freq_table <- table(selected_ranks)
  prop_table <- prop.table(freq_table)
  correct_rate <- mean(selected_ranks == r_RR, na.rm = TRUE)
  
  cat(sprintf("\nCorrect rate: %.1f%%\n", correct_rate * 100))
  print(prop_table)
  
  return(list(
    selected_ranks = selected_ranks,
    freq_table     = freq_table,
    prop_table     = prop_table,
    correct_rate   = correct_rate
  ))
}


# after run main/ load main result data to resuse main sim parameters, then run the rank sim，
set.seed(123)


rank_results_all <- list()
setting_names <- sprintf("me_%s_effect_%s", rep(me_weight_list, each = 2), effect_weight_list)
for (idx in seq_along(setting_names)) {
  sim_name <- setting_names[idx]
  cat("\n====", sim_name, "====\n")
  params <- simulate_result_prediction$parameters_list[[idx]]
  # rank_results_all[[sim_name]] <- run_rank_selection_sim(params, eloop = 1000, bt_loop = 500)
  rank_results_all[[sim_name]] <- run_rank_selection_sim(params, eloop = 1000)
}


save(rank_results_all, file = "results/simulate_result_rank_260506_M0.RData")


# bootstrap rank test - old ver 
# run_rank_selection_sim <- function(parameters, eloop = 300, bt_loop = 500, alpha = 0.05) {
#   n  = 1000
#   py = parameters$py
#   px = parameters$px
#   VX_tilde = parameters$VX_tilde
#   Sigma_X  = parameters$Sigma_X
#   Sigma_Y  = parameters$Sigma_Y
#   C  = parameters$C
#   r_RR = parameters$r_RR
#   W  = parameters$weight.matrix
#   
#   selected_ranks <- rep(NA, eloop)
#   
#   for (i in 1:eloop) {
#     set.seed(i)
#     
#     # 生成数据（和 main sim 相同的 DGP）
#     gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
#     Gamma_j_star = gamma_j_star %*% t(C)
#     
#     x_j_hat = matrix(0, n, px)
#     y_j_hat = matrix(0, n, py)
#     for (j in 1:n) {
#       x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#       y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
#     }
#     
#     # rank test（独立 seed，不影响 main sim）
#     selected_ranks[i] <- rank_test(W, y_j_hat, x_j_hat, Sigma_X,
#                                    print = FALSE, bt_loop = bt_loop,
#                                    alpha = alpha)
#     
#     # 判断正确性
#     label <- ifelse(selected_ranks[i] == r_RR, "eq",
#                     ifelse(selected_ranks[i] < r_RR, "small", "large"))
#     cat(sprintf("Rank sim %d/%d: selected=%d, true=%d (%s)\n",
#                 i, eloop, selected_ranks[i], r_RR, label))
#   }
#   
#   # 汇总
#   freq_table <- table(selected_ranks)
#   prop_table <- prop.table(freq_table)
#   correct_rate <- mean(selected_ranks == r_RR, na.rm = TRUE)
#   
#   cat(sprintf("\nCorrect rate: %.1f%%\n", correct_rate * 100))
#   print(prop_table)
#   
#   return(list(
#     selected_ranks = selected_ranks,
#     freq_table     = freq_table,
#     prop_table     = prop_table,
#     correct_rate   = correct_rate
#   ))
# }


# rank_results_all <- list()
# for (me in me_weight_list) {
#   for (eff in effect_weight_list) {
#     sim_name <- paste0("me_", me, "_effect_", eff)
#     cat("\n====", sim_name, "====\n")
#     params <- .get_parameters(C, me, eff, r_RR = 2)
#     rank_results_all[[sim_name]] <- run_rank_selection_sim(params, eloop = 1000, bt_loop = 500)
#   }
# }
# 
# save(rank_results_all, file = "results/simulate_result_rank_260416.RData")


# ============================================================ #
# 9. simulation - main + prediction - MR-sparse 260409 ####
# ============================================================ #
# run main and this section separately to generate data. Then using the following section to plot the two data together.

set.seed(123)
.simulation_sparse <- function(parameters, sparse_lambda = rep(1e-3, 9)) {
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
      
      # # test
      # me_weight =1
      # effect_weight=1
      # me_weight_index=1
      # effect_weight_index=1
      
      
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
        sparse_lambda = rep(1e-3, px)
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
  sparse_lambda_list = NULL,  
  eloop = 1000
)

save(sparse_results, file = "results/simulate_result_pred_sparse_1000_260712.RData") # sparse 2e-3 --> 1e-3
load("results/simulate_result_pred_sparse_1000_260712.RData")

# ============================================================ #
# 10. simulation - main boxplot (including sparse) ####
# ============================================================ #
# plot main boxplot
plot_simulation_result <- function(sim_result = simulate_result_prediction,
                                   sparse_result = NULL,
                                   me_weight, effect_weight,
                                   true_C,
                                   width = 16, height = 8,
                                   winsor_prob = c(0.99, 0.99, 0.99)) {
  
  sim_name  <- paste0("me_", me_weight, "_effect_", effect_weight)
  today_str <- format(Sys.Date(), "%y%m%d")
  
  py <- nrow(true_C)
  px <- ncol(true_C)
  
  # Generic labels for simulation
  row_labels <- setNames(paste0("Outcome ", 1:py), as.character(0:(py - 1)))
  col_labels <- setNames(paste0("Exposure ", 1:px), as.character(0:(px - 1)))
  
  estimator_names  <- c("C_ivw_list", "C_adivw_list", "AB_list", "AB_d_list",
                        "AB_d_r_list", "C_sparse_list", "MrDAG_list")
  estimator_labels <- c("IVW", "SRIVW", "Naive MR-rr", "MR-rr",
                        "Reg. MR-rr", "Sparse MR-rr", "MrDAG")
  
  pos_ids <- expand.grid(Row = 0:(py - 1), Col = 0:(px - 1))
  
  # ----- Organize data -----
  base_names <- setdiff(estimator_names, "C_sparse_list")
  df_long <- purrr::map_dfr(base_names, function(est) {
    mat <- sim_result[[est]][[sim_name]]
    purrr::map_dfr(1:(px * py), function(i) {
      tibble::tibble(
        Estimator = est,
        Row = pos_ids$Row[i],
        Col = pos_ids$Col[i],
        Value = mat[i, ]
      )
    })
  })
  
  if (!is.null(sparse_result)) {
    mat <- sparse_result$C_sparse_list[[sim_name]]
    df_sparse <- purrr::map_dfr(1:(px * py), function(i) {
      tibble::tibble(
        Estimator = "C_sparse_list",
        Row = pos_ids$Row[i],
        Col = pos_ids$Col[i],
        Value = mat[i, ]
      )
    })
    df_long <- dplyr::bind_rows(df_long, df_sparse)
  } else {
    estimator_names  <- setdiff(estimator_names, "C_sparse_list")
    estimator_labels <- estimator_labels[-which(estimator_labels == "Sparse MR-rr")]
  }
  
  df_long <- df_long %>%
    dplyr::mutate(
      Estimator = factor(
        Estimator,
        levels = estimator_names,
        labels = estimator_labels
      )
    )
  
  df_true <- expand.grid(Row = 0:(py - 1), Col = 0:(px - 1)) %>%
    dplyr::mutate(TrueValue = as.vector(true_C))
  
  # ----- Row-specific y-limits for display only -----
  ylim_by_row <- purrr::map_dfr(0:(py - 1), function(r) {
    df_long %>%
      dplyr::filter(Row == r) %>%
      dplyr::summarise(
        Row = r,
        ylo = quantile(Value, 1 - winsor_prob[r + 1], na.rm = TRUE),
        yhi = quantile(Value,     winsor_prob[r + 1], na.rm = TRUE)
      )
  })
  
  # ----- One plot per outcome row -----
  make_row_plot <- function(row_id) {
    ylo <- ylim_by_row$ylo[ylim_by_row$Row == row_id]
    yhi <- ylim_by_row$yhi[ylim_by_row$Row == row_id]
    
    show_xtext <- row_id == (py - 1)
    show_top   <- row_id == 0
    
    d <- df_long %>%
      dplyr::filter(Row == row_id) %>%
      dplyr::mutate(RowLab = row_labels[as.character(row_id)])
    
    dt <- df_true %>%
      dplyr::filter(Row == row_id) %>%
      dplyr::mutate(RowLab = row_labels[as.character(row_id)])
    
    ggplot(d, aes(x = Estimator, y = Value, fill = Estimator)) +
      geom_boxplot(outlier.size = 0.3) +
      geom_hline(
        data = dt,
        aes(yintercept = TrueValue),
        color = "red",
        linetype = "dashed",
        linewidth = 0.5,
        inherit.aes = FALSE
      ) +
      coord_cartesian(ylim = c(ylo, yhi)) +
      facet_grid(
        RowLab ~ Col,
        labeller = labeller(Col = col_labels)
      ) +
      theme_bw(base_size = 10) +
      theme(
        axis.text.x = if (show_xtext) {
          element_text(angle = 90, vjust = 0.5, hjust = 1)
        } else {
          element_blank()
        },
        axis.ticks.x = if (show_xtext) element_line() else element_blank(),
        strip.text.x = if (show_top) element_text(size = 9) else element_blank(),
        strip.background.x = if (show_top) element_rect() else element_blank(),
        strip.text.y = element_text(size = 9)
      ) +
      xlab("") +
      ylab("")
  }
  
  p_combined <- (make_row_plot(0) / make_row_plot(1) / make_row_plot(2)) +
    patchwork::plot_layout(guides = "collect") &
    theme(legend.position = "right")
  
  # ----- Add common y-axis label -----
  p_final <- cowplot::ggdraw() +
    cowplot::draw_plot(
      p_combined,
      x = 0.025,
      y = 0,
      width = 0.975,
      height = 1
    ) +
    cowplot::draw_label(
      "Estimated C Value",
      x = 0.012,
      y = 0.5,
      angle = 90,
      size = 12
    )
  
  full_path <- file.path(
    getwd(),
    "results",
    sprintf("lessvar_plot_%s_%s_generic_labels.png", sim_name, today_str)
  )
  
  ggsave(full_path, plot = p_final, width = width, height = height, dpi = 300)
  message("Plot saved to ", full_path)
  print(p_final)
}


plot_simulation_result(sim_result = simulate_result_prediction,
                       sparse_result = sparse_results,
                       me_weight = "2.5",
                       effect_weight = "0.25",
                       true_C = C,
                       width = 10, height = 5,
                       winsor_prob = c(0.98, 0.9925, 0.985))


plot_simulation_result(sim_result = simulate_result_prediction,
                       sparse_result = sparse_results,
                       me_weight = "1",
                       effect_weight = "1",
                       true_C = C,
                       width = 10, height = 5)


# ============================================================ #
# 11. simulation - pred boxplot (including sparse)  ####
# ============================================================ #
# TODO: update for new sparse results
# load("results/simulate_result_pred_sparse_1000_260409.RData")

load("results/simulate_result_pred_sparse_1000_260712.RData")
# new
plot_simulation_result_pred <- function(sim_result, 
                                        me_weight, effect_weight,
                                        Y,  # true Y = C %*% X
                                        sparse_result = NULL,
                                        type = c("box", "ci"),
                                        width = 12, height = 6) {
  type <- match.arg(type)
  
  sim_name  <- paste0("me_", me_weight, "_effect_", effect_weight)
  today_str <- format(Sys.Date(), "%y%m%d")
  filename  <- sprintf("results/pred_plot_with_sparse_%s_type_%s_%s_generic_labels.png",
                       sim_name, type, today_str)
  
  estimator_names <- c(
    "Y_pred_C_ivw_list", 
    "Y_pred_C_adivw_list", 
    "Y_pred_AB_list", 
    "Y_pred_AB_d_list", 
    "Y_pred_AB_d_r_list", 
    "Y_pred_MrDAG_list"
  )
  
  estimator_labels <- c(
    "IVW", 
    "SRIVW", 
    "Naive MR-rr", 
    "MR-rr", 
    "Reg. MR-rr", 
    "MrDAG"
  )
  
  # Final display order: put Sparse before MrDAG
  final_levels <- c(
    "IVW", 
    "SRIVW", 
    "Naive MR-rr", 
    "MR-rr",
    "Reg. MR-rr", 
    "Sparse MR-rr", 
    "MrDAG"
  )
  
  Y_vec <- as.vector(Y)
  py <- length(Y_vec)
  pos_ids <- data.frame(Row = 0:(py - 1))
  
  # Generic outcome labels for simulation
  row_labels <- setNames(
    paste0("Outcome ", 1:py),
    as.character(0:(py - 1))
  )
  
  # ===== Organize data =====
  if (type == "box") {
    df_long <- purrr::map2_dfr(estimator_names, estimator_labels, function(est, label) {
      matrix_data <- sim_result[[est]][[sim_name]]
      purrr::map_dfr(1:py, function(i) {
        tibble::tibble(
          Estimator = label,
          Row = pos_ids$Row[i],
          Value = matrix_data[i, ]
        )
      })
    })
    
    if (!is.null(sparse_result)) {
      matrix_data <- sparse_result$Y_pred_C_sparse_list[[sim_name]]
      df_sparse <- purrr::map_dfr(1:py, function(i) {
        tibble::tibble(
          Estimator = "Sparse MR-rr",
          Row = pos_ids$Row[i],
          Value = matrix_data[i, ]
        )
      })
      df_long <- dplyr::bind_rows(df_long, df_sparse)
    }
    
    df_long$Estimator <- factor(df_long$Estimator, levels = final_levels)
    
  } else {
    df_long <- purrr::map2_dfr(estimator_names, estimator_labels, function(est, label) {
      matrix_data <- sim_result[[est]][[sim_name]]
      purrr::map_dfr(1:py, function(i) {
        vals <- matrix_data[i, ]
        tibble::tibble(
          Estimator = label,
          Row = pos_ids$Row[i],
          Median = median(vals, na.rm = TRUE),
          Lower  = quantile(vals, 0.025, na.rm = TRUE),
          Upper  = quantile(vals, 0.975, na.rm = TRUE)
        )
      })
    })
    
    if (!is.null(sparse_result)) {
      matrix_data <- sparse_result$Y_pred_C_sparse_list[[sim_name]]
      df_sparse <- purrr::map_dfr(1:py, function(i) {
        vals <- matrix_data[i, ]
        tibble::tibble(
          Estimator = "Sparse MR-rr",
          Row = pos_ids$Row[i],
          Median = median(vals, na.rm = TRUE),
          Lower  = quantile(vals, 0.025, na.rm = TRUE),
          Upper  = quantile(vals, 0.975, na.rm = TRUE)
        )
      })
      df_long <- dplyr::bind_rows(df_long, df_sparse)
    }
    
    df_long$Estimator <- factor(df_long$Estimator, levels = final_levels)
  }
  
  # ===== True Y horizontal line =====
  df_true <- tibble::tibble(
    Row = 0:(py - 1),
    TrueValue = Y_vec
  )
  
  # ===== Plot layers =====
  plot_layers <- if (type == "box") {
    list(
      geom_boxplot(aes(y = Value, fill = Estimator), outlier.size = 0.3),
      geom_hline(
        data = df_true,
        aes(yintercept = TrueValue),
        color = "red",
        linetype = "dashed",
        linewidth = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    )
  } else {
    list(
      geom_point(
        aes(y = Median, color = Estimator),
        position = position_dodge(width = 0.5)
      ),
      geom_errorbar(
        aes(ymin = Lower, ymax = Upper, color = Estimator),
        width = 0.2,
        position = position_dodge(width = 0.5)
      ),
      geom_hline(
        data = df_true,
        aes(yintercept = TrueValue),
        color = "red",
        linetype = "dashed",
        linewidth = 0.8,
        inherit.aes = FALSE,
        show.legend = FALSE
      )
    )
  }
  
  # ===== Plot =====
  p <- ggplot(df_long, aes(x = Estimator)) +
    plot_layers +
    facet_wrap(
      ~ Row,
      ncol = py,
      labeller = labeller(Row = row_labels)
    ) +
    coord_cartesian(ylim = c(-2, 2)) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = element_text(size = 10),
      legend.position = "right"
    ) +
    ylab("Estimated Risk Score") +
    xlab("")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}


plot_simulation_result_pred(
  sim_result = simulate_result_prediction,
  sparse_result = sparse_results,  # 加入 sparse 结果
  me_weight = 1,
  effect_weight = 1,
  Y = Y,
  type = "box"
)


# ============================================================ #
# 12. simulation - generate main table ####
# ============================================================ #
gen_table <- function(simulate_result_prediction,
                      me_weight_to_plot,
                      effect_weight_to_plot,
                      sparse_result = NULL)
{
  parameters_list        <- simulate_result_prediction$parameters_list
  result_naive_MRrr_list <- simulate_result_prediction$result_AB_list
  result_MRrr_list       <- simulate_result_prediction$result_AB_d_list
  result_r_MRrr_list     <- simulate_result_prediction$result_AB_d_r_list
  result_ivw_list        <- simulate_result_prediction$result_C_ivw_list
  result_adivw_list      <- simulate_result_prediction$result_C_adivw_list
  result_MrDAG_list      <- simulate_result_prediction$result_MrDAG_list
  result_sparse_list     <- if (!is.null(sparse_result)) sparse_result$result_C_sparse_list else NULL
  
  px <- 9
  py <- 3
  sim_name <- paste0("me_", me_weight_to_plot, "_effect_", effect_weight_to_plot)
  
  # 小工具：给一个 result_list，返回 (abs_mean, mean, sd) 三个长度 px*py 的向量
  compute_entry_stats <- function(res_list) {
    mat <- res_list[[sim_name]]
    abs_mean <- mean_b <- sd_b <- rep(NA, px * py)
    for (i in 1:(px * py)) {
      abs_mean[i] <- abs(mean(mat[i, ], na.rm = TRUE))
      mean_b[i]   <- mean(mat[i, ], na.rm = TRUE)
      sd_b[i]     <- sd(mat[i, ], na.rm = TRUE)
    }
    list(abs_mean = abs_mean, mean = mean_b, sd = sd_b)
  }
  
  s_naive  <- compute_entry_stats(result_naive_MRrr_list)
  s_d      <- compute_entry_stats(result_MRrr_list)
  s_d_r    <- compute_entry_stats(result_r_MRrr_list)
  s_ivw    <- compute_entry_stats(result_ivw_list)
  s_adivw  <- compute_entry_stats(result_adivw_list)
  s_MrDAG  <- compute_entry_stats(result_MrDAG_list)
  s_sparse <- if (!is.null(result_sparse_list)) compute_entry_stats(result_sparse_list) else NULL
  
  # 组织成矩阵（py x px）
  to_mat <- function(v) t(matrix(v, nrow = px, ncol = py))
  
  avg_bias        <- to_mat(s_naive$abs_mean)
  avg_bias_d      <- to_mat(s_d$abs_mean)
  avg_bias_d_r    <- to_mat(s_d_r$abs_mean)
  avg_bias_ivw    <- to_mat(s_ivw$abs_mean)
  avg_bias_adivw  <- to_mat(s_adivw$abs_mean)
  avg_bias_MrDAG  <- to_mat(s_MrDAG$abs_mean)
  avg_bias_sparse <- if (!is.null(s_sparse)) to_mat(s_sparse$abs_mean) else NULL
  
  avg_sd        <- to_mat(s_naive$sd)
  avg_sd_d      <- to_mat(s_d$sd)
  avg_sd_d_r    <- to_mat(s_d_r$sd)
  avg_sd_ivw    <- to_mat(s_ivw$sd)
  avg_sd_adivw  <- to_mat(s_adivw$sd)
  avg_sd_MrDAG  <- to_mat(s_MrDAG$sd)
  avg_sd_sparse <- if (!is.null(s_sparse)) to_mat(s_sparse$sd) else NULL
  
  # 真 C_r
  me_weight_index     <- match(me_weight_to_plot, me_weight_list)
  effect_weight_index <- match(effect_weight_to_plot, effect_weight_list)
  parameters <- parameters_list[[.get_sim_index(me_weight_index, effect_weight_index)]]
  C_r <- parameters$C_r
  mean_abs_C_entry <- round(mean(abs(C_r)), 3)
  
  # summary 小工具
  med       <- function(x) round(median(abs(x), na.rm = TRUE), 3)
  med_sd    <- function(x) round(median(x, na.rm = TRUE), 3)
  qiqr_bias <- function(x) round(quantile(abs(x), c(0.25, 0.75), na.rm = TRUE), 3)
  qiqr_sd   <- function(x) round(quantile(x,      c(0.25, 0.75), na.rm = TRUE), 3)
  
  has_sparse <- !is.null(s_sparse)
  
  # ---- 组装输出：顺序为 IVW → adIVW → Naive → MR-rr → Reg. → Sparse → MrDAG ----
  out <- list(
    # median
    median_abs_bias_ivw    = med(avg_bias_ivw),   median_sd_ivw    = med_sd(avg_sd_ivw),
    median_abs_bias_adivw  = med(avg_bias_adivw), median_sd_adivw  = med_sd(avg_sd_adivw),
    median_abs_bias_naive  = med(avg_bias),       median_sd_naive  = med_sd(avg_sd),
    median_abs_bias        = med(avg_bias_d),     median_sd        = med_sd(avg_sd_d),
    median_abs_bias_r      = med(avg_bias_d_r),   median_sd_r      = med_sd(avg_sd_d_r)
  )
  if (has_sparse) {
    out$median_abs_bias_sparse <- med(avg_bias_sparse)
    out$median_sd_sparse       <- med_sd(avg_sd_sparse)
  }
  out$median_abs_bias_MrDAG <- med(avg_bias_MrDAG)
  out$median_sd_MrDAG       <- med_sd(avg_sd_MrDAG)
  
  # IQR
  out$quantile_abs_bias_ivw   <- qiqr_bias(avg_bias_ivw);   out$quantile_sd_ivw   <- qiqr_sd(avg_sd_ivw)
  out$quantile_abs_bias_adivw <- qiqr_bias(avg_bias_adivw); out$quantile_sd_adivw <- qiqr_sd(avg_sd_adivw)
  out$quantile_abs_bias_naive <- qiqr_bias(avg_bias);       out$quantile_sd_naive <- qiqr_sd(avg_sd)
  out$quantile_abs_bias       <- qiqr_bias(avg_bias_d);     out$quantile_sd       <- qiqr_sd(avg_sd_d)
  out$quantile_abs_bias_r     <- qiqr_bias(avg_bias_d_r);   out$quantile_sd_r     <- qiqr_sd(avg_sd_d_r)
  if (has_sparse) {
    out$quantile_abs_bias_sparse <- qiqr_bias(avg_bias_sparse)
    out$quantile_sd_sparse       <- qiqr_sd(avg_sd_sparse)
  }
  out$quantile_abs_bias_MrDAG <- qiqr_bias(avg_bias_MrDAG)
  out$quantile_sd_MrDAG       <- qiqr_sd(avg_sd_MrDAG)
  
  return(out)
}


for (re in effect_weight_list) {
  for (me in me_weight_list) {
    print(me); print(re)
    print(gen_table(simulate_result_prediction,
                    me_weight_to_plot = me,
                    effect_weight_to_plot = re,
                    sparse_result = sparse_results))
  }
}
