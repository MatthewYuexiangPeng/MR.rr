# ============================================================ #
# 1. loading environment and data ####
# ============================================================ #
rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
# library(GRAPPLE)
# library(matrixStats)
# library(MASS)
library(ggplot2)
# library(tidyverse)
library(patchwork)
# library(pheatmap::pheatmap)
# library(readxl)
# library(glmnet)
library(foreach)
library(doParallel)
library(reshape2)
library(dplyr)
library(ggdist)
library(pheatmap)


set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

source("scripts/MR_rr_estimators.R")

lip_data_real = read.csv('data/dat_1e-4.csv')
lip_corr_real = read.csv('data/rho_mat_1e-4.csv')
n_Y = c(1296908,1241207,1245612,1241619)

var_Z = 2 * lip_data_real$ImpMAF * (1 - lip_data_real$ImpMAF)
gamma_exp_j = as.matrix(lip_data_real[,paste0('gamma_exp',1:9)])
se_exp_j = as.matrix(lip_data_real[,paste0('se_exp',1:9)])
gamma_out_j = as.matrix(lip_data_real[,paste0('gamma_out',1:4)])
se_out_j = as.matrix(lip_data_real[,paste0('se_out',1:4)])
gamma_exp_j = gamma_exp_j * sqrt(var_Z)
se_exp_j = se_exp_j * sqrt(var_Z)
gamma_out_j = gamma_out_j * sqrt(var_Z)
se_out_j = se_out_j * sqrt(var_Z)

cor_exp = as.matrix(lip_corr_real[1:9,1:9])
cor_out = as.matrix(lip_corr_real[10:13,10:13])

exp_name = read.csv('data/traits_1e-4.csv')$x
# out_name = c("IS", "LAS", "CES", "SVS")

# 26/2/8 delete first outcome
gamma_out_j = gamma_out_j[,2:4]
se_out_j = se_out_j[,2:4]
cor_out = as.matrix(lip_corr_real[11:13,11:13])

out_name = c("LAS", "CES", "SVS")


# ============================================================ #
# 2. functions ####
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


# Construct the surrogate exposure matrix using exactly the same steps as
# the current/old mr_rr_sparse() implementation in MR_rr_estimators.R.
.construct_gamma_tilde <- function(GAMMA_hat, gamma_hat, Sigma_X) {
  n <- nrow(gamma_hat)
  P_GAMMA <- GAMMA_hat %*% solve(t(GAMMA_hat) %*% GAMMA_hat) %*% t(GAMMA_hat)
  P_GAMMA_prep <- diag(1, n) - P_GAMMA
  Sigma_gammahat <- t(gamma_hat) %*% gamma_hat / n
  matrix_part1 <- Sigma_gammahat - t(P_GAMMA %*% gamma_hat) %*% (P_GAMMA %*% gamma_hat) / n
  matrix_part2 <- matrix_part1 - Sigma_X

  if (!is_psd(matrix_part2)) {
    matrix_part2 <- .nearest_psd(matrix_part2, epsilon = 1e-6)
  }

  R <- chol(matrix_part1)
  Q <- chol(matrix_part2)
  L <- solve(R) %*% Q
  gamma_tilde <- P_GAMMA %*% gamma_hat + (P_GAMMA_prep %*% gamma_hat) %*% L

  return(gamma_tilde)
}


# Post-selection sparse MR-rr refit.
# The support is selected once by the old penalized mr_rr_sparse() estimator.
# Conditional on that support, this function drops the l1 penalty and alternates
# the same A update as the old sparse algorithm with an unpenalized B update.
mr_rr_sparse_refit <- function(GAMMA_hat, gamma_hat, W, Sigma_X, support,
                               A_init = NULL, B_init = NULL,
                               max_iter = 100, tol = 1e-2) {
  px <- ncol(gamma_hat)
  py <- ncol(GAMMA_hat)
  r <- nrow(support)

  if (!all(dim(support) == c(r, px))) {
    stop("support has an incorrect dimension")
  }

  gamma_tilde <- .construct_gamma_tilde(GAMMA_hat, gamma_hat, Sigma_X)
  W_sqrt <- .sqrt_matrix(W)

  if (is.null(A_init) || is.null(B_init)) {
    init_result <- mr_rr(GAMMA_hat, gamma_hat, r = r, W = W, Sigma_X = Sigma_X)
    A_hat <- init_result$A
    B_hat <- init_result$B
  } else {
    A_hat <- A_init
    B_hat <- B_init
  }
  B_hat[!support] <- 0

  for (iter in 1:max_iter) {
    # Fix A, update only the selected entries of B without an l1 penalty.
    B_hat <- matrix(0, nrow = r, ncol = px)
    for (k in 1:r) {
      selected_idx <- which(support[k, ])
      if (length(selected_idx) > 0) {
        latent_response <- as.vector(GAMMA_hat %*% W %*% A_hat[, k])
        B_hat[k, selected_idx] <- qr.solve(
          gamma_tilde[, selected_idx, drop = FALSE],
          latent_response
        )
      }
    }

    # Fix B, optimize A using the same weighted Procrustes update as mr_rr_sparse().
    svd_result <- svd(B_hat %*% t(gamma_tilde) %*% GAMMA_hat %*% W_sqrt)
    A_hat_new <- solve(W_sqrt) %*% svd_result$v %*% t(svd_result$u)

    # Keep the convergence criterion consistent with the old sparse implementation.
    dist <- norm(A_hat_new %*% B_hat - A_hat %*% B_hat, "F") /
      max(1e-8, norm(A_hat %*% B_hat, "F"))
    A_hat <- A_hat_new
    if (dist < tol) break
  }

  C_hat <- A_hat %*% B_hat
  return(list(A = A_hat, B = B_hat, AB = C_hat, support = support,
              iter = iter, dist = dist))
}


# Canonicalize signs only for displaying A and B. This does not alter AB.
.canonicalize_AB <- function(A, B) {
  for (k in 1:nrow(B)) {
    j <- which.max(abs(B[k, ]))
    if (B[k, j] < 0) {
      B[k, ] <- -B[k, ]
      A[, k] <- -A[, k]
    }
  }
  return(list(A = A, B = B, AB = A %*% B))
}


estimate_C <- function(gamma_exp_j,se_exp_j,gamma_out_j,se_out_j,cor_exp,cor_out,
                       r_RR = NULL, sparse_support = NULL,
                       sparse_A_init = NULL, sparse_B_init = NULL,
                       print_rank = TRUE){
  px = 9
  py = 3
  pz = n = dim(lip_data_real)[1]
  #### 1. ivw #### 
  ivw_list <- vector("list", py)
  
  for (i in 1:py) {
    res <- mvmr.ivw(
      beta.exposure = gamma_exp_j,
      se.exposure = se_exp_j,
      beta.outcome = as.vector(gamma_out_j[, i]),
      se.outcome = se_out_j[, i],
      gen_cor = cor_exp
    )$beta.hat
    ivw_list[[i]] <- t(res)  # 1 x px row
  }
  
  # Combine into matrix
  C_ivw <- do.call(rbind, ivw_list)
  
  
  #### 2. adivw ####
  adivw_list <- vector("list", py)
  for (i in 1:py) {
    res <- mvmr.divw(
      beta.exposure = gamma_exp_j,
      se.exposure = se_exp_j,
      beta.outcome = as.vector(gamma_out_j[, i]),
      se.outcome = se_out_j[, i],
      gen_cor = cor_exp
    )$beta.hat
    adivw_list[[i]] <- t(res)  # 1 x px row
  }
  
  # Combine into matrix
  C_adivw <- do.call(rbind, adivw_list)
  
  
  #### 3.4. naive MRrr and MRrr ####
  Sigma_X_sum <- matrix(0, nrow = px, ncol = px)
  Sigma_Y_sum <- matrix(0, nrow = py, ncol = py)
  
  for (j in 1:pz) {
    D_exp_j <- diag(se_exp_j[j, ])
    D_out_j <- diag(se_out_j[j, ])
    
    Sigma_X_sum <- Sigma_X_sum + D_exp_j %*% cor_exp %*% D_exp_j
    Sigma_Y_sum <- Sigma_Y_sum + D_out_j %*% cor_out %*% D_out_j
  }
  
  Sigma_X <- Sigma_X_sum / pz
  Sigma_Y <- Sigma_Y_sum / pz
  
  W = solve(Sigma_Y)
  
  # Select rank only when r_RR is not supplied. For bootstrap, the rank selected
  # from the original data is passed into this function and remains fixed.
  if (is.null(r_RR)) {
    r_RR = rank_test_M0(W, gamma_out_j, gamma_exp_j, Sigma_X,
                        print = print_rank)
  }
  
  MRrr_naive = mr_rr_naive(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W)
  C_naive_MRrr = MRrr_naive$AB
  
  MRrr = mr_rr(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X)
  C_MRrr = MRrr$AB
  A_MRrr = MRrr$A
  B_MRrr = MRrr$B
  
  #### 5. MRrr with regularization ####
  ## choose rate
  SigmaXX = cov(gamma_exp_j)
  VX_tilde = SigmaXX - Sigma_X
  
  # check if PSD -- pass
  # if (any(eigen(VX_tilde)$values < 0)) {
  #   stop("VX_tilde is not positive semi-definite.")
  # }
  
  W_sqrt = .sqrt_matrix(W)
  
  mu_min = min(eigen(solve(.sqrt_matrix(Sigma_X)) %*% VX_tilde %*% solve(.sqrt_matrix(Sigma_X)))$values)
  # iv_strength = mu_min * sqrt(1000)
  
  # 25.10.2
  D_list = seq(0, 15, 1)
  regu_rate_list = c()
  obj_value = c()
  for (D in D_list) {
    sigma_y2 = mean(eigen(Sigma_X)$values)
    
    regu_rate = sigma_y2^2 * exp(0.3*(D - sqrt(n) * mu_min))/n
    
    regu_rate_list = c(regu_rate_list, regu_rate)
    result_d <- mr_rr_regularized(gamma_out_j, gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate)
    C_hat = result_d$AB
    # update: 260716: need to devided by pz. objective function value
    residual_term <- norm(
      (y_j_hat - x_j_hat %*% t(C_hat)) %*% W_sqrt,
      type = "F"
    )^2 / n   # n = pZ
    
    debias_term <- sum(diag(
      W_sqrt %*%
        C_hat %*%
        Sigma_X %*%
        t(C_hat) %*%
        W_sqrt
    ))
    
    obj <- residual_term - debias_term
    obj_value = c(obj_value, obj)
  }
  
  opt_rate = regu_rate_list[order(obj_value)][1]
  
  
  MRrr_regularized = mr_rr_regularized(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=opt_rate)
  C_MRrr_regularized = MRrr_regularized$AB
  A_MRrr_regularized = MRrr_regularized$A
  B_MRrr_regularized = MRrr_regularized$B
  
  
  #### 6. MR-rr sparse: selection followed by fixed-support refit ####
  MRrr_sparse_selection = NULL
  if (is.null(sparse_support)) {
    # Temporary use of the current/old sparse algorithm with lambda = 1e-3.
    MRrr_sparse_selection = mr_rr_sparse(
      GAMMA_hat = gamma_out_j,
      gamma_hat = gamma_exp_j,
      W = W,
      Sigma_X = Sigma_X,
      lambda = rep(1e-3, px),
      r = r_RR,
      max_iter = 100,
      tol = 1e-2
    )
    # mr_rr_sparse() already thresholds abs(B) < 1e-2 to exactly zero.
    sparse_support = MRrr_sparse_selection$B != 0
    sparse_A_init = MRrr_sparse_selection$A
    sparse_B_init = MRrr_sparse_selection$B
  }
  
  MRrr_sparse = mr_rr_sparse_refit(
    GAMMA_hat = gamma_out_j,
    gamma_hat = gamma_exp_j,
    W = W,
    Sigma_X = Sigma_X,
    support = sparse_support,
    A_init = sparse_A_init,
    B_init = sparse_B_init,
    max_iter = 100,
    tol = 1e-2
  )
  C_MRrr_sparse = MRrr_sparse$AB
  
  #### 7. Mr_DAG ####
  C_Mr_DAG = Mr_DAG(Y=gamma_out_j, X=gamma_exp_j, niter = 10000, burnin = 2000)
  
  output = list(
    C_ivw                    = C_ivw,
    C_adivw                  = C_adivw,
    C_naive_MRrr             = C_naive_MRrr,
    C_MRrr                   = C_MRrr,
    C_MRrr_regularized       = C_MRrr_regularized,
    C_MRrr_sparse            = C_MRrr_sparse,
    C_Mr_DAG                 = C_Mr_DAG,
    A_MRrr                   = A_MRrr,
    B_MRrr                   = B_MRrr,
    A_MRrr_regularized       = A_MRrr_regularized,
    B_MRrr_regularized       = B_MRrr_regularized,
    A_MRrr_sparse            = MRrr_sparse$A,
    B_MRrr_sparse            = MRrr_sparse$B,
    A_MRrr_sparse_selection  = if (is.null(MRrr_sparse_selection)) NULL else MRrr_sparse_selection$A,
    B_MRrr_sparse_selection  = if (is.null(MRrr_sparse_selection)) NULL else MRrr_sparse_selection$B,
    sparse_support           = sparse_support,
    sparse_refit_iter        = MRrr_sparse$iter,
    sparse_refit_dist        = MRrr_sparse$dist,
    r_RR                     = r_RR,
    Sigma_X                  = Sigma_X,
    Sigma_Y                  = Sigma_Y,
    W                        = W,
    opt_rate                 = opt_rate
  )
}


.compute_iv_strength_all_directions <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = NULL) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {
    stop("either beta.exposure or se.exposure only has one column; if univariable MR is performed, please use univariable MR methods.")
  }
  
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {
    P <- diag(K)
  } else {
    P <- as.matrix(gen_cor)
  }
  
  if (ncol(P) != ncol(beta.exposure)) {
    stop("The shared correlation matrix has a different number of columns than the input beta.exposure")
  }
  
  if (nrow(beta.exposure) != length(beta.outcome)) {
    stop("The number of SNPs in beta.exposure and beta.outcome is different")
  }
  
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  p <- nrow(beta.exposure)
  
  # 生成标准化所需矩阵
  P_eigen <- eigen(P)
  P_root_inv <- P_eigen$vectors %*% diag(1 / sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
  
  Vj_root_inv <- lapply(1:p, function(j) P_root_inv %*% diag(1 / se.exposure[j, ]))
  
  # 构建 IV strength matrix
  IV_strength_matrix <- Reduce("+", lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j, ]
    beta.exposure.V %*% t(beta.exposure.V)
  })) - p * diag(K)
  
  # 计算所有方向的 IV strength（特征值）
  eig_vals <- eigen(IV_strength_matrix / sqrt(p), symmetric = TRUE)$values
  
  # 两位小数
  eig_vals <- round(eig_vals, 2)
  
  return(eig_vals)  # 输出所有方向的IV strength
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


rank_test_old = function(W, y_j_hat, x_j_hat, Sigma_X, print = TRUE, min_rank = 1, bt_loop = 100){
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
    C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% C_hat_MRrr_vec
    
    # # consider centering c hat
    # mu_hat = colMeans(C_hat_matrix)
    # C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% (C_hat_MRrr_vec - mu_hat)
    
    C_tilde = matrix(C_tilde_vec, nrow = py, ncol = px)
    
    # singular value of C_tilde
    svd_C_tilde = svd(C_tilde)
    lambda_tilde = svd_C_tilde$d
    
    # M1
    M1 = sum(lambda_tilde[(r + 1):py]**2)
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


# ============================================================ #
# 3. Iv strength ####
# ============================================================ #
# IV_strength_real = c()
# for (i in 1:4){
#   IV_strength_real = c(IV_strength_real, mvmr.ivw(gamma_exp_j, se_exp_j, gamma_out_j[,i], se_out_j[,i], gen_cor = cor_exp)$iv_strength_parameter)
# }

IV_strength_real = .compute_iv_strength_all_directions(gamma_exp_j, se_exp_j, gamma_out_j[,1], se_out_j[,1], gen_cor = cor_exp)
min(IV_strength_real)


# ============================================================ #
# 4. Estimate C using different methods ####
# ============================================================ #
real_result <- estimate_C(
  gamma_out_j = gamma_out_j,
  gamma_exp_j = gamma_exp_j,
  se_out_j = se_out_j,
  se_exp_j = se_exp_j,
  cor_out = cor_out,
  cor_exp = cor_exp
)

C_ivw <- real_result$C_ivw
C_adivw <- real_result$C_adivw
C_naive_MRrr <- real_result$C_naive_MRrr
C_MRrr <- real_result$C_MRrr
C_MRrr_regularized <- real_result$C_MRrr_regularized
C_MRrr_sparse <- real_result$C_MRrr_sparse
C_Mr_DAG <- real_result$C_Mr_DAG

# Rank and sparse support selected once on the original data.
r_RR <- real_result$r_RR
sparse_support <- real_result$sparse_support
cat("Selected rank used by all three MR-rr estimators:", r_RR, "\n")
writeLines(
  paste0("Selected rank used by all three MR-rr estimators: ", r_RR),
  "results/selected_rank_260713.txt"
)

# A and B used for interpretation. Sign canonicalization is for display only.
MRrr_AB_display <- .canonicalize_AB(real_result$A_MRrr, real_result$B_MRrr)
Reg_MRrr_AB_display <- .canonicalize_AB(real_result$A_MRrr_regularized,
                                         real_result$B_MRrr_regularized)
Sparse_MRrr_AB_display <- .canonicalize_AB(real_result$A_MRrr_sparse,
                                            real_result$B_MRrr_sparse)

write.csv(MRrr_AB_display$A, "results/A_MRrr_260713.csv", row.names = FALSE)
write.csv(MRrr_AB_display$B, "results/B_MRrr_260713.csv", row.names = FALSE)
write.csv(Reg_MRrr_AB_display$A, "results/A_MRrr_regularized_260713.csv", row.names = FALSE)
write.csv(Reg_MRrr_AB_display$B, "results/B_MRrr_regularized_260713.csv", row.names = FALSE)
write.csv(Sparse_MRrr_AB_display$A, "results/A_MRrr_sparse_refit_260713.csv", row.names = FALSE)
write.csv(Sparse_MRrr_AB_display$B, "results/B_MRrr_sparse_refit_260713.csv", row.names = FALSE)
write.csv(real_result$B_MRrr_sparse_selection,
          "results/B_MRrr_sparse_selection_260713.csv", row.names = FALSE)
write.csv(1 * sparse_support,
          "results/B_MRrr_sparse_support_260713.csv", row.names = FALSE)

norm(C_adivw, "F")
norm(C_MRrr_regularized, "F")


# ============================================================ #
# 5. visualize C hat with heatmap ####
# ============================================================ #
draw_multiple_heatmaps <- function(
    matrix_list,
    titles = NULL,
    row_labels = out_name,
    col_labels = exp_name,
    filename = NULL,
    ncol = NULL,
    global_title = NULL,
    cellwidth = 30,
    cellheight = 20,
    fontsize_row = 8,
    fontsize_col = 8,
    color_palette = colorRampPalette(c("blue", "white", "red"))(50)
) {
  library(pheatmap)
  library(gridExtra)
  library(grid)
  
  # 默认 title
  if (is.null(titles)) {
    titles <- paste0("Matrix ", seq_along(matrix_list))
  }
  
  if (length(titles) != length(matrix_list)) {
    stop("Length of titles must match length of matrix_list.")
  }
  
  # 画出所有 heatmap（存为 grob）
  heatmap_grobs <- lapply(seq_along(matrix_list), function(i) {
    mat <- matrix_list[[i]]
    
    if (!is.null(row_labels)) {
      if (length(row_labels) != nrow(mat)) stop("row_labels length mismatch")
      rownames(mat) <- row_labels
    }
    
    if (!is.null(col_labels)) {
      if (length(col_labels) != ncol(mat)) stop("col_labels length mismatch")
      colnames(mat) <- col_labels
    }
    
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col,
      main = titles[i],
      color = color_palette,
      cellwidth = cellwidth,
      cellheight = cellheight,
      silent = TRUE  # 关键！返回 grob 而不是直接画
    )$gtable
  })
  
  # 设置 ncol 默认值
  if (is.null(ncol)) ncol <- length(matrix_list)
  
  # 合并成一个图
  combined <- gridExtra::grid.arrange(
    grobs = heatmap_grobs,
    ncol = ncol,
    top = if (!is.null(global_title)) grid::textGrob(global_title, gp = grid::gpar(fontsize = 16)) else NULL
  )
  
  # 如果指定 filename 就保存
  if (!is.null(filename)) {
    ggsave(filename, plot = combined, width = 2.5 * length(matrix_list), height = 6)
  }
  
  return(combined)
}

draw_multiple_heatmaps(
  matrix_list = list(C_ivw, C_adivw, C_naive_MRrr, C_MRrr,
                     C_MRrr_regularized, C_MRrr_sparse, C_Mr_DAG),
  titles = c("IVW","adIVW","Naive MR-rr","MR-rr","Reg. MR-rr","Sparse MR-rr","Mr-DAG"),
  filename = "results/Causal_C_estimation_hm_260713_postselection.png",
  ncol = 3
)


# Loading matrices B for the three MR-rr estimators.
pathway_name <- paste0("Pathway ", 1:r_RR)
draw_multiple_heatmaps(
  matrix_list = list(MRrr_AB_display$B,
                     Reg_MRrr_AB_display$B,
                     Sparse_MRrr_AB_display$B),
  titles = c("MR-rr B", "Reg. MR-rr B", "Sparse MR-rr B (refit)"),
  row_labels = pathway_name,
  col_labels = exp_name,
  filename = "results/B_loading_heatmaps_260713.png",
  ncol = 1,
  cellwidth = 30,
  cellheight = 25
)


# ============================================================ #
# 6. bootsrap confidence interval ####
# ============================================================ #
set.seed(123)
# bootstrap
bootstrap_C_list <- vector("list", 7)
names(bootstrap_C_list) <- c("ivw", "adivw", "naive_MRrr", "MRrr",
                             "MRrr_regularized", "MRrr_sparse", "Mr_DAG")
B <- 1000  # bootstrap次数

for (b in 1:B) {
  boot_index <- sample(1:nrow(gamma_exp_j), replace = TRUE)
  boot_result <- estimate_C(
    gamma_out_j = gamma_out_j[boot_index, , drop = FALSE],
    gamma_exp_j = gamma_exp_j[boot_index, , drop = FALSE],
    se_out_j = se_out_j[boot_index, , drop = FALSE],
    se_exp_j = se_exp_j[boot_index, , drop = FALSE],
    cor_out = cor_out,
    cor_exp = cor_exp,
    r_RR = r_RR,
    sparse_support = sparse_support,
    sparse_A_init = real_result$A_MRrr_sparse,
    sparse_B_init = real_result$B_MRrr_sparse,
    print_rank = FALSE
  )
  for (k in names(bootstrap_C_list)) {
    if (is.null(bootstrap_C_list[[k]])) {
      bootstrap_C_list[[k]] <- array(boot_result[[paste0("C_", k)]], dim = c(3, 9, B))
    } else {
      bootstrap_C_list[[k]][,,b] <- boot_result[[paste0("C_", k)]]
    }
  }
  # 每隔10报告一次进度
  if (b %% 10 == 0) {
    cat("Bootstrap iteration:", b, "\n")
  }
}


# save bootstrap_C_list
save(bootstrap_C_list, file = "bootstrap_C_list_1000_260713_postselection.RData")
# load bootstrap_C_list
load("bootstrap_C_list_1000_260713_postselection.RData")


# ============================================================ #
# 7. plot big boxplot - with sparse ####
# ============================================================ #
exp_name <- read.csv('data/traits_1e-4.csv')$x
out_name <- c("LAS", "CES", "SVS")
names(out_name) <- 0:2
names(exp_name) <- 0:8

estimate_C_list <- list(
  ivw              = real_result$C_ivw,
  adivw            = real_result$C_adivw,
  naive_MRrr       = real_result$C_naive_MRrr,
  MRrr             = real_result$C_MRrr,
  MRrr_regularized = real_result$C_MRrr_regularized,
  MRrr_sparse      = real_result$C_MRrr_sparse,
  Mr_DAG           = real_result$C_Mr_DAG
)

est_levels <- c("ivw","adivw","naive_MRrr","MRrr",
                "MRrr_regularized","MRrr_sparse","Mr_DAG")
est_labels <- c("IVW","SRIVW","Naive MR-rr","MR-rr",
                "Reg. MR-rr","Sparse MR-rr","MrDAG")

build_long_and_sig <- function(bootstrap_C_list, estimate_C_list, ci_level = 0.95) {
  a <- (1 - ci_level) / 2
  long <- data.frame(); sig <- data.frame()
  for (i in 1:3) for (j in 1:9) for (k in names(bootstrap_C_list)) {
    v  <- bootstrap_C_list[[k]][i, j, ]; v <- v[!is.na(v)]
    pe <- estimate_C_list[[k]][i, j]
    lo <- quantile(v, a); up <- quantile(v, 1 - a)
    long <- rbind(long, data.frame(Row = i-1, Col = j-1,
                                   estimator = k, value = v))
    sig  <- rbind(sig,  data.frame(Row = i-1, Col = j-1, estimator = k,
                                   point_est = pe, lower = lo, upper = up,
                                   sig = (lo > 0) | (up < 0)))
  }
  long$estimator <- factor(long$estimator, levels = est_levels, labels = est_labels)
  sig$estimator  <- factor(sig$estimator,  levels = est_levels, labels = est_labels)
  list(long = long, sig = sig)
}

dat <- build_long_and_sig(bootstrap_C_list, estimate_C_list, ci_level = 0.95)


# new
plot_realdata_halfeye <- function(long_df, sig_df, out_name, exp_name,
                                  save_path = NULL,
                                  width = 13, height = 6,
                                  ylim = c(-0.5, 0.5)) {
  p <- ggplot(long_df, aes(x = estimator, y = value, fill = estimator)) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
    stat_halfeye(
      adjust = 1, .width = 0.95, justification = -0.2, scale = 0.8,
      point_colour = NA, interval_colour = "black", interval_size = 0.75,
      slab_alpha = 0.7, normalize = "xy"
    ) +
    geom_point(
      data = sig_df, inherit.aes = FALSE,
      aes(x = estimator, y = point_est, color = sig),
      shape = 18, size = 3
    ) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    coord_cartesian(ylim = ylim) +
    facet_grid(rows = vars(Row), cols = vars(Col),
               labeller = labeller(Row = out_name, Col = exp_name)) +
    theme_bw(base_size = 12) +
    labs(x = "", y = "Estimated C Value", fill = "Estimator") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text  = element_text(size = 10),
          legend.position = "right")
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to ", save_path)
  }
  print(p)
}


plot_realdata_halfeye(
  dat$long, dat$sig, out_name, exp_name,
  ylim = c(-0.5, 0.5),
  save_path = "results/real_data_CI_halfeye_260713_postselection.png"
)









