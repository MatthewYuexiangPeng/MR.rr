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
library(scales)


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
.second_moment <- function(X) {
  crossprod(X) / nrow(X)
}


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


.align_AB_to_reference <- function(A, B, B_ref) {
  r <- nrow(B)
  
  # rank 1: only sign ambiguity
  if (r == 1) {
    if (sum((B - B_ref)^2) > sum((-B - B_ref)^2)) {
      B <- -B
      A <- -A
    }
    return(list(A = A, B = B, AB = A %*% B))
  }
  
  # Current real-data sensitivity analysis only uses r = 2
  if (r != 2) {
    stop(".align_AB_to_reference currently supports r = 1 or 2.")
  }
  
  permutations <- list(c(1, 2), c(2, 1))
  signs <- list(
    c( 1,  1),
    c( 1, -1),
    c(-1,  1),
    c(-1, -1)
  )
  
  best_loss <- Inf
  best_A <- NULL
  best_B <- NULL
  
  for (perm in permutations) {
    for (sgn in signs) {
      
      B_try <- B[perm, , drop = FALSE]
      A_try <- A[, perm, drop = FALSE]
      
      B_try <- diag(sgn) %*% B_try
      A_try <- A_try %*% diag(sgn)
      
      loss <- sum((B_try - B_ref)^2)
      
      if (loss < best_loss) {
        best_loss <- loss
        best_A <- A_try
        best_B <- B_try
      }
    }
  }
  
  list(
    A = best_A,
    B = best_B,
    AB = best_A %*% best_B
  )
}


estimate_C <- function(gamma_exp_j,se_exp_j,gamma_out_j,se_out_j,cor_exp,cor_out,
                       r_RR = NULL, 
                       sparse_eta = 1e-3,
                       sparse_support = NULL,
                       sparse_A_init = NULL, 
                       sparse_B_init = NULL,
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
  SigmaXX = .second_moment(gamma_exp_j)
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
      (gamma_out_j - gamma_exp_j %*% t(C_hat)) %*% W_sqrt,
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
      lambda = rep(sparse_eta, px),
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
    opt_rate                 = opt_rate,
    sparse_eta               = sparse_eta
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
# Select eta on the original real data ####
# ============================================================ #
set.seed(123)

px <- ncol(gamma_exp_j)
py <- ncol(gamma_out_j)
pz <- nrow(gamma_exp_j)

Sigma_X_sum <- matrix(0, px, px)
Sigma_Y_sum <- matrix(0, py, py)

for (j in seq_len(pz)) {
  
  D_exp_j <- diag(se_exp_j[j, ])
  D_out_j <- diag(se_out_j[j, ])
  
  Sigma_X_sum <- Sigma_X_sum +
    D_exp_j %*% cor_exp %*% D_exp_j
  
  Sigma_Y_sum <- Sigma_Y_sum +
    D_out_j %*% cor_out %*% D_out_j
}

Sigma_X_real <- Sigma_X_sum / pz
Sigma_Y_real <- Sigma_Y_sum / pz

W_real <- solve(Sigma_Y_real)
W_sqrt_real <- .sqrt_matrix(W_real)

selected_rank <- rank_test_M0(
  W       = W_real,
  y_j_hat = gamma_out_j,
  x_j_hat = gamma_exp_j,
  Sigma_X = Sigma_X_real,
  print   = TRUE
)

# sensitivity test: set rank = 2
# selected_rank = 2

# real-data eta path
eta_grid <- c(
  1e-4,
  3e-4,
  5e-4,
  7e-4,
  1e-3,
  1.2e-3,
  1.5e-3,
  2e-3,
  2.5e-3,
  3e-3,
  4e-3,
  5e-3,
  1e-2
)

zero_tol <- 1e-2

# Selected eta values
selected_eta_rank1 <- 1.2e-3
selected_eta_rank2 <- 1e-3

compute_real_eta_path <- function(rank_use) {
  
  eta_path_list <- vector("list", length(eta_grid))
  eta_fit_list  <- vector("list", length(eta_grid))
  
  for (i in seq_along(eta_grid)) {
    
    eta <- eta_grid[i]
    
    fit_eta <- mr_rr_sparse(
      GAMMA_hat = gamma_out_j,
      gamma_hat = gamma_exp_j,
      W         = W_real,
      Sigma_X   = Sigma_X_real,
      lambda    = rep(eta, px),
      r         = rank_use,
      max_iter  = 100,
      tol       = 1e-2
    )
    
    B_hat <- if (!is.null(fit_eta$B_raw)) {
      fit_eta$B_raw
    } else {
      fit_eta$B
    }
    
    C_hat <- if (!is.null(fit_eta$AB_raw)) {
      fit_eta$AB_raw
    } else {
      fit_eta$A %*% B_hat
    }
    
    residual_term <- norm(
      (gamma_out_j -
         gamma_exp_j %*% t(C_hat)) %*%
        W_sqrt_real,
      type = "F"
    )^2 / pz
    
    debias_term <- sum(diag(
      W_sqrt_real %*%
        C_hat %*%
        Sigma_X_real %*%
        t(C_hat) %*%
        W_sqrt_real
    ))
    
    eta_path_list[[i]] <- data.frame(
      rank      = rank_use,
      eta       = eta,
      objective = residual_term - debias_term,
      zero_prop = mean(abs(B_hat) < zero_tol),
      n_nonzero = sum(abs(B_hat) >= zero_tol)
    )
    
    eta_fit_list[[i]] <- fit_eta
  }
  
  list(
    path = bind_rows(eta_path_list),
    fits = eta_fit_list
  )
}


eta_result_rank1 <- compute_real_eta_path(rank_use = 1)
eta_result_rank2 <- compute_real_eta_path(rank_use = 2)

eta_path_rank1 <- eta_result_rank1$path
eta_path_rank2 <- eta_result_rank2$path

eta_fit_rank1 <- eta_result_rank1$fits
eta_fit_rank2 <- eta_result_rank2$fits


eta_label <- function(x) {
  lab <- formatC(x, format = "e", digits = 1)
  lab <- sub("\\.0e", "e", lab)
  lab <- sub("e-0", "e-", lab)
  lab
}


make_real_eta_plot <- function(
    eta_path_df,
    rank_use,
    selected_eta
) {
  
  p_obj <- ggplot(
    eta_path_df,
    aes(x = eta, y = objective)
  ) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 2) +
    geom_vline(
      xintercept = selected_eta,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    scale_x_log10(
      breaks = eta_grid,
      labels = NULL
    ) +
    labs(
      title = paste0("Rank ", rank_use),
      x = NULL,
      y = "Debiased objective"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold"
      ),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    )
  
  p_sparse <- ggplot(
    eta_path_df,
    aes(x = eta, y = zero_prop)
  ) +
    geom_step(
      direction = "hv",
      linewidth = 0.5
    ) +
    geom_point(size = 2) +
    geom_vline(
      xintercept = selected_eta,
      linetype = "dashed",
      linewidth = 0.6
    ) +
    scale_x_log10(
      breaks = eta_grid,
      labels = eta_label
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = label_percent()
    ) +
    labs(
      x = expression(eta),
      y = "Zero proportion in B"
    ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 8
      )
    )
  
  p_obj / p_sparse +
    plot_layout(heights = c(1, 1))
}


plot_rank1 <- make_real_eta_plot(
  eta_path_df = eta_path_rank1,
  rank_use = 1,
  selected_eta = selected_eta_rank1
)

plot_rank2 <- make_real_eta_plot(
  eta_path_df = eta_path_rank2,
  rank_use = 2,
  selected_eta = selected_eta_rank2
)

eta_selection_plot_both <- wrap_plots(
  plot_rank1,
  plot_rank2,
  ncol = 2
)

eta_selection_plot_both

ggsave(
  filename = "results/real_data_sparse_eta_selection_rank1_rank2_260729.png",
  plot = eta_selection_plot_both,
  width = 14,
  height = 7.5,
  units = "in",
  dpi = 300,
  bg = "white"
)


# ============================================================ #
# 4. Estimate C using different methods ####
# ============================================================ #
# set final eta
selected_rank <- 1
selected_eta  <- 1.2e-3
analysis_tag  <- "rank1"

# sensitivity analysis, when rank =2, corresponding eta:
# selected_rank <- 2
# selected_eta  <- 1e-3
# analysis_tag  <- "rank2"

eta_tag <- formatC(selected_eta, format = "e", digits = 1)
eta_tag <- gsub("\\.", "p", eta_tag)
eta_tag <- gsub("e-0", "e-", eta_tag)
eta_tag <- paste0("eta_", eta_tag)

out_file <- function(stem, ext) {
  file.path(
    "results",
    paste0(stem, "_260820_", analysis_tag, "_", eta_tag, ext)
  )
}


real_result <- estimate_C(
  gamma_out_j = gamma_out_j,
  gamma_exp_j = gamma_exp_j,
  se_out_j = se_out_j,
  se_exp_j = se_exp_j,
  cor_out = cor_out,
  cor_exp = cor_exp,
  r_RR        = selected_rank,
  sparse_eta  = selected_eta
)

save(
  real_result,
  file = out_file("real_result", ".RData")
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
  out_file("selected_rank", ".txt")
)

# A and B used for interpretation. Sign canonicalization is for display only.
MRrr_AB_display <- .canonicalize_AB(
  real_result$A_MRrr,
  real_result$B_MRrr
)

Reg_MRrr_AB_display <- .align_AB_to_reference(
  real_result$A_MRrr_regularized,
  real_result$B_MRrr_regularized,
  B_ref = MRrr_AB_display$B
)

Sparse_MRrr_AB_display <- .align_AB_to_reference(
  real_result$A_MRrr_sparse,
  real_result$B_MRrr_sparse,
  B_ref = MRrr_AB_display$B
)

stopifnot(
  max(abs(MRrr_AB_display$AB -
            real_result$C_MRrr)) < 1e-8,
  
  max(abs(Reg_MRrr_AB_display$AB -
            real_result$C_MRrr_regularized)) < 1e-8,
  
  max(abs(Sparse_MRrr_AB_display$AB -
            real_result$C_MRrr_sparse)) < 1e-8
)

pathway_name <- paste0("Pathway ", seq_len(selected_rank))

rownames(MRrr_AB_display$A) <- out_name
rownames(Reg_MRrr_AB_display$A) <- out_name
rownames(Sparse_MRrr_AB_display$A) <- out_name

colnames(MRrr_AB_display$A) <- pathway_name
colnames(Reg_MRrr_AB_display$A) <- pathway_name
colnames(Sparse_MRrr_AB_display$A) <- pathway_name

rownames(MRrr_AB_display$B) <- pathway_name
rownames(Reg_MRrr_AB_display$B) <- pathway_name
rownames(Sparse_MRrr_AB_display$B) <- pathway_name

colnames(MRrr_AB_display$B) <- exp_name
colnames(Reg_MRrr_AB_display$B) <- exp_name
colnames(Sparse_MRrr_AB_display$B) <- exp_name

# table for A
A_table <- data.frame(Outcome = out_name)

for (k in seq_len(selected_rank)) {
  A_table[[paste0("MR-rr Pathway ", k)]] <-
    MRrr_AB_display$A[, k]
  
  A_table[[paste0("Reg. MR-rr Pathway ", k)]] <-
    Reg_MRrr_AB_display$A[, k]
  
  A_table[[paste0("Sparse MR-rr Pathway ", k)]] <-
    Sparse_MRrr_AB_display$A[, k]
}

A_table[-1] <- lapply(A_table[-1], function(x) round(x, 4))
print(A_table)

write.csv(
  A_table,
  out_file("A_table", ".csv"),
  row.names = FALSE
)

# table for B
B_table <- data.frame(Protein = exp_name)

for (k in seq_len(selected_rank)) {
  B_table[[paste0("MR-rr Pathway ", k)]] <-
    MRrr_AB_display$B[k, ]
  
  B_table[[paste0("Reg. MR-rr Pathway ", k)]] <-
    Reg_MRrr_AB_display$B[k, ]
  
  B_table[[paste0("Sparse MR-rr Pathway ", k)]] <-
    Sparse_MRrr_AB_display$B[k, ]
}

B_table[-1] <- lapply(B_table[-1], function(x) round(x, 3))
print(B_table)

write.csv(
  B_table,
  out_file("B_table", ".csv"),
  row.names = FALSE
)

# write.csv(MRrr_AB_display$A,
#           out_file("A_MRrr", ".csv"), row.names = FALSE)
# 
# write.csv(MRrr_AB_display$B,
#           out_file("B_MRrr", ".csv"), row.names = FALSE)
# 
# write.csv(Reg_MRrr_AB_display$A,
#           out_file("A_MRrr_regularized", ".csv"), row.names = FALSE)
# 
# write.csv(Reg_MRrr_AB_display$B,
#           out_file("B_MRrr_regularized", ".csv"), row.names = FALSE)
# 
# write.csv(Sparse_MRrr_AB_display$A,
#           out_file("A_MRrr_sparse_refit", ".csv"), row.names = FALSE)
# 
# write.csv(Sparse_MRrr_AB_display$B,
#           out_file("B_MRrr_sparse_refit", ".csv"), row.names = FALSE)

write.csv(real_result$B_MRrr_sparse_selection,
          out_file("B_MRrr_sparse_selection", ".csv"), row.names = FALSE)

write.csv(1 * sparse_support,
          out_file("B_MRrr_sparse_support", ".csv"), row.names = FALSE)


# ============================================================ #
# 5. bootsrap confidence interval ####
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
save(
  bootstrap_C_list,
  file = out_file("bootstrap_C_list_1000_postselection", ".RData")
)

# load bootstrap_C_list
# load(out_file("bootstrap_C_list_1000_postselection", ".RData"))


# ============================================================ #
# 6. plot big boxplot - with sparse ####
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
  save_path = out_file("real_data_CI_halfeye_postselection", ".png")
)

