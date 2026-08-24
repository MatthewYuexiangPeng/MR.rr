#### loading environment and data ####
rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap::pheatmap)
# library(readxl)
# library(glmnet)
library(CVXR)
library(ADMM)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(purrr)


set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

lip_data_real = read.csv('data/dat_1e-4.csv')
lip_corr_real = read.csv('data/rho_mat_1e-4.csv')
n_Y = c(1296908,1241207,1245612,1241619)

exp_name = read.csv('data/traits_1e-4.csv')$x
out_name = c("IS", "LAS", "CES", "SVS")
names(out_name) <- 0:3
names(exp_name) <- 0:8

source("scripts/MR_rr_estimators.R")


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


.simulation = function(parameters, regularized = TRUE, regularization_rate = 1e-13) {
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

  # sample gamma_hat, GAMMA_hat
  gamma_hat = matrix(0, n, px)
  GAMMA_hat = matrix(0, n, py)
  for (j in 1:n) {
    gamma_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    GAMMA_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }

  # compute A_hat, B_hat
  result <- mr_rr_naive(GAMMA_hat, gamma_hat, r=r_RR, W=W) # TODO: changed here need check the result
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB # sample level estimator A_hat * B_hat

  # need to subtract the Sigma_xx by Sigma_x. Sigma_x can be estimated by the regression of exp ~ z (have been approximate in the get parameter function)
  # if (regularized == FALSE) {
  #   result_d <- mr_rr(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  # } else if (regularized == TRUE){
  #   result_d <- mr_rr_regularized(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)
  # } else {
  result_d <- mr_rr(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  result_d_r <- mr_rr_regularized(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)

  # result_ivw <- ivw_multiple_outcomes(GAMMA_hat, gamma_hat, Sigma_X, Sigma_Y)

  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB

  A_d_hat_r = result_d_r$A
  B_d_hat_r = result_d_r$B
  AB_d_hat_r = result_d_r$AB
  return(list(A_hat, B_hat, AB_hat,
              A_d_hat, B_d_hat, AB_d_hat,
              A_d_hat_r, B_d_hat_r, AB_d_hat_r,
              gamma_hat, GAMMA_hat #,result_ivw
              ))
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


# is_psd <- function(matrix) {
#   eigenvalues <- eigen(matrix)$values
#   all(eigenvalues >= 0)
# }


# move to estimators.R
# admm_psd_projection <- function(Sigma_hat, mu = 1, epsilon = 1e-6, max_iter = 10, tol = 1e-6) {
#   p <- nrow(Sigma_hat)
#   B <- matrix(0, p, p)
#   Lambda <- matrix(0, p, p)
#   
#   project_to_psd <- function(Z, epsilon) {
#     eig <- eigen(Z, symmetric = TRUE)
#     eigvals <- pmax(eig$values, epsilon)
#     Z_proj <- eig$vectors %*% diag(eigvals) %*% t(eig$vectors)
#     return((Z_proj + t(Z_proj)) / 2)
#   }
#   
#   vecl <- function(M) {
#     idx <- lower.tri(M, diag = TRUE)
#     return(M[idx])
#   }
#   
#   matl <- function(x, p) {
#     M <- matrix(0, p, p)
#     idx <- lower.tri(M, diag = TRUE)
#     M[idx] <- x
#     M <- M + t(M)
#     diag(M) <- diag(M) / 2
#     return(M)
#   }
#   
#   soft_threshold <- function(x, tau) {
#     return(sign(x) * pmax(abs(x) - tau, 0))
#   }
#   
#   for (iter in 1:max_iter) {
#     A <- project_to_psd(B + Sigma_hat + mu * Lambda, epsilon)
#     
#     tmp <- A - Sigma_hat - mu * Lambda
#     vec_tmp <- vecl(tmp)
#     vec_shrink <- soft_threshold(vec_tmp, mu)
#     B <- matl(vec_tmp - vec_shrink, p)
#     
#     Lambda <- Lambda - (A - B - Sigma_hat) / mu
#     
#     r_norm <- norm(A - B - Sigma_hat, type = "F")
#     if (r_norm < tol) {
#       break
#     }
#   }
#   
#   return(A)
# }


.align_AB_to_ref <- function(Ahat, Bhat, Bref) {
  stopifnot(is.matrix(Ahat), is.matrix(Bhat), is.matrix(Bref))
  stopifnot(nrow(Bhat) == nrow(Bref), ncol(Bhat) == ncol(Bref))
  stopifnot(ncol(Ahat) == nrow(Bhat))
  
  r <- nrow(Bhat)
  
  # 所有行置换
  perms <- gtools::permutations(r, r)
  
  # 所有符号组合 (±1)
  sign_grid <- as.matrix(expand.grid(rep(list(c(-1, 1)), r)))
  
  best_err <- Inf
  best_A <- Ahat
  best_B <- Bhat
  
  for (i in seq_len(nrow(perms))) {
    p <- perms[i, ]
    
    # 行置换：B 行交换，A 列交换
    Bhp <- Bhat[p, , drop = FALSE]
    Ahp <- Ahat[, p, drop = FALSE]
    
    for (j in seq_len(nrow(sign_grid))) {
      s <- sign_grid[j, ]
      S <- diag(as.numeric(s), r, r)
      
      Bc <- S %*% Bhp          # 行符号翻转
      Ac <- Ahp %*% S          # 列符号翻转
      
      err <- sum((Bc - Bref)^2)
      
      if (err < best_err) {
        best_err <- err
        best_A <- Ac
        best_B <- Bc
      }
    }
  }
  
  list(A = best_A, B = best_B)
}



.simulation_sparse_compare <- function(parameters, lambda = rep(2e-3, parameters$px), r_RR = 2) {
  n <- 1000
  py <- parameters$py
  px <- parameters$px
  C <- parameters$C
  W <- parameters$weight.matrix
  Sigma_X <- parameters$Sigma_X
  Sigma_Y <- parameters$Sigma_Y
  VX_tilde <- parameters$VX_tilde
  
  # --- Step 1: simulate gamma_star and GAMMA_star ---
  gamma_star <- MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  GAMMA_star <- gamma_star %*% t(C)
  
  gamma_hat <- matrix(0, n, px)
  GAMMA_hat <- matrix(0, n, py)
  for (j in 1:n) {
    gamma_hat[j, ] <- MASS::mvrnorm(1, mu = gamma_star[j, ], Sigma = Sigma_X, tol = 100)
    GAMMA_hat[j, ] <- MASS::mvrnorm(1, mu = GAMMA_star[j, ], Sigma = Sigma_Y, tol = 100)
  }
  
  # --- Step 2: Naive MR-RR ---
  res_naive <- mr_rr_naive(GAMMA_hat, gamma_hat, r = r_RR, W = W)
  
  # --- Step 3: Standard MR-RR ---
  res_standard <- mr_rr(GAMMA_hat, gamma_hat, r = r_RR, W = W, Sigma_X = Sigma_X)
  
  # --- Step 4: Regularized MR-RR ---
  res_regularized <- mr_rr_regularized(GAMMA_hat, gamma_hat, r = r_RR, W = W, Sigma_X = Sigma_X, regularization_rate = 1e-13)
  
  # --- Step 5: Sparse MR-RR ---
  res_sparse <- mr_rr_sparse(GAMMA_hat, gamma_hat, W = W, Sigma_X = Sigma_X, lambda = lambda, r = r_RR, max_iter = 100)
  
  # # --- Align sparse (A,B) to regularized B (reference) to fix sign/permutation ambiguity ---
  # if (r_RR == 2) {
  #   aligned <- .align_AB_to_ref(
  #     Ahat = res_sparse$A,
  #     Bhat = res_sparse$B,
  #     Bref = res_regularized$B
  #   )
  #   res_sparse$A <- aligned$A
  #   res_sparse$B <- aligned$B
  #   res_sparse$AB <- res_sparse$A %*% res_sparse$B
  # } else {
  #   # r>2: 需要更一般的对齐（Hungarian + sign），先不在这里展开
  # }
  
  return(list(
    A_naive = res_naive$A, B_naive = res_naive$B, AB_naive = res_naive$AB,
    A_standard = res_standard$A, B_standard = res_standard$B, AB_standard = res_standard$AB,
    A_regularized = res_regularized$A, B_regularized = res_regularized$B, AB_regularized = res_regularized$AB,
    A_sparse = res_sparse$A, B_sparse = res_sparse$B, AB_sparse = res_sparse$AB
  ))
}



#### generate low rank C with sparse B ####
set.seed(123)

px <- 9
py <- 3
r <- 2

var_Z = 2 * lip_data_real$ImpMAF * (1 - lip_data_real$ImpMAF)

pz_lip_data = length(var_Z)

# get W
sigma_Gamma_j = as.matrix(lip_data_real[,paste0('se_out',2:(py + 1))])
sigma_Gamma_j = sigma_Gamma_j * sqrt(var_Z)
Corr_Y = lip_corr_real[11:13,11:13]
Corr_Y = as.matrix(Corr_Y)

Sigma_Yj = lapply(1:pz_lip_data, function(j)
  diag(sigma_Gamma_j[j,]) %*% Corr_Y %*% diag(sigma_Gamma_j[j,]))
array_3d_Y <- array(unlist(Sigma_Yj), dim = c(3, 3, length(Sigma_Yj)))
Sigma_Y <- apply(array_3d_Y, c(1, 2), mean)

W = solve(Sigma_Y)


# Construct A with W-orthonormal columns: A^T W A = I
# Step 1: generate orthonormal Q (Q^T Q = I)
Q <- qr.Q(qr(matrix(rnorm(py * r), py, r)))  # Q: py x r

# Step 2: make it W-orthonormal via inverse sqrt transformation
W_sqrt_inv <- solve(.sqrt_matrix(W))
A <- W_sqrt_inv %*% Q  # A^T W A = Q^T W_sqrt^T W W_sqrt Q = Q^T Q = I


t(A) %*% W %*% A  # should be close to identity matrix

# Construct sparse B
B_sparse <- matrix(0, r, px)
for (j in 1:px) {
  nonzero_rows <- sample(1:r, size = 1)
  B_sparse[nonzero_rows, j] <- rnorm(length(nonzero_rows), mean = 0, sd = 25)
}

# Construct C = A %*% B
C <- A %*% B_sparse

norm(C)
svd(C)$d

#### 260712 simulation & plots ####
run_simulation_sparse_compare <- function(C, B_sparse, lambda = rep(1e-3, ncol(C)), r_RR = 2, eloop = 100) {
  py <- nrow(C)
  px <- ncol(C)
  
  # allocate storage
  mat_dim <- py * px
  result_AB_naive <- matrix(NA, mat_dim, eloop)
  result_AB_standard <- matrix(NA, mat_dim, eloop)
  result_AB_regularized <- matrix(NA, mat_dim, eloop)
  result_AB_sparse <- matrix(NA, mat_dim, eloop)
  
  result_B_sparse <- matrix(NA, r_RR * px, eloop)
  
  # get shared parameters (fixed once)
  parameters <- .get_parameters(C, me_weight = 1, effect_weight = 1, r_RR = r_RR)
  # B_sparse <- parameters$B  # true sparse B
  B_true_vec <- as.numeric(B_sparse != 0)  # logical vector
  
  for (i in 1:eloop) {
    sim <- .simulation_sparse_compare(parameters = parameters, lambda = lambda, r_RR = r_RR)
    
    result_AB_naive[, i] <- as.vector(sim$AB_naive)
    result_AB_standard[, i] <- as.vector(sim$AB_standard)
    result_AB_regularized[, i] <- as.vector(sim$AB_regularized)
    result_AB_sparse[, i] <- as.vector(sim$AB_sparse)
    
    result_B_sparse[, i] <- as.vector(sim$B_sparse)
    
    cat("[", i, "] Done\n")
  }
  
  # ==== Compute sensitivity & specificity for sparse B ====
  B_hat_bin_mat <- apply(result_B_sparse, 2, function(col) as.numeric(abs(col) > 1e-3))
  B_true_mat <- matrix(B_true_vec, nrow = r_RR * px, ncol = eloop)
  
  TP_total <- sum((B_hat_bin_mat == 1) & (B_true_mat == 1))
  FP_total <- sum((B_hat_bin_mat == 1) & (B_true_mat == 0))
  TN_total <- sum((B_hat_bin_mat == 0) & (B_true_mat == 0))
  FN_total <- sum((B_hat_bin_mat == 0) & (B_true_mat == 1))
  
  sensitivity <- TP_total / (TP_total + FN_total)
  specificity <- TN_total / (TN_total + FP_total)
  
  return(list(
    C_true = C,
    B_sparse = B_sparse,
    result_AB_naive = result_AB_naive,
    result_AB_standard = result_AB_standard,
    result_AB_regularized = result_AB_regularized,
    result_AB_sparse = result_AB_sparse,
    result_B_sparse = result_B_sparse,
    sensitivity = sensitivity,
    specificity = specificity
  ))
}

set.seed(123)
res = run_simulation_sparse_compare(C, B_sparse, lambda = rep(1e-3, 9), eloop = 1000)

# save
save(res, file = "results/simulation_sparse_result_060314_1e-3.RData")
load("results/simulation_sparse_result_060314_1e-3.RData")


canonicalize_B <- function(B) {
  r <- nrow(B)
  for (k in 1:r) {
    j_star <- which.max(abs(B[k, ]))
    if (B[k, j_star] < 0) B[k, ] <- -B[k, ]
  }
  B
}


# problem: the sign of some row maybe wrong. plot the canonicalize_B and B_sparse
plot_sparse_B_estimation_canonicalized <- function(res, 
                                                   filename = sprintf("results/sparse_B_boxplot_canonicalized_%s.png", format(Sys.Date(), "%y%m%d")), 
                                                   width = 12, height = 4,
                                                   trim_quantile = 0) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  B_hat_matrix <- res$result_B_sparse   # (r*px) x N
  B_true <- res$B_sparse                # r x px
  
  r <- nrow(B_true)
  px <- ncol(B_true)
  N <- ncol(B_hat_matrix)
  
  # ---------------------------
  # Canonicalize: Btrue + Bsparse(hat)
  # ---------------------------
  # 1) canonicalize sign by row
  B_true_can <- canonicalize_B(B_true)
  
  # 2) get a fixed row order from B_true_can
  key_j_true <- apply(abs(B_true_can), 1, which.max)
  key_v_true <- apply(abs(B_true_can), 1, max)
  ord_true <- order(key_j_true, -key_v_true)
  
  # 3) apply the same row order to B_true and to every B_hat replicate
  B_true_can <- B_true_can[ord_true, , drop = FALSE]
  
  B_hat_can <- matrix(NA_real_, nrow = r * px, ncol = N)
  for (n in 1:N) {
    Bn <- matrix(B_hat_matrix[, n], nrow = r, ncol = px)
    Bn <- canonicalize_B(Bn)
    Bn <- Bn[ord_true, , drop = FALSE]
    B_hat_can[, n] <- as.vector(Bn)
  }
  B_hat_matrix <- B_hat_can
  B_true <- B_true_can
  
  # ---------------------------
  # Build plotting data
  # ---------------------------
  entry_grid <- expand.grid(
    Row = 0:(r - 1),
    Col = 0:(px - 1)
  )
  entry_grid <- entry_grid[order(entry_grid$Col, entry_grid$Row), ]
  
  df_long <- tibble(
    Entry = rep(1:(r * px), each = N),
    Value = as.vector(t(B_hat_matrix)),
    Row = rep(entry_grid$Row, each = N),
    Col = rep(entry_grid$Col, each = N),
    TrueValue = rep(as.vector(B_true), each = N)
  )
  
  if (trim_quantile > 0) {
    df_long <- df_long %>%
      group_by(Row, Col) %>%
      filter(
        Value >= quantile(Value, trim_quantile / 2, na.rm = TRUE),
        Value <= quantile(Value, 1 - trim_quantile / 2, na.rm = TRUE)
      ) %>%
      ungroup()
  }
  
  df_true <- tibble(
    Row = entry_grid$Row,
    Col = entry_grid$Col,
    TrueValue = as.vector(B_true),
    is_zero = as.vector(B_true == 0)
  )
  
  # generic labels
  row_lab <- setNames(paste0("Pathway ", 1:r), as.character(0:(r - 1)))
  col_lab <- setNames(paste0("Exposure ", 1:px), as.character(0:(px - 1)))
  
  p <- ggplot(df_long, aes(x = "", y = Value)) +
    geom_boxplot(fill = "lightgreen", outlier.size = 0.3) +
    geom_hline(
      data = df_true,
      aes(yintercept = TrueValue, color = is_zero),
      linetype = "dashed", linewidth = 0.8,
      inherit.aes = FALSE, show.legend = FALSE
    ) +
    scale_color_manual(values = c(`TRUE` = "gray60", `FALSE` = "red")) +
    facet_grid(
      rows = vars(Row),
      cols = vars(Col),
      labeller = labeller(Row = row_lab, Col = col_lab),
      scales = "free"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 7),
      strip.text = element_text(size = 8),
      legend.position = "none"
    ) +
    ylab("Estimated B Entry") +
    xlab("") +
    ggtitle("Simulation Boxplots of Estimated B Entries\n(Red = Nonzero True B, Gray = Zero True B)")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
  print(B_true)
}


plot_sparse_B_estimation_canonicalized(res, trim_quantile = 0)


# plot C
# plot_sparse_estimation_comparison <- function(res,
#                                               filename = sprintf("results/sparse_comparison_boxplot_%s.png", format(Sys.Date(), "%y%m%d")),
#                                               width = 12, height = 6) {
#   # Ensure directory exists
#   dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
#   
#   # 标签准备
#   out_name <- c("LAS", "CES", "SVS")
#   names(out_name) <- 0:2
#   exp_name <- read.csv("data/traits_1e-4.csv")$x
#   names(exp_name) <- 0:8
#   
#   C_true <- res$C_true
#   px <- ncol(C_true)
#   py <- nrow(C_true)
#   N <- ncol(res$result_AB_naive)
#   
#   # 构造 entry grid: (Row, Col) 按列优先展开顺序
#   entry_grid <- expand.grid(Row = 0:(py - 1), Col = 0:(px - 1))
#   entry_grid <- entry_grid[order(entry_grid$Col, entry_grid$Row), ]
#   
#   # 拼接所有 estimator 的数据
#   df_all <- purrr::imap_dfr(list(
#     `MR-rr Naive` = res$result_AB_naive,
#     `MR-rr` = res$result_AB_standard,
#     `MR-rr Regularized` = res$result_AB_regularized,
#     `MR-rr Sparse` = res$result_AB_sparse
#   ), function(mat, method) {
#     tibble(
#       Entry = rep(1:(py * px), each = N),
#       Value = as.vector(t(mat)),
#       Row = rep(entry_grid$Row, each = N),
#       Col = rep(entry_grid$Col, each = N),
#       Method = method,
#       TrueValue = rep(as.vector(C_true), each = N)
#     )
#   })
#   
#   # True value for horizontal lines
#   df_true <- tibble(
#     Row = entry_grid$Row,
#     Col = entry_grid$Col,
#     TrueValue = as.vector(C_true)
#   )
#   
#   # 绘图
#   p <- ggplot(df_all, aes(x = Method, y = Value, fill = Method)) +
#     geom_boxplot(outlier.size = 0.3) +
#     geom_hline(data = df_true, aes(yintercept = TrueValue),
#                color = "red", linetype = "dashed", linewidth = 0.5) +
#     facet_grid(
#       rows = vars(Row),
#       cols = vars(Col),
#       labeller = labeller(Row = out_name, Col = exp_name),
#       scales = "free_y"
#     ) +
#     theme_bw(base_size = 10) +
#     theme(axis.text.x = element_text(angle = 90, size = 7),
#           axis.text.y = element_text(size = 7),
#           strip.text = element_text(size = 8),
#           legend.position = "right") +
#     ylab("Estimated Value") +
#     xlab("") +
#     ggtitle("Comparison of Estimated C Entries (Red Dashed Line = True C)")
#   
#   ggsave(filename, plot = p, width = width, height = height, dpi = 300)
#   message("Plot saved to ", filename)
#   print(p)
# }
# 
# 
# plot_sparse_estimation_comparison(res)


#### 260314 recovery rate ####
# recovery rate
compute_sparse_support_metrics <- function(res, threshold = 1e-1) {
  B_true <- res$B_sparse                  # r x px
  B_hat_mat <- res$result_B_sparse        # (r*px) x N
  
  true_vec <- as.numeric(as.vector(B_true) != 0)
  N <- ncol(B_hat_mat)
  p <- length(true_vec)
  
  # estimated support
  est_bin_mat <- apply(B_hat_mat, 2, function(x) as.numeric(abs(x) > threshold))
  est_bin_mat <- matrix(est_bin_mat, nrow = p, ncol = N)
  
  # replicate-level counts
  TP <- colSums(est_bin_mat == 1 & true_vec == 1)
  FP <- colSums(est_bin_mat == 1 & true_vec == 0)
  TN <- colSums(est_bin_mat == 0 & true_vec == 0)
  FN <- colSums(est_bin_mat == 0 & true_vec == 1)
  
  # exact recovery: one replicate gets every entry correct
  exact_recovery <- colSums(est_bin_mat == true_vec) == p
  
  # per-replicate metrics
  sensitivity_rep <- TP / (TP + FN)
  specificity_rep <- TN / (TN + FP)
  precision_rep   <- TP / (TP + FP)
  fdr_rep         <- FP / (TP + FP)
  
  # handle 0/0 cases
  precision_rep[TP + FP == 0] <- NA
  fdr_rep[TP + FP == 0] <- NA
  
  # pooled counts over all replicates
  TP_total <- sum(TP)
  FP_total <- sum(FP)
  TN_total <- sum(TN)
  FN_total <- sum(FN)
  
  sensitivity <- TP_total / (TP_total + FN_total)
  specificity <- TN_total / (TN_total + FP_total)
  precision   <- TP_total / (TP_total + FP_total)
  fdr         <- FP_total / (TP_total + FP_total)
  
  list(
    threshold = threshold,
    exact_recovery_rate = mean(exact_recovery),
    n_exact_recovery = sum(exact_recovery),
    
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    fdr = fdr,
    
    sensitivity_mean = mean(sensitivity_rep, na.rm = TRUE),
    specificity_mean = mean(specificity_rep, na.rm = TRUE),
    precision_mean = mean(precision_rep, na.rm = TRUE),
    fdr_mean = mean(fdr_rep, na.rm = TRUE),
    
    TP_total = TP_total,
    FP_total = FP_total,
    TN_total = TN_total,
    FN_total = FN_total,
    
    per_replicate = data.frame(
      TP = TP, FP = FP, TN = TN, FN = FN,
      exact_recovery = exact_recovery,
      sensitivity = sensitivity_rep,
      specificity = specificity_rep,
      precision = precision_rep,
      fdr = fdr_rep
    )
  )
}

load("results/simulation_sparse_result_060314_1e-3.RData")

metrics <- compute_sparse_support_metrics(res, threshold = 1e-3)

metrics$exact_recovery_rate
metrics$n_exact_recovery
metrics$sensitivity
metrics$fdr
metrics$specificity
metrics$precision


cat(sprintf(
  "Over 1000 replicates, exact support recovery occurred in %d replicates (%.1f%%). Sensitivity = %.3f, specificity = %.3f, precision = %.3f, and FDR = %.3f.\n",
  metrics$n_exact_recovery,
  100 * metrics$exact_recovery_rate,
  metrics$sensitivity,
  metrics$specificity,
  metrics$precision,
  metrics$fdr
))


# try more lambda
run_lambda_grid_sparse <- function(C, B_sparse, lambda_grid,
                                   r_RR = 2,
                                   eloop = 1000,
                                   save_dir = "results/lambda_grid_sparse",
                                   prefix = "simulation_sparse") {
  
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
  
  summary_list <- vector("list", length(lambda_grid))
  
  for (i in seq_along(lambda_grid)) {
    lam <- lambda_grid[i]
    lambda_vec <- rep(lam, ncol(C))
    
    cat(sprintf("\n[%d/%d] Running lambda = %.6g\n", i, length(lambda_grid), lam))
    
    res <- run_simulation_sparse_compare(
      C = C,
      B_sparse = B_sparse,
      lambda = lambda_vec,
      r_RR = r_RR,
      eloop = eloop
    )
    
    # 如果你已经定义了这个函数，就直接复用
    metrics <- compute_sparse_support_metrics(res, threshold = 1e-3)
    
    # 文件名里把科学计数法转成安全字符串
    lam_str <- gsub("\\+", "", format(lam, scientific = TRUE))
    lam_str <- gsub("-", "m", lam_str)
    lam_str <- gsub("\\.", "p", lam_str)
    
    save_file <- file.path(save_dir, sprintf("%s_lambda_%s.RData", prefix, lam_str))
    save(res, metrics, lam, file = save_file)
    
    cat(sprintf("Saved to %s\n", save_file))
    
    summary_list[[i]] <- data.frame(
      lambda = lam,
      exact_recovery_rate = metrics$exact_recovery_rate,
      n_exact_recovery = metrics$n_exact_recovery,
      sensitivity = metrics$sensitivity,
      specificity = metrics$specificity,
      precision = metrics$precision,
      fdr = metrics$fdr,
      TP_total = metrics$TP_total,
      FP_total = metrics$FP_total,
      TN_total = metrics$TN_total,
      FN_total = metrics$FN_total
    )
  }
  
  summary_df <- do.call(rbind, summary_list)
  
  summary_csv <- file.path(save_dir, sprintf("%s_summary.csv", prefix))
  write.csv(summary_df, summary_csv, row.names = FALSE)
  
  cat(sprintf("\nAll done. Summary saved to %s\n", summary_csv))
  
  return(summary_df)
}
# lambda_grid <- c(1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3, 1e-2)
lambda_grid <- 10^seq(-4, -2, by = 0.25)

set.seed(123)
summary_lambda <- run_lambda_grid_sparse(
  C = C,
  B_sparse = B_sparse,
  lambda_grid = lambda_grid,
  r_RR = 2,
  eloop = 1000,
  save_dir = "results/lambda_grid_sparse",
  prefix = "simulation_sparse_060223"
)

print(summary_lambda)
