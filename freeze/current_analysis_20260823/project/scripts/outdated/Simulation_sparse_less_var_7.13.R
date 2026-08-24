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
setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
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

  result_ivw <- ivw_multiple_outcomes(GAMMA_hat, gamma_hat, Sigma_X, Sigma_Y)

  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB

  A_d_hat_r = result_d_r$A
  B_d_hat_r = result_d_r$B
  AB_d_hat_r = result_d_r$AB
  return(list(A_hat, B_hat, AB_hat,
              A_d_hat, B_d_hat, AB_d_hat,
              A_d_hat_r, B_d_hat_r, AB_d_hat_r,
              gamma_hat, GAMMA_hat,
              result_ivw))
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
  py = 4 # py can be changed to any number from 1 to 9
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
  
  sigma_Gamma_j = as.matrix(lip_data_real[,paste0('se_out',1:py)])
  sigma_Gamma_j = sigma_Gamma_j * sqrt(var_Z_raw)
  Corr_Y = lip_corr_real[10:13,10:13]
  Corr_Y = as.matrix(Corr_Y)
  
  Sigma_Yj = lapply(1:pz_lip_data, function(j)
    diag(sigma_Gamma_j[j,]) %*% Corr_Y %*% diag(sigma_Gamma_j[j,]))
  Sigma_Yj_sample = Sigma_Yj[z_index] # Q 2025: same index
  Sigma_Y_temp = lapply(1:pz, function(j) Sigma_Yj_sample[[j]]) # Q 2025: changed from lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d_Y <- array(unlist(Sigma_Y_temp), dim = c(4, 4, length(Sigma_Y_temp)))
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


admm_psd_projection <- function(Sigma_hat, mu = 1, epsilon = 1e-6, max_iter = 10, tol = 1e-6) {
  p <- nrow(Sigma_hat)
  B <- matrix(0, p, p)
  Lambda <- matrix(0, p, p)
  
  project_to_psd <- function(Z, epsilon) {
    eig <- eigen(Z, symmetric = TRUE)
    eigvals <- pmax(eig$values, epsilon)
    Z_proj <- eig$vectors %*% diag(eigvals) %*% t(eig$vectors)
    return((Z_proj + t(Z_proj)) / 2)
  }
  
  vecl <- function(M) {
    idx <- lower.tri(M, diag = TRUE)
    return(M[idx])
  }
  
  matl <- function(x, p) {
    M <- matrix(0, p, p)
    idx <- lower.tri(M, diag = TRUE)
    M[idx] <- x
    M <- M + t(M)
    diag(M) <- diag(M) / 2
    return(M)
  }
  
  soft_threshold <- function(x, tau) {
    return(sign(x) * pmax(abs(x) - tau, 0))
  }
  
  for (iter in 1:max_iter) {
    A <- project_to_psd(B + Sigma_hat + mu * Lambda, epsilon)
    
    tmp <- A - Sigma_hat - mu * Lambda
    vec_tmp <- vecl(tmp)
    vec_shrink <- soft_threshold(vec_tmp, mu)
    B <- matl(vec_tmp - vec_shrink, p)
    
    Lambda <- Lambda - (A - B - Sigma_hat) / mu
    
    r_norm <- norm(A - B - Sigma_hat, type = "F")
    if (r_norm < tol) {
      break
    }
  }
  
  return(A)
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
  res_sparse <- mr_rr_sparse(GAMMA_hat, gamma_hat, W = W, Sigma_X = Sigma_X, lambda = lambda, r = r_RR)
  
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
py <- 4
r <- 2

var_Z = 2 * lip_data_real$ImpMAF * (1 - lip_data_real$ImpMAF)

pz_lip_data = length(var_Z)

# get W
sigma_Gamma_j = as.matrix(lip_data_real[,paste0('se_out',1:py)])
sigma_Gamma_j = sigma_Gamma_j * sqrt(var_Z)
Corr_Y = lip_corr_real[10:13,10:13]
Corr_Y = as.matrix(Corr_Y)

Sigma_Yj = lapply(1:pz_lip_data, function(j)
  diag(sigma_Gamma_j[j,]) %*% Corr_Y %*% diag(sigma_Gamma_j[j,]))
array_3d_Y <- array(unlist(Sigma_Yj), dim = c(4, 4, length(Sigma_Yj)))
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
  B_sparse[nonzero_rows, j] <- rnorm(length(nonzero_rows), mean = 0, sd = 20)
}

# Construct C = A %*% B
C <- A %*% B_sparse

norm(C)


#### 0713 new ver ####
run_simulation_sparse_compare <- function(C, B_sparse, lambda = rep(2e-3, ncol(C)), r_RR = 2, eloop = 100) {
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

res = run_simulation_sparse_compare(C, B_sparse, lambda = rep(2e-3, 9), eloop = 1000)

# save
save(res, file = "results/simulation_sparse_result_250713.RData")
load("results/simulation_sparse_result_250713.RData")


# plot C
plot_sparse_estimation_comparison <- function(res,
                                              filename = sprintf("results/sparse_comparison_boxplot_%s.png", format(Sys.Date(), "%y%m%d")),
                                              width = 12, height = 6) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(readr)
  
  # Ensure directory exists
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  
  # 标签准备
  out_name <- c("IS", "LAS", "CES", "SVS")
  names(out_name) <- 0:3
  exp_name <- read.csv("data/traits_1e-4.csv")$x
  names(exp_name) <- 0:8
  
  C_true <- res$C_true
  px <- ncol(C_true)
  py <- nrow(C_true)
  N <- ncol(res$result_naive_est)
  
  # 构造 entry grid: (Row, Col) 按列优先展开顺序
  entry_grid <- expand.grid(Row = 0:(py - 1), Col = 0:(px - 1))
  entry_grid <- entry_grid[order(entry_grid$Col, entry_grid$Row), ]
  
  # 拼接所有 estimator 的数据
  df_all <- purrr::imap_dfr(list(
    `Naive` = res$result_naive_est,
    `MR-rr` = res$result_mrr_est,
    `Regularized` = res$result_mrr_r_est,
    `Sparse` = res$result_sparse_est
  ), function(mat, method) {
    tibble(
      Entry = rep(1:(py * px), each = N),
      Value = as.vector(t(mat)),
      Row = rep(entry_grid$Row, each = N),
      Col = rep(entry_grid$Col, each = N),
      Method = method,
      TrueValue = rep(as.vector(C_true), each = N)
    )
  })
  
  # True value for horizontal lines
  df_true <- tibble(
    Row = entry_grid$Row,
    Col = entry_grid$Col,
    TrueValue = as.vector(C_true)
  )
  
  # 绘图
  p <- ggplot(df_all, aes(x = Method, y = Value, fill = Method)) +
    geom_boxplot(outlier.size = 0.3) +
    geom_hline(data = df_true, aes(yintercept = TrueValue),
               color = "red", linetype = "dashed", linewidth = 0.1) +
    facet_grid(
      rows = vars(Row),
      cols = vars(Col),
      labeller = labeller(Row = out_name, Col = exp_name),
      scales = "free_y"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, size = 7),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 8),
          legend.position = "right") +
    ylab("Estimated Value") +
    xlab("") +
    ggtitle("Comparison of Estimated C Entries (Red Dashed Line = True C)")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}


plot_sparse_estimation_compare(res)


plot_sparse_B_estimation <- function(res, 
                                     filename = sprintf("results/sparse_B_boxplot_%s.png", format(Sys.Date(), "%y%m%d")), 
                                     width = 12, height = 4) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(purrr)
  
  B_hat_matrix <- res$result_B_sparse   # r * px × N
  B_true <- res$B_sparse                # r × px
  
  r <- nrow(B_true)
  px <- ncol(B_true)
  N <- ncol(B_hat_matrix)
  
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
  
  df_true <- tibble(
    Row = entry_grid$Row,
    Col = entry_grid$Col,
    TrueValue = as.vector(B_true),
    is_zero = as.vector(B_true == 0)
  )
  
  p <- ggplot(df_long, aes(x = "", y = Value)) +
    geom_boxplot(fill = "lightgreen", outlier.size = 0.3) +
    geom_hline(data = df_true, 
               aes(yintercept = TrueValue, color = is_zero),
               linetype = "dashed", linewidth = 0.8,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "gray60", `FALSE` = "red")) +
    facet_grid(
      rows = vars(Row),
      cols = vars(Col),
      labeller = labeller(Row = c("Passway 1","Passway 2"), Col = exp_name),
      scale = "free_y"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 8),
          legend.position = "none") +
    ylab("Estimated B Entry") +
    xlab("") +
    ggtitle("Simulation Boxplots of Estimated B Entries\n(Red = Nonzero True B, Gray = Zero True B)")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}

plot_sparse_B_estimation <- function(res, 
                                     filename = sprintf("results/sparse_B_boxplot_%s.png", format(Sys.Date(), "%y%m%d")), 
                                     width = 12, height = 4,
                                     trim_quantile = 0) {
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(purrr)
  
  B_hat_matrix <- res$result_B_sparse   # r * px × N
  B_true <- res$B_sparse                # r × px
  
  r <- nrow(B_true)
  px <- ncol(B_true)
  N <- ncol(B_hat_matrix)
  
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
  
  # 去除极端值（可选）
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
  
  p <- ggplot(df_long, aes(x = "", y = Value)) +
    geom_boxplot(fill = "lightgreen", outlier.size = 0.3) +
    geom_hline(data = df_true, 
               aes(yintercept = TrueValue, color = is_zero),
               linetype = "dashed", linewidth = 0.8,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_color_manual(values = c(`TRUE` = "gray60", `FALSE` = "red")) +
    facet_grid(
      rows = vars(Row),
      cols = vars(Col),
      labeller = labeller(Row = c("Passway 1","Passway 2"), Col = exp_name),
      scales = "free"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 8),
          legend.position = "none") +
    ylab("Estimated B Entry") +
    xlab("") +
    ggtitle("Simulation Boxplots of Estimated B Entries\n(Red = Nonzero True B, Gray = Zero True B)")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}


plot_sparse_B_estimation(res, trim_quantile = 0)


#### old ver: only sparse methos ####
simulation_sparse_once <- function(C, lambda = rep(2e-3, 9), r_RR = 2, n = 1000) {
  # # test
  # set.seed(123)
  # lambda = rep(2e-3, 9)
  # r_RR = 2
  # n = 1000
  
  px <- ncol(C)
  py <- nrow(C)
  
  # === Step 1: Parameter setup ===
  parameters <- .get_parameters(C, me_weight = 1, effect_weight = 1, r_RR = r_RR)
  W <- parameters$weight.matrix
  W_sqrt <- .sqrt_matrix(W)
  Sigma_X <- parameters$Sigma_X
  Sigma_Y <- parameters$Sigma_Y
  VX_tilde <- parameters$VX_tilde
  
  # === Step 2: Simulate data ===
  gamma_j_star <- MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star <- gamma_j_star %*% t(C)
  
  gamma_hat <- matrix(0, n, px)
  GAMMA_hat <- matrix(0, n, py)
  for (j in 1:n) {
    gamma_hat[j, ] <- MASS::mvrnorm(1, mu = gamma_j_star[j, ], Sigma = Sigma_X, tol = 100)
    GAMMA_hat[j, ] <- MASS::mvrnorm(1, mu = Gamma_j_star[j, ], Sigma = Sigma_Y, tol = 100)
  }
  
  # === Step 3: gamma tilde transformation ===
  P_GAMMA <- GAMMA_hat %*% solve(t(GAMMA_hat) %*% GAMMA_hat) %*% t(GAMMA_hat)
  P_GAMMA_prep <- diag(1, n) - P_GAMMA
  Sigma_gammahat <- t(gamma_hat) %*% gamma_hat / n
  matrix_part1 <- Sigma_gammahat - t(P_GAMMA %*% gamma_hat) %*% (P_GAMMA %*% gamma_hat) / n
  matrix_part2 <- matrix_part1 - Sigma_X
  
  if (!is_psd(matrix_part2)) {
    cat("matrix_part2 not PSD — projecting...\n")
    matrix_part2 <- admm_psd_projection(matrix_part2, epsilon = 1e-6)
  }
  
  R <- chol(matrix_part1)
  Q <- chol(matrix_part2)
  L <- solve(R) %*% Q
  gamma_tilde <- P_GAMMA %*% gamma_hat + (P_GAMMA_prep %*% gamma_hat) %*% L
  
  # === Step 4: Initialization ===
  result_init <- mr_rr(GAMMA_hat, gamma_hat, r = r_RR, W = W, Sigma_X = Sigma_X)
  A_hat <- result_init$A
  B_hat <- result_init$B
  
  iter_num <- 100
  for (iter in 1:iter_num) {
    B_var <- CVXR::Variable(rows = r_RR, cols = px)
    loss <- sum_squares(t(A_hat) %*% W %*% t(GAMMA_hat) - B_var %*% t(gamma_tilde)) / n +
      Reduce(`+`, lapply(1:px, function(k) {
        lambda[k] * sum_entries(norm1(B_var[, k, drop = FALSE]))
      })) + 1e-6 * sum_entries(B_var)
    
    prob <- Problem(Minimize(loss))
    result <- solve(prob, solver = "OSQP")
    B_hat <- result$getValue(B_var)
    
    svd_result <- svd(B_hat %*% t(gamma_tilde) %*% GAMMA_hat %*% W_sqrt)
    A_hat <- solve(W_sqrt) %*% svd_result$v %*% t(svd_result$u)
    
    dist_to_C <- norm(A_hat %*% B_hat - C, "F") / norm(C, "F")
    if (dist_to_C < 0.01) break
  }
  
  B_hat[abs(B_hat) < 1e-2] <- 0
  C_hat <- A_hat %*% B_hat
  
  return(list(C_hat = C_hat, A_hat = A_hat, B_hat = B_hat))
}


run_simulation_sparse_prediction <- function(B_sparse, eloop = 100, lambda = rep(2e-3, 9)) {
  # # test
  # eloop = 5
  # lambda = rep(2e-3, 9)
  
  C_hat_matrix <- matrix(NA, py * px, eloop)
  B_hat_matrix <- matrix(NA, r * px, eloop)
  
  for (i in 1:eloop) {
    sim <- simulation_sparse_once(C = C, lambda = lambda, r_RR = r)
    C_hat_matrix[, i] <- as.vector(sim$C_hat)
    B_hat_matrix[, i] <- as.vector(sim$B_hat)
    cat("[", i, "] Done\n")
  }
  
  # sensitivity and specificity of B
  # Step 1: flatten true B to binary vector
  # True binary B
  B_true_vec <- as.numeric(B_sparse != 0)         # length = r * px
  B_hat_bin_mat <- apply(B_hat_matrix, 2, function(col) as.numeric(abs(col) > 1e-3))  # [r*px, eloop]
  
  # Repeat true B across columns to match dimensions
  B_true_mat <- matrix(B_true_vec, nrow = r * px, ncol = eloop)
  
  # Compute total TP, FP, TN, FN over all entries & sims
  TP_total <- sum((B_hat_bin_mat == 1) & (B_true_mat == 1))
  FP_total <- sum((B_hat_bin_mat == 1) & (B_true_mat == 0))
  TN_total <- sum((B_hat_bin_mat == 0) & (B_true_mat == 0))
  FN_total <- sum((B_hat_bin_mat == 0) & (B_true_mat == 1))
  
  # Compute sensitivity and specificity
  sensitivity <- TP_total / (TP_total + FN_total)
  specificity <- TN_total / (TN_total + FP_total)
  
  
  return(list(
    result_sparse_est = C_hat_matrix,
    result_sparse_B = B_hat_matrix,
    C_true = C,
    B_sparse = B_sparse,
    sensitivity = sensitivity,
    specificity = specificity
  ))
}


res = run_simulation_sparse_prediction(B_sparse, eloop = 100, lambda = rep(2e-3, 9))



#### plot ####
exp_name = read.csv('data/traits_1e-4.csv')$x
out_name = c("IS", "LAS", "CES", "SVS")
names(out_name) <- 0:3
names(exp_name) <- 0:8


plot_sparse_estimation <- function(res, 
                                   filename = sprintf("results/sparse_boxplot_%s.png", format(Sys.Date(), "%y%m%d")), 
                                   width = 12, height = 6, ylim = c(-0.5, 0.5)) {
  est_matrix <- res$result_sparse_est   # 36 x N
  C_true <- res$C_true                  # 4 x 9
  
  py <- nrow(C_true)
  px <- ncol(C_true)
  N <- ncol(est_matrix)
  
  # 标签
  out_name <- c("IS", "LAS", "CES", "SVS")
  names(out_name) <- 0:(py - 1)
  exp_name <- read.csv("data/traits_1e-4.csv")$x
  names(exp_name) <- 0:(px - 1)
  
  # 构造 entry grid: 按列展开顺序生成 (row, col)
  entry_grid <- expand.grid(
    Row = 0:(py - 1),
    Col = 0:(px - 1)
  )
  entry_grid <- entry_grid[order(entry_grid$Col, entry_grid$Row), ]
  
  # reshape data
  df_long <- tibble(
    Entry = rep(1:(py * px), each = N),
    Value = as.vector(t(est_matrix)),
    Row = rep(entry_grid$Row, each = N),
    Col = rep(entry_grid$Col, each = N),
    TrueValue = rep(as.vector(C_true), each = N)
  )
  
  df_true <- tibble(
    Row = entry_grid$Row,
    Col = entry_grid$Col,
    TrueValue = as.vector(C_true)
  )
  
  # plot
  p <- ggplot(df_long, aes(x = "", y = Value)) +
    geom_boxplot(fill = "skyblue", outlier.size = 0.3) +
    geom_hline(data = df_true, aes(yintercept = TrueValue),
               color = "red", linetype = "dashed", linewidth = 0.8,
               inherit.aes = FALSE) +
    facet_grid(
      rows = vars(Row),
      cols = vars(Col),
      labeller = labeller(Row = out_name, Col = exp_name)
    ) +
    coord_cartesian(ylim = ylim) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(size = 7),
          strip.text = element_text(size = 8),
          legend.position = "none") +
    ylab("Estimated Value") +
    xlab("") +
    ggtitle("Simulation Boxplots of Estimated C Entries (Red Dashed Line = True C)")
  
  ggsave(filename, plot = p, width = width, height = height, dpi = 300)
  message("Plot saved to ", filename)
  print(p)
}



plot_sparse_estimation(res)




