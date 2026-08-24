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

# # regress lip_data$gamma_out on lip_data$gamma_exp1 to lip_data$gamma_exp24
# mean(abs(lm(lip_data$gamma_out ~ ., data = lip_data[,paste0('gamma_exp',1:24)])$coefficients))
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


is_psd <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  all(eigenvalues >= 0)
}


#### generate low rank true C (sparse) ####
px = 9
py = 4
r = 2
# r =  min(px,py) # maximum rank of true C
# U = matrix(rnorm(py*r), py, r)
# V = matrix(rnorm(px*r), r, px)
# # assume true C to be rank r=3 for now. C=C^(r)
# # use modified THM 2.1 to define true C^(r), (it is equivalent to define it as A_d*B_d)
# eignvalue_matrix = diag(c(rep(0.1, r_approx),
#                           rep(0, r-r_approx))) # set at 0.001 if want C to be rank > r
# C <- U %*% eignvalue_matrix %*% V

# generate C from A times sparse B 
# Step 1: Construct A with orthonormal columns (py x r)
A <- qr.Q(qr(matrix(rnorm(py * r), py, r)))  # A^T A = I

# Step 2: Construct sparse B (r x px)
B <- matrix(0, nrow = r, ncol = px)

for (j in 1:px) {
  # 每列至少一个非零，在两行中随机选择1或2个位置赋值
  nonzero_rows <- sample(1:r, size = sample(1:2, 1))
  B[nonzero_rows, j] <- rnorm(length(nonzero_rows), mean = 1, sd = 0.5)
}

# Step 3: Construct C = A %*% B
C <- A %*% B


#### set universal parameters ####
me_weight_list = c(2.5, 1)
effect_weight_list = c(0.25, 1)
save_sim_filename = "results/simulate_result_sparse_250603.RData"


#### sparse simulation ####
# 25.3.9 sparse simulation
parameters = .get_parameters(C, me_weight=1, effect_weight=1, r_RR = 2)
pz = n = 1000 #(pz)

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


# sample data
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


P_GAMMA = GAMMA_hat %*% solve(t(GAMMA_hat)%*%GAMMA_hat)%*%t(GAMMA_hat)
P_GAMMA_prep = diag(1, nrow = n, ncol = n) - P_GAMMA

Sigma_gammahat_gammahat_hat = t(gamma_hat) %*% gamma_hat / n

matrix_part1 = Sigma_gammahat_gammahat_hat - 
  t(P_GAMMA %*% gamma_hat) %*% (P_GAMMA %*% gamma_hat)/n

matrix_part2 = matrix_part1 - Sigma_X

# check PSD
is_psd(matrix_part1)
is_psd(matrix_part2)
# true, dont need the ADMM

R = chol(matrix_part1)
Q = chol(matrix_part2)

L = solve(R) %*% Q

gamma_tilde = P_GAMMA %*% gamma_hat + (P_GAMMA_prep %*% gamma_hat) %*% L

# check
norm(t(gamma_tilde) %*% gamma_tilde / n - (Sigma_gammahat_gammahat_hat - Sigma_X))
norm(t(gamma_tilde) %*% GAMMA_hat - t(gamma_hat) %*% GAMMA_hat)


#### method1: CVXR ####
# start to implement algorithm
# step 1: optimizing the objective function
# store iteration values of A and B
lambda = rep(2e-3, parameters$px)

A_hat_list = list()
B_hat_list = list()
# initial values of A_hat and B_hat using dRRR
result <- mr_rr(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
A_hat = result$A
B_hat = result$B
A_hat_list[[1]] = A_hat
B_hat_list[[1]] = B_hat

# start the iteration
iter_num = 100
for (iter in 1:iter_num) {
  # given A minimize B: CVXR ----
  A_hat = A_hat_list[[iter]]
  # 设定优化变量 B_hat (r_RR x px)
  B_hat <- Variable(rows = r_RR, cols = px)  # B_hat 是 r_RR × px 的变量
  # 计算目标函数的第一部分 (二次误差项)
  # first_term <- (1 / px) * Reduce("+", lapply(1:px, function(j) {
  #   Gamma_tilde_j <- tilde_Gamma[, j, drop = FALSE]  # 取 tilde_Gamma 的第 j 列 (py x 1)
  #   gamma_tilde_j <- tilde_gamma[, j, drop = FALSE]  # 取 tilde_gamma 的第 j 列 (px x 1)
  # 
  #   # 计算残差
  #   # residual <- W_sqrt %*% (Gamma_tilde_j - A_hat %*% B_hat %*% gamma_tilde_j)
  #   residual <- t(A_hat) %*% W %*% Gamma_tilde_j - B_hat %*% gamma_tilde_j
  #   
  #   sum_entries(norm2(residual)^2)  # 这里返回 CVXR 表达式
  # }))
  first_term <- sum_squares(t(A_hat) %*% W %*% t(GAMMA_hat) - B_hat %*% t(gamma_tilde)) / pz
  # 计算目标函数的第二部分 (L1 正则化)
  second_term <- Reduce("+", lapply(1:px, function(k) {
    lambda[k] * sum_entries(norm1(B_hat[, k, drop = FALSE]))  # L1 正则化
  }))
  + 1e-6 * sum_entries(B_hat)  # 添加一个小的 L2 正则化，有助于收敛
  # 设定最终目标函数
  objective <- first_term + second_term
  # 设定优化问题
  prob <- Problem(Minimize(objective))
  # 求解优化问题
  result <- solve(prob, solver = "OSQP") # "ECOS"    "ECOS_BB" "SCS"     "OSQP" 
  # 获取优化后的 B_hat
  B_hat_opt <- result$getValue(B_hat)
  B_hat_list = c(B_hat_list, list(B_hat_opt))


  # given B minimize A:  Orthogonal Procrustes ----
  B_hat = B_hat_opt

  svd_result <- svd(B_hat %*% t(gamma_tilde) %*% GAMMA_hat %*% W_sqrt)
  U <- svd_result$u
  V <- svd_result$v

  A_hat_opt <- solve(W_sqrt) %*% V %*% t(U)

  # regularization, make t(A_hat) %*% W %*% A_hat= I

  A_hat_list = c(A_hat_list, list(A_hat_opt))


  # print current iteration and "loss"

  # # 计算 dist to C
  # dist_to_C <- norm(A_hat_opt %*% B_hat_opt - C, type = "F")  # 使用 Frobenius 范数
  dist_to_C <- norm(A_hat_opt %*% B_hat_opt - C, type = "F") / norm(C, "F")
  
  cat("Iteration:", iter,
      "| Dist to C:", round(dist_to_C, 4),
      "| Zeros in B_hat:", sum(B_hat < 1e-2), "\n")

  if (dist_to_C < 0.1) {
    break
  }
}

A_hat = A_hat_opt
B_hat = B_hat_opt



# set small entry to 0
B_hat[abs(B_hat) < 1e-2] = 0

norm(A_hat %*% B_hat - C, type = "F") / norm(C, "F")


norm(t(A_hat) %*% W %*% A_hat - diag(1, r_RR))

norm(A_hat_opt%*%B_hat_opt-C)



result <- mr_rr(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X)

norm(result$A %*% result$B - C, type = "F") / norm(C, "F")



#### visualization ####
plot_C_comparison <- function(true_C, C_MR_regu, C_MR_sparse) {
  stopifnot(all(dim(true_C) == dim(C_MR_regu)),
            all(dim(C_MR_regu) == dim(C_MR_sparse)))
  
  px <- nrow(true_C)
  py <- ncol(true_C)
  
  methods <- c("True C", "MR_regu", "MR_sparse")
  
  df <- expand.grid(Row = 1:px, Col = 1:py, Method = methods)
  df$Value <- c(as.vector(true_C),
                as.vector(C_MR_regu),
                as.vector(C_MR_sparse))
  
  df$Method <- factor(df$Method, levels = c("True C", "MR_regu", "MR_sparse"))
  
  row_labels <- paste0("Row", 1:px)
  col_labels <- paste0("Col", 1:py)
  names(row_labels) <- 1:px
  names(col_labels) <- 1:py
  
  ggplot(df) +
    geom_segment(aes(x = 0.5, xend = 3.5, y = Value, yend = Value,
                     color = Method, linetype = Method),
                 size = 0.8) +
    facet_grid(Row ~ Col, labeller = labeller(Row = row_labels, Col = col_labels)) +
    scale_color_manual(values = c("black", "red", "blue")) +
    scale_linetype_manual(values = c("solid", "dashed", "dashed")) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(size = 6),     # ✅ 显示纵轴刻度
          axis.ticks.y = element_line(color = "black", size = 0.3),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size = 10),
          legend.position = "bottom") +
    labs(title = "C Entry-wise Estimates Comparison for Regularized MR and Sparse MR",
         x = NULL, y = NULL)
  
  ggsave("results/Sparse_C_comparison_plot.png", width = 12, height = 6)
}


result_r <- mr_rr_regularized(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X, )
C_MR_regu = result_r$AB
plot_C_comparison(true_C = C,
                  C_MR_regu = C_MR_regu,
                  C_MR_sparse = A_hat %*% B_hat)


#### TRY IMPLEMENT BY COORDINATE DESCENT ####
lambda = rep(1e-3, px)  # L1 正则化参数
# 设定迭代次数
iter_num = 100
max_cd_iter = 10  # Coordinate Descent 最大迭代次数
B_hat_list = list()
A_hat_list = list()

# 初始化 A_hat 和 B_hat
result <- mr_rr(GAMMA_hat, gamma_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
A_hat = result$A
B_hat = result$B
A_hat_list[[1]] = A_hat
B_hat_list[[1]] = B_hat

# Soft-thresholding function for L1 shrinkage
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}

for (iter in 1:iter_num) {

  # **Step 1: Given A, optimize B using Coordinate Descent**
  A_hat = A_hat_list[[iter]]
  B_hat = B_hat_list[[iter]]  # 初始化 B_hat

  # Coordinate Descent on B_hat
  iter_cd = 1
  max_change = 1000
  while (iter_cd < max_cd_iter && max_change > 1e-2) {
    prev_B_hat <- B_hat  # 记录上次 B_hat 用于收敛判断
    for (k in 1:px) {
      update_k = 0
      for (j in 1:pz) {
        r_j_k = t(A_hat) %*% W %*% GAMMA_hat[j,] - B_hat[, -k] %*% as.matrix(gamma_tilde[j, -k])
        update_k = update_k + r_j_k * gamma_tilde[j, k]
        # sum = rowSums(B_hat * matrix(tilde_gamma[j, ], nrow = r_RR, ncol = px, byrow = TRUE))
        # r_j = t(A_hat) %*% W %*% tilde_Gamma[,j] - sum
        # r_j_times_tilde_gamma = r_j * tilde_gamma[j, k]
        # update_k = update_k + r_j_times_tilde_gamma
      }
      B_hat[, k] = soft_threshold(B_hat[, k] + update_k/pz, lambda[k])  # L1 Shrinkage
    }

    # 计算收敛条件
    max_change = max(abs(B_hat - prev_B_hat))
    iter_cd = iter_cd + 1
  }

  B_hat_list[[iter + 1]] = B_hat  # 存储优化后的 B_hat

  # **Step 2: Given B, optimize A using Orthogonal Procrustes**
  svd_result <- svd(B_hat %*% t(gamma_tilde) %*% GAMMA_hat %*% W_sqrt)
  U <- svd_result$u
  V <- svd_result$v
  A_hat_opt <- solve(W_sqrt) %*% V %*% t(U)

  # 存储优化后的 A_hat
  A_hat_list[[iter + 1]] = A_hat_opt

  # **计算 dist to C**
  # dist_to_C <- norm(A_hat_opt %*% B_hat - C, type = "F")  # 使用 Frobenius 范数
  # mean abs
  dist_to_C <- norm(A_hat_opt %*% B_hat - C, "F") / norm(C, "F")

  cat("Iteration:", iter,
      "| Dist to C:", round(dist_to_C, 4),
      "| Zeros in B_hat:", sum(B_hat == 0), "\n")

  # 如果收敛，则提前终止
  if (dist_to_C < 1e-3) {
    break
  }
}

A_hat = A_hat_opt
B_hat = B_hat_list[[iter + 1]]







