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



set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

source("scripts/MR_rr_estimators.R")

# # regress lip_data$gamma_out on lip_data$gamma_exp1 to lip_data$gamma_exp24
# mean(abs(lm(lip_data$gamma_out ~ ., data = lip_data[,paste0('gamma_exp',1:24)])$coefficients))
# # the scale is similar to the generated C in .simulation


# hidden functions ----------------------------------------------------------------
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
  n = 2000
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

  # compute A_hat, B_hat
  result <- mr_rr_naive(y_j_hat, x_j_hat, r=r_RR, W=W) # TODO: changed here need check the result
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB # sample level estimator A_hat * B_hat

  # need to subtract the Sigma_xx by Sigma_x. Sigma_x can be estimated by the regression of exp ~ z (have been approximate in the get parameter function)
  # if (regularized == FALSE) {
  #   result_d <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  # } else if (regularized == TRUE){
  #   result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)
  # } else {
  result_d <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  result_d_r <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)

  result_ivw <- ivw_multiple_outcomes(y_j_hat, x_j_hat, Sigma_X, Sigma_Y)

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


.get_parameters = function(me_weight, r_RR = 5, px, r_approx = 5){
  # ## test
  # me_weight =0.2
  # r_RR = 5
  # px=24
  # r_approx = 5

  # r_RR is the rank chose by user when performing RRR
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 2000 # pz can be changed to any number
  py = 10 # py can be changed to any number from 1 to 24
  gamma_j = as.matrix(lip_data[,paste0('gamma_exp',1:24)])
  var_Z = 2 * lip_data$eaf * (1 - lip_data$eaf) # each Z is sum of two alleles, so var(Z) = var(Z^1+Z^2) = 2*var(Z^1), where Z1, Z2 ~ Binomial(eaf.outcome)
  var_Z_raw = var_Z

  z_index = sample(1:length(var_Z), pz, replace = TRUE)
  var_Z = var_Z[z_index]
  sqrt_var_Z = sqrt(var_Z)

  # TODO: consider use the same index here?
  lip_data_index_z <- sample(1:114, pz, replace = TRUE) # Question: Is this p_Z too large?
  gamma_j_sample <- gamma_j[lip_data_index_z,]
  gamma_j_star_temp <- gamma_j_sample*sqrt_var_Z

  # sample cov
  VX_tilde <- cov(gamma_j_star_temp)

  # # upscale VX_tilde
  # corr_upper_triangle = cor(gamma_j_star_temp)[upper.tri(cor(gamma_j_star_temp),
  #                                                            diag = FALSE)]
  # diag_VX_tilde = diag(VX_tilde)
  # diag_VX_tilde_sample = sample(diag_VX_tilde, px, replace = TRUE)
  # sample_size = px*(px+1)/2
  # VX_tilde_corr <- matrix(0, nrow = px, ncol = px)
  # for (i in 1:px-1){
  #   for (j in (i+1):px) {
  #     VX_tilde_corr[i,j] = sample(corr_upper_triangle, 1)
  #   }
  # }
  # VX_tilde_corr = VX_tilde_corr + t(VX_tilde_corr) + diag(1, px)
  # VX_tilde_upscale = .sqrt_matrix(diag(diag_VX_tilde_sample)) %*% VX_tilde_corr %*%
  #                    .sqrt_matrix(diag(diag_VX_tilde_sample))
  # VX_tilde_upscale <- VX_tilde_upscale*1/2 # Q: why?


  # Sigma_Xj & Sigma_X ----------------------------------------------------------------
  sigma_gamma_j = as.matrix(lip_data[,paste0('se_exp',1:24)])
  Sigma = lip_corr[1:24,1:24]
  # change "Sigma" to matrix
  Sigma = as.matrix(Sigma)
  sqrt_Sigma = .sqrt_matrix(Sigma)
  pz_lip_data = nrow(gamma_j)
  Sigma_Xj = lapply(1:pz_lip_data, function(j)
    diag(sigma_gamma_j[j,]) %*% Sigma %*% diag(sigma_gamma_j[j,]))
  Sigma_Xj_sample = Sigma_Xj[lip_data_index_z]
  Sigma_X_temp = lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d <- array(unlist(Sigma_X_temp), dim = c(24, 24, length(Sigma_X_temp)))
  Sigma_X <- apply(array_3d, c(1, 2), mean)

  # # upsclae Sigma_X
  # Sigma_X_triangle = Sigma[upper.tri(Sigma, diag = FALSE)]
  # sample_size = px*(px+1)/2
  # Sigma_X_corr <- matrix(0, nrow = px, ncol = px)
  # for (i in 1:px-1){
  #   for (j in (i+1):px) {
  #     Sigma_X_corr[i,j] = sample(Sigma_X_triangle, 1)
  #   }
  # }
  # Sigma_X_corr = Sigma_X_corr + t(Sigma_X_corr) + diag(1, px)
  #
  # sigma_gamma_j_upscale = sigma_gamma_j[,sample(24, px, replace = TRUE)]
  # Sigma_Xj_upscale = lapply(1:pz_lip_data, function(j) diag(sigma_gamma_j_upscale[j,]) %*% Sigma_X_corr %*% diag(sigma_gamma_j_upscale[j,]))
  # Sigma_Xj_upscale_sample = Sigma_Xj_upscale[tmp_index_z]
  # Sigma_X_temp_upscale = lapply(1:pz, function(j) Sigma_Xj_upscale_sample[[j]]*var_Z[j])
  # array_3d_upscale <- array(unlist(Sigma_X_temp_upscale), dim = c(px, px, length(Sigma_X_temp_upscale)))
  # Sigma_X_upscale <- apply(array_3d_upscale, c(1, 2), mean)
  # Sigma_X_upscale <- me_weight * Sigma_X_upscale

  ## old version
  # Sigma_X_triangle = Sigma[upper.tri(Sigma, diag = FALSE)]
  # diag_Sigma_X = diag(sigma_gamma_j)  # Q: problem
  # diag_Sigma_X_sample = sample(diag_Sigma_X, px, replace = TRUE)
  # sample_size = px*(px+1)/2
  # Sigma_X_corr <- matrix(0, nrow = px, ncol = px)
  # for (i in 1:px-1){
  #   for (j in (i+1):px) {
  #     Sigma_X_corr[i,j] = sample(Sigma_X_triangle, 1)
  #   }
  # }
  # Sigma_X_corr = Sigma_X_corr + t(Sigma_X_corr) + diag(1, px)
  # Sigma_X_upscale = .sqrt_matrix(diag(diag_Sigma_X_sample)) %*% Sigma_X_corr %*%
  #   .sqrt_matrix(diag(diag_Sigma_X_sample))
  # Sigma_X_upscale <- me_weight * Sigma_X_upscale

  # lets not upscale px for now.
  Sigma_X = me_weight * Sigma_X

  # Sigma_Y & weight.matrix -----------------------------------------------------------
  # Question: Use previous n_y?
  n_Y = median(lip_samplesize$samplesize) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
  # var_Y = sample(lip_data$se_out1[sample(1:114, 1119, replace = TRUE)]^2 * bmi.cad$N.outcome * var_Z_raw[sample(1:length(var_Z_raw), 1119, replace = TRUE)], py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
  var_Y = sample(lip_data$se_out1^2 * lip_samplesize$samplesize[sample(1:nrow(lip_samplesize), 114, replace = TRUE)]* var_Z_raw[sample(1:length(var_Z_raw), 114, replace = TRUE)],py)

  # var_Y = sample(bmi.cad$se.outcome^2 * bmi.cad$N.outcome * var_Z_raw, py) # Question: or use the se_out for lip_data?
  Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))

  # Sigma_Y = me_weight * Sigma_Y # TODO: should not weight Y?
  weight.matrix = solve(Sigma_Y)


  # SigmaXX, SigmaYX, SigmaXY, C ----------------------------------------------------------------
  # py * r matrix
  # TODO: how to approximate low rank?

  r =  min(px,py) # maximum rank of true C
  U = matrix(rnorm(py*r), py, r)
  V = matrix(rnorm(px*r), r, px)
  # assume true C to be rank r=3 for now. C=C^(r)
  # use modified THM 2.1 to define true C^(r), (it is equivalent to define it as A_d*B_d)
  eignvalue_matrix = diag(c(sample(c(sqrt(0.3), sqrt(0.2), sqrt(0.2)), r_approx, replace = TRUE),
                            rep(0, r-r_approx))) # set at 0.001 if want C to be rank > r
  C <- U %*% eignvalue_matrix %*% V


  # try: change the middle matrix of SVD to get the matrix rank of C
  # temp_matrix = matrix(rnorm(py*px), py, px)
  # # svd
  # svd_result = svd(temp_matrix)
  # U = svd_result$u
  # V = svd_result$v
  # d = svd_result$d
  # D_r = diag(c(d[1:r_approx],rep(0, length(d)-r_approx)))
  # C = U %*% D_r %*% t(V)
  # # rankMatrix(C)


  # # 25.3.9: generate C from A times B ----
  # M <- matrix(rnorm(py * r_approx), nrow = py)
  # svd_result <- svd(M)
  # U <- svd_result$u
  # A <- solve(.sqrt_matrix(weight.matrix)) %*% U[, 1:r_approx]
  # t(A) %*% weight.matrix %*% A
  # 
  # # generate a sparse B
  # num_ones_per_row = round(px * 0.2)
  # B <- matrix(0, nrow = r_approx, ncol = px)
  # 
  # for (i in 1:r_approx) {
  #   B[i, sample(1:px, num_ones_per_row)] <- 500
  # }
  # 
  # # 交换列以增加随机性
  # B <- B[, sample(1:px)]
  # 
  # C <- A %*% B
  # 
  # # norm(C, type = "F")
  # 
  # true.A_sparse = A
  # true.B_sparse = B
  # ####


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
                    VY_tilde = VY_tilde#,
                    #true.A_sparse = true.A_sparse, true.B_sparse = true.B_sparse
                    )
  return(parameters)
}
#####

# 25.3.9 sparse simulation
parameters = .get_parameters(me_weight=0.2, px=24)
pz = n = 2000 #(pz)
lambda = rep(1, parameters$px)

# generate B sparse para
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

# true.A_sparse = parameters$true.A_sparse
# true.B_sparse = parameters$true.B_sparse

# parameters$A_d %*% parameters$B_d - parameters$C
# # true.A_sparse %*% true.B_sparse - parameters$C
# mean(abs(A_hat %*% B_hat - C))
## [1] 0.2400043
# norm(A_hat %*% B_hat - C, type = "F")
## [1] 7.22151
# mean(abs(A_hat %*% B_hat - C))/ mean(abs(C))
# [1] 0.3152263


# sample data
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


# start to implement algorithm
# step 1: debiasing the objective function
# compute sample cov matrix
hat_Sigma_xhat_xhat = t(x_j_hat) %*% x_j_hat / n
hat_Sigma_x_x = hat_Sigma_xhat_xhat - Sigma_X

# check if PSD
is_psd <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  all(eigenvalues >= 0)
}


is_psd(hat_Sigma_x_x)
# chol decomposition
tilde_gamma = t(chol(hat_Sigma_x_x)) * sqrt(n)
# check if chol is correct
all.equal(tilde_gamma %*% t(tilde_gamma) / n, hat_Sigma_x_x)

tilde_Gamma = t(solve(tilde_gamma)%*%t(x_j_hat)%*%y_j_hat)

# check the definition of tilde_Gamma
all.equal(tilde_gamma %*% t(tilde_Gamma), t(x_j_hat)%*%y_j_hat)

# step 2: optimizing the objective function
# store iteration values of A and B
A_hat_list = list()
B_hat_list = list()
# initial values of A_hat and B_hat using dRRR
result <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
A_hat = result$A
B_hat = result$B
A_hat_list[[1]] = A_hat
B_hat_list[[1]] = B_hat

# start the iteration
iter_num = 50
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
  first_term <- (1 / px) * sum_squares(t(A_hat) %*% W %*% tilde_Gamma - B_hat %*% tilde_gamma)
  # 计算目标函数的第二部分 (L1 正则化)
  second_term <- Reduce("+", lapply(1:px, function(k) {
    lambda[k] * sum_entries(norm1(B_hat[, k, drop = FALSE]))  # L1 正则化
  }))
  # + 1e-6 * sum_entries(B_hat)  # 添加一个小的 L2 正则化，有助于收敛
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

  svd_result <- svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)
  U <- svd_result$u
  V <- svd_result$v

  A_hat_opt <- solve(W_sqrt) %*% V %*% t(U)

  # regularization, make t(A_hat) %*% W %*% A_hat= I

  A_hat_list = c(A_hat_list, list(A_hat_opt))


  # print current iteration and "loss"

  # # 计算 dist to C
  # dist_to_C <- norm(A_hat_opt %*% B_hat_opt - C, type = "F")  # 使用 Frobenius 范数
  dist_to_C = mean(abs(A_hat_opt %*% B_hat_opt - C))/ mean(abs(C))
  
  cat("Iteration:", iter,
      "| Dist to C:", round(dist_to_C, 4),
      "| Zeros in B_hat:", sum(B_hat == 0), "\n")

  if (dist_to_C < 0.1) {
    break
  }
}

A_hat = A_hat_opt
B_hat = B_hat_opt



# set small entry to 0
B_hat[abs(B_hat) < 1e-2] = 0

mean(abs(A_hat %*% B_hat - C))/ mean(abs(C))


norm(t(A_hat) %*% W %*% A_hat - diag(1, r_RR))

norm(A_hat_opt%*%B_hat_opt-C)


norm(B_hat - true.B_sparse, "F")




#### TRY IMPLEMENT BY COORDINATE DESCENT ----
lambda = rep(10, px)  # L1 正则化参数
# 设定迭代次数
iter_num = 10
max_cd_iter = 100  # Coordinate Descent 最大迭代次数
B_hat_list = list()
A_hat_list = list()

# 初始化 A_hat 和 B_hat
result <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
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
      for (j in 1:px) {
        r_j_k = t(A_hat) %*% W %*% tilde_Gamma[,j] - B_hat[, -k] %*% as.matrix(tilde_gamma[j, -k])
        update_k = update_k + r_j_k * tilde_gamma[j, k]
        # sum = rowSums(B_hat * matrix(tilde_gamma[j, ], nrow = r_RR, ncol = px, byrow = TRUE))
        # r_j = t(A_hat) %*% W %*% tilde_Gamma[,j] - sum
        # r_j_times_tilde_gamma = r_j * tilde_gamma[j, k]
        # update_k = update_k + r_j_times_tilde_gamma
      }
      B_hat[, k] = soft_threshold(B_hat[, k] + update_k/px, lambda[k])  # L1 Shrinkage
    }

    # 计算收敛条件
    max_change = max(abs(B_hat - prev_B_hat))
    iter_cd = iter_cd + 1
  }

  B_hat_list[[iter + 1]] = B_hat  # 存储优化后的 B_hat

  # **Step 2: Given B, optimize A using Orthogonal Procrustes**
  svd_result <- svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)
  U <- svd_result$u
  V <- svd_result$v
  A_hat_opt <- solve(W_sqrt) %*% V %*% t(U)

  # 存储优化后的 A_hat
  A_hat_list[[iter + 1]] = A_hat_opt

  # **计算 dist to C**
  # dist_to_C <- norm(A_hat_opt %*% B_hat - C, type = "F")  # 使用 Frobenius 范数
  # mean abs
  dist_to_C <- mean(abs(A_hat_opt %*% B_hat - C))

  cat("Iteration:", iter,
      "| Dist to C:", round(dist_to_C, 4),
      "| Zeros in B_hat:", sum(B_hat == 0), "\n")

  # 如果收敛，则提前终止
  if (dist_to_C < 0.1) {
    break
  }
}

A_hat = A_hat_opt
B_hat = B_hat_list[[iter + 1]]







#####
# # 2.22 first try: 2.382579e-05  1.470794e-07  1.786200e-08 4.774000e-108 7.578601e-275
# # 2.22 second try: 8.268869e-11  3.624960e-12  8.996573e-13 2.584655e-25 7.505441e-139
#
# # c=2: [1] 8.212083e-11 1.447780e-11 3.202669e-12 3.463476e-13 2.466910e-57
# # c=1: [1] 8.432722e-11 1.470502e-11 4.104934e-12 1.074628e-12 1.032944e-23
# # c=0.5: 8.295011e-11 1.444003e-11 4.106420e-12 1.282637e-12 1.230937e-13
# regularization_rate_list = c(1.279385e-11, 1.913728e-12, 3.941866e-13, 6.660968e-14, 3.191639e-15)
# # The simulation function for the naive MR-rr estimator and the MR-rr estimator (with spectral regularization)
# run_simulation <- function(regularized = TRUE, regularization_rate_list){
#   sample_weight_list = c(5, 2, 1, 0.5, 0.2)
#   eloop = 500
#   result_AB_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
#   result_AB_d_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
#   result_AB_d_r_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
#   result_C_ivw_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
#   iv_strength_list = c()
#   parameters_list = c()
#
#   for (weight_index in seq_along(sample_weight_list)) {
#     me_weight = sample_weight_list[weight_index]
#     parameters = .get_parameters(me_weight, px = 24, r_RR = 5)
#     C = parameters$C
#     A = parameters$A
#     B = parameters$B
#     A_d = parameters$A_d
#     B_d = parameters$B_d
#     C_r = parameters$C_r
#     py = parameters$py
#     px = parameters$px
#
#     iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
#     iv_strength_list = c(iv_strength_list, iv_strength)
#
#     parameters_list = c(parameters_list, list(c(parameters, iv_strength = iv_strength)))
#
#     bias_AB_matrix = bias_AB_d_matrix = bias_AB_d_r_matrix = bias_C_ivw_matrix = matrix(NA, nrow = py*px, ncol = eloop)
#     B_star = t(B) # in order to compute the norm
#     B_d_star = t(B_d)
#     for (i in 1:eloop) {
#       simulation_result = .simulation(parameters, regularized = regularized, regularization_rate=regularization_rate_list[weight_index])
#       A_hat = simulation_result[[1]]
#       B_hat = simulation_result[[2]]
#       B_hat_star = t(simulation_result[[2]])
#       AB_hat = simulation_result[[3]]
#
#       bias_AB_vectorized = as.vector(AB_hat - C)
#       # bias_AB_vectorized = as.vector(AB_hat - C_r)
#       bias_AB_matrix[,i] = bias_AB_vectorized
#
#       A_d_hat = simulation_result[[4]]
#       B_d_hat = simulation_result[[5]]
#       B_d_hat_star = t(simulation_result[[5]])
#       AB_d_hat = simulation_result[[6]]
#
#       # bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
#       bias_AB_d_vectorized = as.vector(AB_d_hat - C)
#       bias_AB_d_matrix[,i] = bias_AB_d_vectorized
#
#       A_d_r_hat = simulation_result[[7]]
#       B_d_r_hat = simulation_result[[8]]
#       B_d_r_hat_star = t(simulation_result[[8]])
#       AB_d_r_hat = simulation_result[[9]]
#
#       # bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
#       bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C)
#       bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
#
#       C_ivw = simulation_result[[12]]
#       bias_C_ivw_vectorized = as.vector(C_ivw - C)
#       bias_C_ivw_matrix[,i] = bias_C_ivw_vectorized
#     }
#     char <- as.character(me_weight)
#     result_AB_list[[char]] = bias_AB_matrix
#     result_AB_d_list[[char]] = bias_AB_d_matrix
#     result_AB_d_r_list[[char]] = bias_AB_d_r_matrix
#     result_C_ivw_list[[char]] = bias_C_ivw_matrix
#   }
#   return(list(result_AB_list=result_AB_list, result_AB_d_list=result_AB_d_list, result_AB_d_r_list=result_AB_d_r_list,
#               iv_strength_list=iv_strength_list, parameters_list=parameters_list,
#               result_C_ivw_list=result_C_ivw_list))
# }
#
#
# simulate_result = run_simulation(regularization_rate_list = regularization_rate_list)  # TODO: choose the regularization rate for each weight based on objective function
# save(simulate_result, file = "results/simulate_result_250306.RData")
#
# # load("results/simulate_result_250222_2.RData")
# round(simulate_result$iv_strength_list, 3)
#
#
# load("results/simulate_result_250306.RData")
#
#
# # The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
# plot_heatmap <- function(parameters_list,
#                          result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
#                          result_ivw_list,
#                          weight_to_plot)
#   {
#   px <- 24
#   py <- 10
#   abs_mean_entry_bias <- rep(NA, (px*py))
#   mean_entry_bias <- rep(NA, (px*py))
#   sd_entry_bias <- rep(NA, (px*py))
#   for (i in 1:(px*py)) {
#     abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
#     mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
#     sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
#   }
#
#   abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
#   mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
#   sd_entry_bias_d <- rep(NA, (px*py))
#   for (i in 1:(px*py)) {
#     abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
#     mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
#     sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
#   }
#
#   abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
#   mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
#   sd_entry_bias_d_r <- rep(NA, (px*py))
#   for (i in 1:(px*py)) {
#     abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
#     mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
#     sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
#   }
#
#   abs_mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for the heatmap
#   mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for reporting the avg bias
#   sd_entry_bias_ivw <- rep(NA, (px*py))
#   for (i in 1:(px*py)) {
#     abs_mean_entry_bias_ivw[i] <- abs(as.numeric(mean(result_ivw_list[[weight_to_plot]][i,])))
#     mean_entry_bias_ivw[i] <- mean(result_ivw_list[[weight_to_plot]][i,])
#     sd_entry_bias_ivw[i] <- sd(result_ivw_list[[weight_to_plot]][i,])
#   }
#
#   # bias heatmap
#   # RR
#   avg_bias <- matrix(abs_mean_entry_bias, nrow = px, ncol = py)
#   avg_bias <- t(avg_bias)
#   avg_bias <- as.data.frame(avg_bias)
#   colnames(avg_bias) <- 1:px
#   rownames(avg_bias) <- 1:py
#
#   # dRR
#   avg_bias_d <- matrix(abs_mean_entry_bias_d, nrow = px, ncol = py)
#   avg_bias_d <- t(avg_bias_d)
#   avg_bias_d <- as.data.frame(avg_bias_d)
#   colnames(avg_bias_d) <- 1:px
#   rownames(avg_bias_d) <- 1:py
#
#   # dRR_r
#   avg_bias_d_r <- matrix(abs_mean_entry_bias_d_r, nrow = px, ncol = py)
#   avg_bias_d_r <- t(avg_bias_d_r)
#   avg_bias_d_r <- as.data.frame(avg_bias_d_r)
#   colnames(avg_bias_d_r) <- 1:px
#   rownames(avg_bias_d_r) <- 1:py
#
#   # IVW
#   avg_bias_ivw <- matrix(abs_mean_entry_bias_ivw, nrow = px, ncol = py)
#   avg_bias_ivw <- t(avg_bias_ivw)
#   avg_bias_ivw <- as.data.frame(avg_bias_ivw)
#   colnames(avg_bias_ivw) <- 1:px
#   rownames(avg_bias_ivw) <- 1:py
#
#
#   # color breaks
#   min_value <- min(min(avg_bias), min(avg_bias_d), min(avg_bias_d_r))
#   max_value <- max(max(avg_bias),max(avg_bias_d), max(avg_bias_d_r))
#
#   # min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
#   # max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))
#
#   breaks <- seq(min_value, max_value, length.out = 101)
#
#   hm_bias_1 = pheatmap::pheatmap(avg_bias, breaks = breaks,
#                      main = "Naive MR-rr estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_bias_2 = pheatmap::pheatmap(avg_bias_d, breaks = breaks,
#                      main = "MR-rr estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_bias_3 = pheatmap::pheatmap(avg_bias_d_r, breaks = breaks,
#                      main = "MR-rr estimator with regularization", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_bias_4 = pheatmap::pheatmap(avg_bias_ivw, breaks = breaks,
#                      main = "IVW estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   combined_hm_bias = gridExtra::grid.arrange(
#     grobs = list(hm_bias_4, hm_bias_1, hm_bias_2, hm_bias_3),
#     ncol = 4,
#     top = grid::textGrob(sprintf("Absolute average bias by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
#   )
#   combined_hm_bias
#   # sd heatmap
#   # RR
#   avg_sd <- sd_entry_bias
#   avg_sd <- matrix(avg_sd, nrow = px, ncol = py)
#   avg_sd <- t(avg_sd)
#   avg_sd <- as.data.frame(avg_sd)
#   colnames(avg_sd) <- 1:px
#   rownames(avg_sd) <- 1:py
#
#   # dRR
#   avg_sd_d <- sd_entry_bias_d
#   avg_sd_d <- matrix(avg_sd_d, nrow = px, ncol = py)
#   avg_sd_d <- t(avg_sd_d)
#   avg_sd_d <- as.data.frame(avg_sd_d)
#   colnames(avg_sd_d) <- 1:px
#   rownames(avg_sd_d) <- 1:py
#
#   # dRR_r
#   avg_sd_d_r <- sd_entry_bias_d_r
#   avg_sd_d_r <- matrix(avg_sd_d_r, nrow = px, ncol = py)
#   avg_sd_d_r <- t(avg_sd_d_r)
#   avg_sd_d_r <- as.data.frame(avg_sd_d_r)
#   colnames(avg_sd_d_r) <- 1:px
#   rownames(avg_sd_d_r) <- 1:py
#
#   # IVW
#   avg_sd_ivw <- sd_entry_bias_ivw
#   avg_sd_ivw <- matrix(avg_sd_ivw, nrow = px, ncol = py)
#   avg_sd_ivw <- t(avg_sd_ivw)
#   avg_sd_ivw <- as.data.frame(avg_sd_ivw)
#   colnames(avg_sd_ivw) <- 1:px
#   rownames(avg_sd_ivw) <- 1:py
#
#
#   min_value_sd <- min(min(avg_sd), min(avg_sd_d), min(avg_sd_d_r))
#   max_value_Sd <- max(max(avg_sd),max(avg_sd_d), max(avg_sd_d_r))
#   breaks_sd <- seq(min_value_sd, max_value_Sd, length.out = 101)
#
#
#   hm_sd_1 = pheatmap::pheatmap(avg_sd, breaks = breaks_sd,
#                      main = "Naive MR-rr estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_sd_2 = pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
#                      main = "MR-rr estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_sd_3 = pheatmap::pheatmap(avg_sd_d_r, breaks = breaks_sd,
#                      main = "MR-rr estimator with regularization", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   hm_sd_4 = pheatmap::pheatmap(avg_sd_ivw, breaks = breaks_sd,
#                      main = "IVW estimator", fontsize = 8,
#                      cluster_rows = FALSE, cluster_cols = FALSE,
#                      color = colorRampPalette(c("white", "red"))(100),
#                      height = 10,
#                      width = 8, silent = TRUE)$gtable
#
#   combined_hm_sd = gridExtra::grid.arrange(
#     grobs = list(hm_sd_4,hm_sd_1, hm_sd_2, hm_sd_3),
#     ncol = 4,
#     top = grid::textGrob(sprintf("SD by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
#   )
#   combined_hm_sd
#
#   # bias scale compares to true C_r and sd
#   weight_index = match(weight_to_plot, c("5", "2", "1", "0.5", "0.2"))
#   parameters = parameters_list[[weight_index]]
#   C_r = parameters$C_r
#   # mean
#   mean_abs_C_entry =  round(mean(abs(as.matrix(C_r))), 3)
#   mean_abs_bias_naive = round(mean(abs(as.matrix(avg_bias))), 3)
#   mean_abs_bias = round(mean(abs(as.matrix(avg_bias_d))), 3)
#   mean_abs_bias_r = round(mean(abs(as.matrix(avg_bias_d_r))), 3)
#   mean_abs_bias_ivw = round(mean(abs(as.matrix(avg_bias_ivw))), 3)
#
#   mean_sd_naive = round(mean(as.matrix(avg_sd)), 3)
#   mean_sd = round(mean(as.matrix(avg_sd_d)), 3)
#   mean_sd_r = round(mean(as.matrix(avg_sd_d_r)), 3)
#   mean_sd_ivw = round(mean(as.matrix(avg_sd_ivw)), 3)
#
#   # median
#   median_abs_bias_naive = round(median(abs(as.matrix(avg_bias))), 3)
#   median_abs_bias = round(median(abs(as.matrix(avg_bias_d))), 3)
#   median_abs_bias_r = round(median(abs(as.matrix(avg_bias_d_r))), 3)
#   median_abs_bias_ivw = round(median(abs(as.matrix(avg_bias_ivw))), 3)
#
#   median_sd_naive = round(median(as.matrix(avg_sd)), 3)
#   median_sd = round(median(as.matrix(avg_sd_d)), 3)
#   median_sd_r = round(median(as.matrix(avg_sd_d_r)), 3)
#   median_sd_ivw = round(median(as.matrix(avg_sd_ivw)), 3)
#
#   # quantile
#   quantile_abs_bias_naive = round(quantile(abs(as.matrix(avg_bias)), c(0.25, 0.75)), 3)
#   quantile_abs_bias = round(quantile(abs(as.matrix(avg_bias_d)), c(0.25, 0.75)), 3)
#   quantile_abs_bias_r = round(quantile(abs(as.matrix(avg_bias_d_r)), c(0.25, 0.75)), 3)
#   quantile_abs_bias_ivw = round(quantile(abs(as.matrix(avg_bias_ivw)), c(0.25, 0.75)), 3)
#
#   quantile_sd_naive = round(quantile(as.matrix(avg_sd), c(0.25, 0.75)), 3)
#   quantile_sd = round(quantile(as.matrix(avg_sd_d), c(0.25, 0.75)), 3)
#   quantile_sd_r = round(quantile(as.matrix(avg_sd_d_r), c(0.25, 0.75)), 3)
#   quantile_sd_ivw = round(quantile(as.matrix(avg_sd_ivw), c(0.25, 0.75)), 3)
#   #####
#
#   ggplot2::ggsave(sprintf("results/hmap_3.7/bias_hm_weight_%s.png",weight_to_plot), plot = combined_hm_bias, width = 12, height = 4)
#   ggplot2::ggsave(sprintf("results/hmap_3.7/sd_hm_weight_%s.png",weight_to_plot), plot = combined_hm_sd, width = 12, height = 4)
#
#   return(list(mean_abs_C_entry = mean_abs_C_entry,
#               # report mean
#               mean_abs_bias_ivw = mean_abs_bias_ivw, mean_sd_ivw = mean_sd_ivw,
#               mean_abs_bias_naive = mean_abs_bias_naive, mean_sd_naive = mean_sd_naive,
#               mean_abs_bias = mean_abs_bias, mean_sd = mean_sd,
#               mean_abs_bias_r = mean_abs_bias_r, mean_sd_r = mean_sd_r,
#               # report median
#               median_abs_bias_ivw = median_abs_bias_ivw, median_sd_ivw = median_sd_ivw,
#               median_abs_bias_naive = median_abs_bias_naive, median_sd_naive = median_sd_naive,
#               median_abs_bias = median_abs_bias, median_sd = median_sd,
#               median_abs_bias_r = median_abs_bias_r, median_sd_r = median_sd_r,
#               # report quantile
#               quantile_abs_bias_ivw = quantile_abs_bias_ivw, quantile_sd_ivw = quantile_sd_ivw,
#               quantile_abs_bias_naive = quantile_abs_bias_naive, quantile_sd_naive = quantile_sd_naive,
#               quantile_abs_bias = quantile_abs_bias, quantile_sd = quantile_sd,
#               quantile_abs_bias_r = quantile_abs_bias_r, quantile_sd_r = quantile_sd_r
#               ))
# }
#
#
# plot_heatmap(parameters_list = simulate_result$parameters_list,
#              result_naive_MRrr_list = simulate_result$result_AB_list,
#              result_MRrr_list = simulate_result$result_AB_d_list,
#              result_r_MRrr_list = simulate_result$result_AB_d_r_list,
#              result_ivw_list = simulate_result$result_C_ivw_list,
#              weight_to_plot="5")
#
#
#
# # choose regularization rate \lambda
# # ----
# c=1
# me_weight_list = c(5,2,1,0.5,0.2)
# eloop = 100
# # eloop * length(me_weight_list) = 50
# result_matrix = matrix(NA, nrow = eloop, ncol = length(me_weight_list))
# for (i in 1:eloop){
#   result_list = list()
#   for (me_weight in me_weight_list){
#     parameters = .get_parameters(me_weight=me_weight, r_RR = 5, px=24,  r_approx = 5)
#     n = 2000
#     py = parameters$py
#     px = parameters$px
#     var_Z = parameters$var_Z
#     VX_tilde = parameters$VX_tilde
#     Sigma_X = parameters$Sigma_X
#     Sigma_Y = parameters$Sigma_Y
#     SigmaXX = parameters$SigmaXX
#     SigmaYX = parameters$SigmaYX
#     SigmaXY = parameters$SigmaXY
#     C = parameters$C
#     r_RR = parameters$r_RR
#     VY_tilde = parameters$VY_tilde
#     W = parameters$weight.matrix
#     W_sqrt = .sqrt_matrix(parameters$weight.matrix)
#     # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
#     gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
#     Gamma_j_star = gamma_j_star %*% t(C)
#
#     # sample x_j_hat, y_j_hat
#     x_j_hat = matrix(0, n, px)
#     y_j_hat = matrix(0, n, py)
#     for (j in 1:n) {
#       x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#       y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
#     }
#
#     iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
#     q_list = seq(1, 70, 1)
#     regu_rate_list = c()
#     obj_value = c()
#     for (q in q_list) {
#       sigma_y2 = mean(eigen(parameters$Sigma_X)$values)
#
#       regu_rate = 1/10 * sigma_y2^2 * exp(c * (q/sqrt(n)-(iv_strength+1)))/n
#
#       # regu_rate = c * sigma_y2^2 * (q-(iv_strength+1)/n)
#       # regu_rate = exp(c*(q-iv_strength*sqrt(n)))
#       regu_rate_list = c(regu_rate_list, regu_rate)
#       result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate)
#       C_hat = result_d$AB
#       # objective function value
#       obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
#       obj_value = c(obj_value, obj)
#     }
#
#     # sigma_y2 = min(eigen(parameters$Sigma_X)$values)
#     # bound = sigma_y2^2 * (iv_strength+1)/n
#
#
#     # if (min(regu_rate_list) > bound){
#     #   opt_rate = bound
#     # } else {
#     #   # 筛掉大于bound的
#     #   obj_value = obj_value[regu_rate_list < bound]
#     #   # q_list = q_list[regu_rate_list < bound]
#     #   regu_rate_list = regu_rate_list[regu_rate_list < bound]
#     #
#     #   #选出obj_value最小的对应的rate
#     #   opt_rate = regu_rate_list[order(obj_value)][1]
#     # }
#
#     opt_rate = regu_rate_list[order(obj_value)][1]
#
#     result_list= c(result_list, opt_rate)
#
#     # test.rate_list = c(1e-5, 1e-6, 1e-7, 1e-8, 1e-9,
#     #                    1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16, 1e-17, 1e-18)
#     # obj_value = c()
#     # for (test.rate in test.rate_list) {
#     #   result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=test.rate)
#     #   C_hat = result_d$AB
#     #   # objective function value
#     #   obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
#     #   obj_value = c(obj_value, obj)
#     # }
#     # opt_rate = test.rate_list[which.min(obj_value)]
#     # # compare opt_rate with the bound set by theory
#     # # \mu_{n,min} = iv_strength
#     # iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
#     # # rate of sigma_y2 ~ eigenvalue of sigma_X
#     # sigma_y2 = mean(eigen(parameters$Sigma_X)$values)
#     # # bound = sigma_y2^2 * (iv_strength+1)/n
#     #
#     # # 2025.2.22 new bound
#     # c = 1
#     # bound = exp(c*(30-iv_strength*sqrt(n)))
#     #
#     # opt_rate = min(opt_rate, bound)
#     # # return the rate that minimize objective function value
#     # result_list= c(result_list, opt_rate)
#   }
#   result_matrix[i,] = unlist(result_list)
# }
# # average
# result_matrix_mean = apply(result_matrix, 2, mean)
# result_matrix_mean
#
# #----
#
# #3.4: [1] 1.197628e-12 1.758998e-13 3.578033e-14 5.484086e-15 3.258540e-16: 第一个cover 有点大，rate太小
#
# # 2.346760e-06 2.717263e-07 2.689552e-08 2.984442e-09 7.296111e-12
# # after setting bound
# # 1.667258e-10 2.881427e-11 8.177887e-12 2.411417e-12 3.597363e-13
#
# # min来算mu 可能太小 2.438109e-16 4.258233e-17 1.201675e-17 3.650834e-18 9.475329e-19
#
# # 8.386425e-11 1.447977e-11 4.013995e-12 1.188179e-12 2.115718e-13
# # note on abs error: we did not take absolute within each entry. only abs after we have average bias for each entry. it is reasonable
#
#
# # c=5: c(8.268869e-11, 3.624960e-12, 8.996573e-13, 2.584655e-25, 7.505441e-139)
# # c=2: [1] 8.212083e-11 1.447780e-11 3.202669e-12 3.463476e-13 2.466910e-57
# # c=1: [1] 8.432722e-11 1.470502e-11 4.104934e-12 1.074628e-12 1.032944e-23
# # c=0.5: [1] 8.295011e-11 1.444003e-11 4.106420e-12 1.282637e-12 1.230937e-13
#
# # new algorithm: [1] 1.350594e-11 2.204258e-12 5.255102e-13 1.112693e-13 1.279553e-14
#
#
# # c=10 [1] 1.031047e-08 6.757877e-10 3.515137e-11 3.530294e-13 5.290825e-18
# # c=2 [1] 2.123233e-10 2.822416e-11 4.900563e-12 6.246488e-13 9.974627e-15
#
# # c=2, 1/10: [1] 2.148572e-11 2.838746e-12 5.020544e-13 5.689990e-14 1.192533e-15
# # c=1. 1/10: [1] 1.279385e-11 1.913728e-12 3.941866e-13 6.660968e-14 3.191639e-15
#
# #----
#
# # archive code
# #####
#
# # # The function to create the boxplot for the bias of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization) comparing to true C for each entry.
# # plot_boxplot_old <- function(result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
# #                          weight_to_plot, individual_plot = FALSE, rank_by = "bias")
# #   {
# #   px <- 24
# #   py <- 10
# #   # individual boxplot the bias of AB_hat
# #   myplots <- apply(result_naive_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias)
# #   for (i in 1:(px*py)){
# #     names(myplots)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
# #   }
# #
# #   abs_mean_entry_bias <- rep(NA, (px*py))
# #   mean_entry_bias <- rep(NA, (px*py))
# #   sd_entry_bias <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
# #     mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   index1 <- order(abs_mean_entry_bias, decreasing = T)
# #   index2 <- order(sd_entry_bias, decreasing = T)
# #
# #   # individual boxplot the bias of AB_d_hat
# #   myplots_d <- apply(result_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d)
# #   for (i in 1:(px*py)){
# #     names(myplots_d)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
# #   }
# #
# #   abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
# #   mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
# #   sd_entry_bias_d <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
# #     mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
# #   index2_d <- order(sd_entry_bias_d, decreasing = T)
# #
# #   if (individual_plot){
# #     if (rank_by == "bias") {
# #       # individual box plot
# #       return(list(myplots[index1][1:2], myplots_d[index1_d][1:2]))
# #     } else if (rank_by == "sd"){
# #       # individual box plot
# #       return(list(myplots[index2][1:2], myplots_d[index2_d][1:2]))
# #     } else {
# #       stop("rank_by should be either 'bias' or 'sd'")
# #     }
# #   } else {
# #     if (rank_by == "bias"){
# #       # top 5 box plot in the same plot
# #       combined_plot <- myplots_d[[index1_d[5]]] + myplots_d[[index1_d[4]]] +
# #         myplots_d[[index1_d[3]]] + myplots_d[[index1_d[2]]] + myplots_d[[index1_d[1]]] +
# #         myplots[[index1[5]]] + myplots[[index1[4]]] + myplots[[index1[3]]] +
# #         myplots[[index1[2]]] + myplots[[index1[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
# #         patchwork::plot_annotation(
# #           title = "Bias boxplots of the entries with top 5 average bias \n for the MR-rr and naive MR-rr estimator",
# #           subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
# #           # caption = "C_r defined as minimizing modified objective function of
# #           # 2.2 or equivalent objective function of 2.1"
# #           )
# #       return(combined_plot)
# #     } else if (rank_by == "sd"){
# #       # rank by var
# #       combined_plot_var <- myplots_d[[index2_d[5]]] + myplots_d[[index2_d[4]]] +
# #         myplots_d[[index2_d[3]]] + myplots_d[[index2_d[2]]] + myplots_d[[index2_d[1]]] +
# #         myplots[[index2[5]]] + myplots[[index2[4]]] + myplots[[index2[3]]] +
# #         myplots[[index2[2]]] + myplots[[index2[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
# #         patchwork::plot_annotation(
# #           title = "Bias boxplots of the entries with top 5 standard deviation \n for the MR-rr and naive MR-rr estimator",
# #           subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
# #           # caption = "C_r defined as minimizing modified objective function of
# #           # 2.2 or equivalent objective function of 2.1"
# #           )
# #       return(combined_plot_var)
# #     } else {
# #       stop("rank_by should be either 'bias' or 'sd'")}
# #   }
# # }
# #
# #
# # # The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
# # plot_heatmap_old <- function(parameters_list, result_naive_MRrr_list, result_MRrr_list, weight_to_plot,regu) {
# #   px <- 24
# #   py <- 10
# #   abs_mean_entry_bias <- rep(NA, (px*py))
# #   mean_entry_bias <- rep(NA, (px*py))
# #   sd_entry_bias <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
# #     mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
# #   mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
# #   sd_entry_bias_d <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
# #     mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   # bias heatmap
# #   # RR
# #   avg_bias <- matrix(abs_mean_entry_bias, nrow = px, ncol = py)
# #   avg_bias <- t(avg_bias)
# #   avg_bias <- as.data.frame(avg_bias)
# #   colnames(avg_bias) <- 1:px
# #   rownames(avg_bias) <- 1:py
# #
# #   # dRR
# #   avg_bias_d <- matrix(abs_mean_entry_bias_d, nrow = px, ncol = py)
# #   avg_bias_d <- t(avg_bias_d)
# #   avg_bias_d <- as.data.frame(avg_bias_d)
# #   colnames(avg_bias_d) <- 1:px
# #   rownames(avg_bias_d) <- 1:py
# #
# #   # color breaks
# #   min_value <- min(min(avg_bias), min(avg_bias_d)) #TODO:
# #   max_value <- max(max(avg_bias),max(avg_bias_d))
# #
# #   # min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
# #   # max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))
# #
# #   breaks <- seq(min_value, max_value, length.out = 101)
# #
# #   pheatmap::pheatmap(avg_bias, breaks = breaks,
# #                      main = sprintf("Absolute bias of the naive MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                      cluster_rows = FALSE, cluster_cols = FALSE,
# #                      color = colorRampPalette(c("white", "red"))(100),
# #                      height = 10,
# #                      width = 8)
# #   if (regu==TRUE){
# #     pheatmap::pheatmap(avg_bias_d, breaks = breaks,
# #                        main = sprintf("Absolute bias of the MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                        cluster_rows = FALSE, cluster_cols = FALSE,
# #                        color = colorRampPalette(c("white", "red"))(100),
# #                        height = 10,
# #                        width = 8)
# #   } else {
# #     pheatmap::pheatmap(avg_bias_d, breaks = breaks,
# #                        main = sprintf("Absolute bias of the MR-rr estimator with regularization by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                        cluster_rows = FALSE, cluster_cols = FALSE,
# #                        color = colorRampPalette(c("white", "red"))(100),
# #                        height = 10,
# #                        width = 8)
# #   }
# #   # sd heatmap
# #   # RR
# #   avg_sd <- sd_entry_bias
# #   avg_sd <- matrix(avg_sd, nrow = px, ncol = py)
# #   avg_sd <- t(avg_sd)
# #   avg_sd <- as.data.frame(avg_sd)
# #   colnames(avg_sd) <- 1:px
# #   rownames(avg_sd) <- 1:py
# #
# #   # dRR
# #   avg_sd_d <- sd_entry_bias_d
# #   avg_sd_d <- matrix(avg_sd_d, nrow = px, ncol = py)
# #   avg_sd_d <- t(avg_sd_d)
# #   avg_sd_d <- as.data.frame(avg_sd_d)
# #   colnames(avg_sd_d) <- 1:px
# #   rownames(avg_sd_d) <- 1:py
# #
# #
# #   # color breaks
# #   # min_value_sd <- min(min(avg_sd), min(avg_sd_d))
# #   # breaks_sd <- c(seq(min_value_sd, max(avg_sd), length.out = 51),
# #   #                seq(max(avg_sd), max(avg_sd_d), length.out = 51)[-1])
# #
# #
# #   min_value_sd <- min(min(avg_sd), min(avg_sd_d))
# #   max_value_Sd <- max(max(avg_sd),max(avg_sd_d))
# #   breaks_sd <- seq(min_value_sd, max_value_Sd, length.out = 101)
# #
# #
# #   pheatmap::pheatmap(avg_sd, breaks = breaks_sd,
# #                      main = sprintf("SD of the naive MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                      cluster_rows = FALSE, cluster_cols = FALSE,
# #                      color = colorRampPalette(c("white", "red"))(100),
# #                      height = 10,
# #                      width = 8)
# #
# #   if (regu==TRUE){
# #     pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
# #                        main = sprintf("SD of the MR-rr estimator with regularization by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                        cluster_rows = FALSE, cluster_cols = FALSE,
# #                        color = colorRampPalette(c("white", "red"))(100),
# #                        height = 10,
# #                        width = 8)
# #   } else {
# #     pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
# #                        main = sprintf("SD of the MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
# #                        cluster_rows = FALSE, cluster_cols = FALSE,
# #                        color = colorRampPalette(c("white", "red"))(100),
# #                        height = 10,
# #                        width = 8)
# #   }
# #
# #   # # compute average iv strength
# #   # iv_list = list()
# #   # for (i in 1:10){
# #   #   parameters = .get_parameters(as.numeric(weight_to_plot), px = 24, r_RR = 5) # weight: 0.2, 0.5, 1 are comparable with Yinxiang's Paper
# #   #   iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
# #   #   iv_list[[i]] = iv_strength
# #   # }
# #   # iv_strength = mean(as.numeric(unlist(iv_list)))
# #
# #   # bias scale compares to true C_r and sd
# #   weight_index = match(weight_to_plot, c("5", "2", "1", "0.5", "0.2"))
# #   parameters = parameters_list[[weight_index]]
# #   C_r = parameters$C_r
# #   mean_abs_C_entry =  mean(abs(as.matrix(C_r)))
# #   mean_abs_bias_naive = mean(abs(as.matrix(avg_bias)))
# #   mean_abs_bias = mean(abs(as.matrix(avg_bias_d)))
# #   mean_sd_naive = mean(as.matrix(avg_sd))
# #   mean_sd = mean(as.matrix(avg_sd_d))
# #
# #
# #   return(list(mean_abs_C_entry = mean_abs_C_entry,
# #               mean_abs_bias_naive = mean_abs_bias_naive, mean_abs_bias = mean_abs_bias,
# #               mean_sd_naive = mean_sd_naive, mean_sd = mean_sd))
# # }
#
# # plot_boxplot <- function(result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
# #                          weight_to_plot, individual_plot = FALSE, rank_by = "bias")
# # {
# #   px <- 24
# #   py <- 10
# #   # individual boxplot the bias of AB_hat
# #   myplots <- apply(result_naive_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias)
# #   for (i in 1:(px*py)){
# #     names(myplots)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
# #   }
# #
# #   abs_mean_entry_bias <- rep(NA, (px*py))
# #   mean_entry_bias <- rep(NA, (px*py))
# #   sd_entry_bias <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
# #     mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   index1 <- order(abs_mean_entry_bias, decreasing = T)
# #   index2 <- order(sd_entry_bias, decreasing = T)
# #
# #   # individual boxplot the bias of AB_d_hat
# #   myplots_d <- apply(result_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d)
# #   for (i in 1:(px*py)){
# #     names(myplots_d)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
# #   }
# #
# #   abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
# #   mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
# #   sd_entry_bias_d <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
# #     mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
# #   index2_d <- order(sd_entry_bias_d, decreasing = T)
# #
# #   # individual boxplot the bias of AB_d_r_hat
# #   myplots_d_r <- apply(result_r_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d_r)
# #   for (i in 1:(px*py)){
# #     names(myplots_d_r)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
# #   }
# #
# #   abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
# #   mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
# #   sd_entry_bias_d_r <- rep(NA, (px*py))
# #   for (i in 1:(px*py)) {
# #     abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
# #     mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
# #     sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
# #   }
# #
# #   index1_d_r <- order(abs_mean_entry_bias_d_r, decreasing = T)
# #   index2_d_r <- order(sd_entry_bias_d_r, decreasing = T)
# #
# #
# #   if (individual_plot){
# #     if (rank_by == "bias") {
# #       # individual box plot
# #       return(list(myplots[index1][1:2], myplots_d[index1_d][1:2]))
# #     } else if (rank_by == "sd"){
# #       # individual box plot
# #       return(list(myplots[index2][1:2], myplots_d[index2_d][1:2]))
# #     } else {
# #       stop("rank_by should be either 'bias' or 'sd'")
# #     }
# #   } else {
# #     if (rank_by == "bias"){
# #       # top 5 box plot in the same plot
# #
# #       # 定义每一行图形
# #       row1 <- myplots[[index1[3]]] + myplots[[index1[2]]] + myplots[[index1[1]]]
# #       row2 <- myplots_d[[index1_d[3]]] + myplots_d[[index1_d[2]]] + myplots_d[[index1_d[1]]]
# #       row3 <- myplots_d_r[[index1_d_r[3]]] + myplots_d_r[[index1_d_r[2]]] + myplots_d_r[[index1_d_r[1]]]
# #
# #       # 定义左侧标题
# #       title1 <- patchwork::wrap_elements(grid::textGrob("Naive MR-rr", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))
# #       title2 <- patchwork::wrap_elements(grid::textGrob("MR-rr", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))
# #       title3 <- patchwork::wrap_elements(grid::textGrob("MR-rr with regularization", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))
# #
# #       # 合并标题和对应的图形
# #       row1_labeled <- title1 + row1 + patchwork::plot_layout(widths = c(1, 5))
# #       row2_labeled <- title2 + row2 + patchwork::plot_layout(widths = c(1, 5))
# #       row3_labeled <- title3 + row3 + patchwork::plot_layout(widths = c(1, 5))
# #
# #       # 合并所有行，并添加全局标题
# #       combined_plot <- row1_labeled / row2_labeled / row3_labeled +
# #         patchwork::plot_annotation(
# #           title = "Bias boxplots of the entries with top 3 average bias",
# #           subtitle = sprintf("Weight = %s", weight_to_plot)
# #         )
# #
# #       combined_plot
# #
# #       ggplot2::ggsave(filename = sprintf("bias_boxplot_%s.png", weight_to_plot), plot = combined_plot, width = 10, height = 10)
# #     } else if (rank_by == "sd"){
# #       # rank by var
# #       combined_plot_var <- myplots_d[[index2_d[5]]] + myplots_d[[index2_d[4]]] +
# #         myplots_d[[index2_d[3]]] + myplots_d[[index2_d[2]]] + myplots_d[[index2_d[1]]] +
# #         myplots[[index2[5]]] + myplots[[index2[4]]] + myplots[[index2[3]]] +
# #         myplots[[index2[2]]] + myplots[[index2[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
# #         patchwork::plot_annotation(
# #           title = "Bias boxplots of the entries with top 3 standard deviation \n for the MR-rr and naive MR-rr estimator",
# #           subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
# #           # caption = "C_r defined as minimizing modified objective function of
# #           # 2.2 or equivalent objective function of 2.1"
# #         )
# #       return(combined_plot_var)
# #     } else {
# #       stop("rank_by should be either 'bias' or 'sd'")}
# #   }
# # }
# #
# #
# # plot_boxplot(result_naive_MRrr_list = simulate_result$result_AB_list,
# #              result_MRrr_list = simulate_result$result_AB_d_list,
# #              result_r_MRrr_list = simulate_result$result_AB_d_r_list,
# #              weight_to_plot="2")
# #####
