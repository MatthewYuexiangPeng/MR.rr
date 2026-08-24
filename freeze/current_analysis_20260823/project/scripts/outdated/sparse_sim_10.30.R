rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
library(matrixStats)
library(MASS)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(pheatmap)
library(readxl)
library(glmnet)

#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/Yuexiang Peng/UW/Research/Ye Ting/sim_ArBr_bias")
data("bmi.cad")
load('data/multivariate_data_medium.rda')

# load csv data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
# read xlsx data: 'data/Kennetu_2016_download_links_updated.xlsx'
samplesize_data = read_excel('data/Kennetu_2016_download_links_updated.xlsx')

source("scripts/RRR.R")
source("scripts/debiased_RRR.R")
source("scripts/debiased_RRR_modified.R")

# # regress lip_data$gamma_out on lip_data$gamma_exp1 to lip_data$gamma_exp24
# mean(abs(lm(lip_data$gamma_out ~ ., data = lip_data[,paste0('gamma_exp',1:24)])$coefficients))
# # the scale is similar to the generated C in simulation


set.seed(123)
regularization_rate = 1e-13 
# (-12 will cause to much bias no matter in with weight)
# (-13 may be the proper one, present small bias and sd for small weight)

# Functions ----------------------------------------------------------------
sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
}

# 24.11.7: need B to be sparse for testing sparse algorithm
get_parameters = function(sample_weight, r_RR = 5, px,  r_approx = 5){
  # ## test
  # sample_weight =1
  # r_RR = 5
  # px=24
  # r_approx = 5

  # r_RR is the rank chosed bu user when performing RRR
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
  # VX_tilde_upscale = sqrt_matrix(diag(diag_VX_tilde_sample)) %*% VX_tilde_corr %*% 
  #                    sqrt_matrix(diag(diag_VX_tilde_sample))
  # VX_tilde_upscale <- VX_tilde_upscale*1/2 # Q: why?
  
  
  # Sigma_Xj & Sigma_X ----------------------------------------------------------------
  sigma_gamma_j = as.matrix(lip_data[,paste0('se_exp',1:24)])
  Sigma = lip_corr[1:24,1:24]
  # change "Sigma" to matrix
  Sigma = as.matrix(Sigma)
  sqrt_Sigma = sqrt_matrix(Sigma)
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
  # Sigma_X_upscale <- sample_weight * Sigma_X_upscale
  
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
  # Sigma_X_upscale = sqrt_matrix(diag(diag_Sigma_X_sample)) %*% Sigma_X_corr %*% 
  #   sqrt_matrix(diag(diag_Sigma_X_sample))
  # Sigma_X_upscale <- sample_weight * Sigma_X_upscale
  
  # lets not upscale px for now.
  Sigma_X = sample_weight * Sigma_X
  
  # Sigma_Y & weight.matrix -----------------------------------------------------------
  # Question: Use previous n_y?
  n_Y = median(samplesize_data$samplesize) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
  # var_Y = sample(lip_data$se_out1[sample(1:114, 1119, replace = TRUE)]^2 * bmi.cad$N.outcome * var_Z_raw[sample(1:length(var_Z_raw), 1119, replace = TRUE)], py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
  var_Y = sample(lip_data$se_out1^2 * samplesize_data$samplesize[sample(1:nrow(samplesize_data), 114, replace = TRUE)]* var_Z_raw[sample(1:length(var_Z_raw), 114, replace = TRUE)],py)
  
  # var_Y = sample(bmi.cad$se.outcome^2 * bmi.cad$N.outcome * var_Z_raw, py) # Question: or use the se_out for lip_data?
  Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))
  
  # Sigma_Y = sample_weight * Sigma_Y # TODO: should not weight Y?
  weight.matrix = solve(Sigma_Y)
  
  
  # SigmaXX, SigmaYX, SigmaXY, C ----------------------------------------------------------------
  # py * r matrix
  # TODO: how to approximate low rank?
  r =  min(px,py) # true rank of C, not specified reduced rank r
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
  
  # # 24.11.7: generate C from A times B ----
  # M <- matrix(rnorm(py * r_approx), nrow = py)
  # svd_result <- svd(M)
  # U <- svd_result$u
  # A <- solve(sqrt_matrix(weight.matrix)) %*% U[, 1:r_approx]
  # # t(A) %*% weight.matrix %*% A
  # 
  # # generate a sparse B
  # B <- matrix(0, nrow = r_approx, ncol = px)  # 初始化全零矩阵
  # 
  # # 确保每行和每列至少一个非零元素
  # for (i in 1:r_approx) {
  #   B[i, sample(1:px, round(px*0.2))] <- sample(-10:10, 1)
  # }
  # for (j in 1:px) {
  #   B[sample(1:r_approx, 1), j] <- sample(-10:10, 1)
  # }
  # 
  # C <- A %*% B
  # true.A_sparse = A
  # true.B_sparse = B
  # ####
  

  # SigmaXX, SigmaYX and SigmaXY
  SigmaXX <- Sigma_X + VX_tilde
  SigmaYX <- C %*% VX_tilde
  # SigmaXX <- Sigma_X_upscale + VX_tilde_upscale
  # SigmaYX <- C %*% VX_tilde_upscale
  SigmaXY <- t(SigmaYX)
  sqrt_Gamma <- sqrt_matrix(weight.matrix)
  sqrt_Gamma_inv <- solve(sqrt_Gamma)
  
  # original RRR population level estimator
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


simulation = function(parameters, n, penalized = FALSE) {
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
  W_sqrt = sqrt_matrix(parameters$weight.matrix)
  
  # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
  gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  # sample x_j_hat, y_j_hat
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    y_j_hat[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }
  
  # compute A_hat, B_hat
  result <- RRR(y_j_hat, x_j_hat, r=r_RR, W_sqrt) # TODO: changed here need check the result
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB # sample level estimator A_hat * B_hat
  
  # need to subtract the Sigma_xx by Sigma_x. Sigma_x can be estimated by the regression of exp ~ z (have been approximate in the get parameter function)
  if (penalized == FALSE) {
    result_d <- debiased_RRR(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X)
  } else {
    result_d <- debiased_RRR_modified(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X, ph=regularization_rate)
  }
  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB
  
  return(list(A_hat, B_hat, AB_hat, A_d_hat, B_d_hat, AB_d_hat, x_j_hat, y_j_hat))
}
#####


# # choose regularization rate \lambda ----
# parameters = get_parameters(sample_weight=0.2, r_RR = 5, px,  r_approx = 5)
# n = 1000
# py = parameters$py
# px = parameters$px
# var_Z = parameters$var_Z
# VX_tilde = parameters$VX_tilde
# Sigma_X = parameters$Sigma_X
# Sigma_Y = parameters$Sigma_Y
# SigmaXX = parameters$SigmaXX
# SigmaYX = parameters$SigmaYX
# SigmaXY = parameters$SigmaXY
# C = parameters$C
# r_RR = parameters$r_RR
# VY_tilde = parameters$VY_tilde
# W_sqrt = sqrt_matrix(parameters$weight.matrix)
# # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
# gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
# Gamma_j_star = gamma_j_star %*% t(C)
# 
# # sample x_j_hat, y_j_hat
# x_j_hat = matrix(0, n, px)
# y_j_hat = matrix(0, n, py)
# for (j in 1:n) {
#   x_j_hat[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#   y_j_hat[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
# }
# 
# test.rate = 1e-11
# result_d <- debiased_RRR_modified(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X, ph=test.rate)
# C_hat = result_d$AB
# # objective function value
# obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
# print(obj)
# ####


penalized = TRUE
# run the simulation -----------------------------------------------------------
sample_weight_list = c(1, 0.5, 0.2, 0.1, 0.05)
eloop = 100
result_AB_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
result_AB_d_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
for (sample_weight in sample_weight_list) {
  parameters = get_parameters(sample_weight, px = 24, r_RR = 5)
  C = parameters$C
  A = parameters$A
  B = parameters$B
  A_d = parameters$FALSEA_d
  B_d = parameters$B_d
  C_r = parameters$C_r
  py = parameters$py
  px = parameters$px

  # norm_A_list = norm_A_d_list = rep(NA, eloop)
  # norm_B_list = norm_B_d_list = rep(NA, eloop)

  bias_AB_matrix = bias_AB_d_matrix = matrix(NA, nrow = py*px, ncol = eloop)
  B_star = t(B) # in order to compute the norm
  B_d_star = t(B_d)
  for (i in 1:eloop) {
    simulation_result = simulation(parameters, n=1000, penalized = penalized)
    A_hat = simulation_result[[1]]
    B_hat = simulation_result[[2]]
    B_hat_star = t(simulation_result[[2]])
    AB_hat = simulation_result[[3]]

    # temp_A = A_hat %*% solve(t(A_hat) %*% A_hat) %*% t(A_hat) -
    #   A_d %*% solve(t(A_d) %*% A_d) %*% t(A_d)
    # temp_B = B_hat_star %*% solve(t(B_hat_star) %*% B_hat_star) %*% t(B_hat_star) -
    #   B_d_star %*% solve(t(B_d_star) %*% B_d_star) %*% t(B_d_star)
    # norm_A_list[i] = norm(temp_A, type = "F")
    # norm_B_list[i] = norm(temp_B, type = "F")

    # bias_AB_vectorized = as.vector(AB_hat - C)
    bias_AB_vectorized = as.vector(AB_hat - C_r)
    bias_AB_matrix[,i] = bias_AB_vectorized

    A_d_hat = simulation_result[[4]]
    B_d_hat = simulation_result[[5]]
    B_d_hat_star = t(simulation_result[[5]])
    AB_d_hat = simulation_result[[6]]

    # temp_A_d = A_d_hat %*% solve(t(A_d_hat) %*% A_d_hat) %*% t(A_d_hat) -
    #   A_d %*% solve(t(A_d) %*% A_d) %*% t(A_d)
    # temp_B_d = B_d_hat_star %*% solve(t(B_d_hat_star) %*% B_d_hat_star) %*% t(B_d_hat_star) -
    #   B_d_star %*% solve(t(B_d_star) %*% B_d_star) %*% t(B_d_star)
    # norm_A_d_list[i] = norm(temp_A_d, type = "F")
    # norm_B_d_list[i] = norm(temp_B_d, type = "F")

    bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
    # bias_AB_d_vectorized = as.vector(AB_d_hat - C)
    bias_AB_d_matrix[,i] = bias_AB_d_vectorized
  }
  char <- as.character(sample_weight)
  result_AB_list[[char]] = bias_AB_matrix
  result_AB_d_list[[char]] = bias_AB_d_matrix

  # result_A_list[[char]] = norm_A_list
  # result_B_list[[char]] = norm_B_list
  # result_A_d_list[[char]] = norm_A_d_list
  # result_B_d_list[[char]] = norm_B_d_list
}

#####
# choose the weight to plot (1, 0.5, 0.2, 0.1, 0.05)
weight_plot = "0.2"
# plot the results -------------------------------------------------------------
# functions
# a function to get location of a matrix from the vector index
get_index <- function(i) {
  return(c(i - ((i-1) %/% py) * py, (i-1) %/% py + 1))
}


plot_entrybias = function (data) {
  ggplot()+
    geom_boxplot(aes(y = data))+
    labs(title = "AB_hat - C_r entry", y = "AB_hat - C_r entry")
}


plot_entrybias_d = function (data) {
  ggplot()+
    geom_boxplot(aes(y = data))+
    labs(title = "debiased_AB_hat - C_r entry", y = "debiased_AB_hat - C_r entry")
}

# boxplot the bias of AB_hat
weight_to_plot <- weight_plot # choose the weight to plot for AB

myplots <- apply(result_AB_list[[weight_to_plot]], MARGIN = 1, plot_entrybias)
for (i in 1:(px*py)){
  names(myplots)[i] <- paste0("(", paste(get_index(i)[1], sep = ", ", get_index(i)[2]), ")")
}

abs_mean_entry_bias <- rep(NA, (px*py))
mean_entry_bias <- rep(NA, (px*py))
sd_entry_bias <- rep(NA, (px*py))
for (i in 1:(px*py)) {
  abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_AB_list[[weight_to_plot]][i,])))
  mean_entry_bias[i] <- mean(result_AB_list[[weight_to_plot]][i,])
  sd_entry_bias[i] <- sd(result_AB_list[[weight_to_plot]][i,])
}

index1 <- order(abs_mean_entry_bias, decreasing = T)
index2 <- order(sd_entry_bias, decreasing = T)
# myplots[index1[5:1]]

# boxplot the bias of AB_d_hat
weight_to_plot <- weight_plot
myplots_d <- apply(result_AB_d_list[[weight_to_plot]], MARGIN = 1, plot_entrybias_d)
for (i in 1:(px*py)){
  names(myplots_d)[i] <- paste0("(", paste(get_index(i)[1], sep = ", ", get_index(i)[2]), ")")
}

abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
sd_entry_bias_d <- rep(NA, (px*py))
for (i in 1:(px*py)) {
  abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_AB_d_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
  mean_entry_bias_d[i] <- mean(result_AB_d_list[[weight_to_plot]][i,])
  sd_entry_bias_d[i] <- sd(result_AB_d_list[[weight_to_plot]][i,])
}

index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
index2_d <- order(sd_entry_bias_d, decreasing = T)
# myplots_d[index1_d[5:1]]


# create heatmap for all entries
# note: need to run the above "boxplot the bias of AB_hat(AB_d_hat)" first
# bias heatmap
# RR
avg_bias <- matrix(abs_mean_entry_bias, nrow = px, ncol = py)
avg_bias <- t(avg_bias)
avg_bias <- as.data.frame(avg_bias)
colnames(avg_bias) <- 1:px
rownames(avg_bias) <- 1:py

# dRR
avg_bias_d <- matrix(abs_mean_entry_bias_d, nrow = px, ncol = py)
avg_bias_d <- t(avg_bias_d)
avg_bias_d <- as.data.frame(avg_bias_d)
colnames(avg_bias_d) <- 1:px
rownames(avg_bias_d) <- 1:py

# color breaks
min_value <- min(min(avg_bias), min(avg_bias_d)) #TODO:
max_value <- max(max(avg_bias),max(avg_bias_d))

# min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
# max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))

breaks <- seq(min_value, max_value, length.out = 101)

pheatmap(avg_bias, breaks = breaks,
         main = sprintf("AB_hat - C_r entry-wise bias' absolute value, weight = %s, px = %s",weight_plot, px), fontsize = 8,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100))
pheatmap(avg_bias_d, breaks = breaks,
         main = sprintf("AB_d_hat - C_r entry-wise bias' absolute value, weight = %s, px = %s",weight_plot, px), fontsize = 8,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100))


# sd heatmap
# RR
avg_sd <- sd_entry_bias
avg_sd <- matrix(avg_sd, nrow = px, ncol = py)
avg_sd <- t(avg_sd)
avg_sd <- as.data.frame(avg_sd)
colnames(avg_sd) <- 1:px
rownames(avg_sd) <- 1:py

# dRR
avg_sd_d <- sd_entry_bias_d
avg_sd_d <- matrix(avg_sd_d, nrow = px, ncol = py)
avg_sd_d <- t(avg_sd_d)
avg_sd_d <- as.data.frame(avg_sd_d)
colnames(avg_sd_d) <- 1:px
rownames(avg_sd_d) <- 1:py


# color breaks
# min_value_sd <- min(min(avg_sd), min(avg_sd_d))
# breaks_sd <- c(seq(min_value_sd, max(avg_sd), length.out = 51),
#                seq(max(avg_sd), max(avg_sd_d), length.out = 51)[-1])


min_value_sd <- min(min(avg_sd), min(avg_sd_d))
max_value_Sd <- max(max(avg_sd),max(avg_sd_d))
breaks_sd <- seq(min_value_sd, max_value_Sd, length.out = 101)


pheatmap(avg_sd, breaks = breaks_sd,
         main = sprintf("AB_hat - C_r entry-wise sd, weight = %s, px = %s",weight_plot,px), fontsize = 8,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100))
pheatmap(avg_sd_d, breaks = breaks_sd,
         main = sprintf("AB_d_hat - C_r entry-wise sd, weight = %s, px = %s",weight_plot,px), fontsize = 8,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "red"))(100))

# # top 5 box plot in the same plot
# combined_plot <- myplots_d[[index1_d[5]]] + myplots_d[[index1_d[4]]] +
#   myplots_d[[index1_d[3]]] + myplots_d[[index1_d[2]]] + myplots_d[[index1_d[1]]] +
#   myplots[[index1[5]]] + myplots[[index1[4]]] + myplots[[index1[3]]] +
#   myplots[[index1[2]]] + myplots[[index1[1]]] + plot_layout(nrow = 2, ncol = 5) +
#   plot_annotation(
#     title = "Top five entry-wise difference between C_r_hat (= AB_hat) and C_r, ranked by mean",
#     subtitle = sprintf("px = 100, rank C = reduced rank = 5, Sigma_X multiplied by %s", weight_plot),
#     caption = "C_r defined as minimizing modified objective function of
#     2.2 or equivalent objective function of 2.1")
# combined_plot
#
# # rank by var
# combined_plot_var <- myplots_d[[index2_d[5]]] + myplots_d[[index2_d[4]]] +
#   myplots_d[[index2_d[3]]] + myplots_d[[index2_d[2]]] + myplots_d[[index2_d[1]]] +
#   myplots[[index2[5]]] + myplots[[index2[4]]] + myplots[[index2[3]]] +
#   myplots[[index2[2]]] + myplots[[index2[1]]] + plot_layout(nrow = 2, ncol = 5) +
#   plot_annotation(
#     title = "Top five entry-wise difference between C_r_hat (= AB_hat) and C_r, ranked by var",
#     subtitle = sprintf("px = 100, rank C = reduced rank = 5, Sigma_X multiplied by %s", weight_plot),
#     caption = "C_r defined as minimizing modified objective function of
#     2.2 or equivalent objective function of 2.1")
# combined_plot_var

parameters = get_parameters(as.numeric(weight_plot), px = 24, r_RR = 5) # weight: 0.2, 0.5, 1 are comparable with Yinxiang's Paper
min(eigen(solve(sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(sqrt_matrix(parameters$Sigma_X)))$values)
mean(as.matrix(avg_bias_d))
mean(as.matrix(avg_bias))
mean(as.matrix(avg_sd_d))
mean(as.matrix(avg_sd))

# bias scale compares to true C_r
mean(abs(as.matrix(C_r)))
mean(abs(as.matrix(avg_bias)))
mean(abs(as.matrix(avg_bias_d)))
#####


# penalized = TRUE
# # parametric bootstrap method to estimate variance of the Drr estimator ----
# #### (a) get the parameters
# parameters = get_parameters(sample_weight = 1, px = 24, r_RR = 5)
# C = parameters$C
# A = parameters$ATRUE
# B = parameters$B
# A_d = parameters$A_d
# B_d = parameters$B_d
# C_r = parameters$C_r
# py = parameters$py
# px = parameters$px
# Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
# Sigma_Y_hat = parameters$Sigma_Y
# 
# #### (b) in each iteration, gen data
# iteration = 100
# # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
# coverage_rate = matrix(NA, px*py, iteration)
# 
# for (loop in 1:iteration){
#   #### (b1) gen data in the trial
#   simulation_result = simulation(parameters, n=1000, penalized = penalized)
#   A_hat = simulation_result[[1]]
#   B_hat = simulation_result[[2]]
#   AB_hat = simulation_result[[3]]
#   A_d_hat = simulation_result[[4]]
#   B_d_hat = simulation_result[[5]]
#   AB_d_hat = simulation_result[[6]]
#   x_j_hat = simulation_result[[7]]
#   y_j_hat = simulation_result[[8]]
# 
#   #### (b2) estimate V_X, C
#   V_X_tilde_star = t(x_j_hat) %*% x_j_hat / 1000 - Sigma_X_hat
#   C_star = AB_d_hat
#   C_r_vec = as.vector(C_r)
# 
#   #### (b3) bootstrap ----
#   bootstrap_size = 100  # sample size to calculate CI in each bt loop
#   pz = 1000
#   C_drr_hat_entrywise_bt_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_list[[i]] = rep(0, bootstrap_size)
#   }
#   #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
#   for (i in 1:bootstrap_size){
#     # sample true effect x_j_star and y_j_star
#     x_j_star = mvrnorm(n = pz, mu = rep(0, px), Sigma = V_X_tilde_star, tol = 100)
#     y_j_star = x_j_star %*% t(C)
# 
#     # sample x_j_hat, y_j_hat
#     x_j_hat_star = matrix(0, pz, px)
#     y_j_hat_star = matrix(0, pz, py)
#     for (j in 1:pz) {
#       x_j_hat_star[j,] = mvrnorm(n = 1, mu = x_j_star[j,], Sigma = Sigma_X_hat, tol = 100)
#       y_j_hat_star[j,] = mvrnorm(n = 1, mu = y_j_star[j,], Sigma = Sigma_Y_hat, tol = 100)
#     }
# 
#     # estimate C_drr_hat
#     if (penalized == FALSE) {
#       result_d_bt <- debiased_RRR(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat)
#     } else {
#       result_d_bt <- debiased_RRR_modified(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat, ph=regularization_rate)
#     }
#     result_d_bt_vec = as.vector(result_d_bt$AB)
#     for (j in 1:(px*py)){
#       C_drr_hat_entrywise_bt_list[[j]][i] = result_d_bt_vec[j]
#     }
#   }
# 
#   #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
#   C_drr_hat_entrywise_bt_95_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_95_list[[i]] = quantile(C_drr_hat_entrywise_bt_list[[i]], c(0.05, 0.95))
#   }
#   #### (b5) calculate the coverage rate for each entry of C_r
#   C_r_vec = as.vector(C_r)
#   index_in95 = rep(1, px*py)
#   for (i in 1:(px*py)){
#     if (C_r_vec[i] < C_drr_hat_entrywise_bt_95_list[[i]][1] | C_r_vec[i] > C_drr_hat_entrywise_bt_95_list[[i]][2]){
#       index_in95[i] = 0
#     }
#   }
#   coverage_rate[,loop] = index_in95
# }
# 
# #### (c) calculate the average coverage rate
# entry_coverage_rate = rowMeans(coverage_rate)
# average_coverage_rate = mean(entry_coverage_rate)
# average_coverage_rate
# #####
# ## result:
# # weight=1, penalized=FALSE, average_coverage_rate=1 (var too large)
# # weight=1, penalized=TRUE, ph=1e-11 , average_coverage_rate=0.6852917
# # weight=1, penalized=TRUE, ph=1e-12 , average_coverage_rate=0.9374167
# # weight=1, penalized=TRUE, ph=1e-13 , average_coverage_rate=0.9910417
# 
# 
# 
# penalized = TRUE
# # non-parametric bootstrap method to estimate variance of the Drr estimator ----
# #### (a) get the parameters
# parameters = get_parameters(sample_weight = 1, px = 24, r_RR = 5)
# C = parameters$C
# A = parameters$A
# B = parameters$B
# A_d = parameters$A_d
# B_d = parameters$B_d
# C_r = parameters$C_r
# C_r_vec = as.vector(C_r)
# py = parameters$py
# px = parameters$px
# Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
# Sigma_Y_hat = parameters$Sigma_Y
# 
# #### (b) in each iteration, gen data
# iteration = 100
# # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
# coverage_rate = matrix(NA, px*py, iteration)
# 
# for (loop in 1:iteration){
#   #### (b1) gen data in the trial
#   simulation_result = simulation(parameters, n=1000, penalized = penalized)
#   A_hat = simulation_result[[1]]
#   B_hat = simulation_result[[2]]
#   AB_hat = simulation_result[[3]]
#   A_d_hat = simulation_result[[4]]
#   B_d_hat = simulation_result[[5]]
#   AB_d_hat = simulation_result[[6]]
#   x_j_hat = simulation_result[[7]]
#   y_j_hat = simulation_result[[8]]
#   
#   #### (b2) estimate V_X, C
#   
#   #### (b3) bootstrap ----
#   bootstrap_size = 100  # sample size to calculate CI in each bt loop
#   pz = 1000
#   C_drr_hat_entrywise_bt_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_list[[i]] = rep(0, bootstrap_size)
#   }
#   #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
#   for (bt in 1:bootstrap_size){
#     # sample bt_sample_number rows from x_j_hat, y_j_hat
#     index_bt = sample(1:pz, pz, replace = TRUE)
#     x_j_hat_star = x_j_hat[index_bt,]
#     y_j_hat_star = y_j_hat[index_bt,]
#     
#     # estimate C_drr_hat
#     if (penalized == FALSE) {
#       result_d_bt <- debiased_RRR(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat)
#     } else {
#       result_d_bt <- debiased_RRR_modified(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat, ph=regularization_rate)
#     }
#     result_d_bt_vec = as.vector(result_d_bt$AB)
#     for (j in 1:(px*py)){
#       C_drr_hat_entrywise_bt_list[[j]][bt] = result_d_bt_vec[j]
#     }
#   }
#   
#   #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
#   C_drr_hat_entrywise_bt_95_list = list()
#   for (i in 1:(px*py)){
#     C_drr_hat_entrywise_bt_95_list[[i]] = quantile(C_drr_hat_entrywise_bt_list[[i]], c(0.025, 0.975))
#   }
#   #### (b5) calculate the coverage rate for each entry of C_r
#   index_in95 = rep(1, px*py)
#   for (i in 1:(px*py)){
#     if (C_r_vec[i] < C_drr_hat_entrywise_bt_95_list[[i]][1] | C_r_vec[i] > C_drr_hat_entrywise_bt_95_list[[i]][2]){
#       index_in95[i] = 0
#     }
#   }
#   coverage_rate[,loop] = index_in95
# }
# 
# #### (c) calculate the average coverage rate
# entry_coverage_rate = rowMeans(coverage_rate)
# average_coverage_rate = mean(entry_coverage_rate)
# average_coverage_rate
# 
# 
# #####
# 
# ## result:
# # weight=1, penalized=FALSE, average_coverage_rate=0.96625
# # weight=1, penalized=TRUE, ph=1e-11 , average_coverage_rate=0.6511667
# # weight=1, penalized=TRUE, ph=1e-13 , average_coverage_rate=0.8910833
# 
# 
# # archive -----
# # penalized = FALSE
# # m = 800
# # # m-out-of-n bootstrap without replacement method to estimate variance of the Drr estimator
# # #### (a) get the parameters
# # parameters = get_parameters(sample_weight = 1, px = 24, r_RR = 5)
# # C = parameters$C
# # A = parameters$A
# # B = parameters$B
# # A_d = parameters$A_d
# # B_d = parameters$B_d
# # C_r = parameters$C_r
# # C_r_vec = as.vector(C_r)
# # py = parameters$py
# # px = parameters$px
# # Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
# # Sigma_Y_hat = parameters$Sigma_Y
# # 
# # #### (b) in each iteration, gen data
# # iteration = 1000
# # # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
# # coverage_rate = matrix(NA, px*py, iteration)
# # 
# # for (loop in 1:iteration){
# #   #### (b1) gen data in the trial
# #   simulation_result = simulation(parameters, n=1000, penalized = penalized)
# #   A_hat = simulation_result[[1]]
# #   B_hat = simulation_result[[2]]
# #   AB_hat = simulation_result[[3]]
# #   A_d_hat = simulation_result[[4]]
# #   B_d_hat = simulation_result[[5]]
# #   AB_d_hat = simulation_result[[6]]
# #   x_j_hat = simulation_result[[7]]
# #   y_j_hat = simulation_result[[8]]
# #   
# #   #### (b2) estimate V_X, C
# #   
# #   #### (b3) bootstrap
# #   bootstrap_size = 10000  # sample size to calculate CI in each bt loop
# #   pz = 1000
# #   C_drr_hat_entrywise_bt_list = list()
# #   for (i in 1:(px*py)){
# #     C_drr_hat_entrywise_bt_list[[i]] = rep(0, bootstrap_size)
# #   }
# #   #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
# #   for (bt in 1:bootstrap_size){
# #     # sample bt_sample_number rows from x_j_hat, y_j_hat
# #     index_bt = sample(1:pz, m, replace = FALSE)
# #     x_j_hat_star = x_j_hat[index_bt,]
# #     y_j_hat_star = y_j_hat[index_bt,]
# #     
# #     # estimate C_drr_hat
# #     if (penalized == FALSE) {
# #       result_d_bt <- debiased_RRR(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat)
# #     } else {
# #       result_d_bt <- debiased_RRR_modified(y_j_hat_star, x_j_hat_star, r=5, sqrt_Gamma = sqrt_matrix(solve(Sigma_Y_hat)), Sigma_X = Sigma_X_hat, ph=regularization_rate)
# #     }
# #     result_d_bt_vec = as.vector(result_d_bt$AB)
# #     for (j in 1:(px*py)){
# #       C_drr_hat_entrywise_bt_list[[j]][bt] = result_d_bt_vec[j]
# #     }
# #   }
# #   
# #   #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
# #   C_drr_hat_entrywise_bt_95_list = list()
# #   for (i in 1:(px*py)){
# #     C_drr_hat_entrywise_bt_95_list[[i]] = quantile(C_drr_hat_entrywise_bt_list[[i]], c(0.05, 0.95))
# #   }
# #   #### (b5) calculate the coverage rate for each entry of C_r
# #   index_in95 = rep(1, px*py)
# #   for (i in 1:(px*py)){
# #     if (C_r_vec[i] < C_drr_hat_entrywise_bt_95_list[[i]][1] | C_r_vec[i] > C_drr_hat_entrywise_bt_95_list[[i]][2]){
# #       index_in95[i] = 0
# #     }
# #   }
# #   coverage_rate[,loop] = index_in95
# # }
# # 
# # #### (c) calculate the average coverage rate
# # entry_coverage_rate = rowMeans(coverage_rate)
# # average_coverage_rate = mean(entry_coverage_rate)
# # average_coverage_rate
# # 
# # 
# # #####




# new sim for sparse 10.30 ----
soft_threshold <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}


parameters = get_parameters(sample_weight=0.1, px = 24, r_RR = 5)
pz = n = 100 #(pz)
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
W_sqrt = sqrt_matrix(parameters$weight.matrix)
true.A_sparse = parameters$true.A_sparse
true.B_sparse = parameters$true.B_sparse

# # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
# gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
# Gamma_j_star = gamma_j_star %*% t(C)
# 
# # sample x_j_hat, y_j_hat
# x_j_hat = matrix(0, n, px)
# y_j_hat = matrix(0, n, py)
# for (j in 1:n) {
#   x_j_hat[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#   y_j_hat[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
# }

# compute sample cov matrix
hat_Sigma_xhat_xhat = t(x_j_hat) %*% x_j_hat / n
hat_Sigma_x_x = hat_Sigma_xhat_xhat - Sigma_X

# check if PSD
is_psd <- function(matrix) {
  eigenvalues <- eigen(matrix)$values
  all(eigenvalues >= 0)
}


is_psd(hat_Sigma_x_x)

# start to implement algorithm
# step 1: debiasing the objective function
# chol decomposition
tilde_gamma = t(chol(hat_Sigma_x_x)) * sqrt(n)
# check if chol is correct
all.equal(tilde_gamma %*% t(tilde_gamma) / n, hat_Sigma_x_x)

tilde_Gamma = t(solve(tilde_gamma)%*%t(x_j_hat)%*%y_j_hat)

# check the definition of tilde_Gamma
all.equal(tilde_gamma %*% t(tilde_Gamma), t(x_j_hat)%*%y_j_hat)

# step 2: solve the objective function
# store iteration values of A and B
A_hat_list = list()
B_hat_list = list()
# initial values of A_hat and B_hat using dRRR
result <- debiased_RRR(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X)
A_hat = result$A
B_hat = result$B
A_hat_list[[1]] = A_hat
B_hat_list[[1]] = B_hat
t(A_hat) %*% W %*% A_hat


# given B minimize A:  Orthogonal Procrustes ----
U = svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)$u
V = svd(B_hat %*% tilde_gamma %*% t(tilde_Gamma) %*% W_sqrt)$v

A_hat = solve(W_sqrt) %*% V %*% t(U)
# regularization, make t(A_hat) %*% W %*% A_hat= I

A_hat_list = c(A_hat_list, list(A_hat))
#####

# given A minimize B:  Cyclic Coordinate Descent ----

# 计算目标函数的第一部分
first_term <- (1 / n) * sum(sapply(1:px, function(j) {
  W_j <- W[, j, drop = FALSE]       # 取矩阵 W 的第 j 列
  gamma_tilde_j <- tilde_gamma[, j, drop = FALSE] # 取矩阵 tilde_gamma 的第 j 列
  norm(t(A_hat) %*% W_j - B_hat %*% gamma_tilde_j, type = "2")^2
}))

# 计算目标函数的第二部分
second_term <- sum(sapply(1:px, function(k) {
  lambda[k] * sum(abs(B_hat[, k])) # 对矩阵 B_hat 的第 k 列计算 L1 范数
}))

# 最终目标函数的值
objective_function_value <- first_term + second_term



# coordinate decent
# tilde_gamma = x_j_hat
# tilde_Gamma = t(y_j_hat)
# W = diag(1, py)
max_change = 1000
iter = 1
# B_hat = matrix(0, r_RR, px)
B_hat_list[[1]] = B_hat
# while iter too large also stop
while (iter < 1000) {
  for (k in 1:px){
    # for given k, compute all r_j, j=1,...,p_Y
    updata_k = 0
    for (j in 1:px){ # px
      # compute the summation of the ith column vect of B_hat times the ji entry of tilde gamma
      sum = 0
      for (i in 1:px){
        if (i != k){
          sum = sum + B_hat[,i] * tilde_gamma[j,i]
        }
      }
      r_j = t(A_hat) %*% W_sqrt %*% W_sqrt %*% tilde_Gamma[,j] - sum
      r_j_times_tilde_gamma = r_j * tilde_gamma[j,k]
      updata_k = updata_k + r_j_times_tilde_gamma
    }
    updata_k = updata_k/pz
    # update B_hat
    B_hat[,k] = soft_threshold(B_hat[,k] + updata_k, lambda[k])
  }
  mean_change = mean(abs(B_hat-B_hat_list[[iter]]))
  B_hat_list[[iter+1]] = B_hat
  iter = iter +1
  # print iter and max_change
  print(c(iter, mean_change))
}


# TODO: test extreme case




# we should see the change of objective function value



# standard coordinate decent (univariate y)----
# Implementation of Cyclic Coordinate Descent for solving Lasso

# Soft-thresholding operator
soft_threshold <- function(x, lambda) {
  return(sign(x) * pmax(abs(x) - lambda, 0))
}

# Generate some example data
set.seed(123)
X <- matrix(rnorm(100 * 10), 100, 10)  # 100 samples, 10 predictors
y <- rnorm(100)  # response vector
lambda <- 0.1  # regularization parameter

# Initialize parameters
beta <- rep(0, ncol(X))  # initial coefficients set to zero
max_iter <- 1000  # maximum number of iterations
tol <- 1e-6  # convergence tolerance

# Coordinate Descent Algorithm
for (iter in 1:max_iter) {
  beta_old <- beta
  for (j in 1:ncol(X)) {
    # Compute partial residual without the contribution from the j-th predictor
    r_j <- y - X[, -j] %*% beta[-j]
    # Update the j-th coefficient using the soft-thresholding operator
    beta[j] <- soft_threshold(mean(X[, j] * r_j), lambda)
  }
  # Check for convergence
  if (max(abs(beta - beta_old)) < tol) {
    cat("Algorithm converged after", iter, "iterations.\n")
    break
  }
}

# Print the final coefficients
cat("Estimated coefficients:", beta, "\n")

# old ver----

# may not apply, this is group lasso
# AT_W_TildeG = t(A_hat) %*% W %*% tilde_Gamma
# fit <- glmnet(t(tilde_gamma), t(AT_W_TildeG), family = "mgaussian", alpha = 1, intercept = FALSE)
# # output B hat as a matrix
# # empty
# B_hat = matrix(0, px, py) 
# for (i in length(coef(fit, s = 0.1))){
#   if (i == 1){
#     B_hat = coef(fit, s = 0.1)$y[,1][-1]
#   } else {
#     B_hat = cbind(B_hat, coef(fit, s = 0.1)$beta[,i][-1])
#   }
# }
# 
# beta_final <- lapply(fit$beta, function(mat) mat[, ncol(mat)])
# B_hat <- t(do.call(cbind, beta_final))
# B_hat_list = c(B_hat_list, B_hat)

#####