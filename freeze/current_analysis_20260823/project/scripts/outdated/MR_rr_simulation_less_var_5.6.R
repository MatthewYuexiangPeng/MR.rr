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

lip_data_real = read.csv('data/dat_1e-4.csv')
lip_corr_real = read.csv('data/rho_mat_1e-4.csv')
n_Y = c(1296908,1241207,1245612,1241619)

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


.get_parameters = function(me_weight, r_RR = 2, px=9, r_approx = 2){
  # ## test
  # me_weight =1
  # r_RR = 2
  # px=9
  # r_approx = 2
  
  # r_RR is the rank chose by user when performing RRR
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 1000 # pz can be changed to any number
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
  
  Sigma_X = me_weight * Sigma_X
  
  # sample cov
  VX_tilde <- cov(gamma_j_star_temp) - Sigma_X
  
  
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
  
  # # 24.11.7: generate C from A times B ----
  # M <- matrix(rnorm(py * r_approx), nrow = py)
  # svd_result <- svd(M)
  # U <- svd_result$u
  # A <- solve(.sqrt_matrix(weight.matrix)) %*% U[, 1:r_approx]
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
              result_MrDAG#,
              # result_GRAPPLE
              ))
}


#### choose regularization rate \lambda ####
# [1] 5.936180e-11 3.457282e-12 1.229468e-18
c=10
me_weight_list = c(2.5,2,1)
eloop = 100
result_matrix = matrix(NA, nrow = eloop, ncol = length(me_weight_list))
for (i in 1:eloop){
  result_list = list()
  for (me_weight in me_weight_list){
    # me_weight = 2.5
    parameters = .get_parameters(me_weight=me_weight, r_RR = 2, px=9, r_approx = 2)
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
    
    mu_min = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
    # mu_min * sqrt(1000)
    # iv_strength = mu_min * sqrt(1000)
    
    q_list = seq(1, 45, 1)
    regu_rate_list = c()
    obj_value = c()
    for (q in q_list) {
      sigma_y2 = mean(eigen(parameters$Sigma_X)$values)
      
      regu_rate = sigma_y2^2 * exp(c * (q/sqrt(n)-(mu_min+1)))/n
      
      # regu_rate = c * sigma_y2^2 * (q-(mu_min+1)/n)
      # regu_rate = exp(c*(q-mu_min*sqrt(n)))
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
  result_matrix[i,] = unlist(result_list)
}
# average
result_matrix_mean = apply(result_matrix, 2, mean)
result_matrix_mean


#### run simulation ####
# [1] 1.426422e-12 7.579985e-13 4.377867e-19
regularization_rate_list = c(1.426422e-12, 7.579985e-13, 4.377867e-19)
# The simulation function for the naive MR-rr estimator and the MR-rr estimator (with spectral regularization)
run_simulation <- function(regularization_rate_list){
  sample_weight_list = c(2.5, 2, 1)
  eloop = 500
  
  result_AB_list = result_AB_d_list = result_AB_d_r_list = 
    result_C_ivw_list = result_C_adivw_list = result_MrDAG_list = # result_C_GRAPPLE_list =
    list("2.5" = NA, "2" = NA, "1" = NA)
  
  iv_strength_list = c()
  parameters_list = c()
  
  for (weight_index in seq_along(sample_weight_list)) {
    me_weight = sample_weight_list[weight_index]
    parameters = .get_parameters(me_weight, px = 9, r_RR = 2)
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
      matrix(NA, nrow = py*px, ncol = eloop)
    B_star = t(B) # in order to compute the norm
    B_d_star = t(B_d)
    for (i in 1:eloop) {
      simulation_result = .simulation(parameters, 
                                      regularization_rate=regularization_rate_list[weight_index])
      A_hat = simulation_result[[1]]
      B_hat = simulation_result[[2]]
      B_hat_star = t(simulation_result[[2]])
      AB_hat = simulation_result[[3]]
      
      bias_AB_vectorized = as.vector(AB_hat - C)
      # bias_AB_vectorized = as.vector(AB_hat - C_r)
      bias_AB_matrix[,i] = bias_AB_vectorized
      
      A_d_hat = simulation_result[[4]]
      B_d_hat = simulation_result[[5]]
      B_d_hat_star = t(simulation_result[[5]])
      AB_d_hat = simulation_result[[6]]
      
      # bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
      bias_AB_d_vectorized = as.vector(AB_d_hat - C)
      bias_AB_d_matrix[,i] = bias_AB_d_vectorized
      
      A_d_r_hat = simulation_result[[7]]
      B_d_r_hat = simulation_result[[8]]
      B_d_r_hat_star = t(simulation_result[[8]])
      AB_d_r_hat = simulation_result[[9]]
      
      # bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
      bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C)
      bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
      
      C_ivw = simulation_result[[12]]
      bias_C_ivw_vectorized = as.vector(C_ivw - C)
      bias_C_ivw_matrix[,i] = bias_C_ivw_vectorized
      
      C_adivw = simulation_result[[13]]
      bias_C_adivw_vectorized = as.vector(C_adivw - C)
      bias_C_adivw_matrix[,i] = bias_C_adivw_vectorized
      
      C_MrDAG = simulation_result[[14]]
      bias_C_MrDAG_vectorized = as.vector(C_MrDAG - C)
      bias_C_MrDAG_matrix[,i] = bias_C_MrDAG_vectorized
      
      # C_GRAPPLE = simulation_result[[15]]
      # bias_C_GRAPPLE_vectorized = as.vector(C_GRAPPLE - C)
      # bias_C_GRAPPLE_matrix[,i] = bias_C_GRAPPLE_vectorized
      
      # print "eloop" and eloop number in a line
      print(paste("eloop:", i, "weight:", me_weight))
    }
    char <- as.character(me_weight)
    result_AB_list[[char]] = bias_AB_matrix
    result_AB_d_list[[char]] = bias_AB_d_matrix
    result_AB_d_r_list[[char]] = bias_AB_d_r_matrix
    result_C_ivw_list[[char]] = bias_C_ivw_matrix
    result_C_adivw_list[[char]] = bias_C_adivw_matrix
    result_MrDAG_list[[char]] = bias_C_MrDAG_matrix
    # result_C_GRAPPLE_list[[char]] = bias_C_GRAPPLE_matrix
  }
  return(list(result_AB_list=result_AB_list, 
              result_AB_d_list=result_AB_d_list, 
              result_AB_d_r_list=result_AB_d_r_list,
              iv_strength_list=iv_strength_list, parameters_list=parameters_list,
              result_C_ivw_list=result_C_ivw_list, 
              result_C_adivw_list=result_C_adivw_list,
              result_MrDAG_list=result_MrDAG_list#,
              # result_C_GRAPPLE_list=result_C_GRAPPLE_list
  ))
}


simulate_result = run_simulation(regularization_rate_list = regularization_rate_list)  # TODO: choose the regularization rate for each weight based on objective function
save(simulate_result, file = "results/simulate_result_250405_2.RData")
round(simulate_result$iv_strength_list, 3)


#### plot heatmap/ generating table ####

load("results/simulate_result_250405.RData")


# The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
plot_heatmap <- function(parameters_list,
                         result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
                         result_ivw_list, result_adivw_list, result_MrDAG_list,
                         weight_to_plot)
{
  px <- 9
  py <- 4
  # naive MR-rr
  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }
  
  # MR-rr
  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
  }
  
  # MR-rr with spectral regularization
  abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d_r <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
  }
  
  # IVW
  abs_mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_ivw <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_ivw[i] <- abs(as.numeric(mean(result_ivw_list[[weight_to_plot]][i,])))
    mean_entry_bias_ivw[i] <- mean(result_ivw_list[[weight_to_plot]][i,])
    sd_entry_bias_ivw[i] <- sd(result_ivw_list[[weight_to_plot]][i,])
  }
  
  # adIVW
  abs_mean_entry_bias_adivw <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_adivw <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_adivw <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_adivw[i] <- abs(as.numeric(mean(result_adivw_list[[weight_to_plot]][i,])))
    mean_entry_bias_adivw[i] <- mean(result_adivw_list[[weight_to_plot]][i,])
    sd_entry_bias_adivw[i] <- sd(result_adivw_list[[weight_to_plot]][i,])
  }
  
  # MrDAG
  abs_mean_entry_bias_MrDAG <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_MrDAG <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_MrDAG <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_MrDAG[i] <- abs(as.numeric(mean(result_MrDAG_list[[weight_to_plot]][i,])))
    mean_entry_bias_MrDAG[i] <- mean(result_MrDAG_list[[weight_to_plot]][i,])
    sd_entry_bias_MrDAG[i] <- sd(result_MrDAG_list[[weight_to_plot]][i,])
  }
  
  
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
  
  # dRR_r
  avg_bias_d_r <- matrix(abs_mean_entry_bias_d_r, nrow = px, ncol = py)
  avg_bias_d_r <- t(avg_bias_d_r)
  avg_bias_d_r <- as.data.frame(avg_bias_d_r)
  colnames(avg_bias_d_r) <- 1:px
  rownames(avg_bias_d_r) <- 1:py
  
  # IVW
  avg_bias_ivw <- matrix(abs_mean_entry_bias_ivw, nrow = px, ncol = py)
  avg_bias_ivw <- t(avg_bias_ivw)
  avg_bias_ivw <- as.data.frame(avg_bias_ivw)
  colnames(avg_bias_ivw) <- 1:px
  rownames(avg_bias_ivw) <- 1:py
  
  # adIVW
  avg_bias_adivw <- matrix(abs_mean_entry_bias_adivw, nrow = px, ncol = py)
  avg_bias_adivw <- t(avg_bias_adivw)
  avg_bias_adivw <- as.data.frame(avg_bias_adivw)
  colnames(avg_bias_adivw) <- 1:px
  rownames(avg_bias_adivw) <- 1:py
  
  # MrDAG
  avg_bias_MrDAG <- matrix(abs_mean_entry_bias_MrDAG, nrow = px, ncol = py)
  avg_bias_MrDAG <- t(avg_bias_MrDAG)
  avg_bias_MrDAG <- as.data.frame(avg_bias_MrDAG)
  colnames(avg_bias_MrDAG) <- 1:px
  rownames(avg_bias_MrDAG) <- 1:py
  
  
  # color breaks
  min_value <- min(min(avg_bias), min(avg_bias_d), min(avg_bias_d_r))
  max_value <- max(max(avg_bias),max(avg_bias_d), max(avg_bias_d_r))
  
  # min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
  # max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))
  
  breaks <- seq(min_value, max_value, length.out = 101)
  
  hm_bias_1 = pheatmap::pheatmap(avg_bias, breaks = breaks,
                                 main = "Naive MR-rr estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_2 = pheatmap::pheatmap(avg_bias_d, breaks = breaks,
                                 main = "MR-rr estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_3 = pheatmap::pheatmap(avg_bias_d_r, breaks = breaks,
                                 main = "MR-rr estimator with regularization", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_4 = pheatmap::pheatmap(avg_bias_ivw, breaks = breaks,
                                 main = "IVW estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_5 = pheatmap::pheatmap(avg_bias_adivw, breaks = breaks,
                                 main = "adIVW estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_6 = pheatmap::pheatmap(avg_bias_MrDAG, breaks = breaks,
                                 main = "MrDAG estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  
  combined_hm_bias = gridExtra::grid.arrange(
    grobs = list(hm_bias_4, hm_bias_5, hm_bias_1, hm_bias_2, hm_bias_3, hm_bias_6),
    ncol = 6,
    top = grid::textGrob(sprintf("Absolute average bias by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
  )
  combined_hm_bias
  
  
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
  
  # dRR_r
  avg_sd_d_r <- sd_entry_bias_d_r
  avg_sd_d_r <- matrix(avg_sd_d_r, nrow = px, ncol = py)
  avg_sd_d_r <- t(avg_sd_d_r)
  avg_sd_d_r <- as.data.frame(avg_sd_d_r)
  colnames(avg_sd_d_r) <- 1:px
  rownames(avg_sd_d_r) <- 1:py
  
  # IVW
  avg_sd_ivw <- sd_entry_bias_ivw
  avg_sd_ivw <- matrix(avg_sd_ivw, nrow = px, ncol = py)
  avg_sd_ivw <- t(avg_sd_ivw)
  avg_sd_ivw <- as.data.frame(avg_sd_ivw)
  colnames(avg_sd_ivw) <- 1:px
  rownames(avg_sd_ivw) <- 1:py
  
  # adIVW
  avg_sd_adivw <- sd_entry_bias_adivw
  avg_sd_adivw <- matrix(avg_sd_adivw, nrow = px, ncol = py)
  avg_sd_adivw <- t(avg_sd_adivw)
  avg_sd_adivw <- as.data.frame(avg_sd_adivw)
  colnames(avg_sd_adivw) <- 1:px
  rownames(avg_sd_adivw) <- 1:py
  
  # MrDAG
  avg_sd_MrDAG <- sd_entry_bias_MrDAG
  avg_sd_MrDAG <- matrix(avg_sd_MrDAG, nrow = px, ncol = py)
  avg_sd_MrDAG <- t(avg_sd_MrDAG)
  avg_sd_MrDAG <- as.data.frame(avg_sd_MrDAG)
  colnames(avg_sd_MrDAG) <- 1:px
  rownames(avg_sd_MrDAG) <- 1:py
  
  min_value_sd <- min(min(avg_sd), min(avg_sd_d), min(avg_sd_d_r))
  max_value_Sd <- max(max(avg_sd),max(avg_sd_d), max(avg_sd_d_r))
  breaks_sd <- seq(min_value_sd, max_value_Sd, length.out = 101)
  
  
  hm_sd_1 = pheatmap::pheatmap(avg_sd, breaks = breaks_sd,
                               main = "Naive MR-rr estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_2 = pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
                               main = "MR-rr estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_3 = pheatmap::pheatmap(avg_sd_d_r, breaks = breaks_sd,
                               main = "MR-rr estimator with regularization", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_4 = pheatmap::pheatmap(avg_sd_ivw, breaks = breaks_sd,
                               main = "IVW estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_5 = pheatmap::pheatmap(avg_sd_adivw, breaks = breaks_sd,
                               main = "adIVW estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_6 = pheatmap::pheatmap(avg_sd_MrDAG, breaks = breaks_sd,
                               main = "MrDAG estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  
  combined_hm_sd = gridExtra::grid.arrange(
    grobs = list(hm_sd_4, hm_sd_5, hm_sd_1, hm_sd_2, hm_sd_3, hm_sd_6),
    ncol = 6,
    top = grid::textGrob(sprintf("SD by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
  )
  combined_hm_sd
  
  # bias scale compares to true C_r and sd
  weight_index = match(weight_to_plot, c("2.5", "2", "1"))
  parameters = parameters_list[[weight_index]]
  C_r = parameters$C_r
  # mean
  mean_abs_C_entry =  round(mean(abs(as.matrix(C_r))), 3)
  mean_abs_bias_naive = round(mean(abs(as.matrix(avg_bias))), 3)
  mean_abs_bias = round(mean(abs(as.matrix(avg_bias_d))), 3)
  mean_abs_bias_r = round(mean(abs(as.matrix(avg_bias_d_r))), 3)
  mean_abs_bias_ivw = round(mean(abs(as.matrix(avg_bias_ivw))), 3)
  mean_abs_bias_adivw = round(mean(abs(as.matrix(avg_bias_adivw))), 3)
  mean_abs_bias_MrDAG = round(mean(abs(as.matrix(avg_bias_MrDAG))), 3)
  
  mean_sd_naive = round(mean(as.matrix(avg_sd)), 3)
  mean_sd = round(mean(as.matrix(avg_sd_d)), 3)
  mean_sd_r = round(mean(as.matrix(avg_sd_d_r)), 3)
  mean_sd_ivw = round(mean(as.matrix(avg_sd_ivw)), 3)
  mean_sd_adivw = round(mean(as.matrix(avg_sd_adivw)), 3)
  mean_sd_MrDAG = round(mean(as.matrix(avg_sd_MrDAG)), 3)
  
  # median
  median_abs_bias_naive = round(median(abs(as.matrix(avg_bias))), 3)
  median_abs_bias = round(median(abs(as.matrix(avg_bias_d))), 3)
  median_abs_bias_r = round(median(abs(as.matrix(avg_bias_d_r))), 3)
  median_abs_bias_ivw = round(median(abs(as.matrix(avg_bias_ivw))), 3)
  median_abs_bias_adivw = round(median(abs(as.matrix(avg_bias_adivw))), 3)
  median_abs_bias_MrDAG = round(median(abs(as.matrix(avg_bias_MrDAG))), 3)
  
  median_sd_naive = round(median(as.matrix(avg_sd)), 3)
  median_sd = round(median(as.matrix(avg_sd_d)), 3)
  median_sd_r = round(median(as.matrix(avg_sd_d_r)), 3)
  median_sd_ivw = round(median(as.matrix(avg_sd_ivw)), 3)
  median_sd_adivw = round(median(as.matrix(avg_sd_adivw)), 3)
  median_sd_MrDAG = round(median(as.matrix(avg_sd_MrDAG)), 3)
  
  # quantile
  quantile_abs_bias_naive = round(quantile(abs(as.matrix(avg_bias)), c(0.25, 0.75)), 3)
  quantile_abs_bias = round(quantile(abs(as.matrix(avg_bias_d)), c(0.25, 0.75)), 3)
  quantile_abs_bias_r = round(quantile(abs(as.matrix(avg_bias_d_r)), c(0.25, 0.75)), 3)
  quantile_abs_bias_ivw = round(quantile(abs(as.matrix(avg_bias_ivw)), c(0.25, 0.75)), 3)
  quantile_abs_bias_adivw = round(quantile(abs(as.matrix(avg_bias_adivw)), c(0.25, 0.75)), 3)
  quantile_abs_bias_MrDAG = round(quantile(abs(as.matrix(avg_bias_MrDAG)), c(0.25, 0.75)), 3)
  
  quantile_sd_naive = round(quantile(as.matrix(avg_sd), c(0.25, 0.75)), 3)
  quantile_sd = round(quantile(as.matrix(avg_sd_d), c(0.25, 0.75)), 3)
  quantile_sd_r = round(quantile(as.matrix(avg_sd_d_r), c(0.25, 0.75)), 3)
  quantile_sd_ivw = round(quantile(as.matrix(avg_sd_ivw), c(0.25, 0.75)), 3)
  quantile_sd_adivw = round(quantile(as.matrix(avg_sd_adivw), c(0.25, 0.75)), 3)
  quantile_sd_MrDAG = round(quantile(as.matrix(avg_sd_MrDAG), c(0.25, 0.75)), 3)
  #####
  
  ggplot2::ggsave(sprintf("results/hmap_5.10_with_MrDAG/bias_hm_weight_%s.png",weight_to_plot), plot = combined_hm_bias, width = 12, height = 4)
  ggplot2::ggsave(sprintf("results/hmap_5.10_with_MrDAG/sd_hm_weight_%s.png",weight_to_plot), plot = combined_hm_sd, width = 12, height = 4)
  
  return(list(mean_abs_C_entry = mean_abs_C_entry,
              # report mean
              mean_abs_bias_ivw = mean_abs_bias_ivw, mean_sd_ivw = mean_sd_ivw,
              mean_abs_bias_adivw = mean_abs_bias_adivw, mean_sd_adivw = mean_sd_adivw,
              mean_abs_bias_naive = mean_abs_bias_naive, mean_sd_naive = mean_sd_naive,
              mean_abs_bias = mean_abs_bias, mean_sd = mean_sd,
              mean_abs_bias_r = mean_abs_bias_r, mean_sd_r = mean_sd_r,
              mean_abs_bias_MrDAG = mean_abs_bias_MrDAG, mean_sd_MrDAG = mean_sd_MrDAG,
              # report median
              median_abs_bias_ivw = median_abs_bias_ivw, median_sd_ivw = median_sd_ivw,
              median_abs_bias_adivw = median_abs_bias_adivw, median_sd_adivw = median_sd_adivw,
              median_abs_bias_naive = median_abs_bias_naive, median_sd_naive = median_sd_naive,
              median_abs_bias = median_abs_bias, median_sd = median_sd,
              median_abs_bias_r = median_abs_bias_r, median_sd_r = median_sd_r,
              median_abs_bias_MrDAG = median_abs_bias_MrDAG, median_sd_MrDAG = median_sd_MrDAG,
              # report quantile
              quantile_abs_bias_ivw = quantile_abs_bias_ivw, quantile_sd_ivw = quantile_sd_ivw,
              quantile_abs_bias_adivw = quantile_abs_bias_adivw, quantile_sd_adivw = quantile_sd_adivw,
              quantile_abs_bias_naive = quantile_abs_bias_naive, quantile_sd_naive = quantile_sd_naive,
              quantile_abs_bias = quantile_abs_bias, quantile_sd = quantile_sd,
              quantile_abs_bias_r = quantile_abs_bias_r, quantile_sd_r = quantile_sd_r,
              quantile_abs_bias_MrDAG = quantile_abs_bias_MrDAG, quantile_sd_MrDAG = quantile_sd_MrDAG
              )
         )
}

simulate_result = simulate_result_big_boxplot
plot_heatmap(parameters_list = simulate_result$parameters_list,
             result_naive_MRrr_list = simulate_result$result_AB_list,
             result_MRrr_list = simulate_result$result_AB_d_list,
             result_r_MRrr_list = simulate_result$result_AB_d_r_list,
             result_ivw_list = simulate_result$result_C_ivw_list,
             result_adivw_list = simulate_result$result_C_adivw_list,
             result_MrDAG_list = simulate_result$result_MrDAG_list,
             weight_to_plot="2.5")



#### run simulation - matrix error ####
regularization_rate_list = c(4.539335e-11, 3.341317e-12, 1.110624e-15)
run_simulation_matrix_error <- function(regularization_rate_list){
  sample_weight_list = c(1.75, 1, 0.4)
  eloop = 1000
  
  result_AB_list = result_AB_d_list = result_AB_d_r_list = 
    result_C_ivw_list = result_C_adivw_list = result_MrDAG_list = # result_C_GRAPPLE_list =
    Frobenius_AB_list = Frobenius_AB_d_list = Frobenius_AB_d_r_list =
    Frobenius_C_ivw_list = Frobenius_C_adivw_list = Frobenius_MrDAG_list = # Frobenius_C_GRAPPLE_list =
    Nuclear_AB_list = Nuclear_AB_d_list = Nuclear_AB_d_r_list =
    Nuclear_C_ivw_list = Nuclear_C_adivw_list = Nuclear_MrDAG_list = # Nuclear_C_GRAPPLE_list =
    list("1.75" = NA, "1" = NA, "0.4" = NA)
  
  iv_strength_list = c()
  parameters_list = c()
  
  for (weight_index in seq_along(sample_weight_list)) {
    me_weight = sample_weight_list[weight_index]
    parameters = .get_parameters(me_weight, px = 9, r_RR = 5)
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
    
    # storing entry-wise bias for each eloop
    bias_AB_matrix = bias_AB_d_matrix = bias_AB_d_r_matrix = 
      bias_C_ivw_matrix = bias_C_adivw_matrix = bias_C_MrDAG_matrix = # bias_C_GRAPPLE_matrix =
      matrix(NA, nrow = py*px, ncol = eloop)
    
    # storing Frobenius norm and Nuclear norm for each eloop
    Frobenius_AB_vec = Frobenius_AB_d_vec = Frobenius_AB_d_r_vec =
      Frobenius_C_ivw_vec = Frobenius_C_adivw_vec = Frobenius_MrDAG_vec = # Frobenius_C_GRAPPLE_vec =
      Nuclear_AB_vec = Nuclear_AB_d_vec = Nuclear_AB_d_r_vec =
      Nuclear_C_ivw_vec = Nuclear_C_adivw_vec = Nuclear_MrDAG_vec = # Nuclear_C_GRAPPLE_vec =
      numeric(eloop)
    
    B_star = t(B) # in order to compute the norm
    B_d_star = t(B_d)
    for (i in 1:eloop) {
      simulation_result = .simulation(parameters, 
                                      regularization_rate=regularization_rate_list[weight_index])
      # Naive MRrr
      A_hat = simulation_result[[1]]
      B_hat = simulation_result[[2]]
      B_hat_star = t(simulation_result[[2]])
      AB_hat = simulation_result[[3]]
      
      bias_AB_vectorized = as.vector(AB_hat - C)
      # bias_AB_vectorized = as.vector(AB_hat - C_r)
      bias_AB_matrix[,i] = bias_AB_vectorized
      Frobenius_AB_vec[i] = norm(AB_hat - C, type = "F")
      Nuclear_AB_vec[i] = .nuclear_norm(AB_hat - C)
      
      # MRrr
      A_d_hat = simulation_result[[4]]
      B_d_hat = simulation_result[[5]]
      B_d_hat_star = t(simulation_result[[5]])
      AB_d_hat = simulation_result[[6]]
      
      # bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
      bias_AB_d_vectorized = as.vector(AB_d_hat - C)
      bias_AB_d_matrix[,i] = bias_AB_d_vectorized
      Frobenius_AB_d_vec[i] = norm(AB_d_hat - C, type = "F")
      Nuclear_AB_d_vec[i] = .nuclear_norm(AB_d_hat - C)
      
      # MRrr with regularization
      A_d_r_hat = simulation_result[[7]]
      B_d_r_hat = simulation_result[[8]]
      B_d_r_hat_star = t(simulation_result[[8]])
      AB_d_r_hat = simulation_result[[9]]
      
      # bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
      bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C)
      bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
      Frobenius_AB_d_r_vec[i] = norm(AB_d_r_hat - C, type = "F")
      Nuclear_AB_d_r_vec[i] = .nuclear_norm(AB_d_r_hat - C)
      
      # IVW
      C_ivw = simulation_result[[12]]
      bias_C_ivw_vectorized = as.vector(C_ivw - C)
      bias_C_ivw_matrix[,i] = bias_C_ivw_vectorized
      Frobenius_C_ivw_vec[i] = norm(C_ivw - C, type = "F")
      Nuclear_C_ivw_vec[i] = .nuclear_norm(C_ivw - C)
      
      # adIVW
      C_adivw = simulation_result[[13]]
      bias_C_adivw_vectorized = as.vector(C_adivw - C)
      bias_C_adivw_matrix[,i] = bias_C_adivw_vectorized
      Frobenius_C_adivw_vec[i] = norm(C_adivw - C, type = "F")
      Nuclear_C_adivw_vec[i] = .nuclear_norm(C_adivw - C)
      
      # MrDAG
      C_MrDAG = simulation_result[[14]]
      bias_C_MrDAG_vectorized = as.vector(C_MrDAG - C)
      bias_C_MrDAG_matrix[,i] = bias_C_MrDAG_vectorized
      Frobenius_MrDAG_vec[i] = norm(C_MrDAG - C, type = "F")
      Nuclear_MrDAG_vec[i] = .nuclear_norm(C_MrDAG - C)
      
      # C_GRAPPLE = simulation_result[[15]]
      # bias_C_GRAPPLE_vectorized = as.vector(C_GRAPPLE - C)
      # bias_C_GRAPPLE_matrix[,i] = bias_C_GRAPPLE_vectorized
      
      # print "eloop" and eloop number in a line
      print(paste("eloop:", i, "weight:", me_weight))
    }
    char <- as.character(me_weight)
    # entry-wise bias across all eloops
    result_AB_list[[char]] = bias_AB_matrix
    result_AB_d_list[[char]] = bias_AB_d_matrix
    result_AB_d_r_list[[char]] = bias_AB_d_r_matrix
    result_C_ivw_list[[char]] = bias_C_ivw_matrix
    result_C_adivw_list[[char]] = bias_C_adivw_matrix
    result_MrDAG_list[[char]] = bias_C_MrDAG_matrix
    # result_C_GRAPPLE_list[[char]] = bias_C_GRAPPLE_matrix
    
    # Frobenius and Nuclear norm across all eloops
    Frobenius_AB_list[[char]] = Frobenius_AB_vec
    Frobenius_AB_d_list[[char]] = Frobenius_AB_d_vec
    Frobenius_AB_d_r_list[[char]] = Frobenius_AB_d_r_vec
    Frobenius_C_ivw_list[[char]] = Frobenius_C_ivw_vec
    Frobenius_C_adivw_list[[char]] = Frobenius_C_adivw_vec
    Frobenius_MrDAG_list[[char]] = Frobenius_MrDAG_vec
    Nuclear_AB_list[[char]] = Nuclear_AB_vec
    Nuclear_AB_d_list[[char]] = Nuclear_AB_d_vec
    Nuclear_AB_d_r_list[[char]] = Nuclear_AB_d_r_vec
    Nuclear_C_ivw_list[[char]] = Nuclear_C_ivw_vec
    Nuclear_C_adivw_list[[char]] = Nuclear_C_adivw_vec
    Nuclear_MrDAG_list[[char]] = Nuclear_MrDAG_vec
    
  }
  return(list(result_AB_list=result_AB_list, 
              result_AB_d_list=result_AB_d_list, 
              result_AB_d_r_list=result_AB_d_r_list,
              iv_strength_list=iv_strength_list, parameters_list=parameters_list,
              result_C_ivw_list=result_C_ivw_list, 
              result_C_adivw_list=result_C_adivw_list,
              result_MrDAG_list=result_MrDAG_list,
              # result_C_GRAPPLE_list=result_C_GRAPPLE_list,
              Frobenius_AB_list=Frobenius_AB_list,
              Frobenius_AB_d_list=Frobenius_AB_d_list,
              Frobenius_AB_d_r_list=Frobenius_AB_d_r_list,
              Frobenius_C_ivw_list=Frobenius_C_ivw_list,
              Frobenius_C_adivw_list=Frobenius_C_adivw_list,
              Frobenius_MrDAG_list=Frobenius_MrDAG_list,
              Nuclear_AB_list=Nuclear_AB_list,
              Nuclear_AB_d_list=Nuclear_AB_d_list,
              Nuclear_AB_d_r_list=Nuclear_AB_d_r_list,
              Nuclear_C_ivw_list=Nuclear_C_ivw_list,
              Nuclear_C_adivw_list=Nuclear_C_adivw_list,
              Nuclear_MrDAG_list=Nuclear_MrDAG_list
              ))
}


start_time <- Sys.time()

simulate_result_matrix_error = run_simulation_matrix_error(regularization_rate_list = regularization_rate_list) 
save(simulate_result_matrix_error, file = "results/simulate_result_250407.RData")
round(simulate_result_matrix_error$iv_strength_list, 3)

end_time <- Sys.time()
print(end_time - start_time)


load("results/simulate_result_250407.RData")


# The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
plot_heatmap_matrix_error <- function(simulate_result_matrix_error, weight_to_plot)
{
  parameters_list = simulate_result_matrix_error$parameters_list
  result_naive_MRrr_list = simulate_result_matrix_error$result_AB_list
  result_MRrr_list = simulate_result_matrix_error$result_AB_d_list
  result_r_MRrr_list = simulate_result_matrix_error$result_AB_d_r_list
  result_ivw_list = simulate_result_matrix_error$result_C_ivw_list
  result_adivw_list = simulate_result_matrix_error$result_C_adivw_list
  result_MrDAG_list = simulate_result_matrix_error$result_MrDAG_list
  Frobenius_AB_list = simulate_result_matrix_error$Frobenius_AB_list
  Frobenius_AB_d_list = simulate_result_matrix_error$Frobenius_AB_d_list
  Frobenius_AB_d_r_list = simulate_result_matrix_error$Frobenius_AB_d_r_list
  Frobenius_C_ivw_list = simulate_result_matrix_error$Frobenius_C_ivw_list
  Frobenius_C_adivw_list = simulate_result_matrix_error$Frobenius_C_adivw_list
  Frobenius_MrDAG_list = simulate_result_matrix_error$Frobenius_MrDAG_list
  Nuclear_AB_list = simulate_result_matrix_error$Nuclear_AB_list
  Nuclear_AB_d_list = simulate_result_matrix_error$Nuclear_AB_d_list
  Nuclear_AB_d_r_list = simulate_result_matrix_error$Nuclear_AB_d_r_list
  Nuclear_C_ivw_list = simulate_result_matrix_error$Nuclear_C_ivw_list
  Nuclear_C_adivw_list = simulate_result_matrix_error$Nuclear_C_adivw_list
  Nuclear_MrDAG_list = simulate_result_matrix_error$Nuclear_MrDAG_list
  
  px <- 9
  py <- 4
  # naive MR-rr
  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }
  
  # MR-rr
  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
  }
  
  # MR-rr with spectral regularization
  abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d_r <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
  }
  
  # IVW
  abs_mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_ivw <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_ivw <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_ivw[i] <- abs(as.numeric(mean(result_ivw_list[[weight_to_plot]][i,])))
    mean_entry_bias_ivw[i] <- mean(result_ivw_list[[weight_to_plot]][i,])
    sd_entry_bias_ivw[i] <- sd(result_ivw_list[[weight_to_plot]][i,])
  }
  
  # adIVW
  abs_mean_entry_bias_adivw <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_adivw <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_adivw <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_adivw[i] <- abs(as.numeric(mean(result_adivw_list[[weight_to_plot]][i,])))
    mean_entry_bias_adivw[i] <- mean(result_adivw_list[[weight_to_plot]][i,])
    sd_entry_bias_adivw[i] <- sd(result_adivw_list[[weight_to_plot]][i,])
  }
  
  # MrDAG
  abs_mean_entry_bias_MrDAG <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_MrDAG <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_MrDAG <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_MrDAG[i] <- abs(as.numeric(mean(result_MrDAG_list[[weight_to_plot]][i,])))
    mean_entry_bias_MrDAG[i] <- mean(result_MrDAG_list[[weight_to_plot]][i,])
    sd_entry_bias_MrDAG[i] <- sd(result_MrDAG_list[[weight_to_plot]][i,])
  }
  
  
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
  
  # dRR_r
  avg_bias_d_r <- matrix(abs_mean_entry_bias_d_r, nrow = px, ncol = py)
  avg_bias_d_r <- t(avg_bias_d_r)
  avg_bias_d_r <- as.data.frame(avg_bias_d_r)
  colnames(avg_bias_d_r) <- 1:px
  rownames(avg_bias_d_r) <- 1:py
  
  # IVW
  avg_bias_ivw <- matrix(abs_mean_entry_bias_ivw, nrow = px, ncol = py)
  avg_bias_ivw <- t(avg_bias_ivw)
  avg_bias_ivw <- as.data.frame(avg_bias_ivw)
  colnames(avg_bias_ivw) <- 1:px
  rownames(avg_bias_ivw) <- 1:py
  
  # adIVW
  avg_bias_adivw <- matrix(abs_mean_entry_bias_adivw, nrow = px, ncol = py)
  avg_bias_adivw <- t(avg_bias_adivw)
  avg_bias_adivw <- as.data.frame(avg_bias_adivw)
  colnames(avg_bias_adivw) <- 1:px
  rownames(avg_bias_adivw) <- 1:py
  
  # MrDAG
  avg_bias_MrDAG <- matrix(abs_mean_entry_bias_MrDAG, nrow = px, ncol = py)
  avg_bias_MrDAG <- t(avg_bias_MrDAG)
  avg_bias_MrDAG <- as.data.frame(avg_bias_MrDAG)
  colnames(avg_bias_MrDAG) <- 1:px
  rownames(avg_bias_MrDAG) <- 1:py
  
  
  # color breaks
  min_value <- min(min(avg_bias), min(avg_bias_d), min(avg_bias_d_r))
  max_value <- max(max(avg_bias),max(avg_bias_d), max(avg_bias_d_r))
  
  # min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
  # max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))
  
  breaks <- seq(min_value, max_value, length.out = 101)
  
  hm_bias_1 = pheatmap::pheatmap(avg_bias, breaks = breaks,
                                 main = "Naive MR-rr estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_2 = pheatmap::pheatmap(avg_bias_d, breaks = breaks,
                                 main = "MR-rr estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_3 = pheatmap::pheatmap(avg_bias_d_r, breaks = breaks,
                                 main = "MR-rr estimator with regularization", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_4 = pheatmap::pheatmap(avg_bias_ivw, breaks = breaks,
                                 main = "IVW estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_5 = pheatmap::pheatmap(avg_bias_adivw, breaks = breaks,
                                 main = "adIVW estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  hm_bias_6 = pheatmap::pheatmap(avg_bias_MrDAG, breaks = breaks,
                                 main = "MrDAG estimator", fontsize = 8,
                                 cluster_rows = FALSE, cluster_cols = FALSE,
                                 color = colorRampPalette(c("white", "red"))(100),
                                 height = 10,
                                 width = 8, silent = TRUE)$gtable
  
  
  combined_hm_bias = gridExtra::grid.arrange(
    grobs = list(hm_bias_4, hm_bias_5, hm_bias_1, hm_bias_2, hm_bias_3, hm_bias_6),
    ncol = 6,
    top = grid::textGrob(sprintf("Absolute average bias by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
  )
  combined_hm_bias
  
  
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
  
  # dRR_r
  avg_sd_d_r <- sd_entry_bias_d_r
  avg_sd_d_r <- matrix(avg_sd_d_r, nrow = px, ncol = py)
  avg_sd_d_r <- t(avg_sd_d_r)
  avg_sd_d_r <- as.data.frame(avg_sd_d_r)
  colnames(avg_sd_d_r) <- 1:px
  rownames(avg_sd_d_r) <- 1:py
  
  # IVW
  avg_sd_ivw <- sd_entry_bias_ivw
  avg_sd_ivw <- matrix(avg_sd_ivw, nrow = px, ncol = py)
  avg_sd_ivw <- t(avg_sd_ivw)
  avg_sd_ivw <- as.data.frame(avg_sd_ivw)
  colnames(avg_sd_ivw) <- 1:px
  rownames(avg_sd_ivw) <- 1:py
  
  # adIVW
  avg_sd_adivw <- sd_entry_bias_adivw
  avg_sd_adivw <- matrix(avg_sd_adivw, nrow = px, ncol = py)
  avg_sd_adivw <- t(avg_sd_adivw)
  avg_sd_adivw <- as.data.frame(avg_sd_adivw)
  colnames(avg_sd_adivw) <- 1:px
  rownames(avg_sd_adivw) <- 1:py
  
  # MrDAG
  avg_sd_MrDAG <- sd_entry_bias_MrDAG
  avg_sd_MrDAG <- matrix(avg_sd_MrDAG, nrow = px, ncol = py)
  avg_sd_MrDAG <- t(avg_sd_MrDAG)
  avg_sd_MrDAG <- as.data.frame(avg_sd_MrDAG)
  colnames(avg_sd_MrDAG) <- 1:px
  rownames(avg_sd_MrDAG) <- 1:py
  
  min_value_sd <- min(min(avg_sd), min(avg_sd_d), min(avg_sd_d_r))
  max_value_Sd <- max(max(avg_sd),max(avg_sd_d), max(avg_sd_d_r))
  breaks_sd <- seq(min_value_sd, max_value_Sd, length.out = 101)
  
  
  hm_sd_1 = pheatmap::pheatmap(avg_sd, breaks = breaks_sd,
                               main = "Naive MR-rr estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_2 = pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
                               main = "MR-rr estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_3 = pheatmap::pheatmap(avg_sd_d_r, breaks = breaks_sd,
                               main = "MR-rr estimator with regularization", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_4 = pheatmap::pheatmap(avg_sd_ivw, breaks = breaks_sd,
                               main = "IVW estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_5 = pheatmap::pheatmap(avg_sd_adivw, breaks = breaks_sd,
                               main = "adIVW estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  hm_sd_6 = pheatmap::pheatmap(avg_sd_MrDAG, breaks = breaks_sd,
                               main = "MrDAG estimator", fontsize = 8,
                               cluster_rows = FALSE, cluster_cols = FALSE,
                               color = colorRampPalette(c("white", "red"))(100),
                               height = 10,
                               width = 8, silent = TRUE)$gtable
  
  
  combined_hm_sd = gridExtra::grid.arrange(
    grobs = list(hm_sd_4, hm_sd_5, hm_sd_1, hm_sd_2, hm_sd_3, hm_sd_6),
    ncol = 6,
    top = grid::textGrob(sprintf("SD by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
  )
  combined_hm_sd
  
  # bias scale compares to true C_r and sd
  weight_index = match(weight_to_plot, c("1.75", "1", "0.4"))
  parameters = parameters_list[[weight_index]]
  C_r = parameters$C_r
  
  F_norm_C = norm(as.matrix(C_r), type = "F")
  N_norm_C = .nuclear_norm(as.matrix(C_r))
  # mean
  mean_abs_C_entry =  round(mean(abs(as.matrix(C_r))), 3)
  mean_abs_bias_naive = round(mean(abs(as.matrix(avg_bias))), 3)
  mean_abs_bias = round(mean(abs(as.matrix(avg_bias_d))), 3)
  mean_abs_bias_r = round(mean(abs(as.matrix(avg_bias_d_r))), 3)
  mean_abs_bias_ivw = round(mean(abs(as.matrix(avg_bias_ivw))), 3)
  mean_abs_bias_adivw = round(mean(abs(as.matrix(avg_bias_adivw))), 3)
  mean_abs_bias_MrDAG = round(mean(abs(as.matrix(avg_bias_MrDAG))), 3)
  
  mean_sd_naive = round(mean(as.matrix(avg_sd)), 3)
  mean_sd = round(mean(as.matrix(avg_sd_d)), 3)
  mean_sd_r = round(mean(as.matrix(avg_sd_d_r)), 3)
  mean_sd_ivw = round(mean(as.matrix(avg_sd_ivw)), 3)
  mean_sd_adivw = round(mean(as.matrix(avg_sd_adivw)), 3)
  mean_sd_MrDAG = round(mean(as.matrix(avg_sd_MrDAG)), 3)
  
  # median
  median_abs_bias_naive = round(median(abs(as.matrix(avg_bias))), 3)
  median_abs_bias = round(median(abs(as.matrix(avg_bias_d))), 3)
  median_abs_bias_r = round(median(abs(as.matrix(avg_bias_d_r))), 3)
  median_abs_bias_ivw = round(median(abs(as.matrix(avg_bias_ivw))), 3)
  median_abs_bias_adivw = round(median(abs(as.matrix(avg_bias_adivw))), 3)
  median_abs_bias_MrDAG = round(median(abs(as.matrix(avg_bias_MrDAG))), 3)
  
  median_sd_naive = round(median(as.matrix(avg_sd)), 3)
  median_sd = round(median(as.matrix(avg_sd_d)), 3)
  median_sd_r = round(median(as.matrix(avg_sd_d_r)), 3)
  median_sd_ivw = round(median(as.matrix(avg_sd_ivw)), 3)
  median_sd_adivw = round(median(as.matrix(avg_sd_adivw)), 3)
  median_sd_MrDAG = round(median(as.matrix(avg_sd_MrDAG)), 3)
  
  # quantile
  quantile_abs_bias_naive = round(quantile(abs(as.matrix(avg_bias)), c(0.25, 0.75)), 3)
  quantile_abs_bias = round(quantile(abs(as.matrix(avg_bias_d)), c(0.25, 0.75)), 3)
  quantile_abs_bias_r = round(quantile(abs(as.matrix(avg_bias_d_r)), c(0.25, 0.75)), 3)
  quantile_abs_bias_ivw = round(quantile(abs(as.matrix(avg_bias_ivw)), c(0.25, 0.75)), 3)
  quantile_abs_bias_adivw = round(quantile(abs(as.matrix(avg_bias_adivw)), c(0.25, 0.75)), 3)
  quantile_abs_bias_MrDAG = round(quantile(abs(as.matrix(avg_bias_MrDAG)), c(0.25, 0.75)), 3)
  
  quantile_sd_naive = round(quantile(as.matrix(avg_sd), c(0.25, 0.75)), 3)
  quantile_sd = round(quantile(as.matrix(avg_sd_d), c(0.25, 0.75)), 3)
  quantile_sd_r = round(quantile(as.matrix(avg_sd_d_r), c(0.25, 0.75)), 3)
  quantile_sd_ivw = round(quantile(as.matrix(avg_sd_ivw), c(0.25, 0.75)), 3)
  quantile_sd_adivw = round(quantile(as.matrix(avg_sd_adivw), c(0.25, 0.75)), 3)
  quantile_sd_MrDAG = round(quantile(as.matrix(avg_sd_MrDAG), c(0.25, 0.75)), 3)
  
  ggplot2::ggsave(sprintf("results/hmap_4.8/bias_hm_weight_%s.png",weight_to_plot), plot = combined_hm_bias, width = 12, height = 4)
  ggplot2::ggsave(sprintf("results/hmap_4.8/sd_hm_weight_%s.png",weight_to_plot), plot = combined_hm_sd, width = 12, height = 4)
  
  # report Frobenius norm and Nuclear norm
  Average_F_norm_naive_MR = round(mean(Frobenius_AB_list[[weight_to_plot]])/F_norm_C,2)
  Average_F_norm_MR = round(mean(Frobenius_AB_d_list[[weight_to_plot]])/F_norm_C,2)
  Average_F_norm_MR_r = round(mean(Frobenius_AB_d_r_list[[weight_to_plot]])/F_norm_C,2)
  Average_F_norm_ivw = round(mean(Frobenius_C_ivw_list[[weight_to_plot]])/F_norm_C,2)
  Average_F_norm_adivw = round(mean(Frobenius_C_adivw_list[[weight_to_plot]])/F_norm_C,2)
  Average_F_norm_MrDAG = round(mean(Frobenius_MrDAG_list[[weight_to_plot]])/F_norm_C,2)
  
  Average_N_norm_naive_MR = round(mean(Nuclear_AB_list[[weight_to_plot]])/N_norm_C,2)
  Average_N_norm_MR = round(mean(Nuclear_AB_d_list[[weight_to_plot]])/N_norm_C,2)
  Average_N_norm_MR_r = round(mean(Nuclear_AB_d_r_list[[weight_to_plot]])/N_norm_C,2)
  Average_N_norm_ivw = round(mean(Nuclear_C_ivw_list[[weight_to_plot]])/N_norm_C,2)
  Average_N_norm_adivw = round(mean(Nuclear_C_adivw_list[[weight_to_plot]])/N_norm_C,2)
  Average_N_norm_MrDAG = round(mean(Nuclear_MrDAG_list[[weight_to_plot]])/N_norm_C,2)
  
  return(list(mean_abs_C_entry = mean_abs_C_entry,
              # report mean
              mean_abs_bias_ivw = mean_abs_bias_ivw, mean_sd_ivw = mean_sd_ivw,
              mean_abs_bias_adivw = mean_abs_bias_adivw, mean_sd_adivw = mean_sd_adivw,
              mean_abs_bias_naive = mean_abs_bias_naive, mean_sd_naive = mean_sd_naive,
              mean_abs_bias = mean_abs_bias, mean_sd = mean_sd,
              mean_abs_bias_r = mean_abs_bias_r, mean_sd_r = mean_sd_r,
              mean_abs_bias_MrDAG = mean_abs_bias_MrDAG, mean_sd_MrDAG = mean_sd_MrDAG,
              # report median
              median_abs_bias_ivw = median_abs_bias_ivw, median_sd_ivw = median_sd_ivw,
              median_abs_bias_adivw = median_abs_bias_adivw, median_sd_adivw = median_sd_adivw,
              median_abs_bias_naive = median_abs_bias_naive, median_sd_naive = median_sd_naive,
              median_abs_bias = median_abs_bias, median_sd = median_sd,
              median_abs_bias_r = median_abs_bias_r, median_sd_r = median_sd_r,
              median_abs_bias_MrDAG = median_abs_bias_MrDAG, median_sd_MrDAG = median_sd_MrDAG,
              # report quantile
              quantile_abs_bias_ivw = quantile_abs_bias_ivw, quantile_sd_ivw = quantile_sd_ivw,
              quantile_abs_bias_adivw = quantile_abs_bias_adivw, quantile_sd_adivw = quantile_sd_adivw,
              quantile_abs_bias_naive = quantile_abs_bias_naive, quantile_sd_naive = quantile_sd_naive,
              quantile_abs_bias = quantile_abs_bias, quantile_sd = quantile_sd,
              quantile_abs_bias_r = quantile_abs_bias_r, quantile_sd_r = quantile_sd_r,
              quantile_abs_bias_MrDAG = quantile_abs_bias_MrDAG, quantile_sd_MrDAG = quantile_sd_MrDAG,
              # report Frobenius norm and Nuclear norm
              Average_F_norm_ivw = Average_F_norm_ivw,
              Average_F_norm_adivw = Average_F_norm_adivw,
              Average_F_norm_naive_MR = Average_F_norm_naive_MR,
              Average_F_norm_MR = Average_F_norm_MR,
              Average_F_norm_MR_r = Average_F_norm_MR_r,
              Average_F_norm_MrDAG = Average_F_norm_MrDAG,
              Average_N_norm_ivw = Average_N_norm_ivw,
              Average_N_norm_adivw = Average_N_norm_adivw,
              Average_N_norm_naive_MR = Average_N_norm_naive_MR,
              Average_N_norm_MR = Average_N_norm_MR,
              Average_N_norm_MR_r = Average_N_norm_MR_r,
              Average_N_norm_MrDAG = Average_N_norm_MrDAG
              )
         )
}


plot_heatmap_matrix_error(simulate_result_matrix_error = simulate_result_matrix_error,
                          weight_to_plot="0.4")


#### run simulation - plot big box plot ####
# 5.936180e-11 3.457282e-12 1.229468e-18
regularization_rate_list = c(5.936180e-11, 3.457282e-12, 1.229468e-18)
# The simulation function for the naive MR-rr estimator and the MR-rr estimator (with spectral regularization)
run_simulation_big_boxplot <- function(regularization_rate_list){
  sample_weight_list = c(2.5, 2, 1)
  eloop = 200
  
  result_AB_list = result_AB_d_list = result_AB_d_r_list = 
    result_C_ivw_list = result_C_adivw_list = result_MrDAG_list = # result_C_GRAPPLE_list =
    AB_list = AB_d_list = AB_d_r_list = 
    C_ivw_list = C_adivw_list = MrDAG_list =
    list("2.5" = NA, "2" = NA, "1" = NA)
  
  iv_strength_list = c()
  parameters_list = c()
  
  for (weight_index in seq_along(sample_weight_list)) {
    me_weight = sample_weight_list[weight_index]
    parameters = .get_parameters(me_weight, px = 9, r_RR = 2)
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
    B_star = t(B) # in order to compute the norm
    B_d_star = t(B_d)
    for (i in 1:eloop) {
      simulation_result = .simulation(parameters, 
                                      regularization_rate=regularization_rate_list[weight_index])
      A_hat = simulation_result[[1]]
      B_hat = simulation_result[[2]]
      B_hat_star = t(simulation_result[[2]])
      AB_hat = simulation_result[[3]]
      
      bias_AB_vectorized = as.vector(AB_hat - C)
      # bias_AB_vectorized = as.vector(AB_hat - C_r)
      bias_AB_matrix[,i] = bias_AB_vectorized
      AB_matrix[,i] = as.vector(AB_hat)
      
      A_d_hat = simulation_result[[4]]
      B_d_hat = simulation_result[[5]]
      B_d_hat_star = t(simulation_result[[5]])
      AB_d_hat = simulation_result[[6]]
      
      # bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
      bias_AB_d_vectorized = as.vector(AB_d_hat - C)
      bias_AB_d_matrix[,i] = bias_AB_d_vectorized
      AB_d_matrix[,i] = as.vector(AB_d_hat)
      
      A_d_r_hat = simulation_result[[7]]
      B_d_r_hat = simulation_result[[8]]
      B_d_r_hat_star = t(simulation_result[[8]])
      AB_d_r_hat = simulation_result[[9]]
      
      # bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
      bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C)
      bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
      AB_d_r_matrix[,i] = as.vector(AB_d_r_hat)
      
      C_ivw = simulation_result[[12]]
      bias_C_ivw_vectorized = as.vector(C_ivw - C)
      bias_C_ivw_matrix[,i] = bias_C_ivw_vectorized
      C_ivw_matrix[,i] = as.vector(C_ivw)
      
      C_adivw = simulation_result[[13]]
      bias_C_adivw_vectorized = as.vector(C_adivw - C)
      bias_C_adivw_matrix[,i] = bias_C_adivw_vectorized
      C_adivw_matrix[,i] = as.vector(C_adivw)
      
      C_MrDAG = simulation_result[[14]]
      bias_C_MrDAG_vectorized = as.vector(C_MrDAG - C)
      bias_C_MrDAG_matrix[,i] = bias_C_MrDAG_vectorized
      C_MrDAG_matrix[,i] = as.vector(C_MrDAG)
      
      # C_GRAPPLE = simulation_result[[15]]
      # bias_C_GRAPPLE_vectorized = as.vector(C_GRAPPLE - C)
      # bias_C_GRAPPLE_matrix[,i] = bias_C_GRAPPLE_vectorized
      
      # print "eloop" and eloop number in a line
      print(paste("eloop:", i, "weight:", me_weight))
    }
    char <- as.character(me_weight)
    # store estimator bias
    result_AB_list[[char]] = bias_AB_matrix
    result_AB_d_list[[char]] = bias_AB_d_matrix
    result_AB_d_r_list[[char]] = bias_AB_d_r_matrix
    result_C_ivw_list[[char]] = bias_C_ivw_matrix
    result_C_adivw_list[[char]] = bias_C_adivw_matrix
    result_MrDAG_list[[char]] = bias_C_MrDAG_matrix
    # result_C_GRAPPLE_list[[char]] = bias_C_GRAPPLE_matrix
    # store estimator value
    AB_list[[char]] = AB_matrix
    AB_d_list[[char]] = AB_d_matrix
    AB_d_r_list[[char]] = AB_d_r_matrix
    C_ivw_list[[char]] = C_ivw_matrix
    C_adivw_list[[char]] = C_adivw_matrix
    MrDAG_list[[char]] = C_MrDAG_matrix
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
              MrDAG_list=MrDAG_list
  ))
}
simulate_result_big_boxplot = run_simulation_big_boxplot(regularization_rate_list = regularization_rate_list)
save(simulate_result_big_boxplot, file = "results/simulate_result_big_boxplot_250510.RData")
round(simulate_result_big_boxplot$iv_strength_list, 3)

load("results/simulate_result_big_boxplot_250510.RData")


sim_result <- simulate_result_big_boxplot

# # test
# weight = "1"
# width = 16
# height = 8
# filename = NULL

plot_simulation_boxplot <- function(sim_result, weight = "1",
                                    filename = NULL, width = 16, height = 8) {
  estimator_names <- c("AB_list", "AB_d_list", "AB_d_r_list", "C_ivw_list", "C_adivw_list", "MrDAG_list")
  
  # 整理数据
  df_long <- map_dfr(estimator_names, function(est) {
    matrix_data <- sim_result[[est]][[weight]]
    n_sim <- ncol(matrix_data)
    n_pos <- nrow(matrix_data)
    
    tibble(
      Estimator = est,
      Pos = rep(1:n_pos, each = n_sim),
      Row = rep((0:3), each = 9 * n_sim),
      Col = rep((0:8), times = 4 * n_sim),
      Value = as.vector(t(matrix_data))
    )
  })
  
  df_long$Estimator <- factor(df_long$Estimator,
                              levels = estimator_names,
                              labels = c("Naive MR-rr", "MR-rr", "regularized MR-rr", "IVW", "adIVW", "MrDAG"))
  
  # 画图
  p <- ggplot(df_long, aes(x = Estimator, y = Value, fill = Estimator)) +
    geom_boxplot(outlier.size = 0.3) +
    facet_grid(rows = vars(Row), cols = vars(Col)) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text = element_text(size = 8)) +
    ggtitle(paste("Simulation Boxplots by Estimator (Weight =", weight, ")")) +
    ylab("Estimated Value") +
    xlab("")
  
  if (!is.null(filename)) {
    ggsave(filename, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to ", filename)
    print(p)
  } else {
    print(p)
  }
}

# 用法示例：
# 只显示图
# plot_simulation_boxplot(simulate_result_big_boxplot, weight = "2.5")

# 保存图
# plot_simulation_boxplot(simulate_result_big_boxplot, weight = "2.5", filename = "boxplot_weight2.5.png")

plot_simulation_boxplot(simulate_result_big_boxplot, weight = "1", filename = "lessvar_big_boxplot_1.png")

