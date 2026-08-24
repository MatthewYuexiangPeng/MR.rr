rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
library(matrixStats)
library(MASS)
library(ggplot2)
library(tidyverse)
library(patchwork)

#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
data("bmi.cad")
load('data/multivariate_data_medium.rda')
source("scripts/RRR.R")
source("scripts/debiased_RRR.R")


set.seed(123)


sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
}


get_parameters = function(sample_weight, r_RR = 3, px,  r_approx = 3){
  # r_RR is the rank chosed bu user when performing RRR
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 2000 # pz can be changed to any number
  py = 8
  gamma_j = as.matrix(tmp$data[,paste0('gamma_exp',1:9)])
  var_Z = 2 * bmi.cad$eaf.outcome * (1 - bmi.cad$eaf.outcome) # each Z is sum of two alleles, so var(Z) = var(Z^1+Z^2) = 2*var(Z^1), where Z1, Z2 ~ Binomial(eaf.outcome)
  var_Z_raw = var_Z
  
  bmi_index = sample(1:1119, pz, replace = TRUE)
  var_Z = var_Z[bmi_index]
  sqrt_var_Z = sqrt(var_Z)
  
  tmp_index_z <- sample(1:273, pz, replace = TRUE)
  gamma_j_sample <- gamma_j[tmp_index_z,]
  gamma_j_star_temp <- gamma_j_sample*sqrt_var_Z
  
  # sample cov
  VX_tilde <- cov(gamma_j_star_temp)
  
  corr_upper_triangle = cor(gamma_j_star_temp)[upper.tri(cor(gamma_j_star_temp), 
                                                             diag = FALSE)]
  diag_VX_tilde = diag(VX_tilde)
  diag_VX_tilde_sample = sample(diag_VX_tilde, px, replace = TRUE)
  sample_size = px*(px+1)/2
  VX_tilde_corr <- matrix(0, nrow = px, ncol = px)
  for (i in 1:px-1){
    for (j in (i+1):px) {
      VX_tilde_corr[i,j] = sample(corr_upper_triangle, 1)
    }
  }
  VX_tilde_corr = VX_tilde_corr + t(VX_tilde_corr) + diag(1, px)
  VX_tilde_upscale = sqrt_matrix(diag(diag_VX_tilde_sample)) %*% VX_tilde_corr %*% 
                     sqrt_matrix(diag(diag_VX_tilde_sample))
  VX_tilde_upscale <- VX_tilde_upscale*1/2
  
  
  # Sigma_Xj & Sigma_X ----------------------------------------------------------------
  sigma_gamma_j = as.matrix(tmp$data[,paste0('se_exp',1:9)])
  Sigma = tmp$cor.mat[1:9,1:9]
  sqrt_Sigma = sqrt_matrix(Sigma)
  pz_tmp = nrow(gamma_j)
  Sigma_Xj = lapply(1:pz_tmp, function(j) diag(sigma_gamma_j[j,]) %*% Sigma %*% diag(sigma_gamma_j[j,]))
  Sigma_Xj_sample = Sigma_Xj[tmp_index_z]
  Sigma_X_temp = lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d <- array(unlist(Sigma_X_temp), dim = c(9, 9, length(Sigma_X_temp)))
  Sigma_X <- apply(array_3d, c(1, 2), mean)
  
  Sigma_X_triangle = Sigma[upper.tri(Sigma, diag = FALSE)]
  diag_Sigma_X = diag(sigma_gamma_j)
  diag_Sigma_X_sample = sample(diag_Sigma_X, px, replace = TRUE)
  sample_size = px*(px+1)/2
  Sigma_X_corr <- matrix(0, nrow = px, ncol = px)
  for (i in 1:px-1){
    for (j in (i+1):px) {
      Sigma_X_corr[i,j] = sample(Sigma_X_triangle, 1)
    }
  }
  Sigma_X_corr = Sigma_X_corr + t(Sigma_X_corr) + diag(1, px)
  Sigma_X_upscale = sqrt_matrix(diag(diag_Sigma_X_sample)) %*% Sigma_X_corr %*% 
    sqrt_matrix(diag(diag_Sigma_X_sample))
  
  Sigma_X_upscale <- sample_weight * Sigma_X_upscale
  
  # Sigma_Y & weight.matrix -----------------------------------------------------------
  n_Y = median(bmi.cad$N.outcome) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
  var_Y = sample(bmi.cad$se.outcome^2 * bmi.cad$N.outcome * var_Z_raw, py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
  Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))
  
  Sigma_Y = sample_weight * Sigma_Y
  
  weight.matrix = solve(Sigma_Y)
  
  
  # SigmaXX, SigmaYX, SigmaXY, C ----------------------------------------------------------------
  # py * r matrix
  r =  min(px,py) # true rank of C, not specified reduced rank r
  U = matrix(rnorm(py*r), py, r)
  V = matrix(rnorm(px*r), r, px)
  # assume true C to be rank r=3 for now. C=C^(r)
  # use modified THM 2.1 to define true C^(r), (it is equivalent to define it as A_d*B_d)
  eignvalue_matrix = diag(c(sample(c(sqrt(0.3), sqrt(0.2), sqrt(0.2)), r_approx, replace = TRUE), 
                            rep(0, r-r_approx)))
  C <- U %*% eignvalue_matrix %*% V
  
  # SigmaXX, SigmaYX and SigmaXY
  SigmaXX <- Sigma_X_upscale + VX_tilde_upscale
  SigmaYX <- C %*% VX_tilde_upscale
  SigmaXY <- t(SigmaYX)
  sqrt_Gamma <- sqrt_matrix(weight.matrix)
  sqrt_Gamma_inv <- solve(sqrt_Gamma)
  
  # original RRR population level estimator
  M = sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX) %*% SigmaXY %*% sqrt_Gamma
  V = eigen(M)$vec[, 1:r_RR, drop = FALSE]
  A = sqrt_Gamma_inv %*% V
  B = t(V) %*% sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX)
  
  # debiased RRR population level estimator
  M_d = sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde_upscale) %*% SigmaXY %*% sqrt_Gamma
  V_d = eigen(M_d)$vec[, 1:r_RR, drop = FALSE]
  V_d = matrix(as.numeric(V_d),py,r_RR)
  A_d = sqrt_Gamma_inv %*% V_d
  B_d = t(V_d) %*% sqrt_Gamma %*% SigmaYX %*% solve(VX_tilde_upscale)
  
  
  # VY_tilde ----------------------------------------------------------------
  VY_tilde_upscale = C %*% VX_tilde_upscale %*% t(C)
  
  parameters = list(py = py, px = px, var_Z = var_Z, VX_tilde = VX_tilde_upscale, 
                    Sigma_X = Sigma_X_upscale, Sigma_Y = Sigma_Y, 
                    weight.matrix = weight.matrix, SigmaXX = SigmaXX, 
                    SigmaYX = SigmaYX, SigmaXY = SigmaXY, C = C, A = A, B = B, 
                    A_d = A_d, B_d = B_d, C_r = A_d %*% B_d, r_RR = r_RR, 
                    VY_tilde = VY_tilde_upscale)
  return(parameters)
}


simulation = function(parameters, n) {
  # n: number of samples (pz)
  # get parameters
  py = parameters$py
  px = parameters$px
  var_Z = parameters$var_Z
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  weight.matrix = parameters$weight.matrix
  SigmaXX = parameters$SigmaXX
  SigmaYX = parameters$SigmaYX
  SigmaXY = parameters$SigmaXY
  C = parameters$C
  r_RR = parameters$r_RR
  VY_tilde = parameters$VY_tilde
  
  # sample true effect gamma_j_star and Gamma_j_star
  gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  # sample x_j, y_j
  x_j = matrix(0, n, px)
  y_j = matrix(0, n, py)
  for (j in 1:n) {
    x_j[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
    y_j[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
  }
  
  # compute A_hat, B_hat
  result <- RRR(y_j, x_j, r=r_RR, weight.matrix)
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB # sample level estimator A_hat * B_hat
  
  # need to subtract the Sigma_xx by Sigma_x. Sigma_x can be estimated by the regression of exp ~ z (have been approximate in the get parameter function)
  # TODO: now we are using true Sigma_X to debias. should we use Sigma_X_hat?
  result_d <- debiased_RRR(y_j, x_j, r=r_RR, sqrt_Gamma = weight.matrix, Sigma_X = Sigma_X)
  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB
  
  return(list(A_hat, B_hat, AB_hat, A_d_hat, B_d_hat, AB_d_hat))
}


# run the simulation -----------------------------------------------------------
sample_weight_list = c(5, 2, 1, 0.5, 0.2)
eloop = 100
result_AB_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_AB_d_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)

# result of A and B separately ----
# result_A_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
# result_B_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
# result_A_d_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
# result_B_d_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
# ----


for (sample_weight in sample_weight_list) {
  parameters = get_parameters(sample_weight, px = 50, r_RR = 3)
  C = parameters$C
  A = parameters$A
  B = parameters$B
  A_d = parameters$A_d
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
    simulation_result = simulation(parameters, n=1000)
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
  
  mean(bias_AB_matrix)
  mean(bias_AB_d_matrix)
}


# plot the results -------------------------------------------------------------
# functions ----
# a function to get location of a matrix from the vector index
get_index <- function(i) {
  return(c(i - ((i-1) %/% py) * py, (i-1) %/% py + 1))
}


plot_entrybias = function (data) {
  ggplot()+
    geom_boxplot(aes(y = data))+
    labs(title = "Box plot for AB_hat - C entry", y = "AB_hat - C entry")
}


plot_entrybias_d = function (data) {
  ggplot()+
    geom_boxplot(aes(y = data))+
    labs(title = "Box plot for debiased_AB_hat - C entry", y = "debiased_AB_hat - C entry")
}


# boxplot the bias of AB_hat -----------------------------------------------------------
weight_to_plot <- "0.2" # choose the weight to plot for AB

myplots <- apply(result_AB_list[[weight_to_plot]], MARGIN = 1, plot_entrybias)
for (i in 1:(px*py)){
  names(myplots)[i] <- paste0("(", paste(get_index(i)[1], sep = ", ", get_index(i)[2]), ")")
}

abs_mean_entry_bias <- rep(NA, (px*py))
sd_entry_bias <- rep(NA, (px*py))
for (i in 1:(px*py)) {
  abs_mean_entry_bias[i] <- abs(mean(result_AB_list[[weight_to_plot]][i,]))
  sd_entry_bias[i] <- sd(result_AB_list[[weight_to_plot]][i,])
}

index1 <- order(abs_mean_entry_bias, decreasing = T)
index2 <- order(sd_entry_bias, decreasing = T)
myplots[index1[5:1]]

# boxplot the bias of AB_d_hat ------------------------------------------------------------
weight_to_plot <- "0.2"
myplots_d <- apply(result_AB_d_list[[weight_to_plot]], MARGIN = 1, plot_entrybias_d)
for (i in 1:(px*py)){
  names(myplots_d)[i] <- paste0("(", paste(get_index(i)[1], sep = ", ", get_index(i)[2]), ")")
}

abs_mean_entry_bias_d <- rep(NA, (px*py))
sd_entry_bias_d <- rep(NA, (px*py))
for (i in 1:(px*py)) {
  abs_mean_entry_bias_d[i] <- abs(mean(result_AB_d_list[[weight_to_plot]][i,]))
  sd_entry_bias_d[i] <- sd(result_AB_d_list[[weight_to_plot]][i,])
}

index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
index2_d <- order(sd_entry_bias_d, decreasing = T)
myplots_d[index1_d[5:1]]

# comment ----
# myplots[index2[10:1]]

#A_bias_boxplot <- ggplot()+
#  geom_boxplot(aes(y = result_B_list[[weight_to_plot]]))+
#  labs(title = "Box plot for A - A_hat norm", y = "A - A_hat norm")
#B_bias_boxplot <- ggplot()+
#  geom_boxplot(aes(y = result_B_list[[weight_to_plot]]))+
#  labs(title = "Box plot for B - B_hat norm", y = "B - B_hat norm")

# myplots[["A"]] <- A_bias_boxplot
# myplots[["B"]] <- B_bias_boxplot

# # save pdf
# pdf("myplots.pdf")
# 
# print(A_bias_boxplot)
# print(B_bias_boxplot)
# for (plot in myplots[index1]){
#   print(plot)
# }
# 
# dev.off()



# # plot A and B ----
# result_A_df <- as.data.frame(result_A_list)
# 
# long_data <- pivot_longer(result_A_df, cols = everything(), names_to = "Category", values_to = "Values")
# A_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot of the bias of A_hat for Each weight", x = "weight", y = "A_hat-A norm") + theme_minimal()
# A_bias
# 
# result_B_df <- as.data.frame(result_B_list)
# 
# long_data <- pivot_longer(result_B_df, cols = everything(), names_to = "Category", values_to = "Values")
# B_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot of the bias of B_hat for Each weight", x = "weight", y = "B_hat-B norm") + theme_minimal()
# B_bias
# 
# # plot debiased A and B ----
# result_A_d_df <- as.data.frame(result_A_d_list)
# 
# long_data <- pivot_longer(result_A_d_df, cols = everything(), names_to = "Category", values_to = "Values")
# A_d_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot of the bias of debiased_A_hat for Each weight", x = "weight", y = "A_hat-A norm") + theme_minimal()
# A_d_bias
# 
# result_B_d_df <- as.data.frame(result_B_d_list)
# 
# long_data <- pivot_longer(result_B_d_df, cols = everything(), names_to = "Category", values_to = "Values")
# B_d_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot of the bias of debiased_B_hat for Each weight", x = "weight", y = "B_hat-B norm") + theme_minimal()
# B_d_bias




            