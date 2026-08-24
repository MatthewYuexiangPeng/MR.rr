rm(list=ls())
library(mr.divw)
library(matrixStats)
library(MASS)
library(ggplot2)
library(tidyverse)
library(patchwork)

# setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias_4.23")
data("bmi.cad")
load('data/multivariate_data_medium.rda')
source("scripts/RRR.R")

set.seed(2333)


sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
}


get_parameters = function(sample_weight){
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 2000 # pz can be changed to any number
  py = 8
  px = 9
  gamma_j = as.matrix(tmp$data[,paste0('gamma_exp',1:9)])
  var_Z = 2 * bmi.cad$eaf.outcome * (1 - bmi.cad$eaf.outcome) # each Z is sum of two alleles, so var(Z) = var(Z^1+Z^2) = 2*var(Z^1), where Z1, Z2 ~ Binomial(eaf.outcome)
  var_Z_raw = var_Z
  
  bmi_index = sample(1:1119, pz, replace = TRUE)
  var_Z = var_Z[bmi_index]
  sqrt_var_Z = sqrt(var_Z)
  
  tmp_index <- sample(1:273, pz, replace = TRUE)
  gamma_j_sample <- gamma_j[tmp_index,]
  gamma_j_star_temp <- gamma_j_sample*sqrt_var_Z
  # sample cov
  VX_tilde <- cov(gamma_j_star_temp)
  VX_tilde <- VX_tilde*1/2
  
  
  # Sigma_Xj & Sigma_X ----------------------------------------------------------------
  sigma_gamma_j = as.matrix(tmp$data[,paste0('se_exp',1:9)])
  Sigma = tmp$cor.mat[1:9,1:9]
  sqrt_Sigma = sqrt_matrix(Sigma)
  pz_tmp = nrow(gamma_j)
  sqrt_Sigma_gamma_j = lapply(1:pz_tmp, function(j) sqrt_Sigma %*% diag(sigma_gamma_j[j,])) # delete?
  Sigma_Xj = lapply(1:pz_tmp, function(j) diag(sigma_gamma_j[j,]) %*% Sigma %*% diag(sigma_gamma_j[j,]))
  Sigma_Xj_sample = Sigma_Xj[tmp_index]
  Sigma_X_temp = lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
  array_3d <- array(unlist(Sigma_X_temp), dim = c(9, 9, length(Sigma_X_temp)))
  Sigma_X <- apply(array_3d, c(1, 2), mean)
  Sigma_X <- sample_weight * Sigma_X
  
  
  # Sigma_Y & weight.matrix -----------------------------------------------------------
  n_Y = median(bmi.cad$N.outcome) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
  var_Y = sample(bmi.cad$se.outcome^2 * bmi.cad$N.outcome * var_Z_raw, py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
  Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))
  Sigma_Y = sample_weight * Sigma_Y
  weight.matrix = solve(Sigma_Y)
  
  
  # SigmaXX, SigmaYX, SigmaXY, AB ----------------------------------------------------------------
  temp.A = cbind(c(rnorm(3), rep(0, 5)),
                 c(rep(0, 3), rnorm(3), rep(0, 2)),
                 c(rep(0, 6), rnorm(2)))
  temp.B = rbind(c(rnorm(3), rep(0, 6)),
                 c(rep(0, 3), rnorm(3), rep(0, 3)),
                 c(rep(0, 6), rnorm(3)))
  C <- temp.A %*% temp.B
  r = 3
  
  # SigmaXX, SigmaYX and SigmaXY
  SigmaXX <- Sigma_X + VX_tilde
  SigmaYX <- C %*% VX_tilde
  SigmaXY <- t(SigmaYX)
  sqrt_Gamma <- sqrt_matrix(weight.matrix)
  sqrt_Gamma_inv <- solve(sqrt_Gamma)
  
  M = sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX) %*% SigmaXY %*% sqrt_Gamma
  V = eigen(M)$vec[, 1:r, drop = FALSE]
  A = sqrt_Gamma_inv %*% V
  B = t(V) %*% sqrt_Gamma %*% SigmaYX %*% solve(SigmaXX)
  AB = C
  
  # test the condition (2.13)
  t(A) %*% weight.matrix %*% A
  B %*% SigmaXX %*% t(B)
  eigen(M)$val[1:r]
  
  
  # VY_tilde ----------------------------------------------------------------
  VY_tilde = C %*% VX_tilde %*% t(C)
  
  parameters = list(py = py, px = px, var_Z = var_Z, VX_tilde = VX_tilde, Sigma_X = Sigma_X, Sigma_Y = Sigma_Y, weight.matrix = weight.matrix, SigmaXX = SigmaXX, SigmaYX = SigmaYX, SigmaXY = SigmaXY, AB = AB, A = A, B = B, VY_tilde = VY_tilde)
  return(parameters)
}


simulation = function(parameters, n) {
  # n: number of samples
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
  AB = parameters$AB
  VY_tilde = parameters$VY_tilde
  
  # sample true effect gamma_j_star and Gamma_j_star
  gamma_j_star = mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde)
  Gamma_j_star = gamma_j_star %*% t(AB)
  
  # sample x_j, y_j
  x_j = matrix(0, n, px)
  y_j = matrix(0, n, py)
  for (j in 1:n) {
    x_j[j,] = mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X)
    y_j[j,] = mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y)
  }
  
  # compute A_hat, B_hat
  result <- RRR(y_j, x_j, r=3, weight.matrix)
  A_hat = result$A
  B_hat = result$B
  AB_hat = result$AB
  
  # norm
  # norm_AB_bias = norm(AB_hat-AB, type = "F")
  # norm_A_bias = norm(A_hat-A, type = "F")
  # norm_B_bias = norm(B_hat-B, type = "F")
  
  return(list(A_hat, B_hat, AB_hat))
}


# run the simulation -----------------------------------------------------------
sample_weight_list = list(5, 2, 1, 0.5, 0.2)
eloop = 1000
result_AB_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_A_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_B_list = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_A_list_raw = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_B_list_raw = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
result_AB_list_raw = list("5" = NA, "2" = NA, "1" = NA, "0.5" = NA, "0.2" = NA)
for (sample_weight in sample_weight_list) {
  parameters = get_parameters(sample_weight)
  AB = parameters$AB
  A = parameters$A
  B = parameters$B
  norm_A_list = rep(NA, eloop)
  norm_B_list = rep(NA, eloop)
  bias_AB_matrix = matrix(NA, nrow = 72, ncol = eloop)
  weighted_A_list = rep(NA, eloop)
  weighted_B_list = rep(NA, eloop)
  weighted_AB_matrix = matrix(NA, nrow = 72, ncol = eloop)
  B_star = t(B) # in order to compute the norm
  for (i in 1:eloop) {
    simulation_result = simulation(parameters, n=1000)
    A_hat = simulation_result[[1]]
    B_hat = simulation_result[[2]]
    B_hat_star = t(simulation_result[[2]])
    AB_hat = simulation_result[[3]]
    temp_A = A_hat %*% solve(t(A_hat) %*% A_hat) %*% t(A_hat) - A %*% solve(t(A) %*% A) %*% t(A)
    temp_B = B_hat_star %*% solve(t(B_hat_star) %*% B_hat_star) %*% t(B_hat_star) - B_star %*% solve(t(B_star) %*% B_star) %*% t(B_star)
    norm_A_list[i] = norm(temp_A, type = "F")
    norm_B_list[i] = norm(temp_B, type = "F")
    bias_AB_vectorized = as.vector(AB_hat - AB)
    bias_AB_matrix[,i] = bias_AB_vectorized
  }
  
  weighted_A_list = norm_A_list/norm(A %*% solve(t(A) %*% A) %*% t(A), type = "F")
  weighted_B_list = norm_B_list/norm(B_star %*% solve(t(B_star) %*% B_star) %*% t(B_star), type = "F")
  weighted_AB_matrix = bias_AB_matrix/norm(AB, type = "F")
  
  char <- as.character(sample_weight)
  result_AB_list_raw[[char]] = bias_AB_matrix
  result_A_list_raw[[char]] = norm_A_list
  result_B_list_raw[[char]] = norm_B_list
  result_AB_list[[char]] = weighted_AB_matrix
  result_A_list[[char]] = weighted_A_list
  result_B_list[[char]] = weighted_B_list
}

# we actually dont need to weight AB
result_AB_list = result_AB_list_raw


# plot the results -------------------------------------------------------------
# functions ----
# a function to get location of a matrix from the vector index
get_index <- function(i) {
  return(c(i - ((i-1) %/% 8) * 8, (i-1) %/% 8 + 1))
}


plot_entrybias = function (data) {
  ggplot()+
    geom_boxplot(aes(y = data))+
    labs(title = "Box plot for AB_hat - AB entry", y = "AB_hat - AB entry")
}


# boxplot the bias of AB ----
weight_to_plot <- "0.2" # choose the weight to plot for AB

myplots <- apply(result_AB_list[[weight_to_plot]], MARGIN = 1, plot_entrybias)
for (i in 1:72){
  names(myplots)[i] <- paste0("(", paste(get_index(i)[1], sep = ", ", get_index(i)[2]), ")")
}

abs_mean_entry_bias <- rep(NA, 72)
sd_entry_bias <- rep(NA, 72)
for (i in 1:72) {
  abs_mean_entry_bias[i] <- abs(mean(result_AB_list[[weight_to_plot]][i,]))
  sd_entry_bias[i] <- sd(result_AB_list[[weight_to_plot]][i,])
}

index1 <- order(abs_mean_entry_bias, decreasing = T)
index2 <- order(sd_entry_bias, decreasing = T)
myplots[index1[10:1]]
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


# plot A and B ----
result_A_df <- as.data.frame(result_A_list)

long_data <- pivot_longer(result_A_df, cols = everything(), names_to = "Category", values_to = "Values")
A_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot for Each weight", x = "weight", y = "norm(A_hat-A)/norm(A)") + theme_minimal()
A_bias

result_B_df <- as.data.frame(result_B_list)

long_data <- pivot_longer(result_B_df, cols = everything(), names_to = "Category", values_to = "Values")
B_bias <- ggplot(long_data, aes(x = Category, y = Values)) + geom_boxplot() + labs(title = "Boxplot for Each weight", x = "weight", y = "norm(B_hat-B)/norm(B)") + theme_minimal()
B_bias



# ----
mean_matrix <- matrix(NA, nrow = 72, ncol = 5)
sd_matrix <- matrix(NA, nrow = 72, ncol = 5)
for (i in 1:72) {
  for (j in 1:5) {
    mean_matrix[i,j] <- mean(result_AB_list[[j]][i,])
    sd_matrix[i,j] <- sd(result_AB_list[[j]][i,])
  }
}

mean_df <- as.data.frame(mean_matrix)
colnames(mean_df) <- sample_weight_list
mean_long_data <- pivot_longer(mean_df, cols = everything(), names_to = "Category", values_to = "Values")
ggplot(mean_long_data, aes(x = Category, y = Values)) + 
  geom_boxplot() + 
  labs(title = "Boxplot for Entry Mean Bias Each weight", x = "weight", y = "AB_hat - AB entries mean bias") + 
  theme_minimal()

ggplot()+
  geom_boxplot(aes(y = mean_df[, 1]))+
  labs(title = "Boxplot for Entry Mean Bias For 0.2 weight", y = "AB_hat - AB entries mean bias")


sd_df <- as.data.frame(sd_matrix)
colnames(sd_df) <- sample_weight_list
sd_long_data <- pivot_longer(sd_df, cols = everything(), names_to = "Category", values_to = "Values")
ggplot(sd_long_data, aes(x = Category, y = Values)) + 
  geom_boxplot() + 
  labs(title = "Boxplot for Entry SD Bias Each weight", x = "weight", y = "AB_hat - AB entries SD bias") + 
  theme_minimal()



