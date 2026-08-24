rm(list=ls())
library(mr.divw)
library(matrixStats)
data("bmi.cad")
load('data/multivariate_data_medium.rda')
source("scripts/RRR.R")

set.seed(2024)


sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
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
  Gamma_j_star = gamma_j_star %*% t(C)
  
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
  norm_AB_bias = norm(AB_hat-AB, type = "F")
  # norm_A_bias = norm(A_hat-A, type = "F")
  # norm_B_bias = norm(B_hat-B, type = "F")
  
  return(norm_AB_bias)
}


# define parameters
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


# Sigma_Y & weight.matrix -----------------------------------------------------------
n_Y = median(bmi.cad$N.outcome) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
var_Y = sample(bmi.cad$se.outcome^2 * bmi.cad$N.outcome * var_Z_raw, py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
Sigma_Y = diag(sqrt(var_Y / n_Y)) %*% Sigma[1:py, 1:py] %*% diag(sqrt(var_Y / n_Y))
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

parameters = list(pz = pz, py = py, px = px, var_Z = var_Z, VX_tilde = VX_tilde, Sigma_X = Sigma_X, Sigma_Y = Sigma_Y, weight.matrix = weight.matrix, SigmaXX = SigmaXX, SigmaYX = SigmaYX, SigmaXY = SigmaXY, AB = AB, VY_tilde = VY_tilde)
# -----------------------------------------------------------------------------


# run the simulation -----------------------------------------------------------
eloop = 10
norm_AB_bias_mean_list = rep(NA, 10)
for (n in seq(100, 3000, 100)) {
  norm_AB_bias_list = rep(NA, eloop)
  for (i in 1:eloop) {
    norm_AB_bias = simulation(parameters, n=1000)
    norm_AB_bias_list[i] = norm_AB_bias
  }
  norm_AB_bias_mean_list[which(seq(100, 3000, 100) == n)] = mean(norm_AB_bias_list)
}

norm_AB = norm(AB, type = "F")
norm_ratio = norm_AB_bias_mean_list / norm_AB
plot(norm_ratio)












