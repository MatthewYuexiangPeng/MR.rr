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


.get_parameters = function(me_weight, r_RR = 5, px=24, r_approx = 5){
  # ## test
  # me_weight =1
  # r_RR = 5
  # px=24
  # r_approx = 5
  
  # r_RR is the rank chose by user when performing RRR
  # var_Z & VX_tilde ----------------------------------------------------------------
  pz = 1000 # pz can be changed to any number
  py = 10 # py can be changed to any number from 1 to 24
  gamma_j = as.matrix(lip_data[,paste0('gamma_exp',1:24)])
  var_Z = 2 * lip_data$eaf * (1 - lip_data$eaf) # each Z is sum of two alleles, so var(Z) = var(Z^1+Z^2) = 2*var(Z^1), where Z1, Z2 ~ Binomial(eaf.outcome)
  var_Z_raw = var_Z
  
  
  z_index = sample(1:length(var_Z), pz, replace = TRUE)
  var_Z = var_Z[z_index]
  sqrt_var_Z = sqrt(var_Z)
  
  # TODO: consider use the same index here?
  # lip_data_index_z <- sample(1:114, pz, replace = TRUE) # Question: Is this p_Z too large?
  gamma_j_sample <- gamma_j[z_index,]
  gamma_j_star_temp <- gamma_j_sample * sqrt_var_Z # Q 2025: suppose to be devided instead of multiply? 
  
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
  # Q 2025: devided by var z
  sigma_gamma_j = sigma_gamma_j * sqrt(var_Z_raw)
  Sigma = lip_corr[1:24,1:24]
  # change "Sigma" to matrix
  Sigma = as.matrix(Sigma)
  sqrt_Sigma = .sqrt_matrix(Sigma)
  pz_lip_data = nrow(gamma_j)
  Sigma_Xj = lapply(1:pz_lip_data, function(j)
    diag(sigma_gamma_j[j,]) %*% Sigma %*% diag(sigma_gamma_j[j,]))
  Sigma_Xj_sample = Sigma_Xj[z_index] # Q 2025: same index
  Sigma_X_temp = lapply(1:pz, function(j) Sigma_Xj_sample[[j]]) # Q 2025: changed from lapply(1:pz, function(j) Sigma_Xj_sample[[j]]*var_Z[j])
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
  # n_Y = median(lip_samplesize$samplesize) # assume all outcomes are from the same dataset. You can adjust this to be up to 500K
  # var_Y = sample(lip_data$se_out1[sample(1:114, 1119, replace = TRUE)]^2 * bmi.cad$N.outcome * var_Z_raw[sample(1:length(var_Z_raw), 1119, replace = TRUE)], py) # Var(Y_k) can be in the range of this (although looks strange probably because the coef is from logistic model, but I ignore this for now)
  # var_Y = sample(lip_data$se_out1^2 * lip_samplesize$samplesize[sample(1:nrow(lip_samplesize), 114, replace = TRUE)]/ var_Z_raw[sample(1:length(var_Z_raw), 114, replace = TRUE)],py) # Q 2025: delete the var_Z[j]?
  
  # var_Y = sample((lip_data$se_out1^2 /var_Z_raw) * lip_samplesize$samplesize[sample(1:nrow(lip_samplesize), 114, replace = TRUE)],py)
  var_Y = sample((lip_data$se_out1^2 * var_Z_raw) * n_Y,py) # Q 2025: devided by var z?
  
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


#### real data ####
#### Iv strength #### 
# read.csv('data/traits_1e-4.csv')

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

IV_strength_real = c()
for (i in 1:4){
  IV_strength_real = c(IV_strength_real, mvmr.ivw(gamma_exp_j, se_exp_j, gamma_out_j[,i], se_out_j[,i], gen_cor = cor_exp)$iv_strength_parameter)
}

#### Estimate C using different methods #### 
px = 9
py = 4
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
r_RR = 2

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

C_naive_MRrr = mr_rr_naive(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W)$AB

C_MRrr = mr_rr(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X)$AB

#### 5. MRrr with regularization ####
## choose rate
SigmaXX = cov(gamma_exp_j)
VX_tilde = SigmaXX - Sigma_X

# check if PSD -- pass
# if (any(eigen(VX_tilde)$values < 0)) {
#   stop("VX_tilde is not positive semi-definite.")
# }

W_sqrt = .sqrt_matrix(W)

mu_min = min(eigen(solve(.sqrt_matrix(Sigma_X)) %*% VX_tilde %*% solve(.sqrt_matrix(Sigma_X)))$values)
# iv_strength = mu_min * sqrt(1000)

c=10
q_list = seq(1, 40, 1)
regu_rate_list = c()
obj_value = c()
for (q in q_list) {
  sigma_y2 = mean(eigen(Sigma_X)$values)
  
  regu_rate = sigma_y2^2 * exp(c * (q/sqrt(n)-(mu_min+1)))/n
  
  # regu_rate = c * sigma_y2^2 * (q-(mu_min+1)/n)
  # regu_rate = exp(c*(q-mu_min*sqrt(n)))
  regu_rate_list = c(regu_rate_list, regu_rate)
  result_d <- mr_rr_regularized(gamma_out_j, gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate)
  C_hat = result_d$AB
  # objective function value
  obj = norm((gamma_out_j - gamma_exp_j %*% t(C_hat))%*%W_sqrt, type = "F")^2
  # Q 2025: add debias term
  debias_term = sum(diag(W_sqrt %*% C_hat %*% Sigma_X %*% t(C_hat) %*% W_sqrt)) # tarce
  obj = obj - debias_term
  obj_value = c(obj_value, obj)
}

opt_rate = regu_rate_list[order(obj_value)][1]

MRrr_regularized = mr_rr_regularized(Y = gamma_out_j, X = gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=opt_rate)
C_MRrr_regularized = MRrr_regularized$AB
A_MRrr_regularized = MRrr_regularized$A
B_MRrr_regularized = MRrr_regularized$B


#### 6. Mr_DAG ####
C_Mr_DAG = Mr_DAG(Y=gamma_out_j, X=gamma_exp_j, niter = 10000, burnin = 2000)

#### A time B visualizzation ####
library(pheatmap)

rownames(A_MRrr_regularized) <- paste0("out", 1:4)
colnames(A_MRrr_regularized) <- paste0("comp", 1:2)

pheatmap(A_MRrr_regularized,
         main = "Matrix A (Outcome × Component)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100))
rownames(B_MRrr_regularized) <- paste0("comp", 1:2)
colnames(B_MRrr_regularized) <- paste0("exp", 1:9)

pheatmap(B_MRrr_regularized,
         main = "Matrix B (Component × Exposure)",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("blue", "white", "red"))(100))


library(networkD3)
install.packages("webshot2")
webshot::install_phantomjs()
library(webshot)


plot_AB_sankey_signed <- function(A, B) {
  outcomes <- paste0("out", 1:nrow(A))
  components <- paste0("comp", 1:ncol(A))
  exposures <- paste0("exp", 1:ncol(B))
  
  node_labels <- c(outcomes, components, exposures)
  nodes <- data.frame(name = node_labels, stringsAsFactors = FALSE)
  
  # A: outcome -> component
  A_links <- expand.grid(out = 1:nrow(A), comp = 1:ncol(A))
  A_links$value <- as.vector(abs(A))
  A_links$sign <- ifelse(as.vector(A) >= 0, "pos", "neg")
  A_links$source <- A_links$out - 1
  A_links$target <- A_links$comp + length(outcomes) - 1
  
  # B: component -> exposure
  B_links <- expand.grid(comp = 1:nrow(B), exp = 1:ncol(B))
  B_links$value <- as.vector(abs(B))
  B_links$sign <- ifelse(as.vector(B) >= 0, "pos", "neg")
  B_links$source <- B_links$comp + length(outcomes) - 1
  B_links$target <- B_links$exp + length(outcomes) + length(components) - 1
  
  # Combine links
  links <- rbind(
    A_links[, c("source", "target", "value", "sign")],
    B_links[, c("source", "target", "value", "sign")]
  )
  
  # Define LinkGroup as factor to map color
  links$group <- links$sign
  link_color_list <- "d3.scaleOrdinal().domain(['pos', 'neg']).range(['steelblue', 'crimson'])"
  
  sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    LinkGroup = "group",
    fontSize = 13,
    nodeWidth = 20,
    colourScale = link_color_list
  )
}

plot_AB_sankey_signed(A=A_MRrr_regularized, B=B_MRrr_regularized)

saveWidget(sankey, "sankey_temp.html", selfcontained = TRUE)
webshot("sankey_temp.html", file = "sankey_AB_signed.png", vwidth = 1000, vheight = 600)


#### visualization with heatmap ####
draw_multiple_heatmaps <- function(
    matrix_list,
    titles = NULL,
    row_labels = paste0("out", 1:4),
    col_labels = paste0("exp", 1:9),
    filename = NULL,
    ncol = NULL,
    global_title = NULL,
    cellwidth = 30,
    cellheight = 20,
    fontsize_row = 8,
    fontsize_col = 8,
    color_palette = colorRampPalette(c("blue", "white", "red"))(50)
) {
  library(pheatmap)
  library(gridExtra)
  library(grid)
  
  # 默认 title
  if (is.null(titles)) {
    titles <- paste0("Matrix ", seq_along(matrix_list))
  }
  
  if (length(titles) != length(matrix_list)) {
    stop("Length of titles must match length of matrix_list.")
  }
  
  # 画出所有 heatmap（存为 grob）
  heatmap_grobs <- lapply(seq_along(matrix_list), function(i) {
    mat <- matrix_list[[i]]
    
    if (!is.null(row_labels)) {
      if (length(row_labels) != nrow(mat)) stop("row_labels length mismatch")
      rownames(mat) <- row_labels
    }
    
    if (!is.null(col_labels)) {
      if (length(col_labels) != ncol(mat)) stop("col_labels length mismatch")
      colnames(mat) <- col_labels
    }
    
    pheatmap::pheatmap(
      mat,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      fontsize_row = fontsize_row,
      fontsize_col = fontsize_col,
      main = titles[i],
      color = color_palette,
      cellwidth = cellwidth,
      cellheight = cellheight,
      silent = TRUE  # 关键！返回 grob 而不是直接画
    )$gtable
  })
  
  # 设置 ncol 默认值
  if (is.null(ncol)) ncol <- length(matrix_list)
  
  # 合并成一个图
  combined <- gridExtra::grid.arrange(
    grobs = heatmap_grobs,
    ncol = ncol,
    top = if (!is.null(global_title)) grid::textGrob(global_title, gp = grid::gpar(fontsize = 16)) else NULL
  )
  
  # 如果指定 filename 就保存
  if (!is.null(filename)) {
    ggsave(filename, plot = combined, width = 2.5 * length(matrix_list), height = 4.5)
  }
  
  return(combined)
}

draw_multiple_heatmaps(
  matrix_list = list(C_ivw, C_adivw, C_naive_MRrr, C_MRrr, C_MRrr_regularized, C_Mr_DAG),
  titles = c("IVW", "adIVW", "Naive MR-rr", "MR-rr", "Reg-MR-rr", "Mr-DAG"),
  filename = "Causal_C_estimation_hm.png",
  ncol = 3
)






#### bootsrap interval ####
