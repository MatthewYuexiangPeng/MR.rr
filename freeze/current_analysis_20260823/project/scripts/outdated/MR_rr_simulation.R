rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
# library(mr.divw)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap::pheatmap)
# library(readxl)
# library(glmnet)

set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/Yuexiang Peng/UW/Research/Ye Ting/sim_ArBr_bias")
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

  # need to subtract the Sigma_xx by Sigma_x. Sigma_x can be estimated by the regression of exp ~ z (have been approximate in the get parameter function)
  # if (regularized == FALSE) {
  #   result_d <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  # } else if (regularized == TRUE){
  #   result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)
  # } else {
  result_d <- mr_rr(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X)
  result_d_r <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regularization_rate)

  A_d_hat = result_d$A
  B_d_hat = result_d$B
  AB_d_hat = result_d$AB

  A_d_hat_r = result_d_r$A
  B_d_hat_r = result_d_r$B
  AB_d_hat_r = result_d_r$AB
  return(list(A_hat, B_hat, AB_hat,
              A_d_hat, B_d_hat, AB_d_hat,
              A_d_hat_r, B_d_hat_r, AB_d_hat_r,
              x_j_hat, y_j_hat))
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
  # me_weight =1
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
                            rep(0.001, r-r_approx))) # set at 0.001 if want C to be rank > r
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
                    VY_tilde = VY_tilde
                    # ,true.A_sparse = true.A_sparse, true.B_sparse = true.B_sparse
                    )
  return(parameters)
}
#####


# The simulation function for the naive MR-rr estimator and the MR-rr estimator (with spectral regularization)
run_simulation <- function(regularized = TRUE, regularization_rate = 1e-13){
  sample_weight_list = c(1, 0.5, 0.2, 0.1, 0.05)
  eloop = 5000
  result_AB_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
  result_AB_d_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
  result_AB_d_r_list = list("1" = NA, "0.5" = NA, "0.2" = NA, "0.1" = NA, "0.05" = NA)
  iv_strength_list = c()
  parameters_list = c()

  for (me_weight in sample_weight_list) {
    parameters = .get_parameters(me_weight, px = 24, r_RR = 5)
    C = parameters$C
    A = parameters$A
    B = parameters$B
    A_d = parameters$FALSEA_d
    B_d = parameters$B_d
    C_r = parameters$C_r
    py = parameters$py
    px = parameters$px

    iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
    iv_strength_list = c(iv_strength_list, iv_strength)

    parameters_list = c(parameters_list, list(c(parameters, iv_strength = iv_strength)))

    bias_AB_matrix = bias_AB_d_matrix = bias_AB_d_r_matrix = matrix(NA, nrow = py*px, ncol = eloop)
    B_star = t(B) # in order to compute the norm
    B_d_star = t(B_d)
    for (i in 1:eloop) {
      simulation_result = .simulation(parameters, regularized = regularized, regularization_rate)
      A_hat = simulation_result[[1]]
      B_hat = simulation_result[[2]]
      B_hat_star = t(simulation_result[[2]])
      AB_hat = simulation_result[[3]]

      # bias_AB_vectorized = as.vector(AB_hat - C)
      bias_AB_vectorized = as.vector(AB_hat - C_r)
      bias_AB_matrix[,i] = bias_AB_vectorized

      A_d_hat = simulation_result[[4]]
      B_d_hat = simulation_result[[5]]
      B_d_hat_star = t(simulation_result[[5]])
      AB_d_hat = simulation_result[[6]]

      bias_AB_d_vectorized = as.vector(AB_d_hat - C_r)
      # bias_AB_d_vectorized = as.vector(AB_d_hat - C)
      bias_AB_d_matrix[,i] = bias_AB_d_vectorized

      A_d_r_hat = simulation_result[[7]]
      B_d_r_hat = simulation_result[[8]]
      B_d_r_hat_star = t(simulation_result[[8]])
      AB_d_r_hat = simulation_result[[9]]

      bias_AB_d_r_vectorized = as.vector(AB_d_r_hat - C_r)
      bias_AB_d_r_matrix[,i] = bias_AB_d_r_vectorized
    }
    char <- as.character(me_weight)
    result_AB_list[[char]] = bias_AB_matrix
    result_AB_d_list[[char]] = bias_AB_d_matrix
    result_AB_d_r_list[[char]] = bias_AB_d_r_matrix
  }
  return(list(result_AB_list=result_AB_list, result_AB_d_list=result_AB_d_list, result_AB_d_r_list=result_AB_d_r_list,
              iv_strength_list=iv_strength_list, parameters_list=parameters_list))
}


simulate_result = run_simulation()
save(simulate_result, file = "results/simulate_result_241123_5000_1e13.RData")

load("results/simulate_result_241123_5000_1e13.RData")
# simulate_result$iv_strength_list


# The function to create the boxplot for the bias of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization) comparing to true C for each entry.
plot_boxplot_old <- function(result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
                         weight_to_plot, individual_plot = FALSE, rank_by = "bias")
  {
  px <- 24
  py <- 10
  # individual boxplot the bias of AB_hat
  myplots <- apply(result_naive_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias)
  for (i in 1:(px*py)){
    names(myplots)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
  }

  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }

  index1 <- order(abs_mean_entry_bias, decreasing = T)
  index2 <- order(sd_entry_bias, decreasing = T)

  # individual boxplot the bias of AB_d_hat
  myplots_d <- apply(result_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d)
  for (i in 1:(px*py)){
    names(myplots_d)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
  }

  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
  }

  index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
  index2_d <- order(sd_entry_bias_d, decreasing = T)

  if (individual_plot){
    if (rank_by == "bias") {
      # individual box plot
      return(list(myplots[index1][1:2], myplots_d[index1_d][1:2]))
    } else if (rank_by == "sd"){
      # individual box plot
      return(list(myplots[index2][1:2], myplots_d[index2_d][1:2]))
    } else {
      stop("rank_by should be either 'bias' or 'sd'")
    }
  } else {
    if (rank_by == "bias"){
      # top 5 box plot in the same plot
      combined_plot <- myplots_d[[index1_d[5]]] + myplots_d[[index1_d[4]]] +
        myplots_d[[index1_d[3]]] + myplots_d[[index1_d[2]]] + myplots_d[[index1_d[1]]] +
        myplots[[index1[5]]] + myplots[[index1[4]]] + myplots[[index1[3]]] +
        myplots[[index1[2]]] + myplots[[index1[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
        patchwork::plot_annotation(
          title = "Bias boxplots of the entries with top 5 average bias \n for the MR-rr and naive MR-rr estimator",
          subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
          # caption = "C_r defined as minimizing modified objective function of
          # 2.2 or equivalent objective function of 2.1"
          )
      return(combined_plot)
    } else if (rank_by == "sd"){
      # rank by var
      combined_plot_var <- myplots_d[[index2_d[5]]] + myplots_d[[index2_d[4]]] +
        myplots_d[[index2_d[3]]] + myplots_d[[index2_d[2]]] + myplots_d[[index2_d[1]]] +
        myplots[[index2[5]]] + myplots[[index2[4]]] + myplots[[index2[3]]] +
        myplots[[index2[2]]] + myplots[[index2[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
        patchwork::plot_annotation(
          title = "Bias boxplots of the entries with top 5 standard deviation \n for the MR-rr and naive MR-rr estimator",
          subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
          # caption = "C_r defined as minimizing modified objective function of
          # 2.2 or equivalent objective function of 2.1"
          )
      return(combined_plot_var)
    } else {
      stop("rank_by should be either 'bias' or 'sd'")}
  }
}


# The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
plot_heatmap_old <- function(parameters_list, result_naive_MRrr_list, result_MRrr_list, weight_to_plot,regu) {
  px <- 24
  py <- 10
  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }

  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
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

  # color breaks
  min_value <- min(min(avg_bias), min(avg_bias_d)) #TODO:
  max_value <- max(max(avg_bias),max(avg_bias_d))

  # min_value <- min(min(as.numeric(unlist(avg_bias))), min(as.numeric(unlist(avg_bias_d)))) #TODO:
  # max_value <- max(max(as.numeric(unlist(avg_bias))),max(as.numeric(unlist(avg_bias_d))))

  breaks <- seq(min_value, max_value, length.out = 101)

  pheatmap::pheatmap(avg_bias, breaks = breaks,
                     main = sprintf("Absolute bias of the naive MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("white", "red"))(100),
                     height = 10,
                     width = 8)
  if (regu==TRUE){
    pheatmap::pheatmap(avg_bias_d, breaks = breaks,
                       main = sprintf("Absolute bias of the MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("white", "red"))(100),
                       height = 10,
                       width = 8)
  } else {
    pheatmap::pheatmap(avg_bias_d, breaks = breaks,
                       main = sprintf("Absolute bias of the MR-rr estimator with regularization by entry, weight = %s",weight_to_plot), fontsize = 8,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("white", "red"))(100),
                       height = 10,
                       width = 8)
  }
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


  pheatmap::pheatmap(avg_sd, breaks = breaks_sd,
                     main = sprintf("SD of the naive MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
                     cluster_rows = FALSE, cluster_cols = FALSE,
                     color = colorRampPalette(c("white", "red"))(100),
                     height = 10,
                     width = 8)

  if (regu==TRUE){
    pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
                       main = sprintf("SD of the MR-rr estimator with regularization by entry, weight = %s",weight_to_plot), fontsize = 8,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("white", "red"))(100),
                       height = 10,
                       width = 8)
  } else {
    pheatmap::pheatmap(avg_sd_d, breaks = breaks_sd,
                       main = sprintf("SD of the MR-rr estimator by entry, weight = %s",weight_to_plot), fontsize = 8,
                       cluster_rows = FALSE, cluster_cols = FALSE,
                       color = colorRampPalette(c("white", "red"))(100),
                       height = 10,
                       width = 8)
  }

  # # compute average iv strength
  # iv_list = list()
  # for (i in 1:10){
  #   parameters = .get_parameters(as.numeric(weight_to_plot), px = 24, r_RR = 5) # weight: 0.2, 0.5, 1 are comparable with Yinxiang's Paper
  #   iv_strength = min(eigen(solve(.sqrt_matrix(parameters$Sigma_X)) %*% parameters$VX_tilde %*% solve(.sqrt_matrix(parameters$Sigma_X)))$values)
  #   iv_list[[i]] = iv_strength
  # }
  # iv_strength = mean(as.numeric(unlist(iv_list)))

  # bias scale compares to true C_r and sd
  weight_index = match(weight_to_plot, c("1", "0.5", "0.2", "0.1", "0.05"))
  parameters = parameters_list[[weight_index]]
  C_r = parameters$C_r
  mean_abs_C_entry =  mean(abs(as.matrix(C_r)))
  mean_abs_bias_naive = mean(abs(as.matrix(avg_bias)))
  mean_abs_bias = mean(abs(as.matrix(avg_bias_d)))
  mean_sd_naive = mean(as.matrix(avg_sd))
  mean_sd = mean(as.matrix(avg_sd_d))
  #####

  return(list(mean_abs_C_entry = mean_abs_C_entry,
              mean_abs_bias_naive = mean_abs_bias_naive, mean_abs_bias = mean_abs_bias,
              mean_sd_naive = mean_sd_naive, mean_sd = mean_sd))
}


# The function to generate the heatmap of the average absolute bias and the standard deviation of the naive MR-rr estimator and the MR-rr estimator (with spectral regularization).
plot_heatmap <- function(parameters_list,
                         result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
                         weight_to_plot)
  {
  px <- 24
  py <- 10
  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }

  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
  }

  abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d_r <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
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

  combined_hm_bias = gridExtra::grid.arrange(
    grobs = list(hm_bias_1, hm_bias_2, hm_bias_3),
    ncol = 3,
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

  combined_hm_sd = gridExtra::grid.arrange(
    grobs = list(hm_sd_1, hm_sd_2, hm_sd_3),
    ncol = 3,
    top = grid::textGrob(sprintf("SD by entry, weight = %s",weight_to_plot), gp = grid::gpar(fontsize = 16))
  )
  combined_hm_sd

  # bias scale compares to true C_r and sd
  weight_index = match(weight_to_plot, c("1", "0.5", "0.2", "0.1", "0.05"))
  parameters = parameters_list[[weight_index]]
  C_r = parameters$C_r
  mean_abs_C_entry =  round(mean(abs(as.matrix(C_r))), 3)
  mean_abs_bias_naive = round(mean(abs(as.matrix(avg_bias))), 3)
  mean_abs_bias = round(mean(abs(as.matrix(avg_bias_d))), 3)
  mean_abs_bias_r = round(mean(abs(as.matrix(avg_bias_d_r))), 3)
  mean_sd_naive = round(mean(as.matrix(avg_sd)), 3)
  mean_sd = round(mean(as.matrix(avg_sd_d)), 3)
  mean_sd_r = round(mean(as.matrix(avg_sd_d_r)), 3)
  #####

  ggplot2::ggsave(sprintf("results/combined_hm_bias_weight_%s.png",weight_to_plot), plot = combined_hm_bias, width = 12, height = 4)
  ggplot2::ggsave(sprintf("results/combined_hm_sd_weight_%s.png",weight_to_plot), plot = combined_hm_sd, width = 12, height = 4)

  return(list(mean_abs_C_entry = mean_abs_C_entry,
              mean_abs_bias_naive = mean_abs_bias_naive, mean_abs_bias = mean_abs_bias, mean_abs_bias_r = mean_abs_bias_r,
              mean_sd_naive = mean_sd_naive, mean_sd = mean_sd, mean_sd_r = mean_sd_r))
}


plot_heatmap(parameters_list = simulate_result$parameters_list,
             result_naive_MRrr_list = simulate_result$result_AB_list,
             result_MRrr_list = simulate_result$result_AB_d_list,
             result_r_MRrr_list = simulate_result$result_AB_d_r_list,
             weight_to_plot="0.5")


plot_boxplot <- function(result_naive_MRrr_list, result_MRrr_list, result_r_MRrr_list,
                         weight_to_plot, individual_plot = FALSE, rank_by = "bias")
{
  px <- 24
  py <- 10
  # individual boxplot the bias of AB_hat
  myplots <- apply(result_naive_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias)
  for (i in 1:(px*py)){
    names(myplots)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
  }

  abs_mean_entry_bias <- rep(NA, (px*py))
  mean_entry_bias <- rep(NA, (px*py))
  sd_entry_bias <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias[i] <- abs(as.numeric(mean(result_naive_MRrr_list[[weight_to_plot]][i,])))
    mean_entry_bias[i] <- mean(result_naive_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias[i] <- sd(result_naive_MRrr_list[[weight_to_plot]][i,])
  }

  index1 <- order(abs_mean_entry_bias, decreasing = T)
  index2 <- order(sd_entry_bias, decreasing = T)

  # individual boxplot the bias of AB_d_hat
  myplots_d <- apply(result_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d)
  for (i in 1:(px*py)){
    names(myplots_d)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
  }

  abs_mean_entry_bias_d <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d[i] <- abs(as.numeric(mean(result_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d[i] <- mean(result_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d[i] <- sd(result_MRrr_list[[weight_to_plot]][i,])
  }

  index1_d <- order(abs_mean_entry_bias_d, decreasing = T)
  index2_d <- order(sd_entry_bias_d, decreasing = T)

  # individual boxplot the bias of AB_d_r_hat
  myplots_d_r <- apply(result_r_MRrr_list[[weight_to_plot]], MARGIN = 1, .plot_entrybias_d_r)
  for (i in 1:(px*py)){
    names(myplots_d_r)[i] <- paste0("(", paste(.get_index(i, py)[1], sep = ", ", .get_index(i, py)[2]), ")")
  }

  abs_mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for the heatmap
  mean_entry_bias_d_r <- rep(NA, (px*py)) # this is for reporting the avg bias
  sd_entry_bias_d_r <- rep(NA, (px*py))
  for (i in 1:(px*py)) {
    abs_mean_entry_bias_d_r[i] <- abs(as.numeric(mean(result_r_MRrr_list[[weight_to_plot]][i,]))) # TODO: 1. Careful about the absolute. Misleading very much, should consider the sign of bias for all entries. 2. generate complex number here for unknown reason
    mean_entry_bias_d_r[i] <- mean(result_r_MRrr_list[[weight_to_plot]][i,])
    sd_entry_bias_d_r[i] <- sd(result_r_MRrr_list[[weight_to_plot]][i,])
  }

  index1_d_r <- order(abs_mean_entry_bias_d_r, decreasing = T)
  index2_d_r <- order(sd_entry_bias_d_r, decreasing = T)


  if (individual_plot){
    if (rank_by == "bias") {
      # individual box plot
      return(list(myplots[index1][1:2], myplots_d[index1_d][1:2]))
    } else if (rank_by == "sd"){
      # individual box plot
      return(list(myplots[index2][1:2], myplots_d[index2_d][1:2]))
    } else {
      stop("rank_by should be either 'bias' or 'sd'")
    }
  } else {
    if (rank_by == "bias"){
      # top 5 box plot in the same plot

      # 定义每一行图形
      row1 <- myplots[[index1[3]]] + myplots[[index1[2]]] + myplots[[index1[1]]]
      row2 <- myplots_d[[index1_d[3]]] + myplots_d[[index1_d[2]]] + myplots_d[[index1_d[1]]]
      row3 <- myplots_d_r[[index1_d_r[3]]] + myplots_d_r[[index1_d_r[2]]] + myplots_d_r[[index1_d_r[1]]]

      # 定义左侧标题
      title1 <- patchwork::wrap_elements(grid::textGrob("Naive MR-rr", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))
      title2 <- patchwork::wrap_elements(grid::textGrob("MR-rr", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))
      title3 <- patchwork::wrap_elements(grid::textGrob("MR-rr with regularization", rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold")))

      # 合并标题和对应的图形
      row1_labeled <- title1 + row1 + patchwork::plot_layout(widths = c(1, 5))
      row2_labeled <- title2 + row2 + patchwork::plot_layout(widths = c(1, 5))
      row3_labeled <- title3 + row3 + patchwork::plot_layout(widths = c(1, 5))

      # 合并所有行，并添加全局标题
      combined_plot <- row1_labeled / row2_labeled / row3_labeled +
        patchwork::plot_annotation(
          title = "Bias boxplots of the entries with top 3 average bias",
          subtitle = sprintf("Weight = %s", weight_to_plot)
        )

      combined_plot

      ggplot2::ggsave(filename = sprintf("bias_boxplot_%s.png", weight_to_plot), plot = combined_plot, width = 10, height = 10)
    } else if (rank_by == "sd"){
      # rank by var
      combined_plot_var <- myplots_d[[index2_d[5]]] + myplots_d[[index2_d[4]]] +
        myplots_d[[index2_d[3]]] + myplots_d[[index2_d[2]]] + myplots_d[[index2_d[1]]] +
        myplots[[index2[5]]] + myplots[[index2[4]]] + myplots[[index2[3]]] +
        myplots[[index2[2]]] + myplots[[index2[1]]] + patchwork::plot_layout(nrow = 2, ncol = 5) +
        patchwork::plot_annotation(
          title = "Bias boxplots of the entries with top 3 standard deviation \n for the MR-rr and naive MR-rr estimator",
          subtitle = sprintf("rank C = 5, Sigma_X multiplied by %s", weight_to_plot)
          # caption = "C_r defined as minimizing modified objective function of
          # 2.2 or equivalent objective function of 2.1"
        )
      return(combined_plot_var)
    } else {
      stop("rank_by should be either 'bias' or 'sd'")}
  }
}


plot_boxplot(result_naive_MRrr_list = simulate_result$result_AB_list,
             result_MRrr_list = simulate_result$result_AB_d_list,
             result_r_MRrr_list = simulate_result$result_AB_d_r_list,
             weight_to_plot="0.2")


# # choose regularization rate \lambda ----
# parameters = .get_parameters(me_weight=0.2, r_RR = 5, px,  r_approx = 5)
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
# W_sqrt = .sqrt_matrix(parameters$weight.matrix)
# # sample true effect gamma_j_star(xj) and Gamma_j_star(yj)
# gamma_j_star = MASS::mvrnorm(n = n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
# Gamma_j_star = gamma_j_star %*% t(C)
#
# # sample x_j_hat, y_j_hat
# x_j_hat = matrix(0, n, px)
# y_j_hat = matrix(0, n, py)
# for (j in 1:n) {
#   x_j_hat[j,] = MASS::mvrnorm(n = 1, mu = gamma_j_star[j,], Sigma = Sigma_X, tol = 100)
#   y_j_hat[j,] = MASS::mvrnorm(n = 1, mu = Gamma_j_star[j,], Sigma = Sigma_Y, tol = 100)
# }
#
# test.rate = 1e-11
# result_d <- mr_rr_regularized(y_j_hat, x_j_hat, r=r_RR, sqrt_Gamma = W_sqrt, Sigma_X = Sigma_X, regularization_rate=test.rate)
# C_hat = result_d$AB
# # objective function value
# obj = norm((y_j_hat - x_j_hat %*% t(C_hat))%*%W_sqrt, type = "F")^2
# print(obj)
# ####
