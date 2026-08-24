# library(mr.divw)
# library(matrixStats)
# library(MASS)
# library(ggplot2)
# library(tidyverse)
# library(patchwork)
# library(pheatmap)
# library(readxl)
# library(glmnet)


# source("R/MR_rr_simulation.R")
# source("R/MR_rr_estimators.R")
# data("lip_data")
# data("lip_corr")
# data("lip_samplesize")

set.seed(123)
setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
source("scripts/MR_rr_estimators.R")
# source("scripts/MR_rr_simulation_1.27.R")

#####
#' Non-parametric bootstrap method to estimate variance of the Drr estimator
#'
#' @param me_weight a character value of a positive number, recommend to be chosen between 0.05 and 1. The argument indicates the measurement error weight, which means the scaler we multiply to the measurement error of the coefficients of exposure to SNPs, which is the \eqn{\Sigma_X} in the manuscript. Smaller measurement error weight implies larger IV strength.
#' @param regularized can be TRUE or FALSE, indicating whether we use the MR-rr estiomator with the spectral regularization. See the manuscript for more details.
#' @param regularization_rate a small positive numerical value, indicating the regularization rate. It is only used when regularized = TRUE. The larger regularization_rate make the estimator having smaller variance, but the rate can not be too large to ensure consistency holds, see more details in the manuscript.
#'
#' @return The average coverage rate of the entry-wise 95% confidence interval produced by non-parametric bootstrap method.
#' @export
#'
# very slow to run
# nonpara_bootstrap_withIVW <- function(parameters_list, me_weight = "0.2", regularization_rate=1e-13){
#   #### (a) get the parameters
#   weight_index = match(me_weight, c("5", "2", "1", "0.5", "0.2"))
#   parameters = parameters_list[[weight_index]]
#   C = parameters$C
#   A = parameters$A
#   B = parameters$B
#   A_d = parameters$A_d
#   B_d = parameters$B_d
#   C_r = parameters$C_r
#   C_r_vec = as.vector(C_r)
#   py = parameters$py
#   px = parameters$px
#   Sigma_X_hat = parameters$Sigma_X # view Sigma_X as an unbiased estimation, Sigma_X_hat = Sigam_X?
#   Sigma_Y_hat = parameters$Sigma_Y
#
#   #### (b) in each iteration, gen data
#   iteration = 100
#   # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
#   coverage_rate_naive = coverage_rate = coverage_rate_r = coverage_rate_IVW = matrix(NA, px*py, iteration)
#
#   for (loop in 1:iteration){
#     #### (b1) gen data in the trial
#     simulation_result = .simulation(parameters)
#     A_hat = simulation_result[[1]]
#     B_hat = simulation_result[[2]]
#     AB_hat = simulation_result[[3]]
#     A_d_hat = simulation_result[[4]]
#     B_d_hat = simulation_result[[5]]
#     AB_d_hat = simulation_result[[6]]
#     A_d_r_hat = simulation_result[[7]]
#     B_d_r_hat = simulation_result[[8]]
#     AB_d_r_hat = simulation_result[[9]]
#     AB_IVW_hat = simulation_result[[12]]
#
#     x_j_hat = simulation_result[[10]]
#     y_j_hat = simulation_result[[11]]
#
#     #### (b2) estimate V_X, C
#
#     #### (b3) bootstrap
#     bootstrap_size = 100  # sample size to calculate CI in each bt loop
#     pz = 2000
#     result_IVW_list = result_naive_list = result_list = result_r_list = list()
#     for (i in 1:(px*py)){
#       result_naive_list[[i]] = rep(0, bootstrap_size)
#       result_list[[i]] = rep(0, bootstrap_size)
#       result_r_list[[i]] = rep(0, bootstrap_size)
#       result_IVW_list[[i]] = rep(0, bootstrap_size)
#     }
#     #### (b3.1-4) gen x_j_star, y_j_star, x_j_hat_star, y_j_hat_star, estimate C_drr_hat
#     for (bt in 1:bootstrap_size){
#       # sample bt_sample_number rows from x_j_hat, y_j_hat
#       index_bt = sample(1:pz, pz, replace = TRUE)
#       x_j_hat_star = x_j_hat[index_bt,]
#       y_j_hat_star = y_j_hat[index_bt,]
#
#       # estimate C_drr_hat
#
#       result_IVW <- ivw_multiple_outcomes(y_j_hat_star, x_j_hat_star, Sigma_X = Sigma_X_hat, Sigma_Y = Sigma_Y_hat)
#       result_naive <- mr_rr_naive(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat))
#       result <- mr_rr(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
#       result_r <- mr_rr_regularized(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat, regularization_rate=regularization_rate)
#
#       result_IVW_vec = as.vector(result_IVW)
#       result_naive_vec = as.vector(result_naive$AB)
#       result_vec = as.vector(result$AB)
#       result_r_vec = as.vector(result_r$AB)
#
#       for (j in 1:(px*py)){
#         result_IVW_list[[j]][bt] = result_IVW_vec[j]
#         result_naive_list[[j]][bt] = result_naive_vec[j]
#         result_list[[j]][bt] = result_vec[j]
#         result_r_list[[j]][bt] = result_r_vec[j]
#       }
#     }
#
#     #### (b4) calculate the 95% percentile interval for each entry of C_drr_hat
#     result_IVW_95_list = result_naive_95_list = result_95_list = result_r_95_list = list()
#     for (i in 1:(px*py)){
#       result_IVW_95_list[[i]] = quantile(result_IVW_list[[i]], c(0.025, 0.975))
#       result_naive_95_list[[i]] = quantile(result_naive_list[[i]], c(0.025, 0.975))
#       result_95_list[[i]] = quantile(result_list[[i]], c(0.025, 0.975))
#       result_r_95_list[[i]] = quantile(result_r_list[[i]], c(0.025, 0.975))
#     }
#     #### (b5) calculate the coverage rate for each entry of C_r
#     index_in95_IVW = index_in95_naive = index_in95 = index_in95_r = rep(1, px*py)
#     for (i in 1:(px*py)){
#       if (C_r_vec[i] < result_IVW_95_list[[i]][1] | C_r_vec[i] > result_IVW_95_list[[i]][2]){
#         index_in95_IVW[i] = 0
#       }
#       if (C_r_vec[i] < result_naive_95_list[[i]][1] | C_r_vec[i] > result_naive_95_list[[i]][2]){
#         index_in95_naive[i] = 0
#       }
#       if (C_r_vec[i] < result_95_list[[i]][1] | C_r_vec[i] > result_95_list[[i]][2]){
#         index_in95[i] = 0
#       }
#       if (C_r_vec[i] < result_r_95_list[[i]][1] | C_r_vec[i] > result_r_95_list[[i]][2]){
#         index_in95_r[i] = 0
#       }
#     }
#     coverage_rate_IVW[,loop] = index_in95_IVW
#     coverage_rate_naive[,loop] = index_in95_naive
#     coverage_rate[,loop] = index_in95
#     coverage_rate_r[,loop] = index_in95_r
#   }
#
#   #### (c) calculate the average coverage rate
#   entry_coverage_rate_IVW = rowMeans(coverage_rate_IVW)
#   average_coverage_rate_IVW = mean(entry_coverage_rate_IVW)
#
#   entry_coverage_rate_naive = rowMeans(coverage_rate_naive)
#   average_coverage_rate_naive = mean(entry_coverage_rate_naive)
#   entry_coverage_rate = rowMeans(coverage_rate)
#   average_coverage_rate = mean(entry_coverage_rate)
#   entry_coverage_rate_r = rowMeans(coverage_rate_r)
#   average_coverage_rate_r = mean(entry_coverage_rate_r)
#
#   return(list(average_coverage_rate_IVW = average_coverage_rate_IVW,
#               average_coverage_rate_naive = average_coverage_rate_naive,
#               average_coverage_rate = average_coverage_rate,
#               average_coverage_rate_r = average_coverage_rate_r))
# }
#####


nonpara_bootstrap <- function(parameters_list, me_weight = "1", regularization_rate=1e-13){
  #### (a) get the parameters
  weight_index = match(me_weight, c("1.75", "1", "0.4"))
  parameters = parameters_list[[weight_index]]
  C = parameters$C
  A = parameters$A
  B = parameters$B
  A_d = parameters$A_d
  B_d = parameters$B_d
  C_r = parameters$C_r
  C_r_vec = as.vector(C_r)
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X
  Sigma_Y_hat = parameters$Sigma_Y

  #### (b) in each iteration, gen data
  iteration = 100
  # store the entry-wise CI coverage rate in each bt loop with a (px*py), iteration matrix
  coverage_rate_naive = coverage_rate = coverage_rate_r = 
    coverage_rate_IVW = coverage_rate_adIVW
    matrix(NA, px*py, iteration)

  for (loop in 1:iteration){
    #### (b1) gen data in the trial
    simulation_result = .simulation(parameters)
    x_j_hat = simulation_result[[10]]
    y_j_hat = simulation_result[[11]]

    #### (b3) bootstrap
    bootstrap_size = 100  # sample size to calculate CI in each bt loop
    pz = 2000
    result_naive_list = result_list = result_r_list = list()
    for (i in 1:(px*py)){
      result_naive_list[[i]] = rep(0, bootstrap_size)
      result_list[[i]] = rep(0, bootstrap_size)
      result_r_list[[i]] = rep(0, bootstrap_size)
    }

    for (bt in 1:bootstrap_size){
      # sample bt_sample_number rows from x_j_hat, y_j_hat
      index_bt = sample(1:pz, pz, replace = TRUE)
      x_j_hat_star = x_j_hat[index_bt,]
      y_j_hat_star = y_j_hat[index_bt,]

      result_naive <- mr_rr_naive(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat))
      result <- mr_rr(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
      result_r <- mr_rr_regularized(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat, regularization_rate=regularization_rate)

      result_naive_vec = as.vector(result_naive$AB)
      result_vec = as.vector(result$AB)
      result_r_vec = as.vector(result_r$AB)

      for (j in 1:(px*py)){
        result_naive_list[[j]][bt] = result_naive_vec[j]
        result_list[[j]][bt] = result_vec[j]
        result_r_list[[j]][bt] = result_r_vec[j]
      }
    }

    #### (b4) calculate the 95% percentile interval for each entry
    result_naive_95_list = result_95_list = result_r_95_list = list()
    for (i in 1:(px*py)){
      result_naive_95_list[[i]] = quantile(result_naive_list[[i]], c(0.025, 0.975))
      result_95_list[[i]] = quantile(result_list[[i]], c(0.025, 0.975))
      result_r_95_list[[i]] = quantile(result_r_list[[i]], c(0.025, 0.975))
    }

    #### (b5) calculate the coverage rate for each entry
    index_in95_naive = index_in95 = index_in95_r = rep(1, px*py)
    for (i in 1:(px*py)){
      if (C_r_vec[i] < result_naive_95_list[[i]][1] | C_r_vec[i] > result_naive_95_list[[i]][2]){
        index_in95_naive[i] = 0
      }
      if (C_r_vec[i] < result_95_list[[i]][1] | C_r_vec[i] > result_95_list[[i]][2]){
        index_in95[i] = 0
      }
      if (C_r_vec[i] < result_r_95_list[[i]][1] | C_r_vec[i] > result_r_95_list[[i]][2]){
        index_in95_r[i] = 0
      }
    }
    coverage_rate_naive[,loop] = index_in95_naive
    coverage_rate[,loop] = index_in95
    coverage_rate_r[,loop] = index_in95_r
  }

  #### (c) calculate the average coverage rate
  entry_coverage_rate_naive = rowMeans(coverage_rate_naive)
  average_coverage_rate_naive = mean(entry_coverage_rate_naive)
  entry_coverage_rate = rowMeans(coverage_rate)
  average_coverage_rate = mean(entry_coverage_rate)
  entry_coverage_rate_r = rowMeans(coverage_rate_r)
  average_coverage_rate_r = mean(entry_coverage_rate_r)

  return(list(average_coverage_rate_naive = average_coverage_rate_naive,
              average_coverage_rate = average_coverage_rate,
              average_coverage_rate_r = average_coverage_rate_r))
}



load("results/simulate_result_250306.RData")

# 8.386425e-11 1.447977e-11 4.013995e-12 1.188179e-12 2.115718e-13
# regularization_rate_list = c(8.295011e-11, 1.444003e-11, 4.106420e-12, 1.282637e-12, 1.230937e-13)
# sample_weight_list = c(5, 2, 1, 0.5, 0.2)

# [1] 8.212083e-11 1.447780e-11 3.202669e-12 3.463476e-13 2.466910e-57 (c2)
# new algorithm: [1] 1.350594e-11 2.204258e-12 5.255102e-13 1.112693e-13 1.279553e-14
# c=1. 1/10: [1] 1.279385e-11 1.913728e-12 3.941866e-13 6.660968e-14 3.191639e-15
regularization_rate_list = c(1.279385e-11, 1.913728e-12, 3.941866e-13, 6.660968e-14, 3.191639e-15)
sample_weight_list = c(5, 2, 1, 0.5, 0.2)

for (i in 1:length(sample_weight_list)){
  result = nonpara_bootstrap(parameters_list = simulate_result$parameters_list,
                             me_weight = sample_weight_list[i],
                             regularization_rate = regularization_rate_list[i])
  print(paste("weight=", sample_weight_list[i]))
  # print(paste("Coverage rate for IVW=", result$average_coverage_rate_IVW))
  print(paste("Coverage rate for naive=", result$average_coverage_rate_naive))
  print(paste("Coverage rate for MR=", result$average_coverage_rate))
  print(paste("Coverage rate for MR_r=", result$average_coverage_rate_r))
}


#####
## result:
# weight=0.2, regularized=TRUE, regularization_rate=1e-13 , average_coverage_rate=0.9065833
# weight=0.5, regularized=TRUE, regularization_rate=1e-13 , average_coverage_rate=0.9202917


# regularized = TRUE
# # parametric bootstrap method to estimate variance of the Drr estimator ----
# #### (a) get the parameters
# parameters = .get_parameters(me_weight = 0.2, px = 24, r_RR = 5)
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
#   simulation_result = .simulation(parameters, regularized = regularized)
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
#   #### (b3) bootstrap
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
#     if (regularized == FALSE) {
#       result_d_bt <- mr_rr(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
#     } else {
#       result_d_bt <- mr_rr_regularized(y_j_hat_star, x_j_hat_star, r=5, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat, regularization_rate=regularization_rate)
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
