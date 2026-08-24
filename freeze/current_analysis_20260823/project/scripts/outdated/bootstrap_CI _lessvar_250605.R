rm(list=ls())
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


####
.get_sim_index = function(me_weight_index, effect_weight_index, len_me = 2){
  sim_index = (me_weight_index-1)*len_me + effect_weight_index
  return(sim_index)
}

.simulation_data_only <- function(parameters) {
  n = 1000
  py = parameters$py
  px = parameters$px
  VX_tilde = parameters$VX_tilde
  Sigma_X = parameters$Sigma_X
  Sigma_Y = parameters$Sigma_Y
  C = parameters$C
  
  # true latent gamma and Gamma
  gamma_j_star = MASS::mvrnorm(n, mu = rep(0, px), Sigma = VX_tilde, tol = 100)
  Gamma_j_star = gamma_j_star %*% t(C)
  
  x_j_hat = matrix(0, n, px)
  y_j_hat = matrix(0, n, py)
  for (j in 1:n) {
    x_j_hat[j, ] = MASS::mvrnorm(1, gamma_j_star[j, ], Sigma_X, tol = 100)
    y_j_hat[j, ] = MASS::mvrnorm(1, Gamma_j_star[j, ], Sigma_Y, tol = 100)
  }
  
  return(list(x_j_hat = x_j_hat, y_j_hat = y_j_hat))
}


me_weight_list = c("2.5", "1")
effect_weight_list = c("0.25", "1")
regularization_rate_list = c(1.422259e-11, 3.600086e-12, 1.335627e-12, 3.706175e-14)



# MRDAG included
nonpara_bootstrap <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                              regularization_rate = 1e-13, bootstrap_size = 100, iteration = 100, r_rank = 2) {
  me_index = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters = parameters_list[[param_index]]
  
  C = parameters$C
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X
  Sigma_Y_hat = parameters$Sigma_Y
  C_vec = as.vector(C)
  pz = 1000
  
  estimator_order = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG")
  coverage_list = lapply(estimator_order, function(i) matrix(NA, px*py, iteration))
  names(coverage_list) = estimator_order
  
  for (loop in 1:iteration) {
    sim = .simulation(parameters, regularization_rate)
    x_j_hat = sim[[10]]
    y_j_hat = sim[[11]]
    
    result_lists = lapply(estimator_order, function(i) replicate(px*py, numeric(bootstrap_size), simplify = FALSE))
    names(result_lists) = estimator_order
    
    for (bt in 1:bootstrap_size) {
      idx = sample(1:pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]
      
      r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
      r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
      r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat,
                                 regularization_rate = regularization_rate)
      r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
      r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
      r_mrdag = Mr_DAG(y_bt, x_bt)
      
      res_list = list(
        IVW = as.vector(r_ivw),
        adIVW = as.vector(r_adivw),
        Naive = as.vector(r_naive$AB),
        MR = as.vector(r_mr$AB),
        MR_r = as.vector(r_mr_r$AB),
        MrDAG = as.vector(r_mrdag)
      )
      
      for (name in estimator_order) {
        for (j in 1:(px*py)) {
          result_lists[[name]][[j]][bt] = res_list[[name]][j]
        }
      }
    }
    
    for (name in estimator_order) {
      for (j in 1:(px*py)) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        coverage_list[[name]][j, loop] = as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      }
    }
    
    message(sprintf("Bootstrap loop %d done (me = %s, effect = %s)", loop, me_weight, effect_weight))
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  return(avg)
}


# no MRDAG
nonpara_bootstrap <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                              regularization_rate = 1e-13, bootstrap_size = 100, iteration = 100, r_rank = 2) {
  me_index = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters = parameters_list[[param_index]]
  
  C = parameters$C
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X
  Sigma_Y_hat = parameters$Sigma_Y
  C_vec = as.vector(C)
  pz = 1000
  
  estimator_order = c("IVW", "adIVW", "Naive", "MR", "MR_r")
  coverage_list = lapply(estimator_order, function(i) matrix(NA, px*py, iteration))
  names(coverage_list) = estimator_order
  
  for (loop in 1:iteration) {
    sim = .simulation_data_only(parameters)
    x_j_hat = sim$x_j_hat
    y_j_hat = sim$y_j_hat
    
    result_lists = lapply(estimator_order, function(i) replicate(px*py, numeric(bootstrap_size), simplify = FALSE))
    names(result_lists) = estimator_order
    
    for (bt in 1:bootstrap_size) {
      idx = sample(1:pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]
      
      r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
      r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
      r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat,
                                 regularization_rate = regularization_rate)
      r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
      r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
      
      res_list = list(
        IVW = as.vector(r_ivw),
        adIVW = as.vector(r_adivw),
        Naive = as.vector(r_naive$AB),
        MR = as.vector(r_mr$AB),
        MR_r = as.vector(r_mr_r$AB)
      )
      
      for (name in estimator_order) {
        for (j in 1:(px*py)) {
          result_lists[[name]][[j]][bt] = res_list[[name]][j]
        }
      }
    }
    
    for (name in estimator_order) {
      for (j in 1:(px*py)) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        coverage_list[[name]][j, loop] = as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      }
    }
    
    message(sprintf("Bootstrap loop %d done (me = %s, effect = %s)", loop, me_weight, effect_weight))
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  return(avg)
}


# only MR
nonpara_bootstrap <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                              regularization_rate = 1e-13, bootstrap_size = 100, iteration = 100, r_rank = 2) {
  me_index = match(as.character(me_weight), me_weight_list)
  effect_index = match(as.character(effect_weight), effect_weight_list)
  param_index = .get_sim_index(me_index, effect_index, length(effect_weight_list))
  parameters = parameters_list[[param_index]]
  
  C = parameters$C
  py = parameters$py
  px = parameters$px
  Sigma_X_hat = parameters$Sigma_X
  Sigma_Y_hat = parameters$Sigma_Y
  C_vec = as.vector(C)
  pz = 1000
  
  estimator_order = c("Naive", "MR", "MR_r")
  coverage_list = lapply(estimator_order, function(i) matrix(NA, px*py, iteration))
  names(coverage_list) = estimator_order
  
  for (loop in 1:iteration) {
    sim = .simulation_data_only(parameters)
    x_j_hat = sim$x_j_hat
    y_j_hat = sim$y_j_hat
    
    result_lists = lapply(estimator_order, function(i) replicate(px*py, numeric(bootstrap_size), simplify = FALSE))
    names(result_lists) = estimator_order
    
    for (bt in 1:bootstrap_size) {
      idx = sample(1:pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]
      
      r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
      r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
      r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat,
                                 regularization_rate = regularization_rate)
      
      res_list = list(
        Naive = as.vector(r_naive$AB),
        MR = as.vector(r_mr$AB),
        MR_r = as.vector(r_mr_r$AB)
      )
      
      for (name in estimator_order) {
        for (j in 1:(px*py)) {
          result_lists[[name]][[j]][bt] = res_list[[name]][j]
        }
      }
    }
    
    for (name in estimator_order) {
      for (j in 1:(px*py)) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        coverage_list[[name]][j, loop] = as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      }
    }
    
    message(sprintf("Bootstrap loop %d done (me = %s, effect = %s)", loop, me_weight, effect_weight))
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  return(avg)
}



load("results/simulate_result_big_boxplot_250521.RData")


CI_cov_result = list()
for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    result = nonpara_bootstrap(parameters_list = simulate_result_big_boxplot$parameters_list,
                               me_weight = me,
                               effect_weight = eff,
                               regularization_rate = reg,
                               bootstrap_size = 300, 
                               iteration = 1000)
    CI_cov_result[[paste0("me_", me, "_eff_", eff)]] = result
  }
}

print(CI_cov_result)











