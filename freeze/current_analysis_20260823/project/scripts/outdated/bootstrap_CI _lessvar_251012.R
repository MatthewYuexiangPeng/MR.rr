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


.sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
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
regularization_rate_list = c(1.502516e-10, 3.716852e-13, 1.523810e-12, 3.777015e-17)


#### non paralell version ####
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



test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)


CI_cov_result = list()
for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    result = nonpara_bootstrap(parameters_list = simulate_result_prediction$parameters_list,
                               me_weight = me,
                               effect_weight = eff,
                               regularization_rate = reg,
                               bootstrap_size = 50, 
                               iteration = 1000)
    CI_cov_result[[paste0("me_", me, "_eff_", eff)]] = result
  }
}

print(CI_cov_result)

# save data
save(CI_cov_result, file = "results/simulate_CI_cov_251014_no_MrDAG.RData")


#### non parallel with 进度 ####
#### Non-parallel version with progress bar ####
set.seed(123)
# MRDAG included
nonpara_bootstrap <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                              regularization_rate = 1e-13, bootstrap_size = 100, 
                              iteration = 100, r_rank = 2, show_progress = TRUE) {
  
  # 加载进度条包（如果需要）
  if (show_progress) {
    if (!require(pbapply, quietly = TRUE)) {
      message("Installing pbapply package for progress bar...")
      install.packages("pbapply")
      library(pbapply)
    }
  }
  
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
  
  message(sprintf("Starting bootstrap (me = %s, effect = %s)", me_weight, effect_weight))
  
  # 使用 pblapply 或普通 lapply
  loop_func = if (show_progress) pblapply else lapply
  
  results = loop_func(1:iteration, function(loop) {
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
    
    # 计算 coverage
    coverage_vec = lapply(estimator_order, function(name) {
      sapply(1:(px*py), function(j) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      })
    })
    names(coverage_vec) = estimator_order
    
    return(coverage_vec)
  })
  
  # 整理结果
  for (loop in 1:iteration) {
    for (name in estimator_order) {
      coverage_list[[name]][, loop] = results[[loop]][[name]]
    }
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  return(avg)
}


# No MRDAG version
nonpara_bootstrap_no_mrdag <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                                       regularization_rate = 1e-13, bootstrap_size = 100, 
                                       iteration = 100, r_rank = 2, show_progress = TRUE) {
  
  # 加载进度条包（如果需要）
  if (show_progress) {
    if (!require(pbapply, quietly = TRUE)) {
      message("Installing pbapply package for progress bar...")
      install.packages("pbapply")
      library(pbapply)
    }
  }
  
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
  
  message(sprintf("Starting bootstrap without MrDAG (me = %s, effect = %s)", me_weight, effect_weight))
  
  # 使用 pblapply 或普通 lapply
  loop_func = if (show_progress) pblapply else lapply
  
  results = loop_func(1:iteration, function(loop) {
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
    
    # 计算 coverage
    coverage_vec = lapply(estimator_order, function(name) {
      sapply(1:(px*py), function(j) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      })
    })
    names(coverage_vec) = estimator_order
    
    return(coverage_vec)
  })
  
  # 整理结果
  for (loop in 1:iteration) {
    for (name in estimator_order) {
      coverage_list[[name]][, loop] = results[[loop]][[name]]
    }
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  return(avg)
}


# Only MR version
nonpara_bootstrap_mr_only <- function(parameters_list, me_weight = "1", effect_weight = "1", 
                                      regularization_rate = 1e-13, bootstrap_size = 100, 
                                      iteration = 100, r_rank = 2, show_progress = TRUE) {
  
  # 加载进度条包（如果需要）
  if (show_progress) {
    if (!require(pbapply, quietly = TRUE)) {
      message("Installing pbapply package for progress bar...")
      install.packages("pbapply")
      library(pbapply)
    }
  }
  
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
  
  message(sprintf("Starting bootstrap (MR only) (me = %s, effect = %s)", me_weight, effect_weight))
  
  # 使用 pblapply 或普通 lapply
  loop_func = if (show_progress) pblapply else lapply
  
  results = loop_func(1:iteration, function(loop) {
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
    
    # 计算 coverage
    coverage_vec = lapply(estimator_order, function(name) {
      sapply(1:(px*py), function(j) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        as.integer(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      })
    })
    names(coverage_vec) = estimator_order
    
    return(coverage_vec)
  })
  
  # 整理结果
  for (loop in 1:iteration) {
    for (name in estimator_order) {
      coverage_list[[name]][, loop] = results[[loop]][[name]]
    }
  }
  
  avg = lapply(coverage_list, function(mat) mean(rowMeans(mat, na.rm = TRUE)))
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  return(avg)
}


# ============================================================================
# 主程序示例
# ============================================================================

test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)

# 计算总任务数
total_tasks = length(me_weight_list) * length(effect_weight_list)
current_task = 0

CI_cov_result = list()

cat(sprintf("\n=== Starting Bootstrap Analysis ===\n"))
cat(sprintf("Total tasks: %d\n\n", total_tasks))

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    current_task = current_task + 1
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    
    cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))
    
    # 选择要使用的函数版本：
    # 1. nonpara_bootstrap() - 包含所有 estimators (包括 MrDAG)
    # 2. nonpara_bootstrap_no_mrdag() - 不包含 MrDAG (更快)
    # 3. nonpara_bootstrap_mr_only() - 只有 MR 系列 (最快)
    
    result = nonpara_bootstrap_no_mrdag(
      parameters_list = simulate_result_prediction$parameters_list,
      me_weight = me,
      effect_weight = eff,
      regularization_rate = reg,
      bootstrap_size = 500,    # 可根据需要调整
      iteration = 1000,       # 可根据需要调整
      show_progress = TRUE    # 显示进度条
    )
    
    CI_cov_result[[paste0("me_", me, "_eff_", eff)]] = result
    
    # 实时保存结果
    save(CI_cov_result, file = "results/CI_coverage_temp.RData")
  }
}

# 最终保存
save(CI_cov_result, file = "results/CI_coverage_without_MrDAG_251014.RData")

print(CI_cov_result)


#### parallel for all estimator ####
# rm(list=ls())

# run the hiddenfunction section in the main simulation file then proceed

set.seed(123)
setwd("~/Yuexiang_Peng/UW/Research/Ye Ting/sim_ArBr_bias")
source("scripts/MR_rr_estimators.R")

# 加载并行计算包
library(foreach)
library(doParallel)

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
regularization_rate_list = c(1.502516e-10, 3.716852e-13, 1.523810e-12, 3.777015e-17)

# 并行版本的 bootstrap 函数
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = 1e-13, bootstrap_size = 100,
                                       iteration = 100, r_rank = 2, n_cores = NULL) {

  # 设置核心数（默认使用所有核心-1）
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }

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

  # 使用 parLapply
  cl <- makeCluster(n_cores)

  # 方法1: 导出整个全局环境（最简单但不优雅）
  clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)

  # 方法2: 或者明确列出所有需要的函数（推荐）
  # 找出所有以 . 开头的辅助函数
  helper_functions <- c(".simulation", ".sqrt_matrix", ".get_sim_index",
                        "mr_rr_naive", "mr_rr", "mr_rr_regularized",
                        "ivw_multiple_outcomes", "adivw_multiple_outcomes", "Mr_DAG")

  # 导出当前环境的变量
  clusterExport(cl, varlist = c("parameters", "regularization_rate", "r_rank",
                                "estimator_order", "C_vec", "Sigma_X_hat",
                                "Sigma_Y_hat", "pz", "px", "py", "bootstrap_size"),
                envir = environment())

  # 导出所有辅助函数（尝试从全局环境和当前环境）
  for (func in helper_functions) {
    tryCatch({
      if (exists(func, envir = .GlobalEnv)) {
        clusterExport(cl, varlist = func, envir = .GlobalEnv)
      }
    }, error = function(e) {
      warning(paste("Could not export function:", func))
    })
  }

  # 加载必要的包
  clusterEvalQ(cl, {
    library(MASS)
  })

  # 重新 source 脚本文件（最保险的方法）
  clusterEvalQ(cl, {
    source("scripts/MR_rr_estimators.R")
  })

  message(sprintf("Starting parallel bootstrap with %d cores (me = %s, effect = %s)",
                  n_cores, me_weight, effect_weight))

  # 并行执行每个 iteration
  coverage_results <- parLapply(cl, 1:iteration, function(loop) {

    # 每个 iteration 生成一次数据
    sim = .simulation(parameters, regularization_rate)
    x_j_hat = sim[[10]]
    y_j_hat = sim[[11]]

    # 为每个 estimator 创建存储空间
    result_lists = lapply(estimator_order, function(i) {
      replicate(px*py, numeric(bootstrap_size), simplify = FALSE)
    })
    names(result_lists) = estimator_order

    # Bootstrap 循环
    for (bt in 1:bootstrap_size) {
      idx = sample(1:pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]

      # 计算各个估计量
      r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
      r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
      r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat),
                                 Sigma_X = Sigma_X_hat,
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

      # 存储 bootstrap 结果
      for (name in estimator_order) {
        for (j in 1:(px*py)) {
          result_lists[[name]][[j]][bt] = res_list[[name]][j]
        }
      }
    }

    # 计算每个参数的 coverage
    coverage_vec = lapply(estimator_order, function(name) {
      sapply(1:(px*py), function(j) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      })
    })
    names(coverage_vec) = estimator_order

    return(coverage_vec)
  })

  # 停止并行集群
  stopCluster(cl)

  # 整理结果
  coverage_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result[[name]])
    })

    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }

    return(mat)
  })
  names(coverage_list) = estimator_order

  # 计算平均 coverage
  avg = lapply(coverage_list, function(mat) {
    if (is.matrix(mat) || is.data.frame(mat)) {
      mean(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
    } else {
      mean(mat, na.rm = TRUE)
    }
  })

  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))

  return(avg)
}

# 主程序
test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)

CI_cov_result = list()

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    
    result = nonpara_bootstrap_parallel(
      parameters_list = simulate_result_prediction$parameters_list,
      me_weight = me,
      effect_weight = eff,
      regularization_rate = reg,
      bootstrap_size = 100,
      iteration = 100,
      n_cores = 10
    )
    
    CI_cov_result[[paste0("me_", me, "_eff_", eff)]] = result
    
    # 实时保存结果
    save(CI_cov_result, file = "results/CI_coverage_parallel_temp.RData")
  }
}

print(CI_cov_result)

# 最终保存
save(CI_cov_result, file = "results/CI_coverage_parallel_final_100_100.RData")



#### parallel except MrDAG with 进度条 ####
set.seed(123)
# 并行版本的 bootstrap 函数
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = 1e-13, bootstrap_size = 500,
                                       iteration = 1000, r_rank = 2, n_cores = NULL,
                                       include_mrdag = FALSE, show_progress = TRUE) {  # 新增参数
  
  # 设置核心数（默认使用所有核心-1）
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }
  
  # 加载进度条包（如果需要）
  if (show_progress) {
    if (!require(pbapply, quietly = TRUE)) {
      message("Installing pbapply package for progress bar...")
      install.packages("pbapply")
      library(pbapply)
    }
  }
  
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
  
  # 根据参数决定使用哪些 estimator
  if (include_mrdag) {
    estimator_order = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG")
    message("Using all estimators including MrDAG")
  } else {
    estimator_order = c("IVW", "adIVW", "Naive", "MR", "MR_r")
    message("Using all estimators EXCEPT MrDAG (faster)")
  }
  
  # 使用 parLapply
  cl <- makeCluster(n_cores)
  
  # 方法1: 导出整个全局环境（最简单但不优雅）
  clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  
  # 方法2: 或者明确列出所有需要的函数（推荐）
  # 找出所有以 . 开头的辅助函数
  helper_functions <- c(".simulation", ".sqrt_matrix", ".get_sim_index",
                        "mr_rr_naive", "mr_rr", "mr_rr_regularized",
                        "ivw_multiple_outcomes", "adivw_multiple_outcomes", "Mr_DAG")
  
  # 导出当前环境的变量
  clusterExport(cl, varlist = c("parameters", "regularization_rate", "r_rank",
                                "estimator_order", "C_vec", "Sigma_X_hat",
                                "Sigma_Y_hat", "pz", "px", "py", "bootstrap_size",
                                "include_mrdag"),  # 导出新参数
                envir = environment())
  
  # 导出所有辅助函数（尝试从全局环境和当前环境）
  for (func in helper_functions) {
    tryCatch({
      if (exists(func, envir = .GlobalEnv)) {
        clusterExport(cl, varlist = func, envir = .GlobalEnv)
      }
    }, error = function(e) {
      warning(paste("Could not export function:", func))
    })
  }
  
  # 加载必要的包
  clusterEvalQ(cl, {
    library(MASS)
  })
  
  # 重新 source 脚本文件（最保险的方法）
  clusterEvalQ(cl, {
    source("scripts/MR_rr_estimators.R")
  })
  
  message(sprintf("Starting parallel bootstrap with %d cores (me = %s, effect = %s)",
                  n_cores, me_weight, effect_weight))
  
  # 并行执行每个 iteration（带进度条）
  if (show_progress) {
    # 使用 pbapply 的并行版本，带进度条
    coverage_results <- pblapply(1:iteration, function(loop) {
      
      # 每个 iteration 生成一次数据
      sim = .simulation(parameters, regularization_rate)
      x_j_hat = sim[[10]]
      y_j_hat = sim[[11]]
      
      # 为每个 estimator 创建存储空间
      result_lists = lapply(estimator_order, function(i) {
        replicate(px*py, numeric(bootstrap_size), simplify = FALSE)
      })
      names(result_lists) = estimator_order
      
      # Bootstrap 循环
      for (bt in 1:bootstrap_size) {
        idx = sample(1:pz, pz, replace = TRUE)
        x_bt = x_j_hat[idx, ]
        y_bt = y_j_hat[idx, ]
        
        # 计算各个估计量
        r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
        r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
        r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat),
                                   Sigma_X = Sigma_X_hat,
                                   regularization_rate = regularization_rate)
        r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        
        # 只在需要时计算 MrDAG
        if (include_mrdag) {
          r_mrdag = Mr_DAG(y_bt, x_bt)
        }
        
        # 构建结果列表
        res_list = list(
          IVW = as.vector(r_ivw),
          adIVW = as.vector(r_adivw),
          Naive = as.vector(r_naive$AB),
          MR = as.vector(r_mr$AB),
          MR_r = as.vector(r_mr_r$AB)
        )
        
        # 只在需要时添加 MrDAG
        if (include_mrdag) {
          res_list$MrDAG = as.vector(r_mrdag)
        }
        
        # 存储 bootstrap 结果
        for (name in estimator_order) {
          for (j in 1:(px*py)) {
            result_lists[[name]][[j]][bt] = res_list[[name]][j]
          }
        }
      }
      
      # 计算每个参数的 coverage
      coverage_vec = lapply(estimator_order, function(name) {
        sapply(1:(px*py), function(j) {
          ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
          as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
        })
      })
      names(coverage_vec) = estimator_order
      
      return(coverage_vec)
    }, cl = cl)
  } else {
    # 不使用进度条的版本
    coverage_results <- parLapply(cl, 1:iteration, function(loop) {
      
      # 每个 iteration 生成一次数据
      sim = .simulation(parameters, regularization_rate)
      x_j_hat = sim[[10]]
      y_j_hat = sim[[11]]
      
      # 为每个 estimator 创建存储空间
      result_lists = lapply(estimator_order, function(i) {
        replicate(px*py, numeric(bootstrap_size), simplify = FALSE)
      })
      names(result_lists) = estimator_order
      
      # Bootstrap 循环
      for (bt in 1:bootstrap_size) {
        idx = sample(1:pz, pz, replace = TRUE)
        x_bt = x_j_hat[idx, ]
        y_bt = y_j_hat[idx, ]
        
        # 计算各个估计量
        r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
        r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
        r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat),
                                   Sigma_X = Sigma_X_hat,
                                   regularization_rate = regularization_rate)
        r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        
        # 只在需要时计算 MrDAG
        if (include_mrdag) {
          r_mrdag = Mr_DAG(y_bt, x_bt)
        }
        
        # 构建结果列表
        res_list = list(
          IVW = as.vector(r_ivw),
          adIVW = as.vector(r_adivw),
          Naive = as.vector(r_naive$AB),
          MR = as.vector(r_mr$AB),
          MR_r = as.vector(r_mr_r$AB)
        )
        
        # 只在需要时添加 MrDAG
        if (include_mrdag) {
          res_list$MrDAG = as.vector(r_mrdag)
        }
        
        # 存储 bootstrap 结果
        for (name in estimator_order) {
          for (j in 1:(px*py)) {
            result_lists[[name]][[j]][bt] = res_list[[name]][j]
          }
        }
      }
      
      # 计算每个参数的 coverage
      coverage_vec = lapply(estimator_order, function(name) {
        sapply(1:(px*py), function(j) {
          ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
          as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
        })
      })
      names(coverage_vec) = estimator_order
      
      return(coverage_vec)
    })
  }
  
  # 停止并行集群
  stopCluster(cl)
  
  # 整理结果
  coverage_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(coverage_list) = estimator_order
  
  # 计算平均 coverage
  avg = lapply(coverage_list, function(mat) {
    if (is.matrix(mat) || is.data.frame(mat)) {
      mean(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
    } else {
      mean(mat, na.rm = TRUE)
    }
  })
  
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  
  return(avg)
}

# ============================================================================
# 主程序示例
# ============================================================================

test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)

# 计算总任务数，用于显示整体进度
total_tasks = length(me_weight_list) * length(effect_weight_list)
current_task = 0

CI_cov_result_no_mrdag = list()
current_task = 0

cat(sprintf("\n\n=== Running without MrDAG (faster) ===\n"))
cat(sprintf("Total tasks: %d\n\n", total_tasks))

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    current_task = current_task + 1
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    
    cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))
    
    result = nonpara_bootstrap_parallel(
      parameters_list = simulate_result_prediction$parameters_list,
      me_weight = me,
      effect_weight = eff,
      regularization_rate = reg,
      bootstrap_size = 500,
      iteration = 1000,
      n_cores = NULL,
      include_mrdag = FALSE,  # 排除 MrDAG
      show_progress = TRUE    # 显示进度条
    )
    
    CI_cov_result_no_mrdag[[paste0("me_", me, "_eff_", eff)]] = result
    save(CI_cov_result_no_mrdag, file = "results/CI_coverage_parallel_temp_no_mrdag.RData")
  }
}

save(CI_cov_result_no_mrdag, file = "results/CI_coverage_parallel_final_500_1000_no_mrdag.RData")

# # 打印结果
# print("Results with all estimators:")
# print(CI_cov_result_full)

print("\nResults without MrDAG:")
print(CI_cov_result_no_mrdag)



#### parallel MRrr ####
#### 改进版本：更好的进度显示 + 灵活的 estimator 选择 ####

set.seed(123)

# 并行版本的 bootstrap 函数
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = 1e-13, bootstrap_size = 500,
                                       iteration = 1000, r_rank = 2, n_cores = NULL,
                                       estimator_set = "no_mrdag") {  # 新参数：控制使用哪些 estimators
  
  # 设置核心数（默认使用所有核心-1）
  if (is.null(n_cores)) {
    n_cores = parallel::detectCores() - 1
  }
  
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
  
  # 根据参数决定使用哪些 estimator
  estimator_order = switch(estimator_set,
                           "all" = c("IVW", "adIVW", "Naive", "MR", "MR_r", "MrDAG"),
                           "no_mrdag" = c("IVW", "adIVW", "Naive", "MR", "MR_r"),
                           "mr_only" = c("Naive", "MR", "MR_r"),
                           stop("estimator_set must be one of: 'all', 'no_mrdag', 'mr_only'")
  )
  
  set_name = switch(estimator_set,
                    "all" = "All estimators (including MrDAG)",
                    "no_mrdag" = "All estimators EXCEPT MrDAG",
                    "mr_only" = "Only MR family (Naive, MR, MR_r)"
  )
  
  message(sprintf("Using: %s", set_name))
  
  # 创建集群
  cl <- makeCluster(n_cores)
  
  # 导出必要的变量和函数
  clusterExport(cl, varlist = ls(envir = .GlobalEnv), envir = .GlobalEnv)
  
  helper_functions <- c(".simulation", ".sqrt_matrix", ".get_sim_index",
                        "mr_rr_naive", "mr_rr", "mr_rr_regularized",
                        "ivw_multiple_outcomes", "adivw_multiple_outcomes", "Mr_DAG")
  
  clusterExport(cl, varlist = c("parameters", "regularization_rate", "r_rank",
                                "estimator_order", "C_vec", "Sigma_X_hat",
                                "Sigma_Y_hat", "pz", "px", "py", "bootstrap_size",
                                "estimator_set"),
                envir = environment())
  
  for (func in helper_functions) {
    tryCatch({
      if (exists(func, envir = .GlobalEnv)) {
        clusterExport(cl, varlist = func, envir = .GlobalEnv)
      }
    }, error = function(e) {
      warning(paste("Could not export function:", func))
    })
  }
  
  clusterEvalQ(cl, {
    library(MASS)
  })
  
  clusterEvalQ(cl, {
    source("scripts/MR_rr_estimators.R")
  })
  
  message(sprintf("Starting parallel bootstrap with %d cores (me = %s, effect = %s)",
                  n_cores, me_weight, effect_weight))
  
  # 创建进度跟踪文件
  progress_file <- tempfile()
  writeLines("0", progress_file)
  
  # 启动进度监控线程
  progress_thread <- function(total, interval = 2) {
    while(TRUE) {
      Sys.sleep(interval)
      if (file.exists(progress_file)) {
        current <- as.numeric(readLines(progress_file, warn = FALSE)[1])
        if (is.na(current)) current <- 0
        pct <- round(100 * current / total, 1)
        cat(sprintf("\r  Progress: %d/%d iterations (%.1f%%) completed", 
                    current, total, pct))
        flush.console()
        if (current >= total) break
      } else {
        break
      }
    }
  }
  
  # 在后台启动进度显示（仅Unix系统）
  if (.Platform$OS.type == "unix") {
    system(sprintf("Rscript -e 'for(i in 1:1000){Sys.sleep(2); if(file.exists(\"%s\")){current=as.numeric(readLines(\"%s\",warn=F)[1]); cat(sprintf(\"\\r  Progress: %%d/%d (%.1f%%%%)  \", current, 100*current/%d)); flush.console(); if(current>=%d) break}}' &",
                   progress_file, progress_file, iteration, iteration, iteration),
           wait = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
  }
  
  # 并行执行每个 iteration
  coverage_results <- parLapply(cl, 1:iteration, function(loop) {
    
    # 每个 iteration 生成一次数据
    sim = .simulation(parameters, regularization_rate)
    x_j_hat = sim[[10]]
    y_j_hat = sim[[11]]
    
    # 为每个 estimator 创建存储空间
    result_lists = lapply(estimator_order, function(i) {
      replicate(px*py, numeric(bootstrap_size), simplify = FALSE)
    })
    names(result_lists) = estimator_order
    
    # Bootstrap 循环
    for (bt in 1:bootstrap_size) {
      idx = sample(1:pz, pz, replace = TRUE)
      x_bt = x_j_hat[idx, ]
      y_bt = y_j_hat[idx, ]
      
      # 根据 estimator_set 决定计算哪些
      if (estimator_set %in% c("all", "no_mrdag", "mr_only")) {
        r_naive = mr_rr_naive(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat))
        r_mr = mr_rr(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat), Sigma_X = Sigma_X_hat)
        r_mr_r = mr_rr_regularized(y_bt, x_bt, r = r_rank, W = solve(Sigma_Y_hat),
                                   Sigma_X = Sigma_X_hat,
                                   regularization_rate = regularization_rate)
      }
      
      if (estimator_set %in% c("all", "no_mrdag")) {
        r_ivw = ivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
        r_adivw = adivw_multiple_outcomes(y_bt, x_bt, Sigma_X_hat, Sigma_Y_hat)
      }
      
      if (estimator_set == "all") {
        r_mrdag = Mr_DAG(y_bt, x_bt)
      }
      
      # 构建结果列表
      res_list = list()
      if (estimator_set %in% c("all", "no_mrdag")) {
        res_list$IVW = as.vector(r_ivw)
        res_list$adIVW = as.vector(r_adivw)
      }
      res_list$Naive = as.vector(r_naive$AB)
      res_list$MR = as.vector(r_mr$AB)
      res_list$MR_r = as.vector(r_mr_r$AB)
      if (estimator_set == "all") {
        res_list$MrDAG = as.vector(r_mrdag)
      }
      
      # 存储 bootstrap 结果
      for (name in estimator_order) {
        for (j in 1:(px*py)) {
          result_lists[[name]][[j]][bt] = res_list[[name]][j]
        }
      }
    }
    
    # 更新进度
    if (file.exists(progress_file)) {
      current <- as.numeric(readLines(progress_file, warn = FALSE)[1])
      if (is.na(current)) current <- 0
      writeLines(as.character(current + 1), progress_file)
    }
    
    # 计算每个参数的 coverage
    coverage_vec = lapply(estimator_order, function(name) {
      sapply(1:(px*py), function(j) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
      })
    })
    names(coverage_vec) = estimator_order
    
    return(coverage_vec)
  })
  
  # 停止并行集群
  stopCluster(cl)
  
  # 清理进度文件
  if (file.exists(progress_file)) {
    unlink(progress_file)
  }
  
  cat("\n")  # 换行
  
  # 整理结果
  coverage_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(coverage_list) = estimator_order
  
  # 计算平均 coverage
  avg = lapply(coverage_list, function(mat) {
    if (is.matrix(mat) || is.data.frame(mat)) {
      mean(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
    } else {
      mean(mat, na.rm = TRUE)
    }
  })
  
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  
  return(avg)
}

# ============================================================================
# 主程序示例
# ============================================================================

test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)

# 计算总任务数
total_tasks = length(me_weight_list) * length(effect_weight_list)

# ============================================================================
# 选项 1: 所有 estimators (包括 MrDAG) - 最慢但最全面
# ============================================================================
# CI_cov_result_all = list()
# current_task = 0
# 
# cat(sprintf("\n=== Option 1: All estimators (including MrDAG) ===\n"))
# cat(sprintf("Total tasks: %d\n\n", total_tasks))
# 
# for (i in seq_along(me_weight_list)) {
#   for (j in seq_along(effect_weight_list)) {
#     current_task = current_task + 1
#     me = me_weight_list[i]
#     eff = effect_weight_list[j]
#     reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
#     
#     cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))
#     
#     result = nonpara_bootstrap_parallel(
#       parameters_list = simulate_result_prediction$parameters_list,
#       me_weight = me,
#       effect_weight = eff,
#       regularization_rate = reg,
#       bootstrap_size = 500,
#       iteration = 1000,
#       n_cores = NULL,
#       estimator_set = "all"  # 使用所有 estimators
#     )
#     
#     CI_cov_result_all[[paste0("me_", me, "_eff_", eff)]] = result
#     save(CI_cov_result_all, file = "results/CI_coverage_parallel_temp_all.RData")
#   }
# }
# 
# save(CI_cov_result_all, file = "results/CI_coverage_parallel_final_500_1000_all.RData")

# ============================================================================
# 选项 2: 除了 MrDAG 的所有 estimators - 平衡速度和全面性
# ============================================================================
CI_cov_result_no_mrdag = list()
current_task = 0

cat(sprintf("\n=== Option 2: All estimators EXCEPT MrDAG (recommended) ===\n"))
cat(sprintf("Total tasks: %d\n\n", total_tasks))

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    current_task = current_task + 1
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]

    cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))

    result = nonpara_bootstrap_parallel(
      parameters_list = simulate_result_prediction$parameters_list,
      me_weight = me,
      effect_weight = eff,
      regularization_rate = reg,
      bootstrap_size = 500,
      iteration = 1000,
      n_cores = 4,
      estimator_set = "no_mrdag"  # 排除 MrDAG
    )

    CI_cov_result_no_mrdag[[paste0("me_", me, "_eff_", eff)]] = result
    save(CI_cov_result_no_mrdag, file = "results/CI_coverage_parallel_temp_no_mrdag.RData")
  }
}

save(CI_cov_result_no_mrdag, file = "results/CI_coverage_parallel_final_500_1000_no_mrdag.RData")

# ============================================================================
# 选项 3: 只有 MR 家族 (Naive, MR, MR_r) - 最快
# ============================================================================
CI_cov_result_mr_only = list()
current_task = 0

cat(sprintf("\n=== Option 3: Only MR family estimators (fastest) ===\n"))
cat(sprintf("Total tasks: %d\n\n", total_tasks))

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    current_task = current_task + 1
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]

    cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))

    result = nonpara_bootstrap_parallel(
      parameters_list = simulate_result_prediction$parameters_list,
      me_weight = me,
      effect_weight = eff,
      regularization_rate = reg,
      bootstrap_size = 500,
      iteration = 1000,
      n_cores = 4,
      estimator_set = "mr_only"  # 只用 MR 家族
    )

    CI_cov_result_mr_only[[paste0("me_", me, "_eff_", eff)]] = result
    save(CI_cov_result_mr_only, file = "results/CI_coverage_parallel_temp_mr_only.RData")
  }
}

save(CI_cov_result_mr_only, file = "results/CI_coverage_parallel_final_500_1000_mr_only.RData")

# 打印结果
cat("\n\n=== Results Summary ===\n")
print(CI_cov_result_mr_only)