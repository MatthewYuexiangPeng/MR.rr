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


me_weight_list = c("2.5", "1")
effect_weight_list = c("0.25", "1")
regularization_rate_list = c(1.502516e-10, 3.716852e-13, 1.523810e-12, 3.777015e-17)


#### 10.25 ####
#### 改进版本：存储数据 + CI长度计算 ####

# 并行版本的 bootstrap 函数
nonpara_bootstrap_parallel <- function(parameters_list, me_weight = "1", effect_weight = "1",
                                       regularization_rate = 1e-13, bootstrap_size = 500,
                                       iteration = 1000, r_rank = 2, n_cores = NULL,
                                       estimator_set = "mr_only",
                                       save_estimates = TRUE) {  # 新参数：是否保存估计值
  
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
                                "estimator_set", "save_estimates"),
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
    
    # 计算每个参数的 coverage 和 CI 长度
    coverage_vec = list()
    ci_length_vec = list()
    estimates_stored = list()  # 存储所有bootstrap估计值
    
    for (name in estimator_order) {
      cov_vec = numeric(px*py)
      length_vec = numeric(px*py)
      
      if (save_estimates) {
        estimates_stored[[name]] = result_lists[[name]]
      }
      
      for (j in 1:(px*py)) {
        ci = quantile(result_lists[[name]][[j]], c(0.025, 0.975), na.rm = TRUE)
        cov_vec[j] = as.numeric(C_vec[j] >= ci[1] & C_vec[j] <= ci[2])
        length_vec[j] = ci[2] - ci[1]  # CI 长度
      }
      
      coverage_vec[[name]] = cov_vec
      ci_length_vec[[name]] = length_vec
    }
    
    return(list(
      coverage = coverage_vec,
      ci_length = ci_length_vec,
      estimates = if(save_estimates) estimates_stored else NULL
    ))
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
      as.numeric(iter_result$coverage[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(coverage_list) = estimator_order
  
  # 整理 CI 长度结果
  ci_length_list = lapply(estimator_order, function(name) {
    mat = sapply(coverage_results, function(iter_result) {
      as.numeric(iter_result$ci_length[[name]])
    })
    
    if (is.vector(mat)) {
      mat = matrix(mat, nrow = 1)
    } else if (nrow(mat) != px*py) {
      mat = t(mat)
    }
    
    return(mat)
  })
  names(ci_length_list) = estimator_order
  
  # 整理估计值（如果保存）
  estimates_list = NULL
  if (save_estimates) {
    estimates_list = lapply(estimator_order, function(name) {
      lapply(coverage_results, function(iter_result) {
        iter_result$estimates[[name]]
      })
    })
    names(estimates_list) = estimator_order
  }
  
  # 计算平均 coverage
  avg_coverage = lapply(coverage_list, function(mat) {
    mean(rowMeans(mat, na.rm = TRUE), na.rm = TRUE)
  })
  
  # 计算 CI 长度统计
  # 对每个 entry，计算所有 iteration 的平均 CI 长度，然后取中位数
  ci_length_summary = lapply(ci_length_list, function(mat) {
    # mat 是 (px*py) x iteration 矩阵
    # 对每个 entry 计算平均长度
    mean_lengths = rowMeans(mat, na.rm = TRUE)
    
    list(
      median = median(mean_lengths, na.rm = TRUE),  # 所有 entry 的中位数
      mean = mean(mean_lengths, na.rm = TRUE),      # 所有 entry 的平均值
      sd = sd(mean_lengths, na.rm = TRUE),          # 标准差
      min = min(mean_lengths, na.rm = TRUE),
      max = max(mean_lengths, na.rm = TRUE)
    )
  })
  
  message(sprintf("Completed (me = %s, effect = %s)", me_weight, effect_weight))
  
  # 打印 CI 长度摘要
  cat("\nCI Length Summary:\n")
  for (name in estimator_order) {
    cat(sprintf("  %s: Median = %.4f, Mean = %.4f, SD = %.4f\n", 
                name, 
                ci_length_summary[[name]]$median,
                ci_length_summary[[name]]$mean,
                ci_length_summary[[name]]$sd))
  }
  
  return(list(
    avg_coverage = avg_coverage,
    ci_length_summary = ci_length_summary,
    coverage_matrix = coverage_list,
    ci_length_matrix = ci_length_list,
    estimates = estimates_list  # 所有bootstrap估计值（如果保存）
  ))
}

# ============================================================================
# 主程序：MR+IVW 家族
# ============================================================================
set.seed(123)


test_sim_pred_filename = "results/simulate_result_pred_250919.RData"
load(test_sim_pred_filename)

# 计算总任务数
total_tasks = length(me_weight_list) * length(effect_weight_list)

CI_cov_result_mr_only = list()
current_task = 0

cat(sprintf("\n=== Only MR family estimators (Naive, MR, MR_r) ===\n"))
cat(sprintf("Total tasks: %d\n\n", total_tasks))

for (i in seq_along(me_weight_list)) {
  for (j in seq_along(effect_weight_list)) {
    current_task = current_task + 1
    me = me_weight_list[i]
    eff = effect_weight_list[j]
    reg = regularization_rate_list[.get_sim_index(i, j, length(effect_weight_list))]
    
    cat(sprintf("\n[Task %d/%d] me=%s, effect=%s\n", current_task, total_tasks, me, eff))
    
    result = tryCatch({
      nonpara_bootstrap_parallel(
        parameters_list = simulate_result_prediction$parameters_list,
        me_weight = me,
        effect_weight = eff,
        regularization_rate = reg,
        bootstrap_size = 300,
        iteration = 500,
        n_cores = 4,
        estimator_set = "no_mrdag", # no_mrdag mr_only
        save_estimates = TRUE  # 保存所有估计值
      )
    }, error = function(e) {
      cat(sprintf("ERROR in task %d: %s\n", current_task, e$message))
      return(NULL)
    })
    
    CI_cov_result_mr_only[[paste0("me_", me, "_eff_", eff)]] = result
    
    # 每完成一个任务就保存
    save(CI_cov_result_mr_only, file = "results/CI_coverage_parallel_temp_mr_IVW_300_500_10.28.RData")
    
    # 强制垃圾回收
    gc()
  }
}

# 最终保存
save(CI_cov_result_mr_only, file = "results/CI_coverage_MR_IVW_300_500_10.28.RData")

# ============================================================================
# 结果汇总和可视化
# ============================================================================

cat("\n\n=== Overall Results Summary ===\n\n")

# 创建汇总表
summary_table = data.frame(
  Setting = character(),
  Estimator = character(),
  Coverage = numeric(),
  CI_Median_Length = numeric(),
  CI_Mean_Length = numeric(),
  stringsAsFactors = FALSE
)

for (setting_name in names(CI_cov_result_mr_only)) {
  result = CI_cov_result_mr_only[[setting_name]]
  
  if (!is.null(result)) {
    for (est_name in names(result$avg_coverage)) {
      summary_table = rbind(summary_table, data.frame(
        Setting = setting_name,
        Estimator = est_name,
        Coverage = result$avg_coverage[[est_name]],
        CI_Median_Length = result$ci_length_summary[[est_name]]$median,
        CI_Mean_Length = result$ci_length_summary[[est_name]]$mean
      ))
    }
  }
}

print(summary_table)

# 保存汇总表
write.csv(summary_table, file = "results/summary_table_mr_only_251026.csv", row.names = FALSE)

cat("\n\nResults saved to:\n")
cat("  - results/CI_coverage_parallel_final_500_1000_mr_only.RData\n")
cat("  - results/summary_table_mr_only.csv\n")