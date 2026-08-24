# ============================================================ #
# 1. loading environment and data ####
# ============================================================ #
rm(list=ls())
# library(devtools)
# install_github("tye27/mr.divw")
library(mr.divw)
library(GRAPPLE)
# library(matrixStats)
# library(MASS)
library(ggplot2)
# library(tidyverse)
library(patchwork)
# library(pheatmap::pheatmap)
# library(readxl)
# library(glmnet)
library(foreach)
library(doParallel)
library(reshape2)
library(dplyr)
library(ggdist)
library(pheatmap)


set.seed(123)
#setwd("D:/24 Winter UW/Reduced Rank Regression/sim_V_bias")
setwd("~/UW/Research/Ye Ting/sim_ArBr_bias")
# data("bmi.cad")
# load('data/multivariate_data_medium.rda')

# load Lipid data
lip_data = read.csv('data/lipids_total24_5e-08.csv')
lip_corr = read.csv('data/lipids_total24_5e-08_cor_mat.csv')
lip_samplesize = readxl::read_excel('data/Kennetu_2016_download_links_updated.xlsx')

source("scripts/MR_rr_estimators.R")

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

exp_name = read.csv('data/traits_1e-4.csv')$x
# out_name = c("IS", "LAS", "CES", "SVS")

# 26/2/8 delete first outcome
gamma_out_j = gamma_out_j[,2:4]
se_out_j = se_out_j[,2:4]
cor_out = as.matrix(lip_corr_real[11:13,11:13])

out_name = c("LAS", "CES", "SVS")


# ============================================================ #
# 2. functions ####
# ============================================================ #
.sqrt_matrix = function(mat, inv = FALSE) {
  eigen_mat = eigen(mat)
  if (inv) {
    d = 1 / sqrt(eigen_mat$val)
  } else {
    d = sqrt(eigen_mat$val)
  }
  eigen_mat$vec %*% diag(d) %*% t(eigen_mat$vec)
}


estimate_C <- function(gamma_exp_j,se_exp_j,gamma_out_j,se_out_j,cor_exp,cor_out){
  px = 9
  py = 3
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
  
  r_RR = rank_test(W, gamma_out_j, gamma_exp_j, Sigma_X, print = TRUE, min_rank = 1)
  
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
  
  # old version of selecting rate
  # c=10
  # q_list = seq(1, 40, 1)
  # regu_rate_list = c()
  # obj_value = c()
  # for (q in q_list) {
  #   sigma_y2 = mean(eigen(Sigma_X)$values)
  #   
  #   regu_rate = sigma_y2^2 * exp(c * (q/sqrt(n)-(mu_min+1)))/n
  #   
  #   # regu_rate = c * sigma_y2^2 * (q-(mu_min+1)/n)
  #   # regu_rate = exp(c*(q-mu_min*sqrt(n)))
  #   regu_rate_list = c(regu_rate_list, regu_rate)
  #   result_d <- mr_rr_regularized(gamma_out_j, gamma_exp_j, r=r_RR, W = W, Sigma_X = Sigma_X, regularization_rate=regu_rate)
  #   C_hat = result_d$AB
  #   # objective function value
  #   obj = norm((gamma_out_j - gamma_exp_j %*% t(C_hat))%*%W_sqrt, type = "F")^2
  #   # Q 2025: add debias term
  #   debias_term = sum(diag(W_sqrt %*% C_hat %*% Sigma_X %*% t(C_hat) %*% W_sqrt)) # tarce
  #   obj = obj - debias_term
  #   obj_value = c(obj_value, obj)
  # }
  # 
  # opt_rate = regu_rate_list[order(obj_value)][1]
  
  # 25.10.2
  D_list = seq(0, 15, 1)
  regu_rate_list = c()
  obj_value = c()
  for (D in D_list) {
    sigma_y2 = mean(eigen(Sigma_X)$values)
    
    # regu_rate = sigma_y2^2 * exp(c * (q/sqrt(n)-(mu_min+1)))/n
    regu_rate = sigma_y2^2 * exp(0.3*(D - sqrt(n) * mu_min))/n
    
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
  
  
  #### 6. MR-rr sparse ####
  MRrr_sparse = mr_rr_sparse(GAMMA_hat = gamma_out_j, gamma_hat = gamma_exp_j,
                             W = W, Sigma_X = Sigma_X, r = r_RR)
  C_MRrr_sparse = MRrr_sparse$AB
  #### 7. Mr_DAG ####
  C_Mr_DAG = Mr_DAG(Y=gamma_out_j, X=gamma_exp_j, niter = 10000, burnin = 2000)
  
  output = list(
    C_ivw              = C_ivw,
    C_adivw            = C_adivw,
    C_naive_MRrr       = C_naive_MRrr,
    C_MRrr             = C_MRrr,
    C_MRrr_regularized = C_MRrr_regularized,
    C_MRrr_sparse      = C_MRrr_sparse,
    C_Mr_DAG           = C_Mr_DAG
  )
}


.compute_iv_strength_all_directions <- function(beta.exposure, se.exposure, beta.outcome, se.outcome, gen_cor = NULL) {
  if (ncol(beta.exposure) <= 1 | ncol(se.exposure) <= 1) {
    stop("either beta.exposure or se.exposure only has one column; if univariable MR is performed, please use univariable MR methods.")
  }
  
  K <- ncol(beta.exposure)
  if (is.null(gen_cor)) {
    P <- diag(K)
  } else {
    P <- as.matrix(gen_cor)
  }
  
  if (ncol(P) != ncol(beta.exposure)) {
    stop("The shared correlation matrix has a different number of columns than the input beta.exposure")
  }
  
  if (nrow(beta.exposure) != length(beta.outcome)) {
    stop("The number of SNPs in beta.exposure and beta.outcome is different")
  }
  
  beta.exposure <- as.matrix(beta.exposure)
  se.exposure <- as.matrix(se.exposure)
  p <- nrow(beta.exposure)
  
  # 生成标准化所需矩阵
  P_eigen <- eigen(P)
  P_root_inv <- P_eigen$vectors %*% diag(1 / sqrt(P_eigen$values)) %*% t(P_eigen$vectors)
  
  Vj_root_inv <- lapply(1:p, function(j) P_root_inv %*% diag(1 / se.exposure[j, ]))
  
  # 构建 IV strength matrix
  IV_strength_matrix <- Reduce("+", lapply(1:p, function(j) {
    beta.exposure.V <- Vj_root_inv[[j]] %*% beta.exposure[j, ]
    beta.exposure.V %*% t(beta.exposure.V)
  })) - p * diag(K)
  
  # 计算所有方向的 IV strength（特征值）
  eig_vals <- eigen(IV_strength_matrix / sqrt(p), symmetric = TRUE)$values
  
  # 两位小数
  eig_vals <- round(eig_vals, 2)
  
  return(eig_vals)  # 输出所有方向的IV strength
}


# Use this one eventrually in paper - 260603
rank_test_M0 = function(W, y_j_hat, x_j_hat, Sigma_X, print = TRUE, min_rank = 1){
  px = ncol(x_j_hat)
  py = ncol(y_j_hat)
  pz = nrow(y_j_hat)
  sigmaxy = crossprod(x_j_hat, y_j_hat) / pz
  debiased_Sigma_xx = crossprod(x_j_hat) / pz - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  W_sqrt = .sqrt_matrix(W)
  V = W_sqrt %*% t(sigmaxy) %*% debiased_Sigma_xx_inv %*% sigmaxy %*% W_sqrt
  lambda = eigen(V)$values
  
  p_value = c()
  for (r in min_rank:(min(px,py)-1)) {
    log_sum_tail <- sum(log(1 + lambda[(r + 1):length(lambda)]))
    M = (pz-(px+py+1)/2) * log_sum_tail
    # 在卡方分布的百分比
    p_value = c(p_value, 1 - pchisq(M, df = (py - r) * (px - r)))
  }
  # 第一个为False的位置
  res_idx <- which(p_value >= 0.05)[1]
  
  if (is.na(res_idx)) {
    if (print) {
      cat("No rank passes the test; all p-values are < 0.05.\n")
      cat("Defaulting to full rank:", min(px, py), "\n")
    }
    return(min(px, py))
  } else {
    if (print) {
      cat("The lowest rank satisfying the test is:", res_idx+min_rank-1, 
          "; p-value =", round(p_value[res_idx], 3), "\n")
    }
    return(res_idx)
  }
}


rank_test_old = function(W, y_j_hat, x_j_hat, Sigma_X, print = TRUE, min_rank = 1, bt_loop = 100){
  # # test 250826
  # min_rank = 1
  # bt_loop = 100
  
  # use bootstrap to estimate var(\hat C_MRrr)
  # store C hats
  px = ncol(x_j_hat)
  py = ncol(y_j_hat)
  pz = nrow(y_j_hat)
  
  p_value = c()
  for (r in min_rank:(min(px,py)-1)) {
    C_hat_MRrr_vec = as.vector(mr_rr(y_j_hat, x_j_hat, r=r, W = W, Sigma_X = Sigma_X)$AB)
    
    C_hat_matrix = matrix(0, nrow = bt_loop, ncol = px * py)
    for (b in 1:bt_loop) {
      sample_idx = sample(1:nrow(y_j_hat), replace = TRUE)
      y_j_hat_bt = y_j_hat[sample_idx, ]
      x_j_hat_bt = x_j_hat[sample_idx, ]
      C_hat_MRrr_bt_vec = as.vector(mr_rr(y_j_hat_bt, x_j_hat_bt, r=r, W = W, Sigma_X = Sigma_X)$AB)
      C_hat_matrix[b, ] = C_hat_MRrr_bt_vec
    }
    
    # variance of C_hat
    var_C_hat = cov(C_hat_matrix)
    C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% C_hat_MRrr_vec
    
    # # consider centering c hat
    # mu_hat = colMeans(C_hat_matrix)
    # C_tilde_vec = .sqrt_matrix(solve(var_C_hat)) %*% (C_hat_MRrr_vec - mu_hat)
    
    C_tilde = matrix(C_tilde_vec, nrow = py, ncol = px)
    
    # singular value of C_tilde
    svd_C_tilde = svd(C_tilde)
    lambda_tilde = svd_C_tilde$d
    
    # M1
    M1 = sum(lambda_tilde[(r + 1):py]**2)
    # compare to chi-square distribution (py-r,px-r)
    p_value = c(p_value, 1 - pchisq(M1, df = (py - r) * (px - r)))
  }
  
  # 第一个为False的位置
  res_idx <- which(p_value >= 0.05)[1]
  
  if (is.na(res_idx)) {
    if (print) {
      cat("No rank passes the test; all p-values are < 0.05.\n")
      cat("Defaulting to full rank:", min(px, py), "\n")
    }
    return(min(px, py))
  } else {
    if (print) {
      cat("The lowest rank satisfying the test is:", res_idx+min_rank-1, 
          "; p-value =", round(p_value[res_idx], 3), "\n")
    }
    return(res_idx+min_rank-1)
  }
}


# bootstrap new ver 260416
rank_test <- function(W, y_j_hat, x_j_hat, Sigma_X,
                      print = TRUE, min_rank = 1,
                      bt_loop = 500, alpha = 0.05) {
  px <- ncol(x_j_hat)
  py <- ncol(y_j_hat)
  pz <- nrow(y_j_hat)
  k  <- min(px, py)
  
  # 拟合 unconstrained (full-rank) 估计量，只需算一次
  C_hat_full <- mr_rr(y_j_hat, x_j_hat, r = k, W = W, Sigma_X = Sigma_X)$AB
  sv_full    <- svd(C_hat_full, nu = py, nv = px)
  
  out <- data.frame(r      = min_rank:(k - 1),
                    M1_obs = NA_real_,
                    q_crit = NA_real_,
                    reject = NA,
                    stringsAsFactors = FALSE)
  
  for (idx in seq_len(nrow(out))) {
    r <- out$r[idx]
    
    # ---- 1. 观测统计量：full-rank 估计量的尾部奇异值 ----
    M1_obs <- sum(sv_full$d[(r + 1):k]^2)
    
    # U2, V2 从 full-rank SVD 取正交补
    U2 <- sv_full$u[, (r + 1):py, drop = FALSE]
    V2 <- sv_full$v[, (r + 1):px, drop = FALSE]
    
    # ---- 2. Bootstrap ----
    M1_star <- numeric(bt_loop)
    for (b in 1:bt_loop) {
      idx_b <- sample.int(pz, pz, replace = TRUE)
      C_hat_star <- mr_rr(y_j_hat[idx_b, ], x_j_hat[idx_b, ],
                          r = k, W = W, Sigma_X = Sigma_X)$AB
      
      Delta <- C_hat_star - C_hat_full
      Proj  <- t(U2) %*% Delta %*% V2
      M1_star[b] <- sum(Proj^2)
    }
    
    # ---- 3. Critical value ----
    q_crit <- quantile(M1_star, probs = 1 - alpha, type = 8)
    reject <- (M1_obs > q_crit)
    
    out$M1_obs[idx] <- M1_obs
    out$q_crit[idx] <- q_crit
    out$reject[idx] <- reject
    
    if (print) {
      cat(sprintf("r=%d: M1=%.6f, q_crit=%.6f, reject=%s\n",
                  r, M1_obs, q_crit, reject))
    }
  }
  
  # ---- 4. 选最小不被 reject 的 r ----
  idx_keep <- which(out$reject == FALSE)[1]
  if (is.na(idx_keep)) {
    r_sel <- k
    if (print) cat("All rejected; default to full rank:", k, "\n")
  } else {
    r_sel <- out$r[idx_keep]
    if (print) cat("Selected rank =", r_sel, "\n")
  }
  
  return(r_sel)
}


# ============================================================ #
# 3. Iv strength ####
# ============================================================ #
# IV_strength_real = c()
# for (i in 1:4){
#   IV_strength_real = c(IV_strength_real, mvmr.ivw(gamma_exp_j, se_exp_j, gamma_out_j[,i], se_out_j[,i], gen_cor = cor_exp)$iv_strength_parameter)
# }

IV_strength_real = .compute_iv_strength_all_directions(gamma_exp_j, se_exp_j, gamma_out_j[,1], se_out_j[,1], gen_cor = cor_exp)
min(IV_strength_real)


# ============================================================ #
# 4. Estimate C using different methods ####
# ============================================================ #
real_result <- estimate_C(
  gamma_out_j = gamma_out_j,
  gamma_exp_j = gamma_exp_j,
  se_out_j = se_out_j,
  se_exp_j = se_exp_j,
  cor_out = cor_out,
  cor_exp = cor_exp
)

C_ivw <- real_result$C_ivw
C_adivw <- real_result$C_adivw
C_naive_MRrr <- real_result$C_naive_MRrr
C_MRrr <- real_result$C_MRrr
C_MRrr_regularized <- real_result$C_MRrr_regularized
C_MRrr_sparse <- real_result$C_MRrr_sparse
C_Mr_DAG <- real_result$C_Mr_DAG

norm(C_adivw, "F")
norm(C_MRrr_regularized, "F")


# ============================================================ #
# 5. visualize C hat with heatmap ####
# ============================================================ #
draw_multiple_heatmaps <- function(
    matrix_list,
    titles = NULL,
    row_labels = out_name,
    col_labels = exp_name,
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
    ggsave(filename, plot = combined, width = 2.5 * length(matrix_list), height = 6)
  }
  
  return(combined)
}

draw_multiple_heatmaps(
  matrix_list = list(C_ivw, C_adivw, C_naive_MRrr, C_MRrr,
                     C_MRrr_regularized, C_MRrr_sparse, C_Mr_DAG),
  titles = c("IVW","adIVW","Naive MR-rr","MR-rr","Reg. MR-rr","Sparse MR-rr","Mr-DAG"),
  filename = "results/Causal_C_estimation_hm_260417.png",
  ncol = 3
)


# ============================================================ #
# 6. bootsrap confidence interval ####
# ============================================================ #
set.seed(123)
# bootstrap
bootstrap_C_list <- vector("list", 7)
names(bootstrap_C_list) <- c("ivw", "adivw", "naive_MRrr", "MRrr",
                             "MRrr_regularized", "MRrr_sparse", "Mr_DAG")
B <- 1000  # bootstrap次数

for (b in 1:B) {
  boot_index <- sample(1:nrow(gamma_exp_j), replace = TRUE)
  boot_result <- estimate_C(
    gamma_out_j = gamma_out_j[boot_index, , drop = FALSE],
    gamma_exp_j = gamma_exp_j[boot_index, , drop = FALSE],
    se_out_j = se_out_j[boot_index, , drop = FALSE],
    se_exp_j = se_exp_j[boot_index, , drop = FALSE],
    cor_out = cor_out,
    cor_exp = cor_exp
  )
  for (k in names(bootstrap_C_list)) {
    if (is.null(bootstrap_C_list[[k]])) {
      bootstrap_C_list[[k]] <- array(boot_result[[paste0("C_", k)]], dim = c(3, 9, B))
    } else {
      bootstrap_C_list[[k]][,,b] <- boot_result[[paste0("C_", k)]]
    }
  }
  # 每隔10报告一次进度
  if (b %% 10 == 0) {
    cat("Bootstrap iteration:", b, "\n")
  }
}


# save bootstrap_C_list
save(bootstrap_C_list, file = "bootstrap_C_list_1000_260417.RData")
# load bootstrap_C_list
load("bootstrap_C_list_1000_260417.RData")


# ============================================================ #
# 7. plot big boxplot - with sparse ####
# ============================================================ #
exp_name <- read.csv('data/traits_1e-4.csv')$x
out_name <- c("LAS", "CES", "SVS")
names(out_name) <- 0:2
names(exp_name) <- 0:8

estimate_C_list <- list(
  ivw              = real_result$C_ivw,
  adivw            = real_result$C_adivw,
  naive_MRrr       = real_result$C_naive_MRrr,
  MRrr             = real_result$C_MRrr,
  MRrr_regularized = real_result$C_MRrr_regularized,
  MRrr_sparse      = real_result$C_MRrr_sparse,
  Mr_DAG           = real_result$C_Mr_DAG
)

est_levels <- c("ivw","adivw","naive_MRrr","MRrr",
                "MRrr_regularized","MRrr_sparse","Mr_DAG")
est_labels <- c("IVW","adIVW","Naive MR-rr","MR-rr",
                "Reg. MR-rr","Sparse MR-rr","MrDAG")

build_long_and_sig <- function(bootstrap_C_list, estimate_C_list, ci_level = 0.95) {
  a <- (1 - ci_level) / 2
  long <- data.frame(); sig <- data.frame()
  for (i in 1:3) for (j in 1:9) for (k in names(bootstrap_C_list)) {
    v  <- bootstrap_C_list[[k]][i, j, ]; v <- v[!is.na(v)]
    pe <- estimate_C_list[[k]][i, j]
    lo <- quantile(v, a); up <- quantile(v, 1 - a)
    long <- rbind(long, data.frame(Row = i-1, Col = j-1,
                                   estimator = k, value = v))
    sig  <- rbind(sig,  data.frame(Row = i-1, Col = j-1, estimator = k,
                                   point_est = pe, lower = lo, upper = up,
                                   sig = (lo > 0) | (up < 0)))
  }
  long$estimator <- factor(long$estimator, levels = est_levels, labels = est_labels)
  sig$estimator  <- factor(sig$estimator,  levels = est_levels, labels = est_labels)
  list(long = long, sig = sig)
}

dat <- build_long_and_sig(bootstrap_C_list, estimate_C_list, ci_level = 0.95)


plot_realdata_halfeye <- function(long_df, sig_df, out_name, exp_name,
                                  title = "Point estimates and 95% non-parametric bootstrap confidence intervals for real data",
                                  save_path = NULL,
                                  width = 13, height = 6,
                                  ylim = c(-0.5, 0.5)) {
  p <- ggplot(long_df, aes(x = estimator, y = value, fill = estimator)) +
    stat_halfeye(
      adjust = 1, .width = 0.95, justification = -0.2, scale = 0.8,
      point_colour = NA, interval_colour = "black", interval_size = 0.75,
      slab_alpha = 0.7, normalize = "xy"
    ) +
    geom_point(
      data = sig_df, inherit.aes = FALSE,
      aes(x = estimator, y = point_est, color = sig),
      shape = 18, size = 3
    ) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    coord_cartesian(ylim = ylim) +
    facet_grid(rows = vars(Row), cols = vars(Col),
               labeller = labeller(Row = out_name, Col = exp_name)) +
    theme_bw(base_size = 12) +
    labs(title = title, x = "Estimator", y = "Estimated C Value") +
    theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
          strip.text  = element_text(size = 10),
          legend.position = "none")
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to ", save_path)
  }
  print(p)
}


# new
plot_realdata_halfeye <- function(long_df, sig_df, out_name, exp_name,
                                  save_path = NULL,
                                  width = 13, height = 6,
                                  ylim = c(-0.5, 0.5)) {
  p <- ggplot(long_df, aes(x = estimator, y = value, fill = estimator)) +
    stat_halfeye(
      adjust = 1, .width = 0.95, justification = -0.2, scale = 0.8,
      point_colour = NA, interval_colour = "black", interval_size = 0.75,
      slab_alpha = 0.7, normalize = "xy"
    ) +
    geom_point(
      data = sig_df, inherit.aes = FALSE,
      aes(x = estimator, y = point_est, color = sig),
      shape = 18, size = 3
    ) +
    scale_color_manual(values = c(`FALSE` = "black", `TRUE` = "red"), guide = "none") +
    coord_cartesian(ylim = ylim) +
    facet_grid(rows = vars(Row), cols = vars(Col),
               labeller = labeller(Row = out_name, Col = exp_name)) +
    theme_bw(base_size = 12) +
    labs(x = "", y = "Estimated C Value", fill = "Estimator") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.text  = element_text(size = 10),
          legend.position = "right")
  
  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = width, height = height, dpi = 300)
    message("Plot saved to ", save_path)
  }
  print(p)
}


plot_realdata_halfeye(
  dat$long, dat$sig, out_name, exp_name,
  ylim = c(-0.5, 0.5),
  save_path = "results/real_data_CI_halfeye_260417_with_sparse.png"
)


# #### plot big boxplot - old ver ####
# exp_name = read.csv('data/traits_1e-4.csv')$x
# out_name = c("LAS", "CES", "SVS")
# names(out_name) <- 0:2
# names(exp_name) <- 0:8
# 
# 
# # 先构建估计值矩阵
# estimate_C_list <- list(
#   ivw = real_result$C_ivw,
#   adivw = real_result$C_adivw,
#   naive_MRrr = real_result$C_naive_MRrr,
#   MRrr = real_result$C_MRrr,
#   MRrr_regularized = real_result$C_MRrr_regularized,
#   Mr_DAG = real_result$C_Mr_DAG
# )
# 
# # 统计 bootstrap CI + 估计值
# summary_data <- data.frame()
# 
# for (i in 1:3) {
#   for (j in 1:9) {
#     for (k in names(bootstrap_C_list)) {
#       estimates <- bootstrap_C_list[[k]][i, j, ]
#       est_valid <- estimates[!is.na(estimates)]
#       
#       temp_df <- data.frame(
#         Row = i - 1,
#         Col = j - 1,
#         estimator = k,
#         # lower = quantile(est_valid, 0.025, na.rm = TRUE),
#         # upper = quantile(est_valid, 0.975, na.rm = TRUE),
#         lower = quantile(est_valid, 0.05, na.rm = TRUE),
#         upper = quantile(est_valid, 0.95, na.rm = TRUE),
#         true_value = estimate_C_list[[k]][i, j]
#       )
#       summary_data <- rbind(summary_data, temp_df)
#     }
#   }
# }
# 
# 
# # 设置 estimator 为 factor（控制颜色顺序）
# summary_data$estimator <- factor(summary_data$estimator,
#                                  levels = c("ivw", "adivw", "naive_MRrr", "MRrr", "MRrr_regularized", "Mr_DAG"),
#                                  labels = c("IVW", "adIVW", "Naive MR-rr", "MR-rr", "Reg. MR-rr", "MrDAG"))
# 
# 
# # halfeye plot
# build_summary_data_long_with_ci <- function(bootstrap_C_list, estimate_C_list, ci_level = 0.95) {
#   alpha <- (1 - ci_level) / 2
#   lower_q <- alpha
#   upper_q <- 1 - alpha
#   
#   summary_data_long <- data.frame()
#   
#   for (i in 1:3) {
#     for (j in 1:9) {
#       for (k in names(bootstrap_C_list)) {
#         estimates <- bootstrap_C_list[[k]][i, j, ]
#         est_valid <- estimates[!is.na(estimates)]
#         
#         temp_df <- data.frame(
#           Row = i - 1,
#           Col = j - 1,
#           estimator = k,
#           value = est_valid,
#           true_value = estimate_C_list[[k]][i, j],
#           lower = quantile(est_valid, probs = lower_q, na.rm = TRUE),
#           upper = quantile(est_valid, probs = upper_q, na.rm = TRUE)
#         )
#         
#         summary_data_long <- rbind(summary_data_long, temp_df)
#       }
#     }
#   }
#   
#   summary_data_long$estimator <- factor(summary_data_long$estimator,
#                                         levels = c("ivw", "adivw", "naive_MRrr", "MRrr", "MRrr_regularized", "Mr_DAG"),
#                                         labels = c("IVW", "adIVW", "Naive MR-rr", "MR-rr", "Reg. MR-rr", "MrDAG"))
#   
#   return(summary_data_long)
# }
# 
# 
# summary_data_long <- build_summary_data_long_with_ci(bootstrap_C_list, estimate_C_list)
# 
# 
# # 添加显著点
# # 由 summary_data_long 生成每个 entry 的 CI 与显著性标记
# build_sig_df <- function(summary_data_long) {
#   summary_data_long %>%
#     distinct(Row, Col, estimator, true_value, lower, upper) %>%
#     mutate(sig = (lower > 0 & upper > 0) | (lower < 0 & upper < 0))
# }
# 
# sig_df <- build_sig_df(summary_data_long)
# 
# 
# plot_realdata_halfeye <- function(summary_data, out_name, exp_name,
#                                   title = "Point estimates and 95% non-parametric bootstrap confidence intervals for real data",
#                                   save_path = NULL,
#                                   width = 12, height = 6,
#                                   ylim = c(-1, 1),
#                                   sig_mark = "*",
#                                   sig_offset = 0.02) {
#   # 每个 entry 只保留一行 + 显著性
#   sig_df <- summary_data %>%
#     distinct(Row, Col, estimator, true_value, lower, upper) %>%
#     mutate(sig = (lower > 0 & upper > 0) | (lower < 0 & upper < 0),
#            y_mark = pmin(upper + sig_offset, ylim[2] - sig_offset))
#   
#   p <- ggplot(summary_data, aes(x = estimator, y = value, fill = estimator)) +
#     stat_halfeye(
#       adjust = 1,
#       .width = 0.95,
#       justification = -0.2,
#       scale = 0.8,
#       point_colour = NA,
#       interval_colour = "black",
#       interval_size = 0.75,
#       slab_alpha = 0.7,
#       normalize = "xy"
#     ) +
#     # 点估计（所有 estimator 都画）
#     geom_point(
#       data = sig_df,
#       aes(x = estimator, y = true_value),
#       inherit.aes = FALSE,
#       shape = 18, size = 3, color = "black"
#     ) +
#     # 显著标记（仅 sig==TRUE 的 entry 画）
#     geom_point(
#       data = sig_df %>% filter(sig),
#       aes(x = estimator, y = true_value),
#       inherit.aes = FALSE,
#       shape = 18, size = 3, color = "red"
#     ) +
#     # geom_text(
#     #   data = sig_df %>% filter(sig),
#     #   aes(x = estimator, y = y_mark, label = sig_mark),
#     #   inherit.aes = FALSE,
#     #   size = 4, color = "black"
#     # ) +
#     coord_cartesian(ylim = ylim) +
#     facet_grid(
#       rows = vars(Row),
#       cols = vars(Col),
#       labeller = labeller(Row = out_name, Col = exp_name)
#     ) +
#     theme_bw(base_size = 12) +
#     labs(title = title, x = "Estimator", y = "Estimated C Value") +
#     theme(
#       axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
#       strip.text = element_text(size = 10),
#       legend.position = "none"
#     )
#   
#   if (!is.null(save_path)) {
#     ggsave(save_path, plot = p, width = width, height = height, dpi = 300)
#     message("Plot saved to ", save_path)
#   }
#   print(p)
# }
# 
# 
# plot_realdata_halfeye(summary_data_long, out_name, exp_name,
#                       ylim = c(-0.5, 0.5),
#                       save_path = "results/real_data_CI_halfeye_260409_test.png")
# 

#### investigate why extreme estimations?250620 ####
# F = ecdf(bootstrap_C_list$ivw[1,1,])
# F(C_ivw[1,1])
# 
# F = ecdf(bootstrap_C_list$ivw[2,1,])
# F(C_ivw[2,1])
# 
# F = ecdf(bootstrap_C_list$naive_MRrr[4,1,])
# F(C_naive_MRrr[4,1])
# 
# # 准备数据
# x <- bootstrap_C_list$naive_MRrr[1, 1, ]
# target <- C_naive_MRrr[1, 1]
# target_p <- ecdf(x)(target)  # 计算F值
# 
# # 构建ECDF数据框（用于绘图）
# df <- data.frame(x = x)
# 
# # 生成图
# p <- ggplot(df, aes(x = x)) +
#   stat_ecdf(geom = "step", color = "blue", linewidth = 1) +
#   geom_vline(xintercept = target, color = "red", linetype = "dashed") +
#   geom_hline(yintercept = target_p, color = "red", linetype = "dashed") +
#   annotate("point", x = target, y = target_p, color = "red", size = 2) +
#   annotate("text", x = target, y = target_p, 
#            label = sprintf("Real data estimation = %.2f", target_p), 
#            hjust = -0.1, vjust = -0.5, color = "black", size = 3.5) +
#   labs(
#     title = "Empirical CDF of Bootstrap Estimates of C(1,1) for Naive MR-rr",
#     x = "Value",
#     y = "ECDF of the Bootstrap Estimations"
#   ) +
#   theme_minimal()
# 
# 
# # 保存图像
# ggsave("results/ecdf_bootstrap_C_naive_MRrr_1_1.png", plot = p, width = 6, height = 4, bg = "white")
# # the estimation based on real data are indeed in the "middle" of bootstrap estimations. bt estimations are not distributed evenly on the line







