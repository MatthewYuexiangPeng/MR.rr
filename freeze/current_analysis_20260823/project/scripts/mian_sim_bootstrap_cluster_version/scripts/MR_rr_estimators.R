library(foreach)
library(doParallel)
library(mr.divw) # version 0.1.0 remotes::install_github("tye27/mr.divw", ref = "86ec1b8874c83a28799dde339eea2c3a740ea458")
# library(GRAPPLE)
library(MrDAG)
# library(CVXR)


#### helper functions ####
is_psd <- function(A, tol = 1e-8) {
  if (!is.matrix(A)) A <- as.matrix(A)
  if (nrow(A) != ncol(A)) return(FALSE)
  
  A <- (A + t(A)) / 2
  ev <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  min(ev) >= -tol
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


# admm_psd_projection <- function(Sigma_hat, mu = 1, epsilon = 1e-6, max_iter = 10, tol = 1e-6) {
#   p <- nrow(Sigma_hat)
#   B <- matrix(0, p, p)
#   Lambda <- matrix(0, p, p)
#   
#   project_to_psd <- function(Z, epsilon) {
#     eig <- eigen(Z, symmetric = TRUE)
#     eigvals <- pmax(eig$values, epsilon)
#     Z_proj <- eig$vectors %*% diag(eigvals) %*% t(eig$vectors)
#     return((Z_proj + t(Z_proj)) / 2)
#   }
#   
#   vecl <- function(M) {
#     idx <- lower.tri(M, diag = TRUE)
#     return(M[idx])
#   }
#   
#   matl <- function(x, p) {
#     M <- matrix(0, p, p)
#     idx <- lower.tri(M, diag = TRUE)
#     M[idx] <- x
#     M <- M + t(M)
#     diag(M) <- diag(M) / 2
#     return(M)
#   }
#   
#   soft_threshold <- function(x, tau) {
#     return(sign(x) * pmax(abs(x) - tau, 0))
#   }
#   
#   for (iter in 1:max_iter) {
#     A <- project_to_psd(B + Sigma_hat + mu * Lambda, epsilon)
#     
#     tmp <- A - Sigma_hat - mu * Lambda
#     vec_tmp <- vecl(tmp)
#     vec_shrink <- soft_threshold(vec_tmp, mu)
#     B <- matl(vec_tmp - vec_shrink, p)
#     
#     Lambda <- Lambda - (A - B - Sigma_hat) / mu
#     
#     r_norm <- norm(A - B - Sigma_hat, type = "F")
#     if (r_norm < tol) {
#       break
#     }
#   }
#   
#   return(A)
# }


.nearest_psd <- function(M, epsilon = 1e-6) {
  M <- (M + t(M)) / 2                              # 先对称化
  eig <- eigen(M, symmetric = TRUE)
  eigvals <- pmax(eig$values, epsilon)             # clip 负/小特征值
  out <- eig$vectors %*% (eigvals * t(eig$vectors))
  (out + t(out)) / 2                               # 再对称化，消掉数值误差
}


#### Naive MR-rr estimator ####
#' Naive MR-rr estimator
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param W a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param W_inv the inverse of W
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the naive MR-rr estimator, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the naive MR-rr estimator. See more details in the manuscript.
#' @export
mr_rr_naive = function(Y, X, r, W = NULL, W_inv = NULL) {

  XtY = crossprod(X, Y)
  XtX_inv = solve(crossprod(X))

  if (is.null(W)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  } else {
    sqrt_Gamma = .sqrt_matrix(W)
    if (is.null(W_inv)) {
      sqrt_Gamma_inv = solve(sqrt_Gamma)
    } else {
      sqrt_Gamma_inv = .sqrt_matrix(W_inv)
    }
  }

  V = sqrt_Gamma %*% t(XtY) %*% XtX_inv %*% XtY %*% sqrt_Gamma / nrow(Y)
  V = eigen(V)$vec[, 1:r, drop = FALSE]

  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(XtY) %*% XtX_inv)

  return(list(A = A, B = B, AB = A %*% B))

}


#### MR-rr estimator ####
#' MR-rr estimator
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param Sigma_X a px by px numeric matrix, the average conditional covariance matrix of the coefficients of regressing exposures on SNPs, which can be culculated by averaging the regression SEs across all SNPs. It corresponds to the \eqn{\Sigma_X} in the manuscript, see more details in the estimator section of the manuscript.
#' @param W a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param W_inv the inverse of W
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the MR-rr estimator, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the MR-rr estimator. See more details in the manuscript.
#' @export
mr_rr = function(Y, X, r, Sigma_X, W = NULL, W_inv = NULL) {

  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx

  if (is.null(W)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  } else {
    sqrt_Gamma = .sqrt_matrix(W)
    if (is.null(W_inv)) {
      sqrt_Gamma_inv = solve(sqrt_Gamma)
    } else {
      sqrt_Gamma_inv = .sqrt_matrix(W_inv)
    }
  }

  V = sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv %*% sigmaxy %*% sqrt_Gamma
  V = eigen(V)$vec[, 1:r, drop = FALSE]

  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv)

  return(list(A = A, B = B, AB = A %*% B))

}


#### MR-rr estimator with spectral regularization #### 
#' MR-rr estimator with spectral regularization
#'
#' @param Y a n by py numeric matrix, where n is the number of SNPs and py is the number of outcomes. This argument corresponds to the \eqn{\Gamma^T} in the manuscript.
#' @param X a n by px numeric matrix, where n is the number of SNPs and px is the number of exposures. This argument corresponds to the \eqn{\gamma^T} in the manuscript.
#' @param r an integer, indicating the rank of the causal effect matrix C we desire to estimate.
#' @param Sigma_X a px by px numeric matrix, the average conditional covariance matrix of the coefficients of regressing exposures on SNPs, which can be culculated by averaging the regression SEs across all SNPs. It corresponds to the \eqn{\Sigma_X} in the manuscript, see more details in the estimator section of the manuscript.
#' @param regularization_rate a numeric value, the regularization rate used in the MR-rr estimator with spectral regularization. It corresponds to the \eqn{\phi} in the manuscript, see more details in the estimator section of the manuscript.
#' @param W a py by py numeric matrix, the weight matrix used in the reduced rank regression method. It is recommended to be set as the inverse of the covariance matrix of the outcomes. If not provided, it is assumed to be the identity matrix. it corresponds to the \eqn{W^{\frac{1}{2}}} in the manuscript.
#' @param W_inv the inverse of W
#'
#' @return a list of three matrices: A, B, and AB. A is a py by r matrix, B is a r by px matrix, and AB is a py by px matrix. The matrix A and B are the A and B identified by the MR-rr estimator with spectral regularization, and there product AB is the estimated causal effect matrix C from exposures to the outcomes by the MR-rr estimator with spectral regularization. See more details in the manuscript.
#' @export
mr_rr_regularized = function(Y, X, r, Sigma_X, regularization_rate = 1e-13, W = NULL, W_inv = NULL) {

  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx
  stable_Sigma_xx = debiased_Sigma_xx + regularization_rate * debiased_Sigma_xx_inv
  stable_Sigma_xx_inv = solve(stable_Sigma_xx)

  if (is.null(W)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  } else {
    sqrt_Gamma = .sqrt_matrix(W)
    if (is.null(W_inv)) {
      sqrt_Gamma_inv = solve(sqrt_Gamma)
    } else {
        sqrt_Gamma_inv = .sqrt_matrix(W_inv)
        }
    }

  V = sqrt_Gamma %*% t(sigmaxy) %*% stable_Sigma_xx_inv %*% sigmaxy %*% sqrt_Gamma
  V = eigen(V)$vec[, 1:r, drop = FALSE]

  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(sigmaxy) %*% stable_Sigma_xx_inv)

  return(list(A = A, B = B, AB = A %*% B))

}


#### MR-rr sparse estimator ####
# mr_rr_sparse <- function(GAMMA_hat, gamma_hat, W, Sigma_X, lambda = rep(1e-3, ncol(gamma_hat)), r = 2, max_iter = 100, tol = 1e-2) {
#   n <- nrow(gamma_hat)
#   px <- ncol(gamma_hat)
#   py <- ncol(GAMMA_hat)
# 
#   W_sqrt <- .sqrt_matrix(W)
# 
#   # Step 1: gamma_tilde transformation
#   P_GAMMA <- GAMMA_hat %*% solve(t(GAMMA_hat) %*% GAMMA_hat) %*% t(GAMMA_hat)
#   P_GAMMA_prep <- diag(1, n) - P_GAMMA
#   Sigma_gammahat <- t(gamma_hat) %*% gamma_hat / n
#   matrix_part1 <- Sigma_gammahat - t(P_GAMMA %*% gamma_hat) %*% (P_GAMMA %*% gamma_hat) / n
#   matrix_part2 <- matrix_part1 - Sigma_X
# 
#   # if (!is_psd(matrix_part2)) {
#   #   matrix_part2 <- admm_psd_projection(matrix_part2, epsilon = 1e-6)
#   # }
# 
#   if (!is_psd(matrix_part2)) {
#     matrix_part2 <- .nearest_psd(matrix_part2, epsilon = 1e-6)
#   }
# 
#   R <- chol(matrix_part1)
#   Q <- chol(matrix_part2)
#   L <- solve(R) %*% Q
#   gamma_tilde <- P_GAMMA %*% gamma_hat + (P_GAMMA_prep %*% gamma_hat) %*% L
# 
#   # Step 2: Initialization by standard MR-RR
#   init_result <- mr_rr(GAMMA_hat, gamma_hat, r = r, W = W, Sigma_X = Sigma_X)
#   A_hat <- init_result$A
#   px <- ncol(gamma_hat)
# 
#   # Step 3: Alternating minimization
#   for (iter in 1:max_iter) {
#     # Fix A, optimize B
#     B_var <- CVXR::Variable(rows = r, cols = px)
#     loss <- sum_squares(t(A_hat) %*% W %*% t(GAMMA_hat) - B_var %*% t(gamma_tilde)) / n +
#       Reduce(`+`, lapply(1:px, function(k) {
#         lambda[k] * sum_entries(norm1(B_var[, k, drop = FALSE]))
#       })) + 1e-6 * sum_entries(B_var)
# 
#     prob <- Problem(Minimize(loss))
#     result <- solve(prob, solver = "OSQP")
#     B_hat <- result$getValue(B_var)
# 
#     # Fix B, optimize A
#     svd_result <- svd(B_hat %*% t(gamma_tilde) %*% GAMMA_hat %*% W_sqrt)
#     A_hat_new <- solve(W_sqrt) %*% svd_result$v %*% t(svd_result$u)
# 
#     # Check convergence
#     dist <- norm(A_hat_new %*% B_hat - A_hat %*% B_hat, "F") / max(1e-8, norm(A_hat %*% B_hat, "F"))
#     A_hat <- A_hat_new
#     if (dist < tol) break
#   }
# 
#   # Threshold small entries
#   B_hat[abs(B_hat) < 1e-2] <- 0
#   C_hat <- A_hat %*% B_hat
# 
#   return(list(A = A_hat, B = B_hat, AB = C_hat))
# }




#### ivw estimator on each outcome (no parallel) #### 
ivw_multiple_outcomes = function(Y, X, Sigma_X, Sigma_Y) {
  pz = nrow(Y)
  px = ncol(X)
  py = ncol(Y)
  
  beta.exposure = X
  beta.outcome = Y
  
  se.exposure_row = t(as.matrix(sqrt(diag(Sigma_X))))
  se.exposure = matrix(rep(se.exposure_row, each = pz), nrow = pz, byrow = TRUE)
  
  se.outcome_row = sqrt(diag(Sigma_Y))
  se.outcome = matrix(rep(se.outcome_row, each = pz), nrow = pz, byrow = TRUE)
  
  # 直接循环，不用并行
  result_matrix = matrix(NA, nrow = py, ncol = px)
  
  for (i in 1:py) {
    res = mvmr.ivw(
      beta.exposure = beta.exposure,
      se.exposure = se.exposure,
      beta.outcome = as.vector(beta.outcome[, i]),
      se.outcome = se.outcome[, i]
    )
    result_matrix[i, ] = res$beta.hat
  }
  
  return(result_matrix)
}

#### adIVW estimator on each outcome (no parallel) ####
adivw_multiple_outcomes = function(Y, X, Sigma_X, Sigma_Y) {
  pz = nrow(Y)
  px = ncol(X)
  py = ncol(Y)
  
  beta.exposure = X
  beta.outcome = Y
  
  se.exposure_row = t(as.matrix(sqrt(diag(Sigma_X))))
  se.exposure = matrix(rep(se.exposure_row, each = pz), nrow = pz, byrow = TRUE)
  
  se.outcome_row = sqrt(diag(Sigma_Y))
  se.outcome = matrix(rep(se.outcome_row, each = pz), nrow = pz, byrow = TRUE)
  
  # 直接循环，不用并行
  result_matrix = matrix(NA, nrow = py, ncol = px)
  
  for (i in 1:py) {
    res = mvmr.divw(
      beta.exposure = beta.exposure,
      se.exposure = se.exposure,
      beta.outcome = as.vector(beta.outcome[, i]),
      se.outcome = se.outcome[, i],
      phi_cand = NULL
    )
    result_matrix[i, ] = res$beta.hat
  }
  
  return(result_matrix)
}


#### MrDAG ####
# devtools::install_github("lb664/MrDAG")
Mr_DAG = function(Y, X, niter = 1000, burnin = 200) {
  data = data.frame(Y, X)
  px = ncol(X)
  py = ncol(Y)
  
  MrDAGcheck <- NULL
  MrDAGcheck$Y_idx <- 1:py
  MrDAGcheck$X_idx <- (py + 1):(px + py)
  
  # 捕获输出但保留结果
  output <- suppressMessages({
    tmp <- NULL
    capture.output({
      tmp <- MrDAG(
        data = data,
        niter = niter, burnin = burnin, thin = 5,  # 请在真实分析中增大
        MrDAGcheck = MrDAGcheck,
        fileName = NULL
      )
    }, file = NULL)
    tmp  # 正确返回对象
  })
  
  ord <- c((py + 1):(px + py), 1:py)
  result <- t(get_causaleffects(output, ord = ord)$causalEffects)
  
  return(result)
}