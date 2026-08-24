debiased_RRR = function(Y, X, r, Sigma_X, sqrt_Gamma = NULL, sqrt_Gamma_inv = NULL) {

  sigmaxy = crossprod(X, Y) / nrow(Y)
  debiased_Sigma_xx = crossprod(X) / nrow(Y) - Sigma_X
  debiased_Sigma_xx_inv = solve(debiased_Sigma_xx) # use the random effect variance matrix to replace Sigma_xx

  if (is.null(sqrt_Gamma)) {
    sqrt_Gamma = diag(1, ncol(Y), ncol(Y))
    sqrt_Gamma_inv = sqrt_Gamma
  }

  V = sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv %*% sigmaxy %*% sqrt_Gamma
  V = eigen(V)$vec[, 1:r, drop = FALSE]

  if (is.null(sqrt_Gamma_inv)) sqrt_Gamma_inv = solve(sqrt_Gamma)

  A = sqrt_Gamma_inv %*% V
  B = crossprod(V, sqrt_Gamma %*% t(sigmaxy) %*% debiased_Sigma_xx_inv)

  return(list(A = A, B = B, AB = A %*% B))

}

# B = rbind(c(rnorm(3), rep(0, 6)),
#                c(rep(0, 3), rnorm(3), rep(0, 3)),
#                c(rep(0, 6), rnorm(3)))
# A = cbind(c(rnorm(3), rep(0, 5)),
#           c(rep(0, 3), rnorm(3), rep(0, 2)),
#           c(rep(0, 6), rnorm(2)))
#
# n = 1000
# p = 9
# q = 8
#
# Sigma_Y = clusterGeneration::genPositiveDefMat(q)$Sigma
# sqrt_Sigma_Y = sqrt_matrix(Sigma_Y)
# sqrt_Sigma_Y_inv = sqrt_matrix(Sigma_Y, inv = TRUE)
#
# X = matrix(rnorm(n * p), nrow = n)
# E = matrix(rnorm(n * q), nrow = n) %*% sqrt_Sigma_Y
# Y = X %*% t(A %*% B) + E
#
# norm(RRR(Y, X, 3, sqrt_Sigma_Y_inv, sqrt_Sigma_Y)$AB - A %*% B, "F")
# norm(RRR(Y, X, 3)$AB - A %*% B, "F")
# norm(RRR(Y, X, 8, sqrt_Sigma_Y_inv, sqrt_Sigma_Y)$AB - A %*% B, "F")
# norm(RRR(Y, X, 8)$AB - A %*% B, "F")
# norm(t(lsfit(X, Y, intercept = FALSE)$coef) - A %*% B, "F")
#
# library(matrixStats)
# (colVars(Y) - colVars(E)) / colVars(Y)
# desired_R2 = 0.3
# c = sqrt(colVars(E) * desired_R2 / colVars(X %*% t(A %*% B)) / (1 - desired_R2))
# A = diag(c) %*% A
# Y = X %*% t(A %*% B) + E
# (colVars(Y) - colVars(E)) / colVars(Y)
