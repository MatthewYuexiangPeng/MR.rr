#' Measurement-Error-Corrected Reduced-Rank Mendelian Randomization
#'
#' Estimates a rank-constrained causal effect matrix while correcting the
#' empirical second-moment matrix of the SNP-exposure associations for
#' measurement error.
#'
#' @param Y A numeric matrix with one row per genetic instrument and one column
#'   per outcome.
#' @param X A numeric matrix with one row per genetic instrument and one column
#'   per exposure.
#' @param r A positive integer specifying the working rank of the causal effect
#'   matrix.
#' @param Sigma_X A symmetric matrix containing the average measurement-error
#'   covariance of the SNP-exposure association estimates.
#' @param W An optional positive-definite outcome weight matrix. If `NULL`, the
#'   identity matrix is used.
#' @param W_inv An optional precomputed inverse of `W`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{A}{The outcome loading matrix.}
#'   \item{B}{The exposure loading matrix.}
#'   \item{AB}{The estimated causal effect matrix.}
#' }
#'
#' @examples
#' X <- matrix(
#'   c(
#'     0.2, -0.1,
#'     0.5,  0.3,
#'    -0.4,  0.6,
#'     0.7, -0.5
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#'
#' Y <- matrix(
#'   c(
#'     0.1,  0.3,
#'     0.4, -0.2,
#'    -0.3,  0.5,
#'     0.6, -0.4
#'   ),
#'   nrow = 4,
#'   byrow = TRUE
#' )
#'
#' Sigma_X <- diag(c(0.01, 0.01))
#' fit <- mr_rr(Y = Y, X = X, r = 1, Sigma_X = Sigma_X)
#' fit$AB
#'
#' @export
mr_rr <- function(
    Y,
    X,
    r,
    Sigma_X,
    W = NULL,
    W_inv = NULL
) {
  inputs <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = r,
    Sigma_X = Sigma_X,
    W = W,
    W_inv = W_inv,
    require_sigma_x = TRUE
  )

  Y <- inputs$Y
  X <- inputs$X
  r <- inputs$r
  Sigma_X <- inputs$Sigma_X
  W <- inputs$W
  W_inv <- inputs$W_inv

  n <- nrow(Y)

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)

  if (is.null(W)) {
    sqrt_W <- diag(ncol(Y))
    sqrt_W_inv <- sqrt_W
  } else {
    sqrt_W <- .sqrt_matrix(W)

    if (is.null(W_inv)) {
      sqrt_W_inv <- solve(sqrt_W)
    } else {
      sqrt_W_inv <- .sqrt_matrix(W_inv)
    }
  }

  target_matrix <- sqrt_W %*%
    t(Sigma_xy) %*%
    debiased_Sigma_xx_inv %*%
    Sigma_xy %*%
    sqrt_W

  eigenvectors <- eigen(target_matrix)$vectors[
    ,
    seq_len(r),
    drop = FALSE
  ]

  A <- sqrt_W_inv %*% eigenvectors

  B <- crossprod(
    eigenvectors,
    sqrt_W %*%
      t(Sigma_xy) %*%
      debiased_Sigma_xx_inv
  )

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}
