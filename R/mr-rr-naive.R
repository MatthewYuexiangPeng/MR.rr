#' Naive Reduced-Rank Mendelian Randomization
#'
#' Estimates a rank-constrained causal effect matrix without correcting for
#' measurement error in the estimated SNP-exposure associations.
#'
#' @param Y A numeric matrix with one row per genetic instrument and one column
#'   per outcome.
#' @param X A numeric matrix with one row per genetic instrument and one column
#'   per exposure.
#' @param r A positive integer specifying the working rank of the causal effect
#'   matrix.
#' @param W An optional positive-definite outcome weight matrix. If `NULL`, the
#'   identity matrix is used.
#' @param W_inv An optional precomputed inverse of `W`. This argument can reduce
#'   repeated matrix inversion when both matrices are already available.
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
#' fit <- mr_rr_naive(Y = Y, X = X, r = 1)
#' fit$AB
#'
#' @export
mr_rr_naive <- function(Y, X, r, W = NULL, W_inv = NULL) {
  inputs <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = r,
    W = W,
    W_inv = W_inv
  )

  Y <- inputs$Y
  X <- inputs$X
  r <- inputs$r
  W <- inputs$W
  W_inv <- inputs$W_inv

  XtY <- crossprod(X, Y)
  XtX_inv <- solve(crossprod(X))

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
    t(XtY) %*%
    XtX_inv %*%
    XtY %*%
    sqrt_W / nrow(Y)

  eigenvectors <- eigen(target_matrix)$vectors[
    ,
    seq_len(r),
    drop = FALSE
  ]

  A <- sqrt_W_inv %*% eigenvectors

  B <- crossprod(
    eigenvectors,
    sqrt_W %*% t(XtY) %*% XtX_inv
  )

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}
