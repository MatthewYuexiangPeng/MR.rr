.regularized_inverse <- function(S, phi) {
  S <- (S + t(S)) / 2

  if (phi == 0) {
    return(solve(S))
  }

  eig <- eigen(S, symmetric = TRUE)
  filtered <- eig$values / (eig$values^2 + phi)

  eig$vectors %*%
    (filtered * t(eig$vectors))
}


#' Spectrally regularized reduced-rank MR estimator
#'
#' Estimates a low-rank causal effect matrix after correcting the exposure
#' second-moment matrix for measurement error and applying spectral
#' regularization.
#'
#' @inheritParams mr_rr
#' @param regularization_rate A single non-negative number controlling the
#'   amount of spectral regularization.
#' @param implementation Numerical implementation used to evaluate the
#'   regularized inverse. `"spectral"` uses the stable eigenspectral formula
#'   and is the default. `"legacy"` reproduces the implementation used for
#'   the paper analyses.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{An outcome loading matrix with dimensions `ncol(Y)` by `r`.}
#'   \item{B}{An exposure coefficient matrix with dimensions `r` by `ncol(X)`.}
#'   \item{AB}{The estimated causal effect matrix `A %*% B`.}
#' }
#'
#' @export
mr_rr_regularized <- function(
    Y,
    X,
    r,
    Sigma_X,
    regularization_rate = 1e-13,
    W = NULL,
    W_inv = NULL,
    implementation = c("spectral", "legacy")) {
  if (!is.numeric(regularization_rate) ||
      length(regularization_rate) != 1L ||
      is.na(regularization_rate) ||
      !is.finite(regularization_rate) ||
      regularization_rate < 0) {
    stop(
      "`regularization_rate` must be a single finite non-negative number.",
      call. = FALSE
    )
  }

  implementation <- match.arg(implementation)

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

  if (implementation == "legacy") {
    debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)

    stable_Sigma_xx <-
      debiased_Sigma_xx +
      regularization_rate * debiased_Sigma_xx_inv

    stable_Sigma_xx_inv <- solve(stable_Sigma_xx)
  } else {
    stable_Sigma_xx_inv <- .regularized_inverse(
      debiased_Sigma_xx,
      regularization_rate
    )
  }

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

  target_matrix <-
    sqrt_W %*%
    t(Sigma_xy) %*%
    stable_Sigma_xx_inv %*%
    Sigma_xy %*%
    sqrt_W

  eigenvectors <- eigen(target_matrix)$vectors[
    , seq_len(r), drop = FALSE
  ]

  A <- sqrt_W_inv %*% eigenvectors

  B <- crossprod(
    eigenvectors,
    sqrt_W %*%
      t(Sigma_xy) %*%
      stable_Sigma_xx_inv
  )

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}
