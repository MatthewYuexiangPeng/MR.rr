#' Construct the surrogate exposure matrix used by sparse MR-rr
#'
#' This helper intentionally follows the transformation used by the frozen
#' paper implementation.
#'
#' @noRd
.construct_gamma_tilde <- function(Y, X, Sigma_X) {
  n <- nrow(X)

  P_Y <- tryCatch(
    Y %*% solve(t(Y) %*% Y) %*% t(Y),
    error = function(e) {
      stop(
        "Cannot construct the sparse MR-rr projection because `Y` is rank deficient.",
        call. = FALSE
      )
    }
  )

  P_Y_complement <- diag(1, n) - P_Y
  Sigma_X_hat <- t(X) %*% X / n
  projected_X <- P_Y %*% X
  matrix_part1 <-
    Sigma_X_hat - t(projected_X) %*% projected_X / n
  matrix_part2 <- matrix_part1 - Sigma_X

  if (!.is_psd(matrix_part2)) {
    matrix_part2 <- .nearest_psd(matrix_part2, epsilon = 1e-6)
  }

  R <- tryCatch(
    chol(matrix_part1),
    error = function(e) {
      stop(
        "The sparse MR-rr surrogate covariance is not positive definite.",
        call. = FALSE
      )
    }
  )

  Q <- tryCatch(
    chol(matrix_part2),
    error = function(e) {
      stop(
        "The corrected sparse MR-rr covariance is not positive definite.",
        call. = FALSE
      )
    }
  )

  L <- solve(R) %*% Q

  P_Y %*% X + (P_Y_complement %*% X) %*% L
}


#' Sparse reduced-rank Mendelian randomization estimator
#'
#' Estimates a low-rank causal effect matrix while applying column-specific
#' L1 penalties to the exposure coefficient matrix.
#'
#' @inheritParams mr_rr
#' @param lambda A non-negative number or a numeric vector of length
#'   `ncol(X)`. A scalar is recycled across exposures.
#' @param max_iter A positive integer giving the maximum number of alternating
#'   optimization iterations.
#' @param tol A positive number used as the relative convergence tolerance.
#' @param threshold A non-negative number. Entries of the final coefficient
#'   matrix with absolute value below this threshold are set to zero. The
#'   default, `1e-2`, matches the frozen paper implementation.
#' @param solver A single character string naming a CVXR solver. The default is
#'   `"OSQP"`.
#'
#' @details
#' `lambda` controls the optimization penalty, whereas `threshold` is applied
#' only after the alternating optimization has finished. The unthresholded
#' estimates are retained in `B_raw` and `AB_raw`.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{An outcome loading matrix with dimensions `ncol(Y)` by `r`.}
#'   \item{B}{The thresholded exposure coefficient matrix with dimensions
#'     `r` by `ncol(X)`.}
#'   \item{AB}{The thresholded causal effect matrix `A \%*\% B`.}
#'   \item{B_raw}{The exposure coefficient matrix before thresholding.}
#'   \item{AB_raw}{The causal effect matrix before thresholding.}
#'   \item{iter}{The number of alternating optimization iterations performed.}
#'   \item{dist}{The final relative convergence distance.}
#'   \item{converged}{Whether `dist` was smaller than `tol`.}
#' }
#'
#' @export
mr_rr_sparse <- function(
    Y,
    X,
    r,
    Sigma_X,
    lambda = rep(1e-3, ncol(X)),
    W = NULL,
    max_iter = 100L,
    tol = 1e-2,
    threshold = 1e-2,
    solver = "OSQP") {
  inputs <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = r,
    Sigma_X = Sigma_X,
    W = W,
    W_inv = NULL,
    require_sigma_x = TRUE
  )

  Y <- inputs$Y
  X <- inputs$X
  r <- inputs$r
  Sigma_X <- inputs$Sigma_X
  W <- inputs$W

  n <- nrow(X)
  px <- ncol(X)
  py <- ncol(Y)

  if (is.null(W)) {
    W <- diag(py)
  }

  if (!is.numeric(lambda) ||
      anyNA(lambda) ||
      any(!is.finite(lambda)) ||
      any(lambda < 0)) {
    stop(
      "`lambda` must contain finite non-negative numbers.",
      call. = FALSE
    )
  }

  if (length(lambda) == 1L) {
    lambda <- rep(lambda, px)
  }

  if (length(lambda) != px) {
    stop(
      "`lambda` must have length 1 or `ncol(X)`.",
      call. = FALSE
    )
  }

  if (!is.numeric(max_iter) ||
      length(max_iter) != 1L ||
      is.na(max_iter) ||
      !is.finite(max_iter) ||
      max_iter < 1 ||
      max_iter != floor(max_iter) ||
      max_iter > .Machine$integer.max) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }
  max_iter <- as.integer(max_iter)

  if (!is.numeric(tol) ||
      length(tol) != 1L ||
      is.na(tol) ||
      !is.finite(tol) ||
      tol <= 0) {
    stop("`tol` must be a single finite positive number.", call. = FALSE)
  }

  if (!is.numeric(threshold) ||
      length(threshold) != 1L ||
      is.na(threshold) ||
      !is.finite(threshold) ||
      threshold < 0) {
    stop(
      "`threshold` must be a single finite non-negative number.",
      call. = FALSE
    )
  }

  if (!is.character(solver) ||
      length(solver) != 1L ||
      is.na(solver) ||
      !nzchar(solver)) {
    stop("`solver` must be a non-empty character string.", call. = FALSE)
  }

  W_sqrt <- .sqrt_matrix(W)
  gamma_tilde <- .construct_gamma_tilde(Y, X, Sigma_X)

  init_result <- mr_rr(
    Y = Y,
    X = X,
    r = r,
    Sigma_X = Sigma_X,
    W = W
  )
  A_hat <- init_result$A

  dist <- Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    B_var <- .cvxr_matrix_variable(r, px)

    penalty <- Reduce(
      `+`,
      lapply(seq_len(px), function(k) {
        lambda[k] *
          CVXR::sum_entries(
            CVXR::norm1(B_var[, k, drop = FALSE])
          )
      })
    )

    loss <-
      CVXR::sum_squares(
        t(A_hat) %*% W %*% t(Y) - B_var %*% t(gamma_tilde)
      ) / n +
      penalty

    problem <- CVXR::Problem(CVXR::Minimize(loss))
    B_hat <- .solve_cvxr_variable(
      problem = problem,
      variable = B_var,
      solver = solver
    )
    B_hat <- matrix(B_hat, nrow = r, ncol = px)

    svd_result <- svd(
      B_hat %*% t(gamma_tilde) %*% Y %*% W_sqrt
    )

    A_hat_new <-
      solve(W_sqrt) %*%
      svd_result$v %*%
      t(svd_result$u)

    dist <-
      norm(A_hat_new %*% B_hat - A_hat %*% B_hat, "F") /
      max(1e-8, norm(A_hat %*% B_hat, "F"))

    A_hat <- A_hat_new

    if (dist < tol) {
      converged <- TRUE
      break
    }
  }

  B_raw <- B_hat
  AB_raw <- A_hat %*% B_raw

  B_hat[abs(B_hat) < threshold] <- 0
  AB_hat <- A_hat %*% B_hat

  list(
    A = A_hat,
    B = B_hat,
    AB = AB_hat,
    B_raw = B_raw,
    AB_raw = AB_raw,
    iter = iter,
    dist = dist,
    converged = converged
  )
}
