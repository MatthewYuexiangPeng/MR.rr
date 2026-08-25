#' Post-selection refit for sparse reduced-rank MR
#'
#' Refits a sparse reduced-rank Mendelian randomization model on a fixed
#' support without the L1 penalty used during support selection.
#'
#' @inheritParams mr_rr_sparse
#' @param support A logical or binary numeric matrix with dimensions `r` by
#'   `ncol(X)`. Nonzero entries identify the coefficients of `B` that are
#'   included in the refit.
#' @param A_init An optional numeric matrix with dimensions `ncol(Y)` by `r`
#'   used to initialize `A`. `A_init` and `B_init` must either both be supplied
#'   or both be `NULL`.
#' @param B_init An optional numeric matrix with dimensions `r` by `ncol(X)`
#'   used together with `A_init`. Entries outside `support` are set to zero.
#'
#' @details
#' The support is held fixed throughout the refit. Conditional on this support,
#' the function alternates between an unpenalized least-squares update of `B`
#' and the same weighted Procrustes update of `A` used by `mr_rr_sparse()`.
#'
#' For the paper analysis, `support` is selected once on the original data and
#' then held fixed across bootstrap samples. A typical support is
#' `sparse_fit$B != 0`.
#'
#' @return A list containing:
#' \describe{
#'   \item{A}{An outcome loading matrix with dimensions `ncol(Y)` by `r`.}
#'   \item{B}{A refitted exposure coefficient matrix with dimensions `r` by
#'     `ncol(X)`, with exact zeros outside `support`.}
#'   \item{AB}{The refitted causal effect matrix `A \%*\% B`.}
#'   \item{support}{The logical support matrix used in the refit.}
#'   \item{iter}{The number of alternating optimization iterations performed.}
#'   \item{dist}{The final relative convergence distance.}
#'   \item{converged}{Whether `dist` was smaller than `tol`.}
#' }
#'
#' @export
mr_rr_sparse_refit <- function(
    Y,
    X,
    Sigma_X,
    support,
    W = NULL,
    A_init = NULL,
    B_init = NULL,
    max_iter = 100L,
    tol = 1e-2) {
  if (!is.matrix(support) || length(dim(support)) != 2L) {
    stop("`support` must be a matrix.", call. = FALSE)
  }

  if (!is.logical(support) && !is.numeric(support)) {
    stop("`support` must be logical or numeric.", call. = FALSE)
  }

  if (anyNA(support)) {
    stop("`support` must not contain missing values.", call. = FALSE)
  }

  if (is.numeric(support) && any(!support %in% c(0, 1))) {
    stop("A numeric `support` matrix must contain only 0 and 1.", call. = FALSE)
  }

  support <- support != 0
  r <- nrow(support)

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

  px <- ncol(X)
  py <- ncol(Y)

  if (!identical(dim(support), c(r, px))) {
    stop(
      "`support` must have dimensions `r` by `ncol(X)`.",
      call. = FALSE
    )
  }

  if (is.null(W)) {
    W <- diag(py)
  }

  if (xor(is.null(A_init), is.null(B_init))) {
    stop(
      "`A_init` and `B_init` must either both be supplied or both be `NULL`.",
      call. = FALSE
    )
  }

  if (!is.null(A_init)) {
    if (!is.numeric(A_init) ||
        !is.matrix(A_init) ||
        !identical(dim(A_init), c(py, r)) ||
        anyNA(A_init) ||
        any(!is.finite(A_init))) {
      stop(
        "`A_init` must be a finite numeric matrix with dimensions `ncol(Y)` by `r`.",
        call. = FALSE
      )
    }

    if (!is.numeric(B_init) ||
        !is.matrix(B_init) ||
        !identical(dim(B_init), c(r, px)) ||
        anyNA(B_init) ||
        any(!is.finite(B_init))) {
      stop(
        "`B_init` must be a finite numeric matrix with dimensions `r` by `ncol(X)`.",
        call. = FALSE
      )
    }
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

  gamma_tilde <- .construct_gamma_tilde(Y, X, Sigma_X)
  W_sqrt <- .sqrt_matrix(W)

  if (is.null(A_init)) {
    init_result <- mr_rr(
      Y = Y,
      X = X,
      r = r,
      Sigma_X = Sigma_X,
      W = W
    )
    A_hat <- init_result$A
    B_hat <- init_result$B
  } else {
    A_hat <- A_init
    B_hat <- B_init
  }

  B_hat[!support] <- 0

  dist <- Inf
  converged <- FALSE

  for (iter in seq_len(max_iter)) {
    B_hat <- matrix(0, nrow = r, ncol = px)

    for (k in seq_len(r)) {
      selected_idx <- which(support[k, ])

      if (length(selected_idx) > 0L) {
        latent_response <- as.vector(Y %*% W %*% A_hat[, k])

        B_hat[k, selected_idx] <- tryCatch(
          qr.solve(
            gamma_tilde[, selected_idx, drop = FALSE],
            latent_response
          ),
          error = function(e) {
            stop(
              sprintf(
                "The unpenalized B update failed for latent component %d.",
                k
              ),
              call. = FALSE
            )
          }
        )
      }
    }

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

  B_hat[!support] <- 0
  AB_hat <- A_hat %*% B_hat

  list(
    A = A_hat,
    B = B_hat,
    AB = AB_hat,
    support = support,
    iter = iter,
    dist = dist,
    converged = converged
  )
}
