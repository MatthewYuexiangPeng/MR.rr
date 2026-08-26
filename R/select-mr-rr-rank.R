#' Select a working rank for reduced-rank MR
#'
#' Applies the sequential chi-square diagnostic used in the paper to select
#' the smallest working rank that is not rejected.
#'
#' @inheritParams mr_rr
#' @param alpha A single number strictly between zero and one giving the
#'   sequential testing level.
#' @param min_rank A non-negative integer giving the smallest candidate rank.
#'   The default is one because the package estimators require at least one
#'   latent pathway.
#'
#' @details
#' For each candidate rank `r`, the function tests the working null hypothesis
#' that the causal effect matrix has rank at most `r`. It selects the first
#' candidate whose p-value is at least `alpha`; if every candidate is rejected,
#' it returns the full rank `min(ncol(X), ncol(Y))`.
#'
#' The chi-square approximation is a working rank-selection diagnostic. As
#' discussed in the paper, it is not guaranteed to provide a formal test in
#' the weak-instrument regime.
#'
#' @return An object of class `mr_rr_rank_selection`, represented by a list
#'   containing:
#' \describe{
#'   \item{selected_rank}{The selected working rank.}
#'   \item{candidate_ranks}{The ranks tested sequentially.}
#'   \item{statistic}{The test statistic for each candidate rank.}
#'   \item{df}{The chi-square degrees of freedom for each candidate rank.}
#'   \item{p_value}{The p-value for each candidate rank.}
#'   \item{eigenvalues}{The eigenvalues used by the diagnostic.}
#'   \item{alpha}{The testing level supplied by the user.}
#'   \item{min_rank}{The smallest candidate rank.}
#'   \item{full_rank}{The maximum possible rank.}
#' }
#'
#' @export
select_mr_rr_rank <- function(
    Y,
    X,
    Sigma_X,
    W = NULL,
    alpha = 0.05,
    min_rank = 1L) {
  inputs <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = 1L,
    Sigma_X = Sigma_X,
    W = W,
    W_inv = NULL,
    require_sigma_x = TRUE
  )

  Y <- inputs$Y
  X <- inputs$X
  Sigma_X <- inputs$Sigma_X
  W <- inputs$W

  n <- nrow(Y)
  px <- ncol(X)
  py <- ncol(Y)
  full_rank <- min(px, py)

  if (is.null(W)) {
    W <- diag(py)
  }

  if (!is.numeric(alpha) ||
      length(alpha) != 1L ||
      is.na(alpha) ||
      !is.finite(alpha) ||
      alpha <= 0 ||
      alpha >= 1) {
    stop("`alpha` must be a single number strictly between 0 and 1.", call. = FALSE)
  }

  if (!is.numeric(min_rank) ||
      length(min_rank) != 1L ||
      is.na(min_rank) ||
      !is.finite(min_rank) ||
      min_rank < 0 ||
      min_rank != floor(min_rank) ||
      min_rank > full_rank) {
    stop(
      "`min_rank` must be an integer between 0 and `min(ncol(X), ncol(Y))`.",
      call. = FALSE
    )
  }
  min_rank <- as.integer(min_rank)

  scale_factor <- n - (px + py + 1) / 2
  if (scale_factor <= 0) {
    stop(
      "The sample size is too small for the rank-selection diagnostic.",
      call. = FALSE
    )
  }

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx <-
    (debiased_Sigma_xx + t(debiased_Sigma_xx)) / 2

  xx_eigenvalues <- eigen(
    debiased_Sigma_xx,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (min(xx_eigenvalues) <= 0) {
    stop(
      "The measurement-error-corrected exposure matrix must be positive definite for rank selection.",
      call. = FALSE
    )
  }

  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)
  W_sqrt <- .sqrt_matrix(W)

  diagnostic_matrix <-
    W_sqrt %*%
    t(Sigma_xy) %*%
    debiased_Sigma_xx_inv %*%
    Sigma_xy %*%
    W_sqrt
  diagnostic_matrix <-
    (diagnostic_matrix + t(diagnostic_matrix)) / 2

  diagnostic_eigenvalues <- eigen(
    diagnostic_matrix,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  eigen_tolerance <- 1e-8 * max(1, max(abs(diagnostic_eigenvalues)))
  if (min(diagnostic_eigenvalues) < -eigen_tolerance) {
    stop(
      "The rank-selection diagnostic matrix is not positive semidefinite.",
      call. = FALSE
    )
  }

  diagnostic_eigenvalues <- pmax(diagnostic_eigenvalues, 0)
  eigenvalues <- diagnostic_eigenvalues[seq_len(full_rank)]

  if (min_rank < full_rank) {
    candidate_ranks <- seq.int(min_rank, full_rank - 1L)

    statistic <- vapply(candidate_ranks, function(rank) {
      tail_indices <- seq.int(rank + 1L, full_rank)
      scale_factor * sum(log1p(eigenvalues[tail_indices]))
    }, numeric(1))

    df <- (py - candidate_ranks) * (px - candidate_ranks)
    p_value <- stats::pchisq(
      statistic,
      df = df,
      lower.tail = FALSE
    )

    accepted <- which(p_value >= alpha)
    if (length(accepted) == 0L) {
      selected_rank <- full_rank
    } else {
      selected_rank <- candidate_ranks[accepted[1L]]
    }
  } else {
    candidate_ranks <- integer(0)
    statistic <- numeric(0)
    df <- numeric(0)
    p_value <- numeric(0)
    selected_rank <- full_rank
  }

  structure(
    list(
      selected_rank = as.integer(selected_rank),
      candidate_ranks = as.integer(candidate_ranks),
      statistic = statistic,
      df = df,
      p_value = p_value,
      eigenvalues = eigenvalues,
      alpha = alpha,
      min_rank = min_rank,
      full_rank = as.integer(full_rank)
    ),
    class = "mr_rr_rank_selection"
  )
}
