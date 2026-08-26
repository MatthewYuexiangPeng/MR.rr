#' Add informative dimension names to MR-rr estimates
#'
#' @noRd
.label_mr_rr_matrices <- function(fit, data, rank) {
  component_names <- paste0("component_", seq_len(rank))
  exposure_names <- colnames(data$X)
  outcome_names <- colnames(data$Y)

  dimnames(fit$A) <- list(outcome_names, component_names)
  dimnames(fit$B) <- list(component_names, exposure_names)
  dimnames(fit$AB) <- list(outcome_names, exposure_names)

  if (!is.null(fit$B_raw)) {
    dimnames(fit$B_raw) <- list(component_names, exposure_names)
  }

  if (!is.null(fit$AB_raw)) {
    dimnames(fit$AB_raw) <- list(outcome_names, exposure_names)
  }

  if (!is.null(fit$support)) {
    dimnames(fit$support) <- list(component_names, exposure_names)
  }

  fit
}


#' Fit an MR-rr estimator to prepared summary data
#'
#' Provides a common user-facing interface to the naive, measurement-error-
#' corrected, spectrally regularized, and sparse MR-rr estimators.
#'
#' @param data An `mr_rr_data` object created by [prepare_mr_rr_data()].
#' @param method The estimator to fit. One of `"corrected"`, `"naive"`,
#'   `"regularized"`, or `"sparse"`.
#' @param rank A positive integer giving the working rank. When `NULL`, the
#'   rank is selected by [select_mr_rr_rank()].
#' @param rank_alpha The sequential testing level passed to
#'   [select_mr_rr_rank()] when `rank` is `NULL`.
#' @param rank_min The smallest candidate rank considered by
#'   [select_mr_rr_rank()] when `rank` is `NULL`. It must be at least one
#'   because the package estimators require a positive rank.
#' @param regularization_rate A non-negative regularization rate passed to
#'   [mr_rr_regularized()].
#' @param regularized_implementation Either `"spectral"` or `"legacy"`.
#'   The default uses the numerically stable spectral implementation.
#' @param sparse_lambda A non-negative scalar or exposure-specific vector
#'   passed as `lambda` to [mr_rr_sparse()].
#' @param sparse_refit A logical value. If `TRUE`, sparse estimation consists
#'   of penalized support selection followed by [mr_rr_sparse_refit()] on the
#'   selected support. The selection fit is retained in the returned object.
#' @param sparse_threshold A non-negative coefficient threshold passed to
#'   [mr_rr_sparse()].
#' @param sparse_max_iter The maximum number of alternating iterations used by
#'   sparse selection and refitting.
#' @param sparse_tol The convergence tolerance used by sparse selection and
#'   refitting.
#' @param sparse_solver The CVXR solver used by [mr_rr_sparse()].
#'
#' @details
#' `method = "corrected"` calls [mr_rr()] and represents the main MR-rr
#' estimator. `method = "regularized"` uses the stable spectral implementation
#' by default. `method = "sparse"` performs both sparse support selection and
#' the unpenalized fixed-support refit by default.
#'
#' The prepared data are stored in the returned object so that later methods,
#' including bootstrap inference, can reconstruct the same analysis without
#' requiring the user to supply the inputs again.
#'
#' @return An object of class `mr_rr_fit`, represented by a list with the
#'   following components:
#'
#'   - `method`: The fitted estimator.
#'   - `rank`: The working rank used for estimation.
#'   - `rank_selected`: Whether the rank was selected automatically.
#'   - `rank_selection`: The rank-selection result, or `NULL` when the rank
#'     was supplied.
#'   - `A`: The estimated outcome loading matrix.
#'   - `B`: The estimated exposure coefficient matrix.
#'   - `AB`: The estimated causal effect matrix.
#'   - `details`: The complete result returned by the underlying estimator.
#'   - `sparse_selection`: The penalized sparse fit when `method = "sparse"`,
#'     otherwise `NULL`.
#'   - `sparse_refitted`: Whether a sparse post-selection refit was used.
#'   - `data`: The prepared `mr_rr_data` object.
#'   - `control`: The method-specific tuning parameters.
#'   - `call`: The matched function call.
#'
#' @export
fit_mr_rr <- function(
    data,
    method = c("corrected", "naive", "regularized", "sparse"),
    rank = NULL,
    rank_alpha = 0.05,
    rank_min = 1L,
    regularization_rate = 1e-13,
    regularized_implementation = c("spectral", "legacy"),
    sparse_lambda = 1e-3,
    sparse_refit = TRUE,
    sparse_threshold = 1e-2,
    sparse_max_iter = 100L,
    sparse_tol = 1e-2,
    sparse_solver = "OSQP") {
  if (!inherits(data, "mr_rr_data") || !is.list(data)) {
    stop(
      "`data` must be an `mr_rr_data` object created by `prepare_mr_rr_data()`.",
      call. = FALSE
    )
  }

  required_components <- c("X", "Y", "Sigma_X", "W")
  if (!all(required_components %in% names(data))) {
    stop(
      "`data` is missing components required for MR-rr estimation.",
      call. = FALSE
    )
  }

  method <- match.arg(method)
  regularized_implementation <- match.arg(regularized_implementation)

  if (!is.logical(sparse_refit) ||
      length(sparse_refit) != 1L ||
      is.na(sparse_refit)) {
    stop("`sparse_refit` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  full_rank <- min(ncol(data$X), ncol(data$Y))
  rank_selected <- is.null(rank)
  rank_selection <- NULL

  if (rank_selected) {
    if (!is.numeric(rank_min) ||
        length(rank_min) != 1L ||
        is.na(rank_min) ||
        !is.finite(rank_min) ||
        rank_min < 1 ||
        rank_min != floor(rank_min) ||
        rank_min > full_rank) {
      stop(
        "`rank_min` must be an integer between 1 and the maximum possible rank.",
        call. = FALSE
      )
    }

    rank_selection <- select_mr_rr_rank(
      Y = data$Y,
      X = data$X,
      Sigma_X = data$Sigma_X,
      W = data$W,
      alpha = rank_alpha,
      min_rank = as.integer(rank_min)
    )
    rank <- rank_selection$selected_rank
  } else {
    if (!is.numeric(rank) ||
        length(rank) != 1L ||
        is.na(rank) ||
        !is.finite(rank) ||
        rank < 1 ||
        rank != floor(rank) ||
        rank > full_rank) {
      stop(
        "`rank` must be an integer between 1 and the maximum possible rank.",
        call. = FALSE
      )
    }
    rank <- as.integer(rank)
  }

  sparse_selection <- NULL
  sparse_refitted <- NULL

  estimator_fit <- switch(
    method,
    naive = mr_rr_naive(
      Y = data$Y,
      X = data$X,
      r = rank,
      W = data$W
    ),
    corrected = mr_rr(
      Y = data$Y,
      X = data$X,
      r = rank,
      Sigma_X = data$Sigma_X,
      W = data$W
    ),
    regularized = mr_rr_regularized(
      Y = data$Y,
      X = data$X,
      r = rank,
      Sigma_X = data$Sigma_X,
      regularization_rate = regularization_rate,
      W = data$W,
      implementation = regularized_implementation
    ),
    sparse = {
      sparse_selection <- mr_rr_sparse(
        Y = data$Y,
        X = data$X,
        r = rank,
        Sigma_X = data$Sigma_X,
        lambda = sparse_lambda,
        W = data$W,
        max_iter = sparse_max_iter,
        tol = sparse_tol,
        threshold = sparse_threshold,
        solver = sparse_solver
      )
      sparse_selection <- .label_mr_rr_matrices(
        sparse_selection,
        data = data,
        rank = rank
      )

      if (sparse_refit) {
        sparse_refitted <- TRUE
        mr_rr_sparse_refit(
          Y = data$Y,
          X = data$X,
          Sigma_X = data$Sigma_X,
          support = sparse_selection$B != 0,
          W = data$W,
          A_init = sparse_selection$A,
          B_init = sparse_selection$B,
          max_iter = sparse_max_iter,
          tol = sparse_tol
        )
      } else {
        sparse_refitted <- FALSE
        sparse_selection
      }
    }
  )

  estimator_fit <- .label_mr_rr_matrices(
    estimator_fit,
    data = data,
    rank = rank
  )

  control <- switch(
    method,
    naive = list(),
    corrected = list(),
    regularized = list(
      regularization_rate = regularization_rate,
      implementation = regularized_implementation
    ),
    sparse = list(
      lambda = sparse_lambda,
      refit = sparse_refit,
      threshold = sparse_threshold,
      max_iter = sparse_max_iter,
      tol = sparse_tol,
      solver = sparse_solver
    )
  )

  structure(
    list(
      method = method,
      rank = as.integer(rank),
      rank_selected = rank_selected,
      rank_selection = rank_selection,
      A = estimator_fit$A,
      B = estimator_fit$B,
      AB = estimator_fit$AB,
      details = estimator_fit,
      sparse_selection = sparse_selection,
      sparse_refitted = sparse_refitted,
      data = data,
      control = control,
      call = match.call()
    ),
    class = "mr_rr_fit"
  )
}
