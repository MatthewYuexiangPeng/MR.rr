#' Fit multiple MR-rr estimators to the same prepared data
#'
#' Runs a selected collection of MR-rr estimators using one prepared dataset
#' and one common working rank. When the rank is not supplied, rank selection
#' is performed once and the selected rank is reused by every estimator.
#'
#' @inheritParams fit_mr_rr
#' @param methods A non-empty character vector containing any of `"naive"`,
#'   `"corrected"`, `"regularized"`, and `"sparse"`. Each method may appear
#'   at most once.
#'
#' @details
#' This function is a batch interface around [fit_mr_rr()]. It does not
#' implement separate estimator algorithms. Consequently, method-specific
#' results and controls have the same structure as the corresponding
#' [fit_mr_rr()] result.
#'
#' If `rank` is `NULL`, the first requested fit performs rank selection. The
#' resulting rank and rank-selection object are then attached to every fit in
#' the collection, so all estimators are directly comparable at the same rank.
#'
#' @return An object of class `mr_rr_fit_collection`. Requested estimator fits
#'   are available directly by method name, for example `result$corrected` and
#'   `result$sparse`. The object also contains:
#'
#'   - `methods`: The fitted methods in their requested order.
#'   - `rank`: The common working rank.
#'   - `rank_selected`: Whether the rank was selected automatically.
#'   - `rank_selection`: The shared rank-selection result, or `NULL` when rank
#'     was supplied.
#'   - `data`: The prepared `mr_rr_data` object.
#'   - `call`: The matched function call.
#'
#' @examples
#' \dontrun{
#' fits <- fit_all_mr_rr(
#'   data = prepared,
#'   rank = NULL,
#'   regularization_rate = 1e-13,
#'   sparse_lambda = 1e-3
#' )
#'
#' fits$corrected$AB
#' fits$regularized$AB
#' fits$sparse$AB
#' }
#'
#' @export
fit_all_mr_rr <- function(
    data,
    methods = c("naive", "corrected", "regularized", "sparse"),
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
  allowed_methods <- c("naive", "corrected", "regularized", "sparse")

  if (!is.character(methods) ||
      length(methods) < 1L ||
      anyNA(methods) ||
      any(!nzchar(methods))) {
    stop("`methods` must be a non-empty character vector.", call. = FALSE)
  }

  unknown_methods <- setdiff(methods, allowed_methods)
  if (length(unknown_methods) > 0L) {
    stop(
      paste0(
        "Unknown method", if (length(unknown_methods) > 1L) "s" else "",
        ": ", paste(unknown_methods, collapse = ", "), "."
      ),
      call. = FALSE
    )
  }

  if (anyDuplicated(methods)) {
    stop("`methods` must not contain duplicates.", call. = FALSE)
  }

  regularized_implementation <- match.arg(regularized_implementation)

  rank_selected <- is.null(rank)
  common_rank <- rank
  common_rank_selection <- NULL

  estimator_fits <- vector("list", length(methods))
  names(estimator_fits) <- methods

  for (i in seq_along(methods)) {
    rank_argument <- if (rank_selected && i == 1L) NULL else common_rank

    current_fit <- fit_mr_rr(
      data = data,
      method = methods[i],
      rank = rank_argument,
      rank_alpha = rank_alpha,
      rank_min = rank_min,
      regularization_rate = regularization_rate,
      regularized_implementation = regularized_implementation,
      sparse_lambda = sparse_lambda,
      sparse_refit = sparse_refit,
      sparse_threshold = sparse_threshold,
      sparse_max_iter = sparse_max_iter,
      sparse_tol = sparse_tol,
      sparse_solver = sparse_solver
    )

    if (rank_selected && i == 1L) {
      common_rank <- current_fit$rank
      common_rank_selection <- current_fit$rank_selection
    }

    if (rank_selected) {
      current_fit$rank_selected <- TRUE
      current_fit$rank_selection <- common_rank_selection
    }

    estimator_fits[[i]] <- current_fit
  }

  structure(
    c(
      estimator_fits,
      list(
        methods = methods,
        rank = as.integer(common_rank),
        rank_selected = rank_selected,
        rank_selection = common_rank_selection,
        data = data,
        call = match.call()
      )
    ),
    class = "mr_rr_fit_collection"
  )
}
