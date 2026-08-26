#' Print an MR-rr fit
#'
#' Displays the fitted method, working rank, matrix dimensions, relevant
#' method controls, and estimated causal effect matrix without printing the
#' complete prepared dataset stored inside the fit object.
#'
#' @param x An object of class `mr_rr_fit`.
#' @param digits Number of significant digits used to print estimates.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.mr_rr_fit <- function(x, digits = 4L, ...) {
  rank_source <- if (isTRUE(x$rank_selected)) "selected" else "supplied"

  cat("MR-rr fit\n")
  cat("  Method: ", x$method, "\n", sep = "")
  cat(
    "  Rank:   ", x$rank, " (", rank_source, ")\n",
    sep = ""
  )
  cat(
    "  Effect: ", nrow(x$AB), " outcomes x ", ncol(x$AB),
    " exposures\n",
    sep = ""
  )

  if (identical(x$method, "regularized")) {
    cat(
      "  Regularization: rate = ",
      format(x$control$regularization_rate, digits = digits),
      ", implementation = ", x$control$implementation, "\n",
      sep = ""
    )
  }

  if (identical(x$method, "sparse")) {
    cat(
      "  Sparse refit: ", if (isTRUE(x$sparse_refitted)) "yes" else "no",
      "\n",
      sep = ""
    )
    cat(
      "  Nonzero B:    ", sum(x$B != 0), " / ", length(x$B), "\n",
      sep = ""
    )
  }

  cat("\nEstimated causal effect matrix (AB):\n")
  print(x$AB, digits = digits)

  invisible(x)
}


#' Print a collection of MR-rr fits
#'
#' Displays the common rank, data dimensions, and a compact overview of the
#' fitted estimators without printing the complete stored dataset.
#'
#' @param x An object of class `mr_rr_fit_collection`.
#' @param ... Additional arguments. Currently ignored.
#'
#' @return `x`, invisibly.
#'
#' @export
print.mr_rr_fit_collection <- function(x, ...) {
  rank_source <- if (isTRUE(x$rank_selected)) "selected" else "supplied"

  overview <- data.frame(
    method = x$methods,
    rank = rep(x$rank, length(x$methods)),
    nonzero_B = vapply(
      x$methods,
      function(method) sum(x[[method]]$B != 0),
      integer(1)
    ),
    stringsAsFactors = FALSE
  )

  cat("MR-rr fit collection\n")
  cat("  Methods: ", paste(x$methods, collapse = ", "), "\n", sep = "")
  cat(
    "  Rank:    ", x$rank, " (", rank_source, ")\n",
    sep = ""
  )
  cat(
    "  Data:    ", x$data$n_instruments, " instruments, ",
    x$data$n_exposures, " exposures, ",
    x$data$n_outcomes, " outcomes\n",
    sep = ""
  )
  cat("\nEstimator overview:\n")
  print(overview, row.names = FALSE)

  invisible(x)
}
