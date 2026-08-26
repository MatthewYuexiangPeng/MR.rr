#' Validate and standardize a numeric matrix used by MR-rr
#'
#' @noRd
.mr_rr_input_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name), call. = FALSE)
  }

  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(sprintf("`%s` must have at least one row and one column.", name), call. = FALSE)
  }

  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name), call. = FALSE)
  }

  storage.mode(x) <- "double"
  x
}


#' Validate a trait-correlation matrix used by MR-rr
#'
#' @noRd
.mr_rr_correlation_matrix <- function(correlation, dimension, name) {
  if (is.null(correlation)) {
    return(diag(dimension))
  }

  correlation <- .mr_rr_input_matrix(correlation, name)

  if (!identical(dim(correlation), c(dimension, dimension))) {
    stop(
      sprintf("`%s` must be a %d by %d matrix.", name, dimension, dimension),
      call. = FALSE
    )
  }

  symmetry_error <- max(abs(correlation - t(correlation)))
  if (symmetry_error > 1e-8) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }

  if (max(abs(diag(correlation) - 1)) > 1e-8) {
    stop(sprintf("`%s` must have ones on its diagonal.", name), call. = FALSE)
  }

  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1

  correlation_eigenvalues <- eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  if (min(correlation_eigenvalues) < -1e-8) {
    stop(
      sprintf("`%s` must be positive semidefinite.", name),
      call. = FALSE
    )
  }

  correlation
}


#' Validate trait names on a square matrix
#'
#' @noRd
.mr_rr_square_matrix_names <- function(x, trait_names, name) {
  if (!is.null(rownames(x)) && !identical(rownames(x), trait_names)) {
    stop(
      sprintf("Row names of `%s` must match the corresponding trait order.", name),
      call. = FALSE
    )
  }

  if (!is.null(colnames(x)) && !identical(colnames(x), trait_names)) {
    stop(
      sprintf("Column names of `%s` must match the corresponding trait order.", name),
      call. = FALSE
    )
  }

  invisible(TRUE)
}


#' Prepare summary-level data for MR-rr estimators
#'
#' Validates aligned summary-association matrices and constructs the average
#' exposure and outcome measurement-error covariance matrices used by MR-rr.
#'
#' @param beta_exposure A numeric matrix with one row per instrument and one
#'   column per exposure.
#' @param se_exposure A positive numeric matrix with the same dimensions and
#'   ordering as `beta_exposure`.
#' @param beta_outcome A numeric matrix with one row per instrument and one
#'   column per outcome.
#' @param se_outcome A positive numeric matrix with the same dimensions and
#'   ordering as `beta_outcome`.
#' @param cor_exposure An optional exposure-trait correlation matrix. The
#'   identity matrix is used when this argument is `NULL`.
#' @param cor_outcome An optional outcome-trait correlation matrix. The
#'   identity matrix is used when this argument is `NULL`.
#' @param variant_variance An optional positive numeric vector with one entry
#'   per instrument. When supplied, all four association and standard-error
#'   matrices are multiplied rowwise by `sqrt(variant_variance)`. For an
#'   additive genotype coded as 0, 1, and 2, this is commonly
#'   `2 * MAF * (1 - MAF)`.
#' @param W An optional positive-definite outcome weight matrix. When `NULL`,
#'   the function uses the inverse of the constructed outcome measurement-error
#'   covariance matrix.
#'
#' @details
#' All four input matrices must contain the same instruments in the same row
#' order and must already be harmonized to the same effect alleles. This
#' function does not perform allele harmonization, LD clumping, or instrument
#' selection.
#'
#' The average exposure covariance is
#' `Sigma_X = mean_j(D_Xj cor_exposure D_Xj)`, where `D_Xj` is the diagonal
#' matrix of exposure standard errors for instrument `j`. `Sigma_Y` is
#' constructed analogously from the outcome standard errors.
#'
#' @return An object of class `mr_rr_data`, represented by a list containing
#'   `X`, `Y`, `se_X`, `se_Y`, `Sigma_X`, `Sigma_Y`, `W`, `cor_X`, `cor_Y`,
#'   `variant_variance`, and dimension metadata.
#'
#' @export
prepare_mr_rr_data <- function(
    beta_exposure,
    se_exposure,
    beta_outcome,
    se_outcome,
    cor_exposure = NULL,
    cor_outcome = NULL,
    variant_variance = NULL,
    W = NULL) {
  X <- .mr_rr_input_matrix(beta_exposure, "beta_exposure")
  se_X <- .mr_rr_input_matrix(se_exposure, "se_exposure")
  Y <- .mr_rr_input_matrix(beta_outcome, "beta_outcome")
  se_Y <- .mr_rr_input_matrix(se_outcome, "se_outcome")

  if (!identical(dim(X), dim(se_X))) {
    stop(
      "`se_exposure` must have the same dimensions as `beta_exposure`.",
      call. = FALSE
    )
  }

  if (!identical(dim(Y), dim(se_Y))) {
    stop(
      "`se_outcome` must have the same dimensions as `beta_outcome`.",
      call. = FALSE
    )
  }

  if (nrow(X) != nrow(Y)) {
    stop(
      "Exposure and outcome matrices must contain the same number of instruments.",
      call. = FALSE
    )
  }

  if (any(se_X <= 0)) {
    stop("`se_exposure` must contain only positive values.", call. = FALSE)
  }

  if (any(se_Y <= 0)) {
    stop("`se_outcome` must contain only positive values.", call. = FALSE)
  }

  input_row_names <- list(
    rownames(X),
    rownames(se_X),
    rownames(Y),
    rownames(se_Y)
  )
  has_row_names <- !vapply(input_row_names, is.null, logical(1))

  if (any(has_row_names) && !all(has_row_names)) {
    stop(
      "Either all four input matrices must have instrument row names or none may have them.",
      call. = FALSE
    )
  }

  if (all(has_row_names)) {
    reference_row_names <- input_row_names[[1L]]
    aligned <- vapply(
      input_row_names[-1L],
      identical,
      logical(1),
      y = reference_row_names
    )

    if (!all(aligned)) {
      stop(
        "Instrument row names and ordering must agree across all four input matrices.",
        call. = FALSE
      )
    }
  }

  if (!is.null(colnames(X)) &&
      !is.null(colnames(se_X)) &&
      !identical(colnames(X), colnames(se_X))) {
    stop(
      "Exposure column names and ordering must agree between effects and standard errors.",
      call. = FALSE
    )
  }

  if (!is.null(colnames(Y)) &&
      !is.null(colnames(se_Y)) &&
      !identical(colnames(Y), colnames(se_Y))) {
    stop(
      "Outcome column names and ordering must agree between effects and standard errors.",
      call. = FALSE
    )
  }

  n <- nrow(X)
  px <- ncol(X)
  py <- ncol(Y)

  exposure_names <- colnames(X)
  if (is.null(exposure_names)) {
    exposure_names <- paste0("exposure_", seq_len(px))
  }

  outcome_names <- colnames(Y)
  if (is.null(outcome_names)) {
    outcome_names <- paste0("outcome_", seq_len(py))
  }

  colnames(X) <- exposure_names
  colnames(se_X) <- exposure_names
  colnames(Y) <- outcome_names
  colnames(se_Y) <- outcome_names

  cor_X <- .mr_rr_correlation_matrix(
    cor_exposure,
    dimension = px,
    name = "cor_exposure"
  )
  cor_Y <- .mr_rr_correlation_matrix(
    cor_outcome,
    dimension = py,
    name = "cor_outcome"
  )

  .mr_rr_square_matrix_names(cor_X, exposure_names, "cor_exposure")
  .mr_rr_square_matrix_names(cor_Y, outcome_names, "cor_outcome")

  dimnames(cor_X) <- list(exposure_names, exposure_names)
  dimnames(cor_Y) <- list(outcome_names, outcome_names)

  if (!is.null(variant_variance)) {
    if (!is.numeric(variant_variance) ||
        length(variant_variance) != n ||
        anyNA(variant_variance) ||
        any(!is.finite(variant_variance)) ||
        any(variant_variance <= 0)) {
      stop(
        "`variant_variance` must contain one finite positive value per instrument.",
        call. = FALSE
      )
    }

    if (!is.null(names(variant_variance)) &&
        all(has_row_names) &&
        !identical(names(variant_variance), rownames(X))) {
      stop(
        "Names of `variant_variance` must match the instrument row order.",
        call. = FALSE
      )
    }

    variant_variance <- as.numeric(variant_variance)

    row_scale <- sqrt(variant_variance)
    X <- X * row_scale
    se_X <- se_X * row_scale
    Y <- Y * row_scale
    se_Y <- se_Y * row_scale
  }

  Sigma_X_sum <- matrix(0, nrow = px, ncol = px)
  Sigma_Y_sum <- matrix(0, nrow = py, ncol = py)

  for (j in seq_len(n)) {
    D_X_j <- diag(se_X[j, ], nrow = px, ncol = px)
    D_Y_j <- diag(se_Y[j, ], nrow = py, ncol = py)

    Sigma_X_sum <- Sigma_X_sum + D_X_j %*% cor_X %*% D_X_j
    Sigma_Y_sum <- Sigma_Y_sum + D_Y_j %*% cor_Y %*% D_Y_j
  }

  Sigma_X <- Sigma_X_sum / n
  Sigma_Y <- Sigma_Y_sum / n
  Sigma_X <- (Sigma_X + t(Sigma_X)) / 2
  Sigma_Y <- (Sigma_Y + t(Sigma_Y)) / 2

  dimnames(Sigma_X) <- list(exposure_names, exposure_names)
  dimnames(Sigma_Y) <- list(outcome_names, outcome_names)

  if (is.null(W)) {
    W <- tryCatch(
      solve(Sigma_Y),
      error = function(e) {
        stop(
          "The constructed `Sigma_Y` is singular; supply a positive-definite `W` explicitly.",
          call. = FALSE
        )
      }
    )
  } else {
    W <- .mr_rr_input_matrix(W, "W")
    .mr_rr_square_matrix_names(W, outcome_names, "W")
  }

  validated <- .validate_estimator_inputs(
    Y = Y,
    X = X,
    r = 1L,
    Sigma_X = Sigma_X,
    W = W,
    W_inv = NULL,
    require_sigma_x = TRUE
  )

  X <- validated$X
  Y <- validated$Y
  Sigma_X <- validated$Sigma_X
  W <- validated$W
  dimnames(W) <- list(outcome_names, outcome_names)

  structure(
    list(
      X = X,
      Y = Y,
      se_X = se_X,
      se_Y = se_Y,
      Sigma_X = Sigma_X,
      Sigma_Y = Sigma_Y,
      W = W,
      cor_X = cor_X,
      cor_Y = cor_Y,
      variant_variance = variant_variance,
      n_instruments = n,
      n_exposures = px,
      n_outcomes = py
    ),
    class = "mr_rr_data"
  )
}
