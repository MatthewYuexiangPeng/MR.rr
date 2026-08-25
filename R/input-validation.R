# Internal input-validation utilities -------------------------------------

.as_numeric_matrix <- function(x, name) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (!is.matrix(x) || !is.numeric(x)) {
    stop(name, " must be a numeric matrix.", call. = FALSE)
  }

  if (nrow(x) < 1L || ncol(x) < 1L) {
    stop(name, " must have at least one row and one column.", call. = FALSE)
  }

  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }

  storage.mode(x) <- "double"
  x
}


.validate_rank <- function(r, X, Y) {
  if (
    !is.numeric(r) ||
    length(r) != 1L ||
    !is.finite(r) ||
    r != floor(r)
  ) {
    stop("r must be a single integer.", call. = FALSE)
  }

  r <- as.integer(r)
  max_rank <- min(ncol(X), ncol(Y))

  if (r < 1L || r > max_rank) {
    stop(
      "r must be between 1 and min(ncol(X), ncol(Y)).",
      call. = FALSE
    )
  }

  r
}


.validate_square_matrix <- function(
    x,
    name,
    dimension,
    symmetric = TRUE,
    positive_definite = FALSE,
    tol = 1e-10
) {
  x <- .as_numeric_matrix(x, name)

  if (!identical(dim(x), c(dimension, dimension))) {
    stop(
      name,
      " must be a ",
      dimension,
      " by ",
      dimension,
      " matrix.",
      call. = FALSE
    )
  }

  if (
    symmetric &&
    max(abs(x - t(x))) > tol
  ) {
    stop(name, " must be symmetric.", call. = FALSE)
  }

  if (positive_definite) {
    eigenvalues <- eigen(
      (x + t(x)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values

    if (min(eigenvalues) <= 0) {
      stop(name, " must be positive definite.", call. = FALSE)
    }
  }

  x
}


.validate_estimator_inputs <- function(
    Y,
    X,
    r,
    Sigma_X = NULL,
    W = NULL,
    W_inv = NULL,
    require_sigma_x = FALSE
) {
  Y <- .as_numeric_matrix(Y, "Y")
  X <- .as_numeric_matrix(X, "X")

  if (nrow(Y) != nrow(X)) {
    stop(
      "X and Y must have the same number of rows.",
      call. = FALSE
    )
  }

  r <- .validate_rank(r, X, Y)

  px <- ncol(X)
  py <- ncol(Y)

  if (require_sigma_x && is.null(Sigma_X)) {
    stop("Sigma_X must be provided.", call. = FALSE)
  }

  if (!is.null(Sigma_X)) {
    Sigma_X <- .validate_square_matrix(
      Sigma_X,
      name = "Sigma_X",
      dimension = px,
      symmetric = TRUE
    )
  }

  if (is.null(W) && !is.null(W_inv)) {
    stop(
      "W_inv cannot be supplied when W is NULL.",
      call. = FALSE
    )
  }

  if (!is.null(W)) {
    W <- .validate_square_matrix(
      W,
      name = "W",
      dimension = py,
      symmetric = TRUE,
      positive_definite = TRUE
    )
  }

  if (!is.null(W_inv)) {
    W_inv <- .validate_square_matrix(
      W_inv,
      name = "W_inv",
      dimension = py,
      symmetric = TRUE,
      positive_definite = TRUE
    )
  }

  list(
    Y = Y,
    X = X,
    r = r,
    Sigma_X = Sigma_X,
    W = W,
    W_inv = W_inv
  )
}
