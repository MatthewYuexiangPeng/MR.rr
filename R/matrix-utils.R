# Internal matrix utilities ------------------------------------------------

# Compute the uncentered second-moment matrix.
.second_moment <- function(x) {
  crossprod(x) / nrow(x)
}


# Check whether a matrix is positive semidefinite.
.is_psd <- function(x, tol = 1e-8) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != ncol(x)) {
    return(FALSE)
  }

  x <- (x + t(x)) / 2

  eigenvalues <- eigen(
    x,
    symmetric = TRUE,
    only.values = TRUE
  )$values

  min(eigenvalues) >= -tol
}


# Compute a symmetric matrix square root or inverse square root.
.sqrt_matrix <- function(x, inv = FALSE) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != ncol(x)) {
    stop("x must be a square matrix.", call. = FALSE)
  }

  eig <- eigen(x)

  transformed_values <- if (inv) {
    1 / sqrt(eig$values)
  } else {
    sqrt(eig$values)
  }

  eig$vectors %*%
    diag(
      transformed_values,
      nrow = length(transformed_values),
      ncol = length(transformed_values)
    ) %*%
    t(eig$vectors)
}


# Replace eigenvalues below epsilon to obtain a positive-definite matrix.
.nearest_psd <- function(x, epsilon = 1e-6) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }

  if (nrow(x) != ncol(x)) {
    stop("x must be a square matrix.", call. = FALSE)
  }

  x <- (x + t(x)) / 2
  eig <- eigen(x, symmetric = TRUE)
  eigenvalues <- pmax(eig$values, epsilon)

  result <- eig$vectors %*%
    (eigenvalues * t(eig$vectors))

  (result + t(result)) / 2
}
