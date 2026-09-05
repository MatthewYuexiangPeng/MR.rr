#!/usr/bin/env Rscript

options(warn = 1)

repo_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

if (!file.exists(file.path(repo_root, "DESCRIPTION"))) {
  stop(
    "Run this script from the MRrr-rebuild repository root.",
    call. = FALSE
  )
}

legacy_estimator_file <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project",
  "scripts",
  "MR_rr_estimators.R"
)

if (!file.exists(legacy_estimator_file)) {
  stop(
    "Frozen estimator file was not found.",
    call. = FALSE
  )
}

cat("Loading frozen estimator implementation:\n")
cat(legacy_estimator_file, "\n\n")

legacy <- new.env(parent = globalenv())

sys.source(
  legacy_estimator_file,
  envir = legacy
)

set.seed(20260826)

n_instruments <- 120L
n_exposures <- 6L
n_outcomes <- 3L
working_rank <- 2L

A_true <- matrix(
  c(
    0.80, -0.50, 0.30,
    -0.20,  0.40, 0.70
  ),
  nrow = n_outcomes,
  ncol = working_rank
)

B_true <- matrix(
  c(
    0.30, -0.20,  0.00, 0.15, 0.00, 0.00,
    0.00,  0.10, -0.25, 0.00, 0.20, 0.00
  ),
  nrow = working_rank,
  byrow = TRUE
)

C_true <- A_true %*% B_true

X_true <- matrix(
  stats::rnorm(
    n_instruments * n_exposures,
    sd = 0.12
  ),
  nrow = n_instruments,
  ncol = n_exposures
)

exposure_error_sd <- 0.01

Sigma_X <- diag(
  exposure_error_sd^2,
  nrow = n_exposures,
  ncol = n_exposures
)

X <- X_true +
  matrix(
    stats::rnorm(
      n_instruments * n_exposures,
      sd = exposure_error_sd
    ),
    nrow = n_instruments,
    ncol = n_exposures
  )

Y <- X_true %*% t(C_true) +
  matrix(
    stats::rnorm(
      n_instruments * n_outcomes,
      sd = 0.01
    ),
    nrow = n_instruments,
    ncol = n_outcomes
  )

W <- diag(n_outcomes)

fits <- list(
  naive = legacy$mr_rr_naive(
    Y = Y,
    X = X,
    r = working_rank,
    W = W
  ),
  corrected = legacy$mr_rr(
    Y = Y,
    X = X,
    r = working_rank,
    Sigma_X = Sigma_X,
    W = W
  ),
  regularized = legacy$mr_rr_regularized(
    Y = Y,
    X = X,
    r = working_rank,
    Sigma_X = Sigma_X,
    regularization_rate = 1e-6,
    W = W
  ),
  sparse = legacy$mr_rr_sparse(
    GAMMA_hat = Y,
    gamma_hat = X,
    W = W,
    Sigma_X = Sigma_X,
    lambda = rep(2e-3, n_exposures),
    r = working_rank,
    max_iter = 10L,
    tol = 1e-3
  )
)

expected_dimensions <- list(
  A = c(n_outcomes, working_rank),
  B = c(working_rank, n_exposures),
  AB = c(n_outcomes, n_exposures)
)

check_fit <- function(fit, method) {
  for (component in names(expected_dimensions)) {
    actual_dimension <- dim(fit[[component]])
    expected_dimension <- expected_dimensions[[component]]

    if (!identical(actual_dimension, expected_dimension)) {
      stop(
        sprintf(
          "%s returned incorrect dimensions for %s.",
          method,
          component
        ),
        call. = FALSE
      )
    }

    if (!all(is.finite(fit[[component]]))) {
      stop(
        sprintf(
          "%s returned non-finite values in %s.",
          method,
          component
        ),
        call. = FALSE
      )
    }
  }

  reconstruction_error <- max(
    abs(fit$AB - fit$A %*% fit$B)
  )

  if (reconstruction_error > 1e-8) {
    stop(
      sprintf(
        "%s failed the AB = A %%*%% B check.",
        method
      ),
      call. = FALSE
    )
  }

  estimated_rank <- qr(
    fit$AB,
    tol = 1e-8
  )$rank

  if (estimated_rank > working_rank) {
    stop(
      sprintf(
        "%s returned an effect matrix above the working rank.",
        method
      ),
      call. = FALSE
    )
  }

  data.frame(
    method = method,
    rows_A = nrow(fit$A),
    columns_A = ncol(fit$A),
    rows_B = nrow(fit$B),
    columns_B = ncol(fit$B),
    rank_AB = estimated_rank,
    reconstruction_error = reconstruction_error,
    stringsAsFactors = FALSE
  )
}

check_results <- do.call(
  rbind,
  Map(
    check_fit,
    fits,
    names(fits)
  )
)

regularized_zero <- legacy$mr_rr_regularized(
  Y = Y,
  X = X,
  r = working_rank,
  Sigma_X = Sigma_X,
  regularization_rate = 0,
  W = W
)

if (!isTRUE(all.equal(
  fits$corrected$AB,
  regularized_zero$AB,
  tolerance = 1e-8,
  check.attributes = FALSE
))) {
  stop(
    paste(
      "Legacy regularized estimator with",
      "regularization_rate = 0 does not match MR-rr."
    ),
    call. = FALSE
  )
}

print(check_results, row.names = FALSE)

cat("\nLegacy regularized estimator at phi = 0: PASS\n")
cat("Frozen estimator smoke test: PASS\n")
