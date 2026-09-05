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

freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)

data_root <- file.path(freeze_root, "data")
results_root <- file.path(freeze_root, "results")

legacy_estimator_file <- file.path(
  freeze_root,
  "scripts",
  "MR_rr_estimators.R"
)

reference_result_file <- file.path(
  results_root,
  "real_result_260820_rank1_eta_1p2e-3.RData"
)

required_files <- c(
  legacy_estimator_file,
  file.path(data_root, "dat_1e-4.csv"),
  file.path(data_root, "rho_mat_1e-4.csv"),
  reference_result_file
)

missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste(
      "Required files are missing:",
      paste(missing_files, collapse = "\n"),
      sep = "\n"
    ),
    call. = FALSE
  )
}

cat("Loading frozen estimator implementation.\n")

legacy <- new.env(parent = globalenv())

sys.source(
  legacy_estimator_file,
  envir = legacy
)


# -------------------------------------------------------------------------
# Load and transform the frozen real-data inputs
# -------------------------------------------------------------------------

lip_data_real <- utils::read.csv(
  file.path(data_root, "dat_1e-4.csv")
)

lip_corr_real <- utils::read.csv(
  file.path(data_root, "rho_mat_1e-4.csv")
)

variant_variance <-
  2 *
  lip_data_real$ImpMAF *
  (1 - lip_data_real$ImpMAF)

gamma_exp <- as.matrix(
  lip_data_real[, paste0("gamma_exp", seq_len(9L))]
)

se_exp <- as.matrix(
  lip_data_real[, paste0("se_exp", seq_len(9L))]
)

gamma_out <- as.matrix(
  lip_data_real[, paste0("gamma_out", seq_len(4L))]
)

se_out <- as.matrix(
  lip_data_real[, paste0("se_out", seq_len(4L))]
)

scale_factor <- sqrt(variant_variance)

gamma_exp <- gamma_exp * scale_factor
se_exp <- se_exp * scale_factor
gamma_out <- gamma_out * scale_factor
se_out <- se_out * scale_factor

cor_exp <- as.matrix(
  lip_corr_real[seq_len(9L), seq_len(9L)]
)

# The primary paper analysis removes the first outcome.
gamma_out <- gamma_out[, 2:4, drop = FALSE]
se_out <- se_out[, 2:4, drop = FALSE]

cor_out <- as.matrix(
  lip_corr_real[11:13, 11:13]
)

n_instruments <- nrow(gamma_exp)
n_exposures <- ncol(gamma_exp)
n_outcomes <- ncol(gamma_out)

if (!identical(
  c(n_instruments, n_exposures, n_outcomes),
  c(177L, 9L, 3L)
)) {
  stop(
    "Unexpected real-data dimensions.",
    call. = FALSE
  )
}

cat(
  sprintf(
    "Real data: %d instruments, %d exposures, %d outcomes\n",
    n_instruments,
    n_exposures,
    n_outcomes
  )
)


# -------------------------------------------------------------------------
# Reconstruct Sigma_X, Sigma_Y, and W
# -------------------------------------------------------------------------

Sigma_X_sum <- matrix(
  0,
  nrow = n_exposures,
  ncol = n_exposures
)

Sigma_Y_sum <- matrix(
  0,
  nrow = n_outcomes,
  ncol = n_outcomes
)

for (j in seq_len(n_instruments)) {
  D_exp <- diag(se_exp[j, ])
  D_out <- diag(se_out[j, ])

  Sigma_X_sum <-
    Sigma_X_sum +
    D_exp %*% cor_exp %*% D_exp

  Sigma_Y_sum <-
    Sigma_Y_sum +
    D_out %*% cor_out %*% D_out
}

Sigma_X <- Sigma_X_sum / n_instruments
Sigma_Y <- Sigma_Y_sum / n_instruments
W <- solve(Sigma_Y)


# -------------------------------------------------------------------------
# Reproduce the frozen rank-selection calculation
# -------------------------------------------------------------------------

select_rank <- function(
    W,
    Y,
    X,
    Sigma_X,
    alpha = 0.05,
    min_rank = 1L
) {
  px <- ncol(X)
  py <- ncol(Y)
  n <- nrow(Y)

  Sigma_xy <- crossprod(X, Y) / n
  debiased_Sigma_xx <- crossprod(X) / n - Sigma_X
  debiased_Sigma_xx_inv <- solve(debiased_Sigma_xx)

  W_sqrt <- legacy$.sqrt_matrix(W)

  target_matrix <-
    W_sqrt %*%
    t(Sigma_xy) %*%
    debiased_Sigma_xx_inv %*%
    Sigma_xy %*%
    W_sqrt

  eigenvalues <- eigen(target_matrix)$values

  candidate_ranks <- min_rank:(min(px, py) - 1L)

  p_values <- vapply(
    candidate_ranks,
    function(rank) {
      log_sum_tail <- sum(
        log(
          1 + eigenvalues[
            (rank + 1L):length(eigenvalues)
          ]
        )
      )

      statistic <-
        (n - (px + py + 1) / 2) *
        log_sum_tail

      1 - stats::pchisq(
        statistic,
        df = (py - rank) * (px - rank)
      )
    },
    numeric(1)
  )

  selected_position <- which(p_values >= alpha)[1L]

  if (is.na(selected_position)) {
    return(min(px, py))
  }

  candidate_ranks[selected_position]
}

selected_rank <- select_rank(
  W = W,
  Y = gamma_out,
  X = gamma_exp,
  Sigma_X = Sigma_X
)

cat("Selected rank:", selected_rank, "\n")


# -------------------------------------------------------------------------
# Load the frozen reference analysis
# -------------------------------------------------------------------------

reference_environment <- new.env(parent = emptyenv())

load(
  reference_result_file,
  envir = reference_environment
)

if (!exists(
  "real_result",
  envir = reference_environment,
  inherits = FALSE
)) {
  stop(
    "The reference RData file does not contain `real_result`.",
    call. = FALSE
  )
}

reference <- get(
  "real_result",
  envir = reference_environment,
  inherits = FALSE
)

if (!identical(
  as.integer(selected_rank),
  as.integer(reference$r_RR)
)) {
  stop(
    sprintf(
      "Selected rank differs from the frozen result: %d versus %d.",
      selected_rank,
      reference$r_RR
    ),
    call. = FALSE
  )
}

if (
  length(reference$opt_rate) != 1L ||
  !is.finite(reference$opt_rate) ||
  reference$opt_rate < 0
) {
  stop(
    "The frozen regularization rate is invalid.",
    call. = FALSE
  )
}

cat(
  "Frozen regularization rate:",
  format(reference$opt_rate, scientific = TRUE),
  "\n"
)


# -------------------------------------------------------------------------
# Refit the three closed-form MR-rr estimators
# -------------------------------------------------------------------------

fit_naive <- legacy$mr_rr_naive(
  Y = gamma_out,
  X = gamma_exp,
  r = selected_rank,
  W = W
)

fit_corrected <- legacy$mr_rr(
  Y = gamma_out,
  X = gamma_exp,
  r = selected_rank,
  Sigma_X = Sigma_X,
  W = W
)

fit_regularized <- legacy$mr_rr_regularized(
  Y = gamma_out,
  X = gamma_exp,
  r = selected_rank,
  Sigma_X = Sigma_X,
  regularization_rate = reference$opt_rate,
  W = W
)


# -------------------------------------------------------------------------
# Compare regenerated quantities with the frozen reference
# -------------------------------------------------------------------------

max_abs_difference <- function(actual, expected) {
  max(
    abs(actual - expected),
    na.rm = TRUE
  )
}

comparison <- data.frame(
  quantity = c(
    "Sigma_X",
    "Sigma_Y",
    "W",
    "C_naive_MRrr",
    "C_MRrr",
    "C_MRrr_regularized"
  ),
  max_absolute_difference = c(
    max_abs_difference(Sigma_X, reference$Sigma_X),
    max_abs_difference(Sigma_Y, reference$Sigma_Y),
    max_abs_difference(W, reference$W),
    max_abs_difference(
      fit_naive$AB,
      reference$C_naive_MRrr
    ),
    max_abs_difference(
      fit_corrected$AB,
      reference$C_MRrr
    ),
    max_abs_difference(
      fit_regularized$AB,
      reference$C_MRrr_regularized
    )
  ),
  stringsAsFactors = FALSE
)

print(
  comparison,
  row.names = FALSE,
  digits = 6
)

tolerance <- 1e-8

if (any(
  !is.finite(comparison$max_absolute_difference)
)) {
  stop(
    "A real-data comparison produced a non-finite difference.",
    call. = FALSE
  )
}

if (any(
  comparison$max_absolute_difference > tolerance
)) {
  stop(
    paste(
      "The regenerated real-data estimates differ from",
      "the frozen reference beyond tolerance."
    ),
    call. = FALSE
  )
}

cat("\nFrozen real-data core reproduction: PASS\n")
