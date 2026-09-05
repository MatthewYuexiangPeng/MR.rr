#!/usr/bin/env Rscript

options(warn = 1)

repo_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

core_test_file <- file.path(
  repo_root,
  "paper",
  "scripts",
  "02_smoke_test_real_data_core.R"
)

if (!file.exists(core_test_file)) {
  stop(
    "The real-data core smoke test was not found.",
    call. = FALSE
  )
}

# This reconstructs the real-data inputs and loads:
# legacy, gamma_out, gamma_exp, Sigma_X, W, and reference.
source(
  core_test_file,
  local = FALSE
)

if (!identical(as.integer(reference$r_RR), 1L)) {
  stop(
    "The primary frozen real-data analysis is not rank one.",
    call. = FALSE
  )
}

if (
  length(reference$sparse_eta) != 1L ||
  !is.finite(reference$sparse_eta) ||
  reference$sparse_eta <= 0
) {
  stop(
    "The frozen sparse tuning parameter is invalid.",
    call. = FALSE
  )
}

cat(
  "\nFrozen sparse eta:",
  format(reference$sparse_eta, scientific = TRUE),
  "\n"
)


# -------------------------------------------------------------------------
# Sparse support selection using the frozen penalized implementation
# -------------------------------------------------------------------------

sparse_selection <- legacy$mr_rr_sparse(
  GAMMA_hat = gamma_out,
  gamma_hat = gamma_exp,
  W = W,
  Sigma_X = Sigma_X,
  lambda = rep(
    reference$sparse_eta,
    ncol(gamma_exp)
  ),
  r = reference$r_RR,
  max_iter = 100L,
  tol = 1e-2
)

selected_support <- sparse_selection$B != 0

if (!all(
  dim(selected_support) ==
  dim(reference$sparse_support)
)) {
  stop(
    "The regenerated sparse support has incorrect dimensions.",
    call. = FALSE
  )
}

if (!all(
  selected_support ==
  reference$sparse_support
)) {
  stop(
    "The regenerated sparse support differs from the frozen support.",
    call. = FALSE
  )
}

selected_indices <- which(
  selected_support[1L, ]
)

cat(
  "Selected exposure indices:",
  paste(selected_indices, collapse = ", "),
  "\n"
)

cat(
  "Number of selected coefficients:",
  sum(selected_support),
  "\n"
)


# -------------------------------------------------------------------------
# Construct gamma_tilde exactly as in the frozen implementation
# -------------------------------------------------------------------------

construct_gamma_tilde <- function(
    GAMMA_hat,
    gamma_hat,
    Sigma_X
) {
  n <- nrow(gamma_hat)

  P_GAMMA <-
    GAMMA_hat %*%
    solve(crossprod(GAMMA_hat)) %*%
    t(GAMMA_hat)

  P_GAMMA_complement <- diag(n) - P_GAMMA

  Sigma_gamma_hat <- crossprod(gamma_hat) / n

  matrix_part_1 <-
    Sigma_gamma_hat -
    crossprod(
      P_GAMMA %*% gamma_hat
    ) / n

  matrix_part_2 <- matrix_part_1 - Sigma_X

  if (!legacy$is_psd(matrix_part_2)) {
    matrix_part_2 <- legacy$.nearest_psd(
      matrix_part_2,
      epsilon = 1e-6
    )
  }

  R <- chol(matrix_part_1)
  Q <- chol(matrix_part_2)
  L <- solve(R) %*% Q

  P_GAMMA %*% gamma_hat +
    (P_GAMMA_complement %*% gamma_hat) %*% L
}


# -------------------------------------------------------------------------
# Frozen post-selection refit
# -------------------------------------------------------------------------

sparse_refit <- function(
    GAMMA_hat,
    gamma_hat,
    W,
    Sigma_X,
    support,
    A_init,
    B_init,
    max_iter = 100L,
    tol = 1e-2
) {
  n_exposures <- ncol(gamma_hat)
  rank <- nrow(support)

  if (!all(
    dim(support) ==
    c(rank, n_exposures)
  )) {
    stop(
      "`support` has incorrect dimensions.",
      call. = FALSE
    )
  }

  gamma_tilde <- construct_gamma_tilde(
    GAMMA_hat = GAMMA_hat,
    gamma_hat = gamma_hat,
    Sigma_X = Sigma_X
  )

  W_sqrt <- legacy$.sqrt_matrix(W)

  A_hat <- A_init
  B_hat <- B_init
  B_hat[!support] <- 0

  distance <- Inf

  for (iteration in seq_len(max_iter)) {
    B_hat <- matrix(
      0,
      nrow = rank,
      ncol = n_exposures
    )

    for (component in seq_len(rank)) {
      selected <- which(
        support[component, ]
      )

      if (length(selected) > 0L) {
        latent_response <- as.vector(
          GAMMA_hat %*%
            W %*%
            A_hat[, component]
        )

        B_hat[component, selected] <- qr.solve(
          gamma_tilde[, selected, drop = FALSE],
          latent_response
        )
      }
    }

    svd_result <- svd(
      B_hat %*%
        t(gamma_tilde) %*%
        GAMMA_hat %*%
        W_sqrt
    )

    A_new <-
      solve(W_sqrt) %*%
      svd_result$v %*%
      t(svd_result$u)

    distance <-
      norm(
        A_new %*% B_hat -
          A_hat %*% B_hat,
        type = "F"
      ) /
      max(
        1e-8,
        norm(
          A_hat %*% B_hat,
          type = "F"
        )
      )

    A_hat <- A_new

    if (distance < tol) {
      break
    }
  }

  list(
    A = A_hat,
    B = B_hat,
    AB = A_hat %*% B_hat,
    support = support,
    iteration = iteration,
    distance = distance
  )
}

refitted <- sparse_refit(
  GAMMA_hat = gamma_out,
  gamma_hat = gamma_exp,
  W = W,
  Sigma_X = Sigma_X,
  support = selected_support,
  A_init = sparse_selection$A,
  B_init = sparse_selection$B,
  max_iter = 100L,
  tol = 1e-2
)


# -------------------------------------------------------------------------
# Sign canonicalization for comparing rank-one A and B
# -------------------------------------------------------------------------

canonicalize_AB <- function(A, B) {
  for (component in seq_len(nrow(B))) {
    largest_index <- which.max(
      abs(B[component, ])
    )

    if (B[component, largest_index] < 0) {
      B[component, ] <- -B[component, ]
      A[, component] <- -A[, component]
    }
  }

  list(
    A = A,
    B = B,
    AB = A %*% B
  )
}

actual_display <- canonicalize_AB(
  refitted$A,
  refitted$B
)

reference_display <- canonicalize_AB(
  reference$A_MRrr_sparse,
  reference$B_MRrr_sparse
)

reference_selection_AB <-
  reference$A_MRrr_sparse_selection %*%
  reference$B_MRrr_sparse_selection


# -------------------------------------------------------------------------
# Numerical comparisons
# -------------------------------------------------------------------------

max_abs_difference <- function(actual, expected) {
  max(
    abs(actual - expected),
    na.rm = TRUE
  )
}

comparison <- data.frame(
  quantity = c(
    "penalized sparse AB",
    "post-selection AB",
    "post-selection A",
    "post-selection B"
  ),
  max_absolute_difference = c(
    max_abs_difference(
      sparse_selection$AB,
      reference_selection_AB
    ),
    max_abs_difference(
      refitted$AB,
      reference$C_MRrr_sparse
    ),
    max_abs_difference(
      actual_display$A,
      reference_display$A
    ),
    max_abs_difference(
      actual_display$B,
      reference_display$B
    )
  ),
  stringsAsFactors = FALSE
)

print(
  comparison,
  row.names = FALSE,
  digits = 6
)

cat(
  "Refit iterations:",
  refitted$iteration,
  "\n"
)

cat(
  "Final refit distance:",
  format(
    refitted$distance,
    scientific = TRUE
  ),
  "\n"
)

tolerance <- 1e-5

if (any(
  !is.finite(comparison$max_absolute_difference)
)) {
  stop(
    "A sparse comparison produced a non-finite value.",
    call. = FALSE
  )
}

if (any(
  comparison$max_absolute_difference > tolerance
)) {
  stop(
    paste(
      "The regenerated sparse analysis differs from",
      "the frozen reference beyond tolerance."
    ),
    call. = FALSE
  )
}

cat("\nFrozen sparse support selection: PASS\n")
cat("Frozen sparse post-selection refit: PASS\n")
