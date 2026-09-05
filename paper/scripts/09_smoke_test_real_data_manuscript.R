#!/usr/bin/env Rscript

# Validate the frozen real-data results used in manuscript Table 8, Table S3,
# Figure 2, and Figure S12. This script only reads the frozen project (and an
# optional external real-data result directory); it never modifies either.

options(warn = 1)


locate_repo_root <- function(start = getwd()) {
  current <- normalizePath(
    start,
    winslash = "/",
    mustWork = TRUE
  )

  repeat {
    freeze_root <- file.path(
      current,
      "freeze",
      "current_analysis_20260823",
      "project"
    )

    if (dir.exists(freeze_root)) {
      return(current)
    }

    parent <- dirname(current)

    if (identical(parent, current)) {
      stop(
        "Could not locate the MR.rr repository root.",
        call. = FALSE
      )
    }

    current <- parent
  }
}


load_required_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)

  if (!object_name %in% loaded) {
    stop(
      basename(path),
      " does not contain `",
      object_name,
      "`.",
      call. = FALSE
    )
  }

  environment[[object_name]]
}


assert_close <- function(
    actual,
    expected,
    label,
    tolerance = 1e-12) {
  if (!identical(dim(actual), dim(expected)) ||
      length(actual) != length(expected)) {
    stop(
      label,
      " has unexpected dimensions.",
      call. = FALSE
    )
  }

  difference <- max(abs(actual - expected))

  if (!is.finite(difference) || difference > tolerance) {
    stop(
      label,
      " differs by ",
      format(difference, scientific = TRUE),
      "; tolerance = ",
      format(tolerance, scientific = TRUE),
      ".",
      call. = FALSE
    )
  }

  invisible(difference)
}


resolve_result_file <- function(
    directories,
    preferred_names,
    fallback_pattern,
    label) {
  for (name in preferred_names) {
    for (directory in directories) {
      candidate <- file.path(directory, name)

      if (file.exists(candidate)) {
        return(normalizePath(
          candidate,
          winslash = "/",
          mustWork = TRUE
        ))
      }
    }
  }

  fallback_hits <- unique(unlist(
    lapply(
      directories,
      list.files,
      pattern = fallback_pattern,
      full.names = TRUE
    ),
    use.names = FALSE
  ))

  if (length(fallback_hits) == 1L) {
    return(normalizePath(
      fallback_hits[[1L]],
      winslash = "/",
      mustWork = TRUE
    ))
  }

  if (length(fallback_hits) > 1L) {
    stop(
      paste(
        paste0("Multiple candidate files were found for ", label, ":"),
        paste(fallback_hits, collapse = "\n  "),
        sep = "\n  "
      ),
      call. = FALSE
    )
  }

  stop(
    paste(
      paste0("Required result file is missing for ", label, "."),
      "Expected one of:",
      paste(preferred_names, collapse = "\n  "),
      sep = "\n  "
    ),
    call. = FALSE
  )
}


canonicalize_AB <- function(A, B) {
  for (component in seq_len(nrow(B))) {
    largest_index <- which.max(abs(B[component, ]))

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


align_AB_to_reference <- function(A, B, B_reference) {
  rank <- nrow(B)

  if (!identical(dim(B), dim(B_reference))) {
    stop(
      "The loading matrix and its reference have different dimensions.",
      call. = FALSE
    )
  }

  if (rank == 1L) {
    if (sum((B - B_reference)^2) >
        sum((-B - B_reference)^2)) {
      B <- -B
      A <- -A
    }

    return(list(
      A = A,
      B = B,
      AB = A %*% B
    ))
  }

  if (rank != 2L) {
    stop(
      "Display alignment currently supports ranks one and two.",
      call. = FALSE
    )
  }

  permutations <- list(
    c(1L, 2L),
    c(2L, 1L)
  )
  sign_patterns <- list(
    c( 1,  1),
    c( 1, -1),
    c(-1,  1),
    c(-1, -1)
  )

  best_loss <- Inf
  best_A <- NULL
  best_B <- NULL

  for (permutation in permutations) {
    for (sign_pattern in sign_patterns) {
      B_candidate <- B[permutation, , drop = FALSE]
      A_candidate <- A[, permutation, drop = FALSE]

      B_candidate <-
        diag(sign_pattern) %*% B_candidate
      A_candidate <-
        A_candidate %*% diag(sign_pattern)

      loss <- sum((B_candidate - B_reference)^2)

      if (loss < best_loss) {
        best_loss <- loss
        best_A <- A_candidate
        best_B <- B_candidate
      }
    }
  }

  list(
    A = best_A,
    B = best_B,
    AB = best_A %*% best_B
  )
}


validate_real_result <- function(
    result,
    expected_rank,
    expected_eta,
    label) {
  required_components <- c(
    "C_ivw",
    "C_adivw",
    "C_naive_MRrr",
    "C_MRrr",
    "C_MRrr_regularized",
    "C_MRrr_sparse",
    "C_Mr_DAG",
    "A_MRrr",
    "B_MRrr",
    "A_MRrr_regularized",
    "B_MRrr_regularized",
    "A_MRrr_sparse",
    "B_MRrr_sparse",
    "sparse_support",
    "r_RR",
    "opt_rate",
    "sparse_eta"
  )

  missing_components <- setdiff(
    required_components,
    names(result)
  )

  if (length(missing_components) > 0L) {
    stop(
      paste(
        paste0(label, " is missing components:"),
        paste(missing_components, collapse = ", "),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  if (!identical(as.integer(result$r_RR), expected_rank)) {
    stop(
      label,
      " has rank ",
      result$r_RR,
      "; expected ",
      expected_rank,
      ".",
      call. = FALSE
    )
  }

  if (length(result$sparse_eta) != 1L ||
      !is.finite(result$sparse_eta) ||
      abs(result$sparse_eta - expected_eta) > 1e-15) {
    stop(
      label,
      " has an unexpected sparse tuning parameter.",
      call. = FALSE
    )
  }

  if (length(result$opt_rate) != 1L ||
      !is.finite(result$opt_rate) ||
      result$opt_rate < 0) {
    stop(
      label,
      " has an invalid regularization rate.",
      call. = FALSE
    )
  }

  effect_components <- required_components[
    grepl("^C_", required_components)
  ]

  for (component in effect_components) {
    value <- result[[component]]

    if (!identical(dim(value), c(3L, 9L)) ||
        any(!is.finite(value))) {
      stop(
        label,
        " has an invalid `",
        component,
        "` matrix.",
        call. = FALSE
      )
    }
  }

  factor_components <- list(
    A_MRrr = c(3L, expected_rank),
    B_MRrr = c(expected_rank, 9L),
    A_MRrr_regularized = c(3L, expected_rank),
    B_MRrr_regularized = c(expected_rank, 9L),
    A_MRrr_sparse = c(3L, expected_rank),
    B_MRrr_sparse = c(expected_rank, 9L),
    sparse_support = c(expected_rank, 9L)
  )

  for (component in names(factor_components)) {
    if (!identical(
      dim(result[[component]]),
      factor_components[[component]]
    )) {
      stop(
        label,
        " has an invalid `",
        component,
        "` matrix.",
        call. = FALSE
      )
    }
  }

  if (any(!is.finite(result$A_MRrr)) ||
      any(!is.finite(result$B_MRrr)) ||
      any(!is.finite(result$A_MRrr_regularized)) ||
      any(!is.finite(result$B_MRrr_regularized)) ||
      any(!is.finite(result$A_MRrr_sparse)) ||
      any(!is.finite(result$B_MRrr_sparse))) {
    stop(
      label,
      " contains a non-finite factor loading.",
      call. = FALSE
    )
  }

  invisible(result)
}


make_display_factorizations <- function(result, label) {
  mr <- canonicalize_AB(
    result$A_MRrr,
    result$B_MRrr
  )
  regularized <- align_AB_to_reference(
    result$A_MRrr_regularized,
    result$B_MRrr_regularized,
    mr$B
  )
  sparse <- align_AB_to_reference(
    result$A_MRrr_sparse,
    result$B_MRrr_sparse,
    mr$B
  )

  assert_close(
    mr$AB,
    result$C_MRrr,
    paste0(label, " MR-rr factorization"),
    tolerance = 1e-8
  )
  assert_close(
    regularized$AB,
    result$C_MRrr_regularized,
    paste0(label, " regularized factorization"),
    tolerance = 1e-8
  )
  assert_close(
    sparse$AB,
    result$C_MRrr_sparse,
    paste0(label, " sparse factorization"),
    tolerance = 1e-8
  )

  list(
    mr = mr,
    regularized = regularized,
    sparse = sparse
  )
}


compute_real_siv <- function(data_path, correlation_path) {
  real_data <- utils::read.csv(data_path)
  correlation_data <- utils::read.csv(correlation_path)

  variant_variance <-
    2 * real_data$ImpMAF * (1 - real_data$ImpMAF)
  scale_factor <- sqrt(variant_variance)

  gamma_exp <- sweep(
    as.matrix(real_data[, paste0("gamma_exp", 1:9)]),
    1L,
    scale_factor,
    "*"
  )
  se_exp <- sweep(
    as.matrix(real_data[, paste0("se_exp", 1:9)]),
    1L,
    scale_factor,
    "*"
  )
  correlation_exp <- as.matrix(
    correlation_data[1:9, 1:9]
  )

  correlation_decomposition <- eigen(
    (correlation_exp + t(correlation_exp)) / 2,
    symmetric = TRUE
  )

  if (any(correlation_decomposition$values <= 0)) {
    stop(
      "The exposure correlation matrix is not positive definite.",
      call. = FALSE
    )
  }

  correlation_inverse_sqrt <-
    correlation_decomposition$vectors %*%
    diag(1 / sqrt(correlation_decomposition$values)) %*%
    t(correlation_decomposition$vectors)

  n_instruments <- nrow(gamma_exp)
  n_exposures <- ncol(gamma_exp)
  strength_matrix <- matrix(
    0,
    nrow = n_exposures,
    ncol = n_exposures
  )

  for (instrument in seq_len(n_instruments)) {
    standardized_effect <- as.vector(
      correlation_inverse_sqrt %*%
      (gamma_exp[instrument, ] / se_exp[instrument, ])
    )
    strength_matrix <-
      strength_matrix + tcrossprod(standardized_effect)
  }

  strength_matrix <-
    strength_matrix - n_instruments * diag(n_exposures)
  strength_matrix <-
    (strength_matrix + t(strength_matrix)) / 2

  min(eigen(
    strength_matrix / sqrt(n_instruments),
    symmetric = TRUE,
    only.values = TRUE
  )$values)
}


validate_bootstrap <- function(bootstrap, label) {
  required_methods <- c(
    "ivw",
    "adivw",
    "naive_MRrr",
    "MRrr",
    "MRrr_regularized",
    "MRrr_sparse",
    "Mr_DAG"
  )

  missing_methods <- setdiff(required_methods, names(bootstrap))

  if (length(missing_methods) > 0L) {
    stop(
      paste0(
        label,
        " is missing bootstrap methods: ",
        paste(missing_methods, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  for (method in required_methods) {
    result_array <- bootstrap[[method]]

    if (!identical(dim(result_array), c(3L, 9L, 1000L))) {
      stop(
        label,
        " has unexpected dimensions for `",
        method,
        "`.",
        call. = FALSE
      )
    }

    finite_counts <- apply(
      result_array,
      c(1L, 2L),
      function(x) sum(is.finite(x))
    )

    if (any(finite_counts == 0L)) {
      stop(
        label,
        " has an empty bootstrap cell for `",
        method,
        "`.",
        call. = FALSE
      )
    }
  }

  invisible(bootstrap)
}


make_interval_inventory <- function(bootstrap) {
  method_labels <- c(
    ivw = "IVW",
    adivw = "SRIVW",
    naive_MRrr = "Naive MR-rr",
    MRrr = "MR-rr",
    MRrr_regularized = "Reg. MR-rr",
    MRrr_sparse = "Sparse MR-rr",
    Mr_DAG = "MrDAG"
  )
  outcome_names <- c("LAS", "CES", "SVS")
  exposure_names <- c(
    "MMP12",
    "CNTN1",
    "FGL1",
    "MXRA8",
    "CNTFR",
    "SCG3",
    "HTRA1",
    "CLEC3B",
    "ANTXR2"
  )

  rows <- list()
  row_index <- 1L

  for (method in names(method_labels)) {
    for (outcome_index in seq_along(outcome_names)) {
      for (exposure_index in seq_along(exposure_names)) {
        values <- bootstrap[[method]][
          outcome_index,
          exposure_index,
        ]
        values <- values[is.finite(values)]
        limits <- stats::quantile(
          values,
          probs = c(0.025, 0.975),
          names = FALSE,
          type = 7
        )

        if (limits[[1L]] > 0 || limits[[2L]] < 0) {
          rows[[row_index]] <- data.frame(
            Outcome = outcome_names[[outcome_index]],
            Exposure = exposure_names[[exposure_index]],
            Estimator = unname(method_labels[[method]]),
            Lower = limits[[1L]],
            Upper = limits[[2L]],
            stringsAsFactors = FALSE,
            row.names = NULL
          )
          row_index <- row_index + 1L
        }
      }
    }
  }

  if (length(rows) == 0L) {
    return(data.frame(
      Outcome = character(0),
      Exposure = character(0),
      Estimator = character(0),
      Lower = numeric(0),
      Upper = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}


inventory_keys <- function(inventory) {
  sort(paste(
    inventory$Outcome,
    inventory$Exposure,
    inventory$Estimator,
    sep = "|"
  ))
}


compare_interval_inventory <- function(
    actual,
    expected_keys,
    label) {
  actual_keys <- inventory_keys(actual)
  expected_keys <- sort(expected_keys)

  if (!identical(actual_keys, expected_keys)) {
    missing <- setdiff(expected_keys, actual_keys)
    unexpected <- setdiff(actual_keys, expected_keys)

    stop(
      paste(
        paste0(label, " interval-exclusion inventory differs."),
        paste0(
          "Missing: ",
          if (length(missing) == 0L) {
            "none"
          } else {
            paste(missing, collapse = ", ")
          }
        ),
        paste0(
          "Unexpected: ",
          if (length(unexpected) == 0L) {
            "none"
          } else {
            paste(unexpected, collapse = ", ")
          }
        ),
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  invisible(actual_keys)
}


repo_root <- locate_repo_root()
freeze_root <- file.path(
  repo_root,
  "freeze",
  "current_analysis_20260823",
  "project"
)
freeze_results_root <- file.path(freeze_root, "results")

result_directories <- freeze_results_root
external_results_root <- Sys.getenv(
  "MRRR_REAL_DATA_RESULTS_DIR",
  unset = ""
)

if (nzchar(external_results_root)) {
  external_results_root <- normalizePath(
    external_results_root,
    winslash = "/",
    mustWork = TRUE
  )
  result_directories <- unique(c(
    result_directories,
    external_results_root
  ))
}

result_files <- list(
  rank1_result = resolve_result_file(
    result_directories,
    preferred_names = c(
      "real_result_260820_rank1_eta_1p2e-3.RData",
      "real_result_260722_rank1_eta_1p2e-3.RData"
    ),
    fallback_pattern =
      "^real_result_[0-9]+_rank1_eta_1p2e-3\\.RData$",
    label = "rank-one point estimates"
  ),
  rank2_result = resolve_result_file(
    result_directories,
    preferred_names = c(
      "real_result_260820_rank2_eta_1p0e-3.RData",
      "real_result_260722_rank2_eta_1p0e-3.RData"
    ),
    fallback_pattern =
      "^real_result_[0-9]+_rank2_eta_1p0e-3\\.RData$",
    label = "rank-two point estimates"
  ),
  rank1_bootstrap = resolve_result_file(
    result_directories,
    preferred_names = c(
      paste0(
        "bootstrap_C_list_1000_postselection_",
        "260717_eta_1p2e-3.RData"
      ),
      paste0(
        "bootstrap_C_list_1000_postselection_",
        "260820_rank1_eta_1p2e-3.RData"
      ),
      paste0(
        "bootstrap_C_list_1000_postselection_",
        "260722_rank1_eta_1p2e-3.RData"
      )
    ),
    fallback_pattern = paste0(
      "^bootstrap_C_list_1000_postselection_[0-9]+_",
      "(rank1_)?eta_1p2e-3\\.RData$"
    ),
    label = "rank-one bootstrap estimates"
  ),
  rank2_bootstrap = resolve_result_file(
    result_directories,
    preferred_names = c(
      paste0(
        "bootstrap_C_list_1000_postselection_",
        "260820_rank2_eta_1p0e-3.RData"
      ),
      paste0(
        "bootstrap_C_list_1000_postselection_",
        "260722_rank2_eta_1p0e-3.RData"
      )
    ),
    fallback_pattern = paste0(
      "^bootstrap_C_list_1000_postselection_[0-9]+_",
      "rank2_eta_1p0e-3\\.RData$"
    ),
    label = "rank-two bootstrap estimates"
  )
)

cat("Real-data result files:\n")
for (file_label in names(result_files)) {
  cat(
    "  ",
    file_label,
    ": ",
    basename(result_files[[file_label]]),
    "\n",
    sep = ""
  )
}
cat("Real-data point/bootstrap result inventory: PASS\n")


rank1_result <- load_required_object(
  result_files$rank1_result,
  "real_result"
)
rank2_result <- load_required_object(
  result_files$rank2_result,
  "real_result"
)
rank1_bootstrap <- load_required_object(
  result_files$rank1_bootstrap,
  "bootstrap_C_list"
)
rank2_bootstrap <- load_required_object(
  result_files$rank2_bootstrap,
  "bootstrap_C_list"
)

validate_real_result(
  rank1_result,
  expected_rank = 1L,
  expected_eta = 1.2e-3,
  label = "Rank-one real-data result"
)
validate_real_result(
  rank2_result,
  expected_rank = 2L,
  expected_eta = 1e-3,
  label = "Rank-two real-data result"
)
validate_bootstrap(
  rank1_bootstrap,
  "Rank-one real-data bootstrap"
)
validate_bootstrap(
  rank2_bootstrap,
  "Rank-two real-data bootstrap"
)

cat("Real-data result object structure: PASS\n")


# -------------------------------------------------------------------------
# Estimated scaled instrument strength reported in Section 8
# -------------------------------------------------------------------------

data_path <- file.path(
  freeze_root,
  "data",
  "dat_1e-4.csv"
)
correlation_path <- file.path(
  freeze_root,
  "data",
  "rho_mat_1e-4.csv"
)

if (!file.exists(data_path) || !file.exists(correlation_path)) {
  stop(
    "The frozen real-data inputs required for the SIV check are missing.",
    call. = FALSE
  )
}

estimated_siv <- compute_real_siv(
  data_path,
  correlation_path
)

assert_close(
  round(estimated_siv, 2L),
  36.86,
  "Section 8 estimated scaled instrument strength"
)

cat("Section 8 estimated scaled instrument strength: PASS\n")


# -------------------------------------------------------------------------
# Table 8: rank-one A and B after display sign alignment
# -------------------------------------------------------------------------

rank1_display <- make_display_factorizations(
  rank1_result,
  "Rank-one"
)

rank1_A_actual <- cbind(
  rank1_display$mr$A[, 1L],
  rank1_display$regularized$A[, 1L],
  rank1_display$sparse$A[, 1L]
)
rank1_B_actual <- cbind(
  rank1_display$mr$B[1L, ],
  rank1_display$regularized$B[1L, ],
  rank1_display$sparse$B[1L, ]
)

rank1_A_target <- matrix(
  c(
    -0.0122, -0.0122, -0.0123,
    -0.0042, -0.0042, -0.0044,
    -0.0105, -0.0105, -0.0103
  ),
  nrow = 3L,
  ncol = 3L,
  byrow = TRUE
)
rank1_B_target <- matrix(
  c(
    11.523,  11.523,  11.905,
     1.066,   1.066,   0.000,
    -1.704,  -1.704,  -1.812,
    21.527,  21.527,  22.270,
     3.298,   3.298,   0.000,
    -9.300,  -9.300,  -8.690,
   -12.453, -12.453, -11.802,
    -6.435,  -6.435,  -5.432,
    -8.231,  -8.231,  -8.024
  ),
  nrow = 9L,
  ncol = 3L,
  byrow = TRUE
)

rank1_A_difference <- assert_close(
  round(rank1_A_actual, 4L),
  rank1_A_target,
  "Manuscript Table 8 Panel A"
)
rank1_B_difference <- assert_close(
  round(rank1_B_actual, 3L),
  rank1_B_target,
  "Manuscript Table 8 Panel B"
)

cat("Manuscript Table 8 rank-one factorization: PASS\n")


# -------------------------------------------------------------------------
# Table S3: rank-two A and B after pathway/sign alignment to MR-rr
# -------------------------------------------------------------------------

rank2_display <- make_display_factorizations(
  rank2_result,
  "Rank-two"
)

rank2_A_actual <- cbind(
  rank2_display$mr$A,
  rank2_display$regularized$A,
  rank2_display$sparse$A
)
rank2_B_actual <- cbind(
  t(rank2_display$mr$B),
  t(rank2_display$regularized$B),
  t(rank2_display$sparse$B)
)

rank2_A_target <- matrix(
  c(
    -0.0122,  0.0033, -0.0122,  0.0033, -0.0121, -0.0034,
    -0.0042,  0.0096, -0.0042,  0.0096, -0.0086,  0.0061,
    -0.0105, -0.0071, -0.0105, -0.0071, -0.0054, -0.0115
  ),
  nrow = 3L,
  ncol = 6L,
  byrow = TRUE
)
rank2_B_target <- matrix(
  c(
    11.523, -5.981,  11.523, -5.981,  12.981,   0.000,
     1.066,  1.669,   1.066,  1.669,   0.000,   1.908,
    -1.704, -0.678,  -1.704, -0.678,  -1.063,  -1.918,
    21.527, 15.316,  21.527, 15.316,   7.852,  23.097,
     3.298,  0.602,   3.298,  0.602,   0.000,   0.000,
    -9.300,  1.305,  -9.300,  1.305,  -8.208,  -3.388,
   -12.453,  3.456, -12.453,  3.456, -12.335,   0.000,
    -6.435, -1.675,  -6.435, -1.675,   0.000,   0.000,
    -8.231, -8.862,  -8.231, -8.862,   0.000, -12.072
  ),
  nrow = 9L,
  ncol = 6L,
  byrow = TRUE
)

rank2_A_difference <- assert_close(
  round(rank2_A_actual, 4L),
  rank2_A_target,
  "Supplementary Table S3 Panel A"
)
rank2_B_difference <- assert_close(
  round(rank2_B_actual, 3L),
  rank2_B_target,
  "Supplementary Table S3 Panel B"
)

cat("Supplementary Table S3 rank-two factorization: PASS\n")


# -------------------------------------------------------------------------
# Figure 2 and Figure S12: exact interval-exclusion inventories
# -------------------------------------------------------------------------

rank1_inventory <- make_interval_inventory(rank1_bootstrap)
rank2_inventory <- make_interval_inventory(rank2_bootstrap)

rank1_expected_keys <- c(
  "CES|MMP12|IVW",
  "CES|MMP12|SRIVW",
  "CES|MMP12|MrDAG",
  "SVS|MXRA8|MrDAG",
  "SVS|SCG3|SRIVW"
)
rank2_expected_keys <- c(
  paste(
    "CES",
    "MMP12",
    c(
      "IVW",
      "SRIVW",
      "Naive MR-rr",
      "MR-rr",
      "Reg. MR-rr",
      "Sparse MR-rr",
      "MrDAG"
    ),
    sep = "|"
  ),
  "SVS|MXRA8|MrDAG",
  "SVS|SCG3|SRIVW"
)

compare_interval_inventory(
  rank1_inventory,
  rank1_expected_keys,
  "Manuscript Figure 2"
)
compare_interval_inventory(
  rank2_inventory,
  rank2_expected_keys,
  "Supplementary Figure S12"
)

cat("Manuscript Figure 2 interval exclusions: PASS\n")
cat("Supplementary Figure S12 interval exclusions: PASS\n")
cat(
  "Interval exclusions: rank one=",
  nrow(rank1_inventory),
  ", rank two=",
  nrow(rank2_inventory),
  "\n",
  sep = ""
)
cat(
  "Maximum displayed-table differences: Table 8[A=",
  format(rank1_A_difference, digits = 3),
  ", B=",
  format(rank1_B_difference, digits = 3),
  "], Table S3[A=",
  format(rank2_A_difference, digits = 3),
  ", B=",
  format(rank2_B_difference, digits = 3),
  "]\n",
  sep = ""
)
cat("Frozen real-data manuscript results: PASS\n")
