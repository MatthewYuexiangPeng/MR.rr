# Manual end-to-end smoke test for the MR.rr package
#
# Run this script from the root of the MR.rr repository. It reads the frozen
# real-data inputs but uses only the current package functions for data
# preparation and estimation. It does not modify any frozen files.

if (!file.exists("DESCRIPTION")) {
  stop(
    "Run this script from the root of the MR.rr repository.",
    call. = FALSE
  )
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install `devtools` before running this script.", call. = FALSE)
}

devtools::load_all(".")


# -----------------------------------------------------------------------------
# 1. Import the frozen real-data files
# -----------------------------------------------------------------------------

frozen_data_dir <- file.path(
  "freeze",
  "current_analysis_20260823",
  "project",
  "data"
)

association_file <- file.path(frozen_data_dir, "dat_1e-4.csv")
correlation_file <- file.path(frozen_data_dir, "rho_mat_1e-4.csv")
trait_file <- file.path(frozen_data_dir, "traits_1e-4.csv")

required_files <- c(association_file, correlation_file, trait_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0L) {
  stop(
    paste(
      "The following frozen input files were not found:",
      paste(missing_files, collapse = "\n")
    ),
    call. = FALSE
  )
}

association_data <- utils::read.csv(
  association_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
correlation_data <- utils::read.csv(
  correlation_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
trait_data <- utils::read.csv(
  trait_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 2. Extract and label the summary-association matrices
# -----------------------------------------------------------------------------

exposure_indices <- seq_len(9L)

# The original dataset contains four outcomes. The frozen paper analysis drops
# the first outcome and retains outcomes 2, 3, and 4: LAS, CES, and SVS.
outcome_indices <- 2:4

beta_exposure_columns <- paste0("gamma_exp", exposure_indices)
se_exposure_columns <- paste0("se_exp", exposure_indices)
beta_outcome_columns <- paste0("gamma_out", outcome_indices)
se_outcome_columns <- paste0("se_out", outcome_indices)

required_columns <- c(
  "ImpMAF",
  beta_exposure_columns,
  se_exposure_columns,
  beta_outcome_columns,
  se_outcome_columns
)
missing_columns <- setdiff(required_columns, names(association_data))

if (length(missing_columns) > 0L) {
  stop(
    paste0(
      "The association file is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      "."
    ),
    call. = FALSE
  )
}

if (!"x" %in% names(trait_data) || nrow(trait_data) < 9L) {
  stop(
    "The trait file must contain at least nine exposure names in column `x`.",
    call. = FALSE
  )
}

exposure_names <- as.character(trait_data$x[exposure_indices])
outcome_names <- c("LAS", "CES", "SVS")

if (anyNA(exposure_names) ||
    any(!nzchar(exposure_names)) ||
    anyDuplicated(exposure_names)) {
  stop("Exposure names in the trait file must be non-missing and unique.", call. = FALSE)
}

variant_id_candidates <- c(
  "SNP",
  "snp",
  "rsid",
  "RSID",
  "variant",
  "variant_id"
)
variant_id_column <- intersect(variant_id_candidates, names(association_data))

if (length(variant_id_column) > 0L) {
  instrument_names <- as.character(association_data[[variant_id_column[1L]]])
} else {
  instrument_names <- paste0("variant_", seq_len(nrow(association_data)))
}

if (anyNA(instrument_names) ||
    any(!nzchar(instrument_names)) ||
    anyDuplicated(instrument_names)) {
  stop("Instrument identifiers must be non-missing and unique.", call. = FALSE)
}

beta_exposure <- as.matrix(
  association_data[, beta_exposure_columns, drop = FALSE]
)
se_exposure <- as.matrix(
  association_data[, se_exposure_columns, drop = FALSE]
)
beta_outcome <- as.matrix(
  association_data[, beta_outcome_columns, drop = FALSE]
)
se_outcome <- as.matrix(
  association_data[, se_outcome_columns, drop = FALSE]
)

storage.mode(beta_exposure) <- "double"
storage.mode(se_exposure) <- "double"
storage.mode(beta_outcome) <- "double"
storage.mode(se_outcome) <- "double"

dimnames(beta_exposure) <- list(instrument_names, exposure_names)
dimnames(se_exposure) <- list(instrument_names, exposure_names)
dimnames(beta_outcome) <- list(instrument_names, outcome_names)
dimnames(se_outcome) <- list(instrument_names, outcome_names)

allele_frequency <- as.numeric(association_data$ImpMAF)
names(allele_frequency) <- instrument_names


# -----------------------------------------------------------------------------
# 3. Extract the trait-correlation matrices
# -----------------------------------------------------------------------------

if (!all(vapply(correlation_data, is.numeric, logical(1)))) {
  stop(
    "All columns of `rho_mat_1e-4.csv` must be numeric.",
    call. = FALSE
  )
}

correlation_matrix <- as.matrix(correlation_data)
storage.mode(correlation_matrix) <- "double"

if (nrow(correlation_matrix) < 13L || ncol(correlation_matrix) < 13L) {
  stop(
    "The correlation file must contain at least a 13 by 13 matrix.",
    call. = FALSE
  )
}

cor_exposure <- correlation_matrix[
  exposure_indices,
  exposure_indices,
  drop = FALSE
]

# Rows and columns 10:13 correspond to the original four outcomes. Keeping
# indices 11:13 matches the retained outcomes 2:4.
outcome_correlation_indices <- 10L + outcome_indices - 1L
cor_outcome <- correlation_matrix[
  outcome_correlation_indices,
  outcome_correlation_indices,
  drop = FALSE
]

dimnames(cor_exposure) <- list(exposure_names, exposure_names)
dimnames(cor_outcome) <- list(outcome_names, outcome_names)


# -----------------------------------------------------------------------------
# 4. Prepare the data using the package interface
# -----------------------------------------------------------------------------

prepared <- prepare_mr_rr_data(
  beta_exposure = beta_exposure,
  se_exposure = se_exposure,
  beta_outcome = beta_outcome,
  se_outcome = se_outcome,
  cor_exposure = cor_exposure,
  cor_outcome = cor_outcome,
  allele_frequency = allele_frequency
)

stopifnot(
  inherits(prepared, "mr_rr_data"),
  identical(dim(prepared$X), c(nrow(association_data), 9L)),
  identical(dim(prepared$Y), c(nrow(association_data), 3L)),
  identical(dim(prepared$Sigma_X), c(9L, 9L)),
  identical(dim(prepared$Sigma_Y), c(3L, 3L)),
  identical(dim(prepared$W), c(3L, 3L)),
  all(is.finite(prepared$X)),
  all(is.finite(prepared$Y)),
  all(is.finite(prepared$Sigma_X)),
  all(is.finite(prepared$W))
)

# Verify that the package performs the same MAF scaling as the frozen script.
expected_row_scale <- sqrt(
  2 * allele_frequency * (1 - allele_frequency)
)

stopifnot(
  isTRUE(all.equal(
    unname(prepared$X),
    unname(beta_exposure * expected_row_scale),
    tolerance = 1e-12
  )),
  isTRUE(all.equal(
    unname(prepared$Y),
    unname(beta_outcome * expected_row_scale),
    tolerance = 1e-12
  ))
)

cat(
  "Prepared data:",
  prepared$n_instruments, "instruments,",
  prepared$n_exposures, "exposures, and",
  prepared$n_outcomes, "outcomes.\n"
)


# -----------------------------------------------------------------------------
# 5. Fit the package estimators
# -----------------------------------------------------------------------------

fit_corrected <- fit_mr_rr(
  data = prepared,
  method = "corrected",
  rank = NULL
)

selected_rank <- fit_corrected$rank

fit_naive <- fit_mr_rr(
  data = prepared,
  method = "naive",
  rank = selected_rank
)

fit_regularized <- fit_mr_rr(
  data = prepared,
  method = "regularized",
  rank = selected_rank,
  regularization_rate = 1e-13,
  regularized_implementation = "spectral"
)

# The frozen primary real-data analysis selected eta = 1.2e-3 for rank one.
# For a different selected rank, replace this value with the corresponding
# preselected tuning parameter before interpreting the sparse estimate.
sparse_lambda <- if (selected_rank == 1L) 1.2e-3 else 1e-3

fit_sparse <- fit_mr_rr(
  data = prepared,
  method = "sparse",
  rank = selected_rank,
  sparse_lambda = sparse_lambda,
  sparse_refit = TRUE,
  sparse_threshold = 1e-2,
  sparse_max_iter = 100L,
  sparse_tol = 1e-2,
  sparse_solver = "OSQP"
)

fits <- list(
  naive = fit_naive,
  corrected = fit_corrected,
  regularized = fit_regularized,
  sparse = fit_sparse
)

expected_effect_dimension <- c(3L, 9L)

stopifnot(
  all(vapply(
    fits,
    function(fit) identical(dim(fit$AB), expected_effect_dimension),
    logical(1)
  )),
  all(vapply(
    fits,
    function(fit) all(is.finite(fit$AB)),
    logical(1)
  )),
  identical(
    fit_sparse$details$support,
    fit_sparse$sparse_selection$B != 0
  )
)


# -----------------------------------------------------------------------------
# 6. Inspect the results
# -----------------------------------------------------------------------------

cat("Selected rank:", selected_rank, "\n")
print(fit_corrected$rank_selection)

effect_estimates <- lapply(fits, function(fit) fit$AB)

cat("\nNaive estimate:\n")
print(round(effect_estimates$naive, 4))

cat("\nMeasurement-error-corrected estimate:\n")
print(round(effect_estimates$corrected, 4))

cat("\nSpectrally regularized estimate:\n")
print(round(effect_estimates$regularized, 4))

cat("\nSparse post-selection refit estimate:\n")
print(round(effect_estimates$sparse, 4))

cat("\nSparse selected support in B:\n")
print(fit_sparse$details$support)

cat(
  "\nSparse selection converged:",
  fit_sparse$sparse_selection$converged,
  "\nSparse refit converged:",
  fit_sparse$details$converged,
  "\n"
)

cat("\nManual real-data smoke test passed.\n")
