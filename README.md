# MR.rr

`MR.rr` implements reduced-rank Mendelian randomization methods for jointly
analyzing multiple exposures and multiple outcomes from summary-level genetic
association data.

The package provides:

- validation and preparation of aligned GWAS summary data;
- data-driven working-rank selection;
- naive, measurement-error-corrected, spectrally regularized, and sparse
  reduced-rank estimators;
- an optional unpenalized refit after sparse support selection; and
- a unified interface for fitting one or several estimators.

This README documents MR.rr version 0.99.0, a pre-release research
implementation.

## Installation

Install the current development version from GitHub:

```r
install.packages("remotes")

remotes::install_github(
  "MatthewYuexiangPeng/MR.rr",
  ref = "rebuild/paper-release"
)
```

The sparse estimator uses `CVXR`. Install it separately if sparse estimation
is required:

```r
install.packages("CVXR")
```

## Required input

Suppose there are `m` genetic instruments, `p` exposures, and `q` outcomes.
The main data-preparation function accepts:

| Argument | Dimension | Description |
|---|---:|---|
| `beta_exposure` | `m` by `p` | SNP-exposure association estimates |
| `se_exposure` | `m` by `p` | Standard errors of SNP-exposure estimates |
| `beta_outcome` | `m` by `q` | SNP-outcome association estimates |
| `se_outcome` | `m` by `q` | Standard errors of SNP-outcome estimates |
| `cor_exposure` | `p` by `p` | Optional exposure-trait correlation matrix |
| `cor_outcome` | `q` by `q` | Optional outcome-trait correlation matrix |
| `allele_frequency` | length `m` | Optional effect-allele or minor-allele frequencies |
| `variant_variance` | length `m` | Optional genotype variances supplied instead of allele frequencies |
| `W` | `q` by `q` | Optional positive-definite outcome weight matrix |

Rows must refer to the same instruments in the same order. Association
estimates must already be harmonized to the same effect allele. The package
does not perform allele harmonization, LD clumping, or instrument selection.

If neither correlation matrix is supplied, the corresponding identity matrix
is used. If `allele_frequency` is supplied, genotype variance is calculated as
`2 * allele_frequency * (1 - allele_frequency)`.

## Basic workflow

The following reproducible example creates summary data, prepares the inputs,
and fits three MR-rr estimators.

```r
library(MR.rr)

set.seed(20260825)

m <- 300L
p <- 8L
q <- 3L

instrument_names <- paste0("rs", seq_len(m))
exposure_names <- paste0("exposure_", seq_len(p))
outcome_names <- paste0("outcome_", seq_len(q))

beta_exposure <- matrix(
  rnorm(m * p, sd = 0.10),
  nrow = m,
  ncol = p,
  dimnames = list(instrument_names, exposure_names)
)

true_outcome_loading <- c(0.80, -0.50, 0.30)
true_exposure_coefficients <- c(0.30, -0.20, 0, 0.15, 0, 0, 0, 0)

true_effect <- outer(
  true_outcome_loading,
  true_exposure_coefficients
)
dimnames(true_effect) <- list(outcome_names, exposure_names)

beta_outcome <- beta_exposure %*% t(true_effect) +
  matrix(rnorm(m * q, sd = 0.01), nrow = m, ncol = q)
dimnames(beta_outcome) <- list(instrument_names, outcome_names)

se_exposure <- matrix(
  0.01,
  nrow = m,
  ncol = p,
  dimnames = list(instrument_names, exposure_names)
)

se_outcome <- matrix(
  0.015,
  nrow = m,
  ncol = q,
  dimnames = list(instrument_names, outcome_names)
)

allele_frequency <- setNames(
  runif(m, min = 0.10, max = 0.50),
  instrument_names
)

prepared <- prepare_mr_rr_data(
  beta_exposure = beta_exposure,
  se_exposure = se_exposure,
  beta_outcome = beta_outcome,
  se_outcome = se_outcome,
  allele_frequency = allele_frequency
)

fits <- fit_all_mr_rr(
  data = prepared,
  methods = c("naive", "corrected", "regularized"),
  rank = 1,
  regularization_rate = 1e-6
)

fits
round(fits$corrected$AB, digits = 3)
```

Each fitted estimator contains:

- `A`: the outcome loading matrix;
- `B`: the exposure coefficient matrix; and
- `AB`: the estimated causal effect matrix.

## Automatic rank selection

Set `rank = NULL` to select the working rank once and use it for every
requested estimator:

```r
fits_selected <- fit_all_mr_rr(
  data = prepared,
  methods = c("naive", "corrected", "regularized"),
  rank = NULL,
  rank_alpha = 0.05,
  regularization_rate = 1e-6
)

fits_selected$rank

rank_diagnostics <- with(
  fits_selected$rank_selection,
  data.frame(
    candidate_rank = candidate_ranks,
    statistic = round(statistic, digits = 3),
    df = df,
    p_value = signif(p_value, digits = 3)
  )
)

rank_diagnostics
```

The rank-selection procedure is a working diagnostic and should be interpreted
with caution when instruments are weak.

## Sparse estimation and refitting

Sparse MR-rr performs penalized support selection and, by default, an
unpenalized refit on the selected support:

```r
sparse_fit <- fit_mr_rr(
  data = prepared,
  method = "sparse",
  rank = 1,
  sparse_lambda = 1e-2,
  sparse_refit = TRUE,
  sparse_threshold = 1e-2
)

sparse_comparison <- data.frame(
  exposure = exposure_names,
  true_support = true_exposure_coefficients != 0,
  selected_support = as.vector(sparse_fit$sparse_selection$B != 0),
  penalized_B = round(
    as.vector(sparse_fit$sparse_selection$B),
    digits = 3
  ),
  refitted_B = round(as.vector(sparse_fit$B), digits = 3)
)

sparse_comparison
round(sparse_fit$AB, digits = 3)
```

`sparse_fit$B` and `sparse_fit$AB` contain the post-selection refit when
`sparse_refit = TRUE`. The original penalized estimates are retained in
`sparse_fit$sparse_selection`. Exact zeros in
`sparse_fit$sparse_selection$B` define the selected support, while
`sparse_fit$sparse_selection$B_raw` contains the coefficients before the
final numerical threshold is applied and is retained for diagnostics rather
than routine reporting.

## Regularized implementation

The package uses the numerically stable spectral implementation by default:

```r
regularized_fit <- fit_mr_rr(
  data = prepared,
  method = "regularized",
  rank = 1,
  regularization_rate = 1e-6,
  regularized_implementation = "spectral"
)
```

The historical implementation remains available through
`regularized_implementation = "legacy"` for controlled comparisons. Exact
manuscript replication should use the frozen analysis scripts and recorded
software environment in the reproducibility directory.

## Package scope

`MR.rr` is the user-facing estimation package. Manuscript-specific simulation,
bootstrap, cluster-submission, figure-generation, and result-aggregation code
is maintained separately in the frozen reproducibility materials.
