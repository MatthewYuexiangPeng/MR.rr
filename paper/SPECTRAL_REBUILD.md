# Spectral paper rebuild

This branch starts from `archive/pre-spectral-20260906` (supplied baseline
`c9348b7`). The historical analyses remain in `freeze/current_analysis_20260823`.

## Current milestone: shared MR-rr computation core

`paper/lib/paper_engine.R` loads the estimator sources directly from this
checkout's `R/` directory into an isolated, locked environment. It does not use
an installed copy of MR.rr, frozen estimators, historical results, or functions
left in the interactive workspace. The same `R/` sources are included when the
algorithm package is built. This milestone does not yet replace the old
simulation/Slurm dispatchers or provide a complete paper rerun command.

Run the native core check from the repository root:

```bash
Rscript --vanilla paper/scripts/31_validate_spectral_engine.R
```

The default check requires CVXR with OSQP and executes actual sparse selection,
post-selection refit, and the PSD/Cholesky regression example. It writes
`checks.csv`, `engine_sources.csv`, `session_info.txt`, and `STATUS.txt` under
`paper/output/spectral_rebuild/core_validation/`. `SPECTRAL CORE VALIDATION:
PASS` identifies completion of this core check; it does not identify completion
of the paper simulations. A failed check returns a nonzero exit status.

`--base-only` runs the numerical checks that do not require CVXR. Success in
that mode is explicitly `PARTIAL_PASS`; use the default mode before committing
this milestone. Base checks include a genuine fixed-support least-squares refit.
The command also accepts `--repo=PATH` and `--output-dir=PATH`.

## Computational contract

All MR-rr estimates in the new pipeline use `paper_engine_fit()`. Each call must
specify the method, analysis stage, and working rank. Regularized fits must
receive an explicit phi; the engine calls `implementation = "spectral"`.
Phi selection itself belongs to the analysis workflow, and every tuning fit
must use this same entry. The core does not choose new tuning parameters.

The estimator uses uncentered second moments, as in the manuscript. Its
spectral inverse retains signed eigenvalues, including the continuous value
zero for a zero eigenvalue at positive phi. This is independent of the sparse
covariance projection described below.

| `stage` | Sparse behavior |
| --- | --- |
| `simulation` | Penalized estimate with thresholding; no refit |
| `tuning` | Penalized estimate with thresholding; no refit |
| `real_point` | Penalized support selection, then fixed-support unpenalized refit |
| `real_bootstrap` | Refit on the supplied original-data support; no repeated selection |

Sparse defaults are explicitly passed as threshold `1e-2`, maximum iterations
`100`, tolerance `1e-2`, and solver `OSQP`. Lambda must be supplied for selection.
Analysis-specific overrides must be recorded by the future run configuration;
the separate support-recovery reporting threshold is not an estimator default.
The high-level package function `fit_mr_rr()` still has its existing public
defaults; the paper engine calls low-level functions explicitly.

For `real_point`, the returned `bootstrap_state` stores the support and
initial A/B from the **selection** fit. Pass these unchanged to
`real_bootstrap`, matching the current manuscript's fixed-support procedure.
The bootstrap rank must agree with the support matrix's number of rows.

```r
source("paper/lib/paper_engine.R")
engine <- paper_engine_load(".")

# X, Y, Sigma_X, W and phi come from the analysis configuration/tuning step.
regularized <- paper_engine_fit(
  engine, method = "regularized", stage = "simulation",
  Y = Y, X = X, rank = 2L, Sigma_X = Sigma_X, W = W,
  regularization_rate = phi
)
```

## Sparse Cholesky repair

The previous rule projected the corrected surrogate covariance when its
smallest eigenvalue was below `-1e-8`. A covariance with a small negative
eigenvalue could therefore pass that check while `chol()` failed. This caused
the recorded failures at replicate 840 in the earlier additional run.

The package now applies the **existing absolute eigenvalue floor `1e-6`** in
either case: the PSD check fails, or it passes but Cholesky fails. Successful
Cholesky factors are retained. The uncorrected surrogate covariance still
raises an error when it cannot be factored; it is not automatically projected.

Sparse selection and sparse refit return `numerical_diagnostics` with:

- `corrected_covariance_projected`;
- `projection_reason`: `none`, `not_psd`, or `chol_failed`;
- `psd_check_passed`;
- `eigenvalue_floor` and, for projected matrices, `min_eigenvalue_before`.

The value is a numerical repair record, not evidence of convergence or a
change to the regularized estimator. `iter`, `dist`, and `converged` remain
separate. The engine retains warnings in `fit$paper$warnings` and leaves them
visible in the job log. The production worker must record convergence and
projection counts; it must not silently discard failed fits.

## RNG and source provenance

`paper_with_seed()` restores the caller's RNG kind and seed after both success
and error, including when no seed previously existed. Use scoped
`Mersenne-Twister` for historical truth construction and explicit
`L'Ecuyer-CMRG` streams for simulation data, bootstrap resampling and stochastic
methods. The next milestone defines the complete seed allocation and manifest;
this helper alone is not that specification.

`engine$source_manifest` records actual source-file MD5 values.
`paper_engine_check_sources(engine)` detects files changed after loading. The
production runner must additionally record the Git commit, configuration,
input hashes, package/solver versions, and per-task provenance. A code commit
alone does not archive ignored output directories.

## Remaining milestones

1. Build the calibrated design/configuration and per-replicate RNG layer.
   Generic main analysis and rank sensitivity share data and the same rank-two
   result objects, including bootstrap summaries.
2. Add workers/manifests for generic, sparse-loading, approximate-low-rank,
   rank selection, support recovery, tuning paths and real-data inference.
   Add IVW/SRIVW/MrDAG adapters with their specified controls.
3. Run native environment checks, then the full configured calculations from
   a fixed commit. Merge with coverage, finite-value and failure checks.
4. Generate manuscript tables/figures from the new results and update the
   affected text. Validate the extracted candidate release in a clean directory.

The final paper bundle and the algorithm source package are separate parts of
one release. `.Rbuildignore` excludes `paper/`, so `R CMD build` alone is not a
paper reproducibility bundle. Full recomputation and fast table/figure
rebuilding will have distinct documented commands.
