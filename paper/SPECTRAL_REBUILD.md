# Spectral paper rebuild

This branch starts from `archive/pre-spectral-20260906` (supplied baseline
`c9348b7`). The historical analyses remain in `freeze/current_analysis_20260823`.

## Completed core milestone

`paper/lib/paper_engine.R` loads the estimator sources directly from this
checkout's `R/` directory into an isolated, locked environment. It does not use
an installed copy of MR.rr, frozen estimators, historical results, or functions
left in the interactive workspace. The same `R/` sources are included when the
algorithm package is built. The simulation worker milestone below uses this
entry. A complete paper rerun additionally needs the remaining analyses listed
at the end of this document.

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
methods. The design milestone below defines the seed allocation and manifest.

`engine$source_manifest` records actual source-file MD5 values.
`paper_engine_check_sources(engine)` detects files changed after loading. The
production runner must additionally record the Git commit, configuration,
input hashes, package/solver versions, and per-task provenance. A code commit
alone does not archive ignored output directories.

## Completed design milestone

`paper/lib/paper_simulation.R` calibrates the design from the committed
`paper/input/dat_1e-4.csv` and `paper/input/rho_mat_1e-4.csv`. These copy the
supplied calibration inputs with line endings normalized to LF; values are
unchanged and `.gitattributes` preserves LF across platforms. It requires
177 SNPs, nine exposures and three outcomes. Runtime
preparation reads neither frozen estimator code nor archived simulation results.

`paper/config/spectral_simulations.csv` records all twelve design/setting cells.
The four `(measurement-error multiplier, genetic-effect multiplier)` pairs are
`(2.5, 0.25)`, `(1, 0.25)`, `(2.5, 1)`, and `(1, 1)`, in manuscript order.
Defaults remain 1,000 Monte Carlo replicates and 300 SNP bootstrap resamples.
Phi and lambda retain the audited, design/setting-specific values. MrDAG's
simulation controls are explicitly 1,000 iterations and 200 burn-in iterations;
these are separate from the later real-data controls.

| Design | True effect | Working ranks for Reg./Sparse MR-rr | Table views |
| --- | --- | --- | --- |
| `generic` | Generic rank two, singular values `(1,1,0)` | 1, 2, 3 | Main generic and rank misspecification |
| `sparse_loading` | Rank two, one nonzero loading per exposure | 2 | Main sparse-loading |
| `approximate_low_rank` | Generic singular vectors, singular values `(1,1,0.1)` | 2 | Approximate low rank |

The generic and sparse-loading true effects and calibrated covariance matrices
are checked against a small reference extracted from the old archived design
metadata. The generic population SIV values remain approximately 3.601872,
9.004681, 14.407489 and 36.018723. Approximate-low-rank truth uses the same
generic singular vectors and the configured third singular value.

The new generator samples the same independent Gaussian latent exposure and
measurement-error distributions using `rnorm()` and Cholesky factors of the
positive-definite calibrated covariances. This vectorized generator and the
new stream allocation define a **new Monte Carlo run**: historical replicate
draws and numerical table entries are not reproduced. Both new table views
will refer to the same newly computed rank-two results.

### Shared data, result keys and table views

One `(design, setting, replicate)` identifies one data set. Working rank is an
estimator control, never a data-generation control. The main generic analysis
and rank-sensitivity analysis therefore share all data and SNP resample indices.

`catalog.csv` contains 100 unique result keys: 44 generic, 28 sparse-loading,
and 28 approximate-low-rank. `table_map.csv` contains 108 displayed rows because
the eight generic rank-two Reg./Sparse rows each appear in two tables.
For example, both tables refer to
`generic__setting-1__regularized_mr_rr__rank-2`. They do not request a second
fit or a second bootstrap calculation. Sparse MR-rr has no bootstrap entry:
SE and CP remain absent in all of these simulation tables. The `working_rank`
value zero for IVW, SRIVW and MrDAG means "not applicable", not a fitted rank.

The merger summarizes each result key once, then lets both table views select
that same summary. It checks equality of shared rank-two rows. The design
check verifies the mapping and a real spectral point/bootstrapped example;
final LaTeX export is a later milestone.

### Random-stream protocol

The bundle stores a fixed `L'Ecuyer-CMRG` registry with 24 streams per
design/setting/replicate. Its layout reserves 1,000 replicate slots in each of
the twelve cells, regardless of a smaller development run. Separate streams
serve data generation, SNP resampling, and each method/working-rank pair.
Within a method stream, the first state belongs to the point estimate and
successive `nextRNGSubStream()` states belong to bootstrap draws 1 through B.
The resampling stream similarly assigns one substream per bootstrap draw.
Truth generation remains scoped `Mersenne-Twister` with seed 123.

Workers must use the stored state for every stochastic call. They must not
derive seeds from task number, worker number, array order or wall-clock time.
Splitting chunks or changing core counts therefore does not redefine data,
resamples or method states. Increasing B preserves the existing resample prefix.
Version 1 explicitly supports at most 1,000 replicates and 1,000 bootstrap draws.
Cross-platform floating-point/solver differences are still possible; generate
the production bundle on the cluster and record that runtime's source hashes
and package versions. A local Windows bundle is useful for inspecting the
design and must not replace the cluster-generated production bundle.

### Validate and prepare

From the repository root, after applying this milestone:

```bash
Rscript --vanilla paper/scripts/33_validate_spectral_design.R
```

The 14 checks cover archived calibration/truth agreement, stream separation
and replay, shared rank-two result keys, actual spectral point/bootstrap fits,
manifest coverage, invalid task rejection and bundle serialization. This check
uses base/recommended R packages and does not require CVXR or run a cluster
pilot. Reports go to `paper/output/spectral_rebuild/design_validation/`.

After committing the source, prepare the full-size configuration and manifest:

```bash
Rscript --vanilla paper/scripts/32_prepare_spectral_simulations.R
```

The default output is `paper/output/spectral_rebuild/simulations_v1/`. It
contains `simulation_bundle.rds`, configuration, catalog, table map, task and
source manifests, session information and a status file. The command refuses
an existing output directory. For an additional preparation use an explicit
fresh `--output=PATH`; `--replicates=N` and `--bootstrap-size=B` are available
for development. These overrides are retained in the bundle configuration.

The default task grouping is:

| Phase | Resource class | Replicates/chunk | Task records |
| --- | --- | --- | --- |
| Point | Standard | 100 | 120 |
| Point | MrDAG | 25 | 480 |
| Bootstrap | Standard | 20 | 600 |
| Bootstrap | MrDAG | 20 | 600 |

All 1,800 task records are validated for exact coverage before writing. The
preparer performs no production estimation and submits no Slurm jobs. The
worker milestone below consumes this manifest. The old
`submit_additional_full_run.sh` does not consume it.

## Current milestone: workers, checkpoints and cluster execution

See [SPECTRAL_CLUSTER.md](SPECTRAL_CLUSTER.md) for the complete installation,
submission, status and recovery commands. `34_run_spectral_tasks.R` provides
`seal`, `verify`, `task`, `inventory` and `merge` actions.

`paper/lib/paper_external.R` calls the audited comparator interfaces:
`mr.divw::mvmr.ivw`, `mr.divw::mvmr.divw` with `phi_cand = NULL`, and
`MrDAG::MrDAG` with configured iterations/burn-in and thinning 5. The MrDAG
output order is exposures followed by outcomes when extracting causal
effects. The installed packages must expose these interfaces; an API mismatch
fails explicitly. The last working cluster environment used `mr.divw 0.1.0`,
`MrDAG 0.1.1`, and `CVXR 1.0.15` with OSQP. No package is automatically upgraded.

### Comparator standard-error correction

The old `ivw_multiple_outcomes()` and `adivw_multiple_outcomes()` wrappers used
`matrix(rep(sd, each = n_snps), nrow = n_snps, byrow = TRUE)`. For unequal
trait-specific standard errors this mixes the traits within each column.
The new wrapper fills the matrix **by columns**, giving each exposure/outcome
its own calibrated standard error at every SNP. The validator checks this
layout explicitly, and its full native mode checks IVW against uncentered
least squares when outcome weights are constant across SNPs.

This is a correction to the comparator inputs in addition to the spectral
MR-rr change. It can change IVW/SRIVW results, depending on which standard
errors their implementations use. Comparisons with historical tables must
record this correction and the new Monte Carlo streams. The original
wrappers remain in the frozen archive. Other audited comparator arguments
are retained; the new wrapper does not silently switch to a different
package function or tune a different set of controls.

### Execution and result contract

The run seal records the prepared bundle/task checksums, computation/input
source checksums, preparation Git provenance, R version/platform/BLAS, installed
dependency versions, and the direct comparator/solver package code checksums.
Workers verify the seal against their runtime and checkout. Resume requires
the recorded computation commit. A changed algorithm or runtime starts a
new run; it is not mixed into existing completed chunks.

Point tasks use one CPU. Bootstrap tasks maintain a PSOCK worker pool across
replicates and distribute bootstrap draws to that pool. Each draw receives
its explicit per-method RNG state, and BLAS/OpenMP threads are set to one.
The native validator compares serial and PSOCK outputs using exact R-object
equality, including stochastic MrDAG outputs in full mode.

Checkpoints are saved after **each replicate**, including fit errors and
warnings. A completed checkpoint is validated and reused. An interrupted
replicate is recalculated with the same data and streams. Incomplete versions
are preserved before replacement. Checkpoint files and completed chunk files
are installed through a same-directory temporary file and rename. The Slurm
wrapper uses an OS lock for each array task to prevent concurrent duplicates.

Each inferential result requires all B successful bootstrap draws. Its SE is
the sample SD of those B estimates; percentile CI bounds are the type-7 2.5%
and 97.5% quantiles. Coverage compares the true effect with these bounds.
Failed draws never produce partial SE/CP. Per-draw errors and warnings are
retained; raw successful bootstrap coefficient matrices are summarized rather
than all retained, and can be regenerated from the sealed inputs and states.
Sparse point estimates retain convergence and projection diagnostics.
Finite nonconverged sparse estimates are retained and counted, not hidden;
their counts must be assessed before manuscript release.

The merger requires every chunk, every replicate and every configured fit.
It rejects overlapping or non-finite estimates, mismatched source/runtime
identity, inconsistent point/bootstrap data, and inconsistent SNP resamples
between standard and MrDAG tasks. It summarizes the 27 effect entries using
absolute mean bias, empirical SD, mean bootstrap SE, coverage percentage and
RMSE, then their median and 25th/75th percentiles across entries. Both generic
table views select the same rank-two summary. The merged directory appears
only after all validation and all result writes succeed.

`SPECTRAL SIMULATION MERGE: PASS` establishes computational completeness of
these simulation designs. It does not establish completion of real-data
analysis, the rest of the paper, or scientific approval of all diagnostics.

## Remaining milestones

1. Run the native worker checks and these full simulation arrays from a fixed
   cluster checkout, then inspect merged results and diagnostic counts.
2. Extend the workflow for rank selection, support recovery, tuning paths and
   real-data inference, using the same source engine and recorded controls.
3. Generate all remaining configured results and validate their merges.
4. Generate manuscript tables/figures from the new results and update the
   affected text. Validate the extracted candidate release in a clean directory.

The final paper bundle and the algorithm source package are separate parts of
one release. `.Rbuildignore` excludes `paper/`, so `R CMD build` alone is not a
paper reproducibility bundle. Full recomputation and fast table/figure
rebuilding will have distinct documented commands.
