# Remaining numerical analyses for the spectral paper rebuild

This is an additive computation stage on `paper/spectral-rebuild`. It reads the
current `R/` estimators through `paper/lib/paper_engine.R`. The completed main
simulation, its recovered chunks, and its source seal are not modified.

## Scope and fixed configuration

| Analysis | Production design | Tasks | Result for the paper |
|---|---|---:|---|
| Working-rank selection | Generic and sparse-loading designs, four settings, 1,000 datasets per cell | 80 | Working-rank table and accompanying numerical statements |
| Simulation sparse tuning | Generic design, four settings, three independent pilot datasets per setting, 32 eta values | 12 | Simulation objective/sparsity paths |
| Real-data sparse tuning | Working ranks 1 and 2, 13 eta values each | 2 | Real-data objective/sparsity paths |
| Real bootstrap, standard methods | 1,000 shared SNP resamples, 20 per task | 50 | Rank-one/rank-two C estimates and percentile intervals |
| Real bootstrap, MrDAG | The same 1,000 SNP resamples, 20 per task | 50 | MrDAG estimates and intervals shared between both rank analyses |

The controller also computes the original-data estimates. There are 11 unique
method/rank combinations: four MR-rr methods at each of ranks 1 and 2, plus IVW,
SRIVW and MrDAG once each. The latter three estimators have no working-rank
argument. Their same saved objects must be used in both real-data displays.

All values are in `paper/config/spectral_remaining.R`. Main simulation eta and
phi values are unchanged. The simulation eta grid adds the selected `1e-4` to
the archived logarithmic grid and selected `1e-3`, giving 32 candidates. The
archived executed pilot call used **three** pilot replicates, despite a function
default of 20; this stage preserves three. These pilots have independent,
recorded streams. Their curves are diagnostics, not a new selection of the
already completed main-simulation tuning constants.

## Statistical details preserved and issues made explicit

### Working-rank diagnostic

The paper-specific diagnostic reproduces the signed inverse and eigenvalues of
the archived `rank_test_M0`: it evaluates ranks 1 and 2 using the Bartlett factor
`n - (px + py + 1)/2`, then falls back to rank 3. It does not substitute the
public rank selector's positive-definiteness guard, project the inverse, take
absolute eigenvalues, or regularize the statistic.

The old code also assigned rank 3 when every p-value was undefined. The new
output preserves that answer in `archived_rank` and `archived_rank*_percent`.
It additionally records `valid`, an error reason, and `undefined_percent`.
The diagnostic `rank` is NA in that situation; its rank percentages use **all
1,000 datasets as the denominator**, with undefined cases a separate category.
Thus the distinction is visible without dropping or relabeling observations.
If undefined cases occur, the manuscript must explain the fallback when using
the archived-rule percentages, or show the separate undefined category.

The archived numerical search starts at **rank 1**, while the supplied general
method description starts at rank 0. This stage preserves the numerical design;
the manuscript's implementation description should state its numerical minimum
rank. This is different from the completed experiment that deliberately fits
misspecified ranks to a common rank-two truth.

### Real data and comparator SEs

The inputs are the existing `paper/input/dat_1e-4.csv` and
`paper/input/rho_mat_1e-4.csv`: 177 harmonized SNPs, nine proteins, and outcomes
2--4 (LAS, CES, SVS). Each beta **and its matching SE** is multiplied by
`sqrt(2 * MAF * (1 - MAF))`. SNP and trait order is retained.

IVW/SRIVW use the actual SNP-specific exposure and outcome SEs and the exposure
correlation matrix `gen_cor`. The main simulation's constant calibrated SE
adapter is not used for real data. The validation compares these calls with
direct native package calls for each outcome. The paper engine's already
corrected main-simulation SE-column code is reused unchanged for that earlier
stage.

Each real bootstrap resamples beta and SE rows together and recalculates
`Sigma_X`, `Sigma_Y`, and `W`. The working rank stays fixed. Spectral regularized
MR-rr reselects phi over the archived D grid `0:15` in every bootstrap sample;
`D=0` still represents its positive phi value. It is not an extra unregularized
candidate.

Sparse real point estimation selects support with eta `1.2e-3` for rank 1 and
`1e-3` for rank 2, then performs the unpenalized refit. Every bootstrap holds
that support fixed and starts from the original **refitted** A and B, matching
the archived real-data call. Tuning paths use the unthresholded penalized
`B_raw`/`AB_raw`, without the post-selection refit. The plotted objective omits
the L1 penalty. Its zero-proportion threshold is `1e-2`.

MrDAG uses 10,000 iterations, burn-in 2,000, thinning 5 for real data. These
settings differ from the earlier main simulation's shorter sampler.

The original heteroskedastic IV-strength diagnostic is saved as
`real_iv_strength`; its smallest value rounds to 36.86. The homogeneous
covariance diagnostic used to construct the phi grid is separately recorded
as `homogeneous_scaled_siv`. The two definitions should not be interchanged in
the manuscript.

## Numerical consistency and reproducibility

The controller creates every rank/pilot dataset and the shared 177-by-1,000
real bootstrap index matrix once, before array submission. Workers read those
saved inputs. With the recovered main result supplied, **all 8,000 rank-study
datasets must exactly match its recorded X/Y hashes**. These checks happen
before the production arrays are released.

The runtime helper sets `OPENBLAS_CORETYPE=HASWELL` and all BLAS thread counts to
one before starting R. Linux workers check AVX2/FMA support. The run also seals
R/package versions, installed comparator code hashes, source/input hashes,
and a numerical probe. Every worker checks them. A kernel/runtime mismatch
stops the task. It is not resolved by rounding data hashes or ignoring them.
The OpenBLAS override and verbose kernel reporting are documented in the
[official runtime-variable documentation](https://www.openmathlib.org/OpenBLAS/docs/runtime_variables/).

There is **no default pin to u183**. `MRRR_REMAIN_CONSTRAINT` or
`MRRR_REMAIN_NODE` can select compatible available hardware if necessary. The
default module is the successful run's `r/4.5.1-5zezbzn`. The code does not
install or upgrade packages on the cluster.

All bootstrap draws are retained. Merge requires every expected task and
finite estimates for every expected fit. It checks the common SNP indices,
data and recalculated covariances across resource classes, and fixed sparse
support. Warnings, convergence flags, and covariance projections are recorded.
An execution PASS does not replace review of those diagnostics or substantive
manuscript claims.

## Cluster commands

Commit the added files locally and push `paper/spectral-rebuild`, then update
the cluster checkout with an explicit branch refspec:

```bash
(
  set -e
  cd /home/peng.1276/MRrr-spectral-cluster
  git diff --quiet
  git diff --cached --quiet
  git fetch origin refs/heads/paper/spectral-rebuild:refs/remotes/origin/paper/spectral-rebuild
  git switch paper/spectral-rebuild
  git merge --ff-only origin/paper/spectral-rebuild
  bash paper/slurm/submit_spectral_remaining.sh submit
)
```

Default output directory:
`paper/output/spectral_rebuild/cluster_remaining_v1`.

Default reference:
`paper/output/spectral_rebuild/cluster_simulations_v1_blas_recovery/run/merged/spectral_simulation_results.rds`.
For this completed main run, its ID is `08c2f8e7e5d4119ac1e18264a64b711e` and
its file MD5 is `b13e2682325c0861b56338bb8bf50be5`.

The controller runs a complete native development preflight (actual
CVXR/OSQP, mr.divw, MrDAG, serial/PSOCK equality, checkpoint checks, and a strict
development merge), then materializes the full inputs and computes real point
fits. Only after success are the four production arrays released. The final
merger depends on successful completion of all arrays. This follows Slurm's
[array and dependency semantics](https://slurm.schedmd.com/sbatch.html).

```bash
# Read-only progress; rerun as needed.
bash paper/slurm/submit_spectral_remaining.sh status

# Automatic terminal refresh; Ctrl-C stops monitoring, not jobs.
watch -n 60 'bash paper/slurm/submit_spectral_remaining.sh status'

# After failed jobs have ended, retry missing work using saved checkpoints.
bash paper/slurm/submit_spectral_remaining.sh resume
```

`resume` refuses duplicate submissions while recorded jobs are active. It
requires the same computation commit and validates saved inputs/results.
Completed chunks are reused; successful draws within interrupted chunks are
also reused. If a chunk is corrupt or the source changes, inspect the error;
the script deliberately does not silently replace a completed result.

Site options may be exported before submission: `MRRR_REMAIN_PARTITION`,
`MRRR_REMAIN_ACCOUNT`, `MRRR_REMAIN_CONSTRAINT`, and `MRRR_REMAIN_NODE`.
`MRRR_REMAIN_RUN` and `MRRR_REMAIN_REFERENCE` accept absolute paths. The default
concurrent-task limits are 8 rank, 4 tuning, 8 standard-bootstrap and 8 MrDAG
bootstrap tasks; each bootstrap task uses up to eight R processes with
single-threaded BLAS. The default combined ceiling is 140 CPUs. Limits can be
lowered with `MRRR_REMAIN_RANK_CONCURRENT`, `MRRR_REMAIN_ETA_CONCURRENT`,
`MRRR_REMAIN_STD_CONCURRENT`, and `MRRR_REMAIN_DAG_CONCURRENT`.

## Local commands and submission reuse

Local base validation needs the existing package source and base R only:

```bash
Rscript --vanilla paper/scripts/38_validate_spectral_remaining.R --base-only --cores=1
```

To rerun without Slurm, use the same R functions and configuration:

```bash
Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=prepare --profile=full --run=YOUR_NEW_RUN
Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=task --run=YOUR_NEW_RUN --task=1 --cores=1
# Execute all task IDs in tasks.csv, then:
Rscript --vanilla paper/scripts/37_spectral_remaining.R --action=merge --run=YOUR_NEW_RUN
```

For a new raw-data reproduction, the CLI's `--reference` is optional. The
cluster submitter for the current rebuild requires it so that the current rank
study is paired with the completed main simulation. A newly prepared run seals
its own environment; transporting a sealed computation to a different runtime
and resuming there is intentionally rejected. Merged results can later be
read by the presentation-only R exporter on another machine.

## Saved outputs

`merged/spectral_remaining_results.rds` is the single input for the next
presentation/export stage. It contains the source/runtime seal, rank outcomes,
tuning paths, raw original-data point fits including A/B, the common SNP index
matrix, all 1,000 bootstrap C estimates per unique method/rank, interval
summaries, phi selections and numerical diagnostics.

The merger also writes:

- `rank_replicates.csv`, `rank_summary.csv`;
- `simulation_eta.csv`, `simulation_eta_summary.csv`, `real_eta.csv`;
- `real_summary.csv` (297 rows: 11 unique fits times 27 entries);
- `real_phi_point.csv`, `real_phi_bootstrap.csv`, `diagnostics.csv`;
- `STATUS.txt` with a production flag and completeness checks.

Individual tuning fits and full bootstrap fit diagnostics remain in the
checksummed chunks/checkpoints. Preserve the run directory on the cluster.
After a production strict merge PASS, the merged results and run metadata can
be collected with:

```bash
(
  set -e
  cd /home/peng.1276/MRrr-spectral-cluster
  MRRR_REMAIN_REL=paper/output/spectral_rebuild/cluster_remaining_v1
  grep -Fx 'SPECTRAL REMAINING STRICT MERGE: PASS' "$MRRR_REMAIN_REL/merged/STATUS.txt"
  grep -Fx 'Production: TRUE' "$MRRR_REMAIN_REL/merged/STATUS.txt"
  MRRR_REMAIN_ARCHIVE="spectral_remaining_results_$(date -u +%Y%m%dT%H%M%S).tar.gz"
  test ! -e "$MRRR_REMAIN_ARCHIVE"
  tar -czf "$MRRR_REMAIN_ARCHIVE" \
    "$MRRR_REMAIN_REL/merged" "$MRRR_REMAIN_REL/seal.rds" \
    "$MRRR_REMAIN_REL/tasks.csv" "$MRRR_REMAIN_REL/PREPARATION.txt" \
    "$MRRR_REMAIN_REL/job_history.tsv" "$MRRR_REMAIN_REL/submission_commit.txt"
  sha256sum "$MRRR_REMAIN_ARCHIVE"
)
```

This compact export is sufficient for the next figures/tables. It is not a
replacement for preserving all computation checkpoints for the final research
archive. No expensive model fitting is needed in the future presentation-only
exporter.

## Validation supplied with this change

The delivery includes a base-R numerical/checkpoint validation report and
shell integration checks against fake Slurm executables, including execution
from a Slurm spool path and a repository path containing spaces. The available
local validation runtime is WebR; it cannot validate native CVXR/MrDAG or
socket parallelism. Those are required by the automatic cluster preflight.
The supplied local status therefore says `PARTIAL_PASS`, not a production
computation PASS. No production tasks have been submitted from this workspace.
