# Additional MR-rr simulations

## Immutable baseline

- Development branch: `simulation/additional-robustness`
- Baseline commit: `a06d988`
- Historical analysis tag: `analysis-freeze-2026-08-23`
- Files under `freeze/current_analysis_20260823/` are read-only.
- The additional simulations may read the frozen raw calibration data and estimator implementation, but must never modify the frozen snapshot.
- Scripts `paper/scripts/00` through `19` remain the baseline reproduction workflow and should initially remain unchanged.

## Shared simulation settings

The additional simulations use the same four settings as the original generic-C simulation:

| Setting | Measurement-error weight | Genetic-effect weight |
|---|---:|---:|
| 1 | 2.5 | 0.25 |
| 2 | 1.0 | 0.25 |
| 3 | 2.5 | 1.00 |
| 4 | 1.0 | 1.00 |

Additional shared rules:

- Use the same raw calibration data as the original simulations.
- Use the frozen manuscript estimator implementation.
- Use 1,000 Monte Carlo replicates for the final analysis.
- Use deterministic replicate-level random seeds.
- Hold the original scenario-specific regularization rates and sparse penalties fixed unless a separate retuning analysis is explicitly requested.
- Calculate bias relative to the complete true causal-effect matrix.
- Large replicate and bootstrap files are not committed to Git.

## Simulation A: working-rank misspecification

### Data-generating model

- Use the original generic causal-effect matrix.
- True rank: 2.
- True singular values: `(1, 1, 0)`.
- Fit the same generated dataset using working ranks `r = 1, 2, 3`.
- Working rank 2 is included as the correctly specified reference.

### Methods

- Regularized MR-rr.
- Sparse MR-rr.

### Reported summaries

For regularized MR-rr:

- absolute bias;
- empirical standard deviation;
- bootstrap standard error;
- bootstrap coverage probability.

For sparse MR-rr:

- absolute bias;
- empirical standard deviation;
- no bootstrap standard error or coverage probability.

The summaries are calculated for the estimated causal-effect matrix `C`. Sparse-loading support recovery is not reported because the loading factorization is not directly comparable across different working ranks.

### Random-number control

Within each setting and replicate, working ranks 1, 2, and 3 must use exactly the same simulated summary data. The working rank must not alter the data-generation seed.

## Simulation B: approximate low rank

### Data-generating model

Let the singular value decomposition used by the original generic design be

`C = U diag(1, 1, 0) V^T`.

Construct the approximate-low-rank truth as

`C_delta = U diag(1, 1, delta) V^T`,

where `delta > 0` is a single prespecified small singular value.

- Working rank: 2.
- The left and right singular vectors are unchanged.
- `delta` must be fixed before the full cluster run.
- No delta grid will be added unless explicitly requested.

### Methods and table

Use the same methods and column definitions as the original generic-C simulation table:

- IVW;
- SRIVW;
- Naive MR-rr;
- MR-rr;
- Regularized MR-rr;
- Sparse MR-rr;
- MrDAG.

Report absolute bias and empirical standard deviation for all methods. Report bootstrap standard error and coverage probability according to the same rules used by the original manuscript table; Sparse MR-rr has no bootstrap SE or CP.

Bias is calculated relative to the complete `C_delta`, not its rank-two truncation.

## Execution stages

1. Validate the truth matrices and common-seed behavior.
2. Run 2-5 local smoke-test replicates.
3. Run a 50-100 replicate pilot.
4. Lock `delta`, tuning parameters, seed rules, and output schema.
5. Run 1,000 replicates on the cluster.
6. Merge chunks and generate the two manuscript tables.
7. Rerun the original reproducibility checks to confirm no regression.
8. Create a new final analysis tag without modifying the historical tag.
