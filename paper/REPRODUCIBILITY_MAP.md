# MR-rr paper reproducibility map

## Reference snapshot

- Freeze tag: `analysis-freeze-2026-08-23`
- Frozen project:
  `freeze/current_analysis_20260823/project`
- The frozen directory is immutable.
- Exact paper reproduction uses the legacy estimator implementation stored in
  the frozen project.
- The user-facing `MR.rr` package may use numerically improved implementations
  and is not the computational source for exact paper reproduction.

## Canonical workflows

| Workflow | Canonical frozen scripts | Manuscript outputs |
|---|---|---|
| Generic low-rank simulation | `scripts/MR_rr_simulation_main_260717_cov.R` and `scripts/MR_rr_estimators.R` | Tables 1–3; Figures S1, S3–S4, S8–S9 |
| Sparse-loading simulation | `scripts/MR_rr_simulation_main_260717_cov.R`, `scripts/Simulation_sparse_260719.R`, and `scripts/MR_rr_estimators.R` | Tables 4–7; Figures S5–S7, S10–S11 |
| Primary real-data analysis | `scripts/real_data_260820.R` and `scripts/MR_rr_estimators.R` | Table 8; Figure 2; Figure S2 |
| Rank-two sensitivity analysis | `scripts/real_data_260820.R` and `scripts/MR_rr_estimators.R` | Table S3; Figure S12 |
| Simulation bootstrap | `scripts/mian_sim_bootstrap_cluster_version/scripts/run_bootstrap_one.R`, `bootstrap_core.R`, and `make_summary.R` | Bootstrap SE and coverage in Tables 3 and 7 |
| Cluster submission | `scripts/mian_sim_bootstrap_cluster_version/run_no_mrdag.sbatch` and `run_mrdag.sbatch` | Distributed bootstrap execution |

## Excluded scripts

Files under `scripts/outdated/` are historical records. They are retained in
the frozen snapshot but are not part of the canonical reproduction workflow.

## Reproduction rules

1. Never edit files under `freeze/current_analysis_20260823/`.
2. Put portable orchestration scripts under `paper/scripts/`.
3. Preserve the frozen statistical estimator implementation.
4. Path handling, execution control, and output locations may be corrected in
   the orchestration layer.
5. Write newly generated results outside the frozen directory.
6. Compare regenerated numerical outputs against the frozen reference results.
7. Use reduced Monte Carlo and bootstrap counts for smoke tests before running
   the full 1,000-replicate analyses.