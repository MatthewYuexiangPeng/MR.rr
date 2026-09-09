# Manuscript numerical rebuild status

Snapshot: 2026-09-09 (comparator validation fix v3).
The v2 controller `11081460` passed saved-design replay and completed development
preparation (26 tasks), then stopped on a vector-vs-column-matrix comparator
validation error. No production arrays were submitted. V3 corrects only this
comparison and uses a fresh `cluster_remaining_v3` directory. Native v3 preflight
and production completion remain pending. The v2 saved-design fix and completed
main-simulation artifacts remain unchanged.
“Generated” describes numerical artifacts. It does not
mean the user's current manuscript text has already been replaced and checked.

| Paper item | Computational status | Manuscript status / next action |
|---|---|---|
| Main generic low-rank performance table (`tab:main_regular_C`) | Production spectral run and recovered strict merge PASS; updated table generated | User updating text/table |
| Main sparse-loading performance table (`tab:main_sparse_loading_C`) | Same completed run; updated table generated | User updating text/table |
| Rank-misspecification and approximate-low-rank tables | Same completed run; eight shared main/rank rows verified identical | User updating text/table |
| C-entry, aligned sparse B and prediction simulation plots (nine figures) | Generated from the completed RDS; latest no-letter/single-line-table presentation accepted | User updating figures/captions |
| Population truth matrices and fixed prediction exposure profile | Design validation already passed; no new Monte Carlo fits needed | Retain validated values; check label/caption correspondence |
| Rank-selection table (`tab:rank_selection`) and percentages in text | New cluster code supplied; production computation pending | Update after reviewing valid/undefined and archived-fallback counts |
| Simulation eta paths (`fig:sparse_eta_simulation`) | New cluster code supplied; production computation pending | Replace path plot and review visual-selection description |
| Real eta paths (`fig:sparse_eta_realdata`) | New cluster code supplied; production computation pending | Replace both rank panels after diagnostic review |
| Real rank-one A/B table and C interval figure (`tab:realdata_AB_rank1`, `fig:realdata_C_rank1`) | New point/refit/bootstrap code supplied; production computation pending | Export from new merged RDS; verify sign convention and substantive claims |
| Real rank-two A/B table and C interval figure (`tab:realdata_AB_rank2`, `fig:realdata_C_rank2`) | Same remaining run, rank two; production computation pending | Export from the same RDS, sharing rank-independent methods |
| Abstract/results/discussion numerical claims | Depend on the above artifacts | Reconcile after both numerical batches are final |
| Final submission code archive | Main stage and remaining stage use the shared spectral engine; final integration still pending | Record environment/dependencies, recovery provenance, figure/table mapping and end-to-end commands; perform final reproduction checks |

The remaining cluster run covers the numerical inputs needed for the remaining
rank-selection, tuning and real-data tables/figures. The next local stage is
their presentation-only R exporter, using the same accepted style as the
completed simulation figures and tables.

Two wording discrepancies in the supplied manuscript require reconciliation:

1. The general rank-selection description includes rank 0, while the archived
   numerical searches start at rank 1. The remaining workflow records the
   numerical rule explicitly, including its historical full-rank fallback.
2. The general sparse method text describes cross-validation, while the
   numerical sections and frozen executed code use visual inspection of eta
   objective/sparsity paths. This workflow reproduces the latter; it does not
   silently introduce a different tuning experiment.

The main and real-data SE adapters are appropriate to their distinct input
models. The new real-data adapter retains actual heteroskedastic SEs and the
exposure correlation matrix. The remaining cluster workflow materializes common
inputs and pins the BLAS environment. The old computation/recovery provenance
remains part of the final archive; it is not rewritten as though the recovery
had never occurred.
