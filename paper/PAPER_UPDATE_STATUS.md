# Paper rebuild status

Status is based on the supplied manuscript and execution logs, as of
2026-09-09. Source branch: `paper/spectral-rebuild`; last user-confirmed commit
before the figure-style update: `094e466`. Generating new results does not by
itself confirm that the manuscript/Overleaf files have been replaced.

| Component | Computation | Reproducible outputs | Manuscript replacement |
|---|---|---|---|
| Main generic simulation | Complete; strict recovered merge passed | Table and C/prediction figures generated | Ready; replacement not confirmed |
| Main sparse-loading simulation | Complete; strict recovered merge passed | Table, C/B/prediction figures and support statistics generated | Ready; replacement not confirmed |
| Working-rank sensitivity | Complete; shared rank-two results reused | Table generated; 8 shared rows exactly match main generic table | Ready; replacement not confirmed |
| Approximate low rank, delta 0.1 | Complete | Table generated | Ready; replacement not confirmed |
| Simulation figure presentation | Uses the same completed results | Layout v4: nine figures without letters, four compact single-line tables; numerical summaries unchanged | Replace the 9 figures and 4 tables; retain labels and use the accompanying captions |
| Simulation rank-selection experiments | Updated rerun not yet confirmed | Updated rank-selection tables pending | Pending |
| Simulation sparse tuning paths | Updated rerun pending | Updated tuning figure pending | Pending |
| Real-data rank selection and tuning paths | Updated rerun pending | Rank test and rank-specific tuning figures pending | Pending |
| Real-data ranks 1/2 point estimates | Updated spectral rerun pending | Updated A/B and C outputs pending | Pending |
| Real-data ranks 1/2 bootstrap | Updated rerun pending | Updated uncertainty tables/figures pending | Pending |
| Formal from-scratch BLAS policy | Targeted recovery succeeded; permanent runner integration still pending | Native compatible-node verification pending | Reproducibility section pending |
| Submission code archive | In development | Clean-directory reproduction and final archive pending | Data/code availability details pending |

Current result identity: `08c2f8e7e5d4119ac1e18264a64b711e`. The targeted
HASWELL recovery retained 1,523 chunks and recomputed 277 bootstrap chunks;
the original strict merger passed. All original point estimates were retained.
The IVW/SRIVW SE-column correction was already present in the source used for
this run; these saved results do not need a second adjustment.

The next computation stage remains rank selection, tuning paths, and the
real-data analyses. The current user-requested figure revision does not mark
those stages as complete or launch them. Before final packaging, integrate the
BLAS policy into the formal runners and verify the new pipeline in a clean
environment. Preserve the existing sealed computation and its provenance.
