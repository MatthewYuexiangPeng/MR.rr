# Spectral simulation outputs

This step exports the completed unified simulation using one base-R program:
`paper/scripts/36_make_spectral_simulation_outputs.R`.
The same program is intended for the paper's reproducibility package. It reads
the merged numerical results, validates their internal consistency, and writes
LaTeX tables, vector PDF figures, 600-dpi PNG copies and machine-readable statistics. It does not
load archived estimators, run new fits, or change the input RDS.

## Run from the repository root

Use the `paper/spectral-rebuild` branch. Copy the supplied R file into
`paper/scripts/` and the supplied documentation into `paper/`. The generated outputs belong
under `paper/output/spectral_rebuild/`.

If the recovered archive has not been extracted locally, first run:

```bash
mkdir -p paper/output/spectral_rebuild/recovered_input_v1
tar -xzf spectral_simulation_results_v1_recovered.tar.gz \
  -C paper/output/spectral_rebuild/recovered_input_v1
```

The archive contains the computation source snapshot as well as data. Extracting
it into this separate directory preserves that snapshot and avoids overwriting
the working branch's sources. The archive must be the recovered **spectral**
run, not the earlier `additional_results_delta01_b300.tar.gz`.

Then run, in R 4.3 or newer:

```bash
Rscript --vanilla paper/scripts/36_make_spectral_simulation_outputs.R \
  --input=paper/output/spectral_rebuild/recovered_input_v1/paper/output/spectral_rebuild/cluster_simulations_v1_blas_recovery/run/merged/spectral_simulation_results.rds \
  --output=paper/output/spectral_rebuild/manuscript_simulations_layout_v4
```

On the cluster, the input already exists, so use:

```bash
Rscript --vanilla paper/scripts/36_make_spectral_simulation_outputs.R \
  --input=paper/output/spectral_rebuild/cluster_simulations_v1_blas_recovery/run/merged/spectral_simulation_results.rds \
  --output=paper/output/spectral_rebuild/manuscript_simulations_layout_v4
```

Only base/recommended R components distributed with R are used; CVXR, OSQP,
MrDAG and mr.divw are not needed for this export. Supply a new output-directory
name on another run. Existing directories are protected. A failed export leaves
a staging directory for diagnosis and does not install a completed output.

R writes four compact single-line `.tex` tables, nine unlettered figure PDFs
and nine corresponding PNGs directly. Preview documents
can optionally be compiled using LaTeX; this adds no numerical computation:

```bash
(
  set -e
  cd paper/output/spectral_rebuild/manuscript_simulations_layout_v4/tables
  pdflatex -interaction=nonstopmode -halt-on-error tables_preview.tex
  pdflatex -interaction=nonstopmode -halt-on-error tables_preview.tex
)
(
  set -e
  cd paper/output/spectral_rebuild/manuscript_simulations_layout_v4/figures
  pdflatex -interaction=nonstopmode -halt-on-error figures_preview.tex
)
```

For Overleaf, copy the `.tex` tables and figure PDFs into the manuscript project.
The tables require `booktabs`; the figures require `graphicx`. Their numbering
is assigned by the manuscript. Replace the original complete table environments
with `\input{table_main_generic}` and `\input{table_main_sparse_loading}` to
avoid duplicate captions and labels. The existing supplementary inputs for
`table_rank_misspecification` and `table_approximate_low_rank` can retain those
file names.

## Output mapping

| Output | Purpose |
|---|---|
| `tables/table_main_generic.tex` | Generic main simulation table |
| `tables/table_main_sparse_loading.tex` | Sparse-loading main simulation table |
| `tables/table_rank_misspecification.tex` | Generic design, working ranks 1, 2 and 3 |
| `tables/table_approximate_low_rank.tex` | Approximate rank-two simulation table |
| `figures/figure_C_generic_setting1.pdf`, `figure_C_generic_setting4.pdf` | Generic C boxplots: weak and strong instruments |
| `figures/figure_C_sparse_loading_setting1.pdf`, `figure_C_sparse_loading_setting4.pdf` | Sparse-loading C boxplots |
| `figures/figure_prediction_generic_setting1.pdf`, `figure_prediction_generic_setting4.pdf` | Generic risk-score prediction boxplots |
| `figures/figure_prediction_sparse_loading_setting1.pdf`, `figure_prediction_sparse_loading_setting4.pdf` | Sparse-loading prediction boxplots |
| `figures/figure_B_sparse_loading_setting4.pdf` | Aligned sparse loading estimates under strong instruments |
| `statistics/simulation_entrywise.csv` | Numerical checks for each of the 27 C entries |
| `statistics/simulation_summary.csv`, `simulation_table_rows.csv` | 100 unique results and their 108 table rows |
| `statistics/shared_rank2_rows.csv` | The eight shared main/rank-sensitivity rows |
| `statistics/prediction_summary.csv`, `prediction_exposure.csv` | Prediction statistics and the stored fixed exposure profile |
| `statistics/sparse_support_summary.csv`, `sparse_support_replicates.csv`, `sparse_B_entrywise.csv` | Loading/support summaries, all four settings |
| `statistics/sparse_B_aligned.rds` | Aligned B arrays and their true target |
| `statistics/figure_boxplot_statistics.csv` | Boxplot hinges, whiskers, extremes and counts outside each viewport |
| `statistics/figure_viewport_audit.csv` | Per-method lower/upper tail counts, outside percentages and box/median/truth clipping flags |
| `provenance/output_provenance.rds` | Original computation seal, recovery provenance and exporter identity |
| `VALIDATION.txt`, `OUTPUT_MD5.csv` | Export checks and fingerprints of the R-generated files |

## Statistical definitions

For each C entry, Bias is the absolute difference between the mean point
estimate and truth; it is not the mean absolute error. SD uses all 1,000 point
estimates. SE averages the saved bootstrap standard errors. CP recomputes
coverage from the stored lower and upper percentile bounds and is reported in
percentage units. The table shows the median and 25th/75th percentiles across
the 27 entry-specific statistics. Quantiles use R's type 7 convention.

The exporter recomputes the statistics and compares them with the original
strict merger (tolerance `1e-10`). It verifies that the main generic table and
rank-sensitivity table reference the **same eight rank-two results**; their
displayed values must be exactly equal. Sparse simulation SE and CP stay
unavailable. Loading support and risk-score summaries are point-estimation
outputs and do not create bootstrap inference for Sparse MR-rr.

Prediction uses the exposure vector stored in the merged object, and sums
`C[, j] * x[j]` in exposure order. It does not draw another exposure profile.

Loading summaries implement the latest supplied manuscript's alignment rule:
minimize squared Frobenius distance to the stored true B over all two-row
permutations and sign changes. Apply the corresponding transformation to A and
check that AB remains unchanged. Do not use a free rotation or change C.
Estimated support is `abs(aligned B) > 0.001`; the threshold is explicit and can
be changed using `--support-threshold=...`. This reporting threshold is separate
from the estimator's recorded `sparse_threshold=0.01`. Pooled sensitivity,
specificity, precision and FDR, as well as exact recovery counts, are exported.

All observations enter table and boxplot calculations. In layout v4, the two
**Setting 1 C figures** default to a fixed `[-1,1]` magnified viewport, common to
all methods, outcomes and both designs. Triangles mark box hinges outside this
viewport; the caption explicitly identifies the zoom. `--weak-c-range=inherit`
restores version 2's pooled 1st/99th-percentile viewport plus truth and a 6%
margin. `--figure-range=full` overrides the weak-IV zoom and shows the complete
range for every figure. Use a new output directory for either alternative.

The other seven figures keep their version-2 scales: strong-instrument C panels
share the scale within each outcome; prediction panels share one scale across
outcomes within a figure; B panels share one scale across both pathways. Their
central view uses pooled 1st/99th percentiles, truth and a 6% margin. No values
are winsorized or filtered before computing boxes, whiskers or numerical
summaries. The boxplot CSV retains extrema and viewport counts; the new
`figure_viewport_audit.csv` also reports lower/upper tail counts, outside
percentages and whether a hinge, median or true value lies outside the view.

The facet layout introduced in style v2 restores the archived Exposure/Outcome/Pathway facet layout, light
grid lines and full estimator labels. A shared legend sits below the plots.
The default `--figure-palette=paper` uses a fixed, accessible palette;
`--figure-palette=legacy` restores the archived seven-hue palette and green B
boxes. These palettes do not change statistics. The figure contains no overall
title or descriptive caption; copy the generated caption into the manuscript.

Stored B estimates and their summaries keep the original validated alignment.
Only the B **figure** additionally uses one fixed signed row permutation chosen
from the true B: the largest absolute true loading is positive in each pathway,
and pathways are ordered by that loading's exposure. This matches the archived
plot's target orientation. `statistics/figure_B_display_map.csv` records the
mapping; in the current result it leaves pathway 1 unchanged and flips pathway
2. Applying the same signed permutation to the columns of A preserves C. The
export never overwrites saved factors or changes selected support.

For an automatic comparison with your previous local output, add:

```bash
--compare-with=paper/output/spectral_rebuild/manuscript_simulations_local_v1
```

This is an **argument to the Rscript command**, not a separate shell command.
The check requires the same input MD5; identical table captions, setting lines,
method/rank ordering and displayed numbers after removing layout syntax;
numerically unchanged values in nine statistical CSVs (tolerance 1e-10); and
unchanged stored aligned B (tolerance zero). The comparison accepts the old
stacked cells and the new inline cells, while still rejecting changed numbers. Plot axis limits and the B display
orientation are presentation metadata and may differ. A successful comparison
writes `STYLE_COMPARISON.txt`.

PDF and PNG are drawn by the same R function. Cairo PDF embeds fonts when Cairo
is available; otherwise the device uses standard Helvetica and records that
fallback. The default is `--formats=pdf,png --png-dpi=600`; use `--formats=pdf`
for PDF only or 300/1200 for an explicit alternate PNG resolution. The
recommended manuscript image is the vector PDF. PNG is provided for convenient
viewing and existing image workflows. `figures/figure_manifest.csv`,
`method_colors.csv`, `captions.txt` and `alt_text.txt` record the presentation.
See `paper/SPECTRAL_FIGURE_STYLE.md` for the old-to-new figure mapping and
JRSSB guidance.

## What is already fixed, and what remains for the submission release

The source snapshot used by this completed run already includes the
IVW/SRIVW input-SE correction in `paper/lib/paper_external.R`, introduced with
worker commit `30d2926`. Each exposure/outcome's calibrated SE fills its own
column. Both point fits and bootstrap fits use that adapter. This is an input
matrix correction, not a relabelling of the final SE table column. Keep that
adapter in the submission code; do not apply another correction to these
completed results or route the paper back through the old wrappers.

The cross-node discrepancy came from CPU-dependent OpenBLAS numerical kernels.
Recovery fixed `HASWELL`, verified all 12,000 point-data replays, retained 1,523
chunks and recomputed 277 bootstrap chunks. The unchanged strict merger passed.
Those numbers and the recovery policy are preserved in the merged RDS. No
tolerance-based acceptance or hash replacement is part of this export.

The original formal Slurm scripts still set thread counts without fixing the
OpenBLAS kernel. **Their final release integration remains a separate pending
step.** Before distributing a from-scratch submission runner:

1. Enforce a documented compatible BLAS/CPU policy before every R process
   starts, including preparation, controller, array jobs and socket workers.
   Record R/package/BLAS versions, kernel, thread settings and CPU compatibility.
2. Validate the policy in a small native cluster run, including more than one
   compatible node and exact common-data/resample checks. A container alone
   does not prevent dynamic CPU-kernel dispatch.
3. Keep strict source/run checks and use a new run identity after changing the
   formal execution code. Preserve this completed original run and its recovery
   provenance as the source of the present tables.
4. Do not make a hard-coded campus hostname such as `u183` a universal
   submission requirement. Use a configurable, verified compatible environment;
   the recovery's node restriction was specific to that recovery.

The new R export is usable now and does not require a further simulation rerun.
It is not a claim that a fresh full calculation has been validated on every
platform. The final package should expose two documented routes: regenerate
tables/figures from saved results with this exporter, and regenerate numerical
results under the documented computation environment before calling the same
exporter.

This simulation export does not cover rank-selection experiments, tuning-path
experiments, or real-data point fits and bootstrap. Those remaining paper
components, the BLAS integration above, manuscript numerical text updates, and
a clean-directory reproduction check are required before calling the entire
paper package submission-ready.

## Figure and table presentation in layout v4

All nine figures now default to `--panel-labels=none`, as requested by the
authors. Exposure, Outcome and Pathway facet strips remain. The two weak-IV C
figures retain the fixed `[-1,1]` zoom and clipped-box markers. No additional
option is needed to generate the unlettered version. `--panel-labels=letters`
is an explicit optional override.

The four tables use one line per cell: `median (q25, q75)`. The R exporter writes
inline content directly and does not redefine LaTeX's `\shortstack`. Tables
use a 10/12 pt body, fixed local line spacing, a 9/11 pt note and full-width
`tabular*` with compact minimum column padding. See `SPECTRAL_TABLE_LAYOUT.md`.
The R-generated table preview uses A4 paper with 1-inch margins, matching the
supplied manuscript geometry.

## Record the exporter in Git

After placing the five source/documentation files in the existing repository:

```bash
git add paper/scripts/36_make_spectral_simulation_outputs.R paper/SPECTRAL_OUTPUTS.md \
  paper/SPECTRAL_FIGURE_STYLE.md paper/PAPER_UPDATE_STATUS.md paper/SPECTRAL_TABLE_LAYOUT.md &&
git diff --cached --check &&
git commit -m "Compact simulation tables and remove default figure panel letters" &&
git push
```

Keep the original archives and large numerical results outside this source
commit. Layout v4 is tested using R 4.6.0 through WebR against the uploaded
native cluster result and the previous v3 export. All numerical CSVs and the
full-data boxplot statistics remain unchanged. The table files change their
layout while preserving every displayed number and its method/rank placement.
Run the Rscript command on native R to generate the local copy; this is export
validation, not estimator refitting.
