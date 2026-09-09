# Simulation figure presentation, layout v4

Layout v4 uses the version-3 facet design and weak-IV zoom, with **panel
letters off by default in all nine figures**, as explicitly requested by the
authors. Named Exposure/Outcome/Pathway strips, axes, colors and statistics are
unchanged. This document also retains the plotting history for provenance. It reads the same
strictly merged spectral RDS as version 1 and does not fit any estimator or run
bootstrap. It was prepared against the exporter committed as `094e466` on
`paper/spectral-rebuild`.

## What the comparison showed

The archived plotting functions in
`freeze/current_analysis_20260823/project/scripts/MR_rr_simulation_main_260717_cov.R`
(`plot_simulation_result` and `plot_simulation_result_pred`) used labelled
facets, full estimator names, `theme_bw` and PNG output. The archived
`Simulation_sparse_260719.R` (`plot_sparse_B_estimation_canonicalized`) used
Exposure/Pathway facets and a true-loading sign convention. Version 1 of the
new exporter used repeated base-graphics axes and numeric method labels. Those
layout choices made the dense C figure harder to read.

The user's comparison also mixed **Setting 4** (old C/prediction screenshots)
with **Setting 1** (new uploaded C/prediction PDFs). Compare the same setting
before attributing different dispersion to plotting. The completed spectral
rerun also has new numerical results; old and new simulation values need not
agree even when displayed identically.

The second pathway in the stored new B has its largest loading negative.
The archived plot displayed it positive. Version 2 restores that global
display convention **after** the validated replicate-specific alignment to
truth. It flips the display of pathway 2 in all replicates and in the true B,
with the corresponding transformation of A preserving C. It does not repeat
the archived replicate-wise heuristic instead of the current alignment rule.
The mapping is exported as `figure_B_display_map.csv`; original B summaries
remain in their stored orientation.

## Version 3: weak-instrument C detail views

The two Setting 1 C figures (S4 and S6 in the supplied supplement) now use
`[-1,1]` on every outcome row of both designs. This is a fixed, shared window,
not an estimator-specific axis and not a claim that all estimates lie inside
it. The completed data have all true effects and all 378 boxplot medians
inside the window. Each figure has nine boxes with at least one hinge beyond
the boundary; triangles identify those cut boxes. Entire tails, including
extreme weak-IV estimates, still enter the boxes and numerical summaries.

This choice magnifies differences in central estimates while preserving the
reported instability in the full-data SDs and the per-method clipping audit.
Do not infer a method's total variability solely from the visible height of
a clipped box. Include the generated plotting-window note in the figure caption.

The previous view used one pooled 1st/99th percentile interval per outcome.
The archived C plotting call used outcome-specific probabilities
`c(0.98, 0.9925, 0.985)` for its weak-IV example, through `coord_cartesian`.
Consequently, restoring that code's probabilities would not impose a uniformly
smaller viewport (the second outcome uses a wider probability interval).
Version 3 therefore uses an explicit, common zoom instead of treating an old
probability as an invariant numerical y limit.

| Argument | Effect |
|---|---|
| `--weak-c-range=zoom` (default) | Fixed `[-1,1]` in the two Setting 1 C figures |
| `--weak-c-range=inherit` | Restore version-2 central viewport rules for those figures |
| `--figure-range=full` | Full range in all figures; overrides the weak-IV zoom |
| `--panel-labels=none` (default) | No panel letters in all nine figures |
| `--panel-labels=letters` | Explicit override to restore upper-left letters |

`statistics/figure_viewport_audit.csv` records lower/upper tail counts and
percentages, hinge clipping, and checks for medians and true values outside
the window, for each method and facet. `figures/figure_manifest.csv` records
the effective window rule and panel-label option. Existing table calculations,
boxplot definitions, prediction summaries and B alignment are unchanged.

## Presentation choices

| Feature | Layout v4 |
|---|---|
| C | 3 outcome rows by 9 exposure columns; shared y axis within each row |
| B | 2 pathway rows by 9 exposure columns; common scale for both pathways |
| Prediction | 3 outcome panels; common y scale within the figure |
| Labels | Grey facet strips and full method names; panel letters off |
| Legend | One shared legend below C and prediction figures |
| Truth | Red dashed lines; gray dashed lines for zero true B |
| Palette | Default fixed palette; optional archived palette via `--figure-palette=legacy` |
| Files | Individual vector PDFs and matching 600-dpi PNGs |
| Dimensions | C: 7.2 by 5.35 in; B: 7.2 by 3.3 in; prediction: 7.2 by 3.9 in |
| Captions | Written separately in `captions.txt`; draft alt text in `alt_text.txt` |
| Display window | Weak C: fixed `[-1,1]`; other figures: central 98% pooled viewport plus truth and 6% margin; optional full range |

All 1,000 observations enter Tukey box statistics. The viewport may hide tail
points, but never filters observations before computing boxes, whiskers,
support rates, or table statistics. The old C code also used viewport limits,
including manually selected probabilities for some rows. Version 3 documents
the weak-C zoom separately and records all out-of-view counts. It does not
restore the archived prediction plot's hard-coded `[-2,2]` window.

Use `--figure-range=full` in a separate output directory to inspect all tails.
Use the generated caption for the chosen setting; do not reuse a strong-
instrument caption for a weak-instrument plot.

## Manuscript replacement map

The figure numbers below follow the supplied supplementary manuscript; LaTeX
labels are more stable if numbering changes. Files have both `.pdf` and `.png`
versions. Prefer `.pdf` in `\includegraphics` and retain the existing labels.

| Supplied figure | LaTeX label | New basename |
|---|---|---|
| S3 | `fig:sim_strong` | `figure_C_generic_setting4` |
| S4 | `fig:sim_weak` | `figure_C_generic_setting1` |
| S5 | `fig:sim_strong_sparseB_C` | `figure_C_sparse_loading_setting4` |
| S6 | `fig:sim_weak_sparseB_C` | `figure_C_sparse_loading_setting1` |
| S7 | `fig:sparse_B` | `figure_B_sparse_loading_setting4` |
| S8 | `fig:pred_strong` | `figure_prediction_generic_setting4` |
| S9 | `fig:pred_weak` | `figure_prediction_generic_setting1` |
| S10 | `fig:pred_strong_sparseB_C` | `figure_prediction_sparse_loading_setting4` |
| S11 | `fig:pred_weak_sparseB_C` | `figure_prediction_sparse_loading_setting1` |

For example, copy the PDF into the manuscript's `image/` folder, then replace
the existing image reference with:

```latex
\includegraphics[width=\linewidth]{image/figure_C_generic_setting4.pdf}
```

Keep one caption and one label for each figure. Put the plotting-window
explanation in the caption or an explicit common plotting note. For B, also
mention the fixed display sign convention, so its orientation is clear when
compared with the stored B matrix. Replace numerical claims in the surrounding
text from the regenerated statistics rather than reading off a plotted box.

## Letter setting

The current author-selected default is no letters. The complete R script,
individual PDFs, PNGs and combined figure preview all use this same default.
Exposure/Outcome/Pathway labels locate each facet. The optional
`--panel-labels=letters` switch remains available but is never needed for the
default export. The v4 delivery contains only the default unlettered figures.
