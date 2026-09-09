# Simulation figure presentation, version 2

This update modifies the presentation layer of script 36. It reads the same
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

## Presentation choices

| Feature | Version 2 |
|---|---|
| C | 3 outcome rows by 9 exposure columns; shared y axis within each row |
| B | 2 pathway rows by 9 exposure columns; common scale for both pathways |
| Prediction | 3 outcome panels; common y scale within the figure |
| Labels | Grey facet strips, full method names, small panel letters |
| Legend | One shared legend below C and prediction figures |
| Truth | Red dashed lines; gray dashed lines for zero true B |
| Palette | Default fixed palette; optional archived palette via `--figure-palette=legacy` |
| Files | Individual vector PDFs and matching 600-dpi PNGs |
| Dimensions | C: 7.2 by 5.35 in; B: 7.2 by 3.3 in; prediction: 7.2 by 3.9 in |
| Captions | Written separately in `captions.txt`; draft alt text in `alt_text.txt` |
| Display window | Central 98% pooled viewport plus truth and 6% margin; optional full range |

All 1,000 observations enter Tukey box statistics. The viewport may hide tail
points, but never filters observations before computing boxes, whiskers,
support rates, or table statistics. The old C code also used viewport limits,
including manually selected probabilities for some rows. Version 2 uses one
documented probability rule across figures, and records the number of values
outside the viewport. It does not restore the archived prediction plot's
hard-coded `[-2,2]` window.

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

## JRSSB guidance checked on 2026-09-09

The journal recommends vector formats for charts and specifies that figure
titles/captions belong in the manuscript, with each multipanel figure supplied
as one file. It also asks for clear labels and accessible descriptions.
These requirements support keeping PDF as the manuscript figure and providing
PNG as a convenient companion. They do not prescribe a ggplot theme, this
palette, these dimensions, or percentile-based display limits; those are
presentation choices for this paper.

Source: [JRSSB, General Instructions, Figures and figure accessibility](https://academic.oup.com/jrsssb/pages/general-instructions).

Inspect the figures at the actual size used in the compiled manuscript. The
27-panel C figure is suited to a full-width supplementary page; avoid shrinking
it to a half-width figure. If a smaller format is required later, split its
outcome rows into larger panels rather than making the labels smaller.
