# Single-line simulation tables, layout v4

The four simulation tables are exported by the same
`paper/scripts/36_make_spectral_simulation_outputs.R` used for the figures.
Their captions, labels, setting/method/rank ordering, numbers, rounding and
missing-inference dashes are unchanged.

## Format

- One cell is `median (q25, q75)` on one line, emitted directly by R.
- Table body: 10 pt text with 12 pt baseline, a local 12 pt strut and
  `\arraystretch=1.1`. Extra 1.5 pt spacing after every method row is removed.
- Caption: local single spacing with 6 pt separation from the table.
- Note: 9 pt text with 11 pt baseline inside a full-width minipage. End the
  note's paragraph before closing the local font group.
- `tabular*` spans the current `\linewidth`. Minimum `\tabcolsep` is 1.5 pt;
  `\extracolsep{\fill}` distributes remaining space. Thin spaces separate the
  median and parentheses and the two quartiles. There are no outer cell pads.
- `[htbp]` permits ordinary float placement. Each table fits one page at the
  tested geometry; the float can move to a later page if the current page has
  insufficient space.

The proposed 10 pt / 4 pt-column-padding layout was tested. Under A4 paper with
1-inch margins it exceeded the text width by approximately 14--42 pt across
the four tables. Compact cell spacing and `tabular*` retain 10 pt text while
resolving that horizontal overflow; no graphic scaling or reduced number of
decimal places is used.

The global LaTeX `\shortstack` command is not redefined. Other tables and
figures therefore keep their own meaning of `\\` and font/spacing settings.

## Manuscript use

Copy these four files into the manuscript project and include each once:

| File | Existing label |
|---|---|
| `table_main_generic.tex` | `tab:main_regular_C` |
| `table_main_sparse_loading.tex` | `tab:main_sparse_loading_C` |
| `table_rank_misspecification.tex` | `tab:rank_misspecification` |
| `table_approximate_low_rank.tex` | `tab:approximate_low_rank` |

For example:

```latex
\input{table_main_generic.tex}
```

Each file already contains the `table` environment, caption and label. Replace
the old complete table rather than wrapping the input in another `table` or
adding another caption. Remove any earlier manual `\shortstack` override from
those old table blocks. Keep `\usepackage{booktabs}` in the preamble.

Standalone preview numbering is 1--4. The manuscript's own counters assign the
final table numbers; the rank-table cross-reference uses the existing label.

## Reproduction and checks

Run script 36 using the existing recovered input and a new output directory.
`--compare-with=PREVIOUS_OUTPUT` checks numerical CSVs and a layout-independent
table payload containing captions, setting lines, method/rank positions and
every displayed number. Old stacked tables and new single-line tables can be
compared without ignoring numerical changes.

The four-page preview is compiled from the R-generated table files. The tested
page geometry is A4 with 1-inch margins, the geometry in the supplied `main.tex`
and `Supp.tex`. A separate check also uses a 12 pt, double-spaced manuscript
with the caption package loaded; table and note spacing remains local.

All nine figures use `--panel-labels=none` by default in this same exporter.
No simulation or bootstrap recomputation is performed.
