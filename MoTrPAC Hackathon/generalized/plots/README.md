# Query result plots

```sh
python plot_results.py --results-dir QUERY_OUTPUT_DIR --out-dir PLOT_DIR --top 25
```

The renderer reads `matches.csv.gz`, `rank_correlation.csv`, and `summary_by_group.csv` from `query_core`. It writes one subfolder per disease study/contrast, tissue, reference study/species and reference contrast category, plus a root `plot_manifest.json`. Each subfolder contains PNG and SVG figures, the plotted data CSVs, a compressed bubble-selection CSV with row-level inclusion reasons, and `plot_metadata.json` with the exact selection rule and counts. `--top` affects display only; it does not change query results or q values.

The feature bubble heatmap shows one disease feature per row and MoTrPAC reference layer/time per column. It uses **only unique one-to-one mappings** with finite MoTrPAC `log2_fc`. If more than one source or reference identity lands in the same plotted cell, that cell is excluded and recorded in `bubble_selection.csv.gz`. Features are chosen by the magnitude of their disease signed statistic, falling back to the magnitude of disease log2 fold change. This is a display choice, not a significance test. Color encodes MoTrPAC log2 fold change. Dot area increases as the **original MoTrPAC q value** decreases. A missing q is a fixed small, grey-edged dot and is never drawn as significant. When `disease_gene_symbol` is supplied by the query, labels show the gene symbol with the original feature ID in parentheses; the plotted CSV retains both fields.

Pathway rows have a **separate** `pathway_bubble_heatmap.png`/SVG and `pathway_bubble_plot_data.csv`. It colors dots by finite `reference_statistic`, not an invented fold change. For the BioCarta example the legend reads “MoTrPAC CAMERA Z,” while the footer identifies the disease PAGE statistic and the reference CAMERA statistic as different methods. The reference scope is also printed: the bundled BioCarta table contains only six PAH-selected sets, so its plot is not a complete pathway survey. The pathway selection audit is `pathway_bubble_selection.csv.gz`.

Time points use `reference_timepoint_order` when query results provide it. Otherwise recognized during/post times are sorted by their numeric minute/hour value, and labels without a recognizable time retain their source order. Layer grouping follows RNA, protein, phosphosite, metabolite, pathway. The chart labels show source layer and time, and add platform when multiple assays share a layer/time cell.

The rank time course reads the query core's descriptive Spearman rho by reference layer. Rows with no rho or fewer than three usable pairs are saved with their reason but not plotted. No inferential P value is shown. The coverage chart uses the per-layer/time counts of unique, ambiguous and unmatched disease rows; ambiguous mappings are not silently counted as matches.

These are comparisons of molecular contrasts from separate cohorts. The plots do not establish an exercise response in patients or a treatment effect. Test the renderer with `python -m unittest discover -s tests -v`; the tests build synthetic query-core outputs, including an ambiguous mapping, a missing q value, and PAGE versus CAMERA pathway scores.
