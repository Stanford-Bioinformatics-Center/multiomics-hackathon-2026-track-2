# Human MoTrPAC within-tissue discordance pilot

This module builds a descriptive RNA–total-protein event catalog and tests a small predictive model on MoTrPAC **human summary-level exercise effects**. It accepts either the canonical reference CSV used by the generalized query core or the native `MotrpacHumanPreSuspensionAnalysis` loader CSV. No disease name, gene list, tissue, exercise mode, or time is embedded in the analysis.

The package's primary contrasts compare exercise change with time-matched control change: `(exercise post - exercise pre) - (control post - control pre)`. The package distributes summary-level results; it does **not** include individual participant measurements. See the [official package overview](https://motrpac.github.io/MotrpacHumanPreSuspensionAnalysis/articles/package_overview.html).

## Input contract

Required columns:

| Column | Meaning |
| --- | --- |
| `tissue` | Exact tissue label, such as `muscle` or `adipose` |
| `contrast_category` | `EE-CON` or `RE-CON` for the control-adjusted human analysis |
| `timepoint` | MoTrPAC time label, such as `post_24_hr` |
| `layer` | `rna`, `protein`, or `phosphosite` |
| `gene_symbol` | Mapped human gene symbol |
| `log2_fc` | Signed numerator-minus-denominator effect from the same contrast |

Recommended columns: `feature_id`, `uniprot` (parent accession for protein and phosphosite), `ci_lower`, `ci_upper`, `q_value`, `ave_expr`, `species`, `contrast_type`, `contrast_short`, `source_package_version`. The R exporter also writes query-core fields `study_id`, `biospecimen`, `id_namespace`, `entity_id`, `contrast`, `effect_scale`, and `statistic` (the package's `z.std`). The module reads package-native aliases: `assay`, `logFC`, `Timepoint`, `CI.L_calculated`/`CI.R_calculated` (or `CI.L`/`CI.R`), `adj_p_value`, and `AveExpr`. `species` defaults to `human` only when absent, for compatibility with the native human package export. A canonical adapter should provide it explicitly. If `feature_id` is absent, stable gene/accession IDs are synthesized for already-collapsed RNA/protein records; phosphosite identity cannot be established across time without an ID.

Rows with a feature ID mapping to multiple genes or parent UniProt IDs are excluded and counted in `run_summary.json`. Within a gene/UniProt group, one RNA or total-protein feature is selected by coverage across the requested time points, then by average abundance, then by feature ID. Response size and significance are never used to select it. `catalog.csv` includes the selected ID, count of candidate features, accession count, and mapping status. Genes with multiple protein accessions remain in the catalog but are excluded from model fitting.

The R export makes `feature_id` consistent with `id_namespace`: an HGNC symbol for mapped RNA, a UniProt accession for mapped protein, and the native MoTrPAC site ID for phosphosites. It also preserves the original assay identifier in `source_feature_id`, which this module uses for duplicate resolution. Ambiguous or absent gene/accession annotations are blanked in the canonical match fields and retain the original values in `source_gene_symbol`/`source_uniprot`; their namespace falls back to `MoTrPAC` so the query core does not mistake an assay ID for an HGNC symbol or UniProt accession.

## Run

With `MotrpacHumanPreSuspensionAnalysis` 2.0.8 installed in an R library visible to `Rscript`:

```powershell
Rscript .\export_motrpac_human.R --tissue muscle --contrast-category EE-CON --out .\reference_effects.csv.gz
python .\build_discordance.py --input .\reference_effects.csv.gz --out-dir .\muscle_ee_24h --tissue muscle --contrast-category EE-CON --target-time post_24_hr --earlier-times post_15_30_45_min post_3.5_4_hr
```

For another compatible reference CSV, pass its path to `--input`. `--target-time`, `--earlier-times`, `--tissue`, `--contrast-category`, `--alpha`, `--min-effect`, `--equivalence-margin`, `--folds`, and `--seed` are parameters. The target and earlier times must be distinct. The script checks that each selected time has only one `contrast_short` when that field is supplied. If `biospecimen` is present, target RNA and protein must share it. In human MoTrPAC blood, RNA is whole blood while circulating protein is plasma, so the within-compartment model rejects that pairing.

For a broad query-core reference export, use `--tissue all --contrast-category all --assays all`; `--contrast-type all` additionally includes within-group, baseline, modality, and control-only contrasts. Comma-separated choices are also accepted. The default export is the narrower RNA/protein/phosphosite, control-adjusted selection used here. The loader includes only assay/tissue combinations actually present; it does not fill nonexistent layers.

Exporter smoke checks with package 2.0.8: `--tissue all --contrast-category all --assays rna` produced 424,854 primary-contrast rows across three tissues and both exercise modes; `--tissue adipose --contrast-category EE-CON --assays all` produced 78,040 rows spanning RNA, protein, phosphosite, and metabolite assays. Both files passed the query core's reference schema loader. The full all-tissue/all-assay/all-contrast combination was not exported for this module.

## Outputs and rules

- `catalog.csv`: one target-time entry per protein gene/UniProt, with RNA-only genes retained as unmatched. `supported_concordant` or `supported_opposite` requires both layers to have a signed effect of at least `--min-effect` (default absolute log2 FC 0.2), assay `q_value <= --alpha` when available, and a confidence interval excluding zero when available. An effect without either q-value or CI is **indeterminate**.
- `rna_response_protein_equivalent` requires supported RNA change and a **protein confidence interval wholly within** the prespecified `± --equivalence-margin` (default ±0.2 log2 FC). This is a conservative equivalence criterion. A large protein p-value alone never means “unchanged.” The analogous class is available for protein response with equivalently small RNA.
- `model_predictions.csv`: out-of-fold predictions for unambiguous, target-time matched genes. It includes observed `RNA log2 FC − protein log2 FC` and the corresponding residual calculated from each protein prediction. These residual columns are an **algebraic transform** of the protein predictions, not a separate validation target. `model_metrics.csv`: MAE, RMSE, and R² for zero-change, contemporaneous RNA-only ridge regression, and temporal ridge regression. Metrics cover all mapped pairs and, when at least 20 eligible genes are available, the RNA-responsive subset. The subset uses the **same all-gene out-of-fold predictions**; no model is refit or tuned on that subset. The temporal model adds earlier RNA and, where observed, the largest absolute early `phosphosite log2 FC − parent protein log2 FC` per gene. Missing predictors are imputed within training folds. All rows for a gene stay in the same fold. Ridge penalty is fixed at 10; there is no tuning on held-out responses.
- `discordance_overview.png`: class counts, same-time RNA/protein effects, and observed versus out-of-fold predicted protein effects. `run_summary.json`: settings, coverage, mapping exclusions, and interpretation limits.

The site-minus-parent difference is **descriptive**, not measured phosphorylation occupancy: covariance between the source estimates is unavailable. Metabolites are not forced into gene-level RNA–protein pairs. They belong in pathway or reaction context with explicit mapping confidence. The small model tests whether earlier RNA/PTM summaries improve prediction of later protein effects; its coefficients do not estimate synthesis/degradation rates or prove mechanisms.

## Reproducibility check

```powershell
python -m unittest -v .\test_discordance.py
```

The synthetic test checks supported opposite signs, conservative equivalence, indeterminate nonsignificant protein, duplicate protein-feature selection, multi-accession model exclusion, native aliases, and end-to-end output files. The exporter and Python pipeline were also run on the local c2.0 / package 2.0.8 muscle EE-CON release; see `demo_muscle_ee/` for a transparently exploratory run. Its weak out-of-fold prediction is a result, not a success claim.

For an exploratory exercise-mode transfer check after running both `EE-CON` and `RE-CON` with the same tissue and times:

```powershell
python .\transfer_model.py --train-dir .\muscle_ee_24h --test-dir .\muscle_re_24h --out-dir .\ee_to_re_transfer
```

This trains on the source effects and evaluates on shared gene/UniProt pairs in the second exercise mode. It is stricter than fitting and testing within one mode, but both contrasts may share the same control participants and release. It is not external replication.

## Phosphosite versus parent-protein audit

The early PTM response offers a separate, potentially better populated Track 2 question. Run:

```powershell
python .\audit_ptm_parent.py --input .\reference_effects.csv.gz --out-dir .\ptm_parent_ee --tissue muscle --contrast-category EE-CON --times post_15_30_45_min post_3.5_4_hr post_24_hr
```

`ptm_parent_summary.csv` gives raw site-feature rows, mapped site rows, exact gene+UniProt matches to parents with one protein feature and one accession per gene, supported site responses, equivalently small parent responses, and their intersection. `ptm_parent_candidates.csv` lists the intersection; `ptm_parent_overview.png` shows its timing and first-time site-versus-parent effects. A supported site needs BH q ≤ 0.05, absolute log2 FC ≥ 0.2, and a sign-consistent CI excluding zero. An equivalently small parent needs its CI wholly inside ±0.2. These are **phosphosite feature rows**, which can include several sites per gene and are not independent mechanisms. The difference between phosphosite and total-protein response is not a measured phosphorylation occupancy or direct kinase-activity estimate.

The local c2.0 / 2.0.8 EE-CON exploratory audit found 735 early site-feature rows across 350 genes satisfying both rules, from 10,938 sites with an exact unique-parent match. At 3.5 hours the corresponding count was 13 rows across eight genes; at 24 hours, one row for one gene. The complete denominators are in `demo_ptm_parent_ee/ptm_parent_summary.csv`.
