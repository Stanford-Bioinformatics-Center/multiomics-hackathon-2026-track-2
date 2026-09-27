# motrpac_probe store: schema

Built once by `mprobe store build` (`src/motrpac_probe/store.py`, about 75 s) from read-only hackathon inputs.
Everything here is derived, so it can be deleted and rebuilt. Only public MoTrPAC summary statistics are used.

| file | rows | what |
|---|---|---|
| `contrasts.parquet` | 4,342,327 | **the store**: one gene per comparison column (dataset × tissue × layer × contrast) |
| `columns.parquet` | 422 | one row per comparison column: labels, kind, time, sex |
| `annotations.parquet` | 21,762 | per human gene: GTEx v8 TPM in matched tissues and GO:CC / MitoCarta class flags |
| `genesets.parquet` | 113,205 | gene sets (gs_name, gene_symbol, source) |
| `id_map.parquet` | 71,918 | UniProt / Ensembl / rat symbol → human symbol, for signature input |
| `provenance.parquet` | 40 | key/value: source package versions, input sha256, build cross-check |
| `library/` | | `mprobe library build` output (GMT + JSON index) |
| `raw_metab/`, `metab*.parquet`, `pathway_map.csv` | | metabolomics layer (see "METAB layer" below and `raw_metab/SOURCES.md`) |

## contrasts.parquet

| column | type | meaning |
|---|---|---|
| column_id | category | `dataset|tissue|layer|contrast`, e.g. `rat_train|SKM-GN|PROT|F_8w` |
| dataset | category | `human_acute` (MotrpacHumanPreSuspensionAnalysis 2.0.8) or `rat_train` (MotrpacRatTraining6moData 2.0.0) |
| species | category | `human` or `rat` |
| tissue | category | human: `VL` (vastus lateralis biopsy), `ADIPOSE`, `BLOOD`; rat: the 19 package tissue codes (SKM-GN, SKM-VL, HEART, WAT-SC, …) |
| layer | category | `RNA`, `PROT`, `PHOSPHO` (gene-collapsed sites), `PROT_OLINK` (human blood Olink) |
| contrast | category | human: `EE_vs_CON_24h` (delta-delta vs control, the 03_join names), `EE_post_vs_pre_24h`, `CON_post_vs_pre_24h`, `EE_vs_RE_24h`, `EE_vs_CON_pre`, …; rat: `F_1w` … `M_8w` (trained vs sex-matched sedentary) |
| gene_symbol_human | category | human gene symbol (rat via FEATURE_TO_GENE → RGD ortholog) |
| feature_id | category | the feature kept (max \|stat\| among features mapping to the gene) |
| logFC | float64 | log2 fold change of that feature, as published |
| stat | float64 | human: limma/dream moderated t; rat RNA: DESeq2 Wald z; rat PROT/PHOSPHO: limma t |
| stat_type | category | `t_moderated`, `zscore`, `t` |
| p | float64 | published p |
| fdr | float64 | published adjusted p (human BH within contrast; rat BY across tissues) |
| fdr_bh | float64 | BH recomputed within dataset × tissue × layer × contrast over all tested features, before collapsing |
| n_collapsed | int16 | number of features that mapped to the gene in this column |
| training_fdr | float32 | rat only: package selection_fdr (overall training test) |
| baseline_expr | float32 | human: limma AveExpr (log2 CPM for RNA; package scale for PROT/PHOSPHO; muscle only); rat RNA: log2(sedentary mean normalized count + 1) |
| prot_n_missing | float32 | rat PROT/PHOSPHO: numNAs of the feature; NA elsewhere |
| source | category | `join_table_v2` (copied verbatim from data/join_table_v2.csv) or `alltissue_v1` (rebuilt from data/raw/rat_alltissue_*, data/raw/human_muscle_*, everywhere/raw/*) |

Rules, identical to `scripts/03_join.py` / `scripts/10_everywhere.py` (the functions `bh` and `collapse` are
executed from `10_everywhere.py` through `legacy.py`): rat RNA rows with min(comparison, reference mean count) < 1
or |logFC| > 5 are dropped before collapsing; one row per gene = the feature with max |stat|.

**Build check.** All 52 join-table columns (542,526 rows) were also rebuilt from the raw dumps: 542,526 / 542,526
rows have identical `stat` and `fdr_bh`. Per-column row counts of the `join_table_v2` rows equal the CSV exactly.

## columns.parquet

column_id, dataset, species, tissue, layer, contrast, source, n_genes, n_fdrbh05 (genes at fdr_bh < 0.05),
stat_type, contrast_raw (package `contrast_short`, human), label (e.g. `EE−CON 24h`, `F 8w`), kind (`exercise vs
control`, `exercise post−pre`, `reference (non-exercise)`, `training vs sedentary`; from `10_everywhere.parse_human`),
group (EE/RE/CON/trained), sex (F/M, rat), time, time_order, early_human_muscle (human VL 15–45 min; NARRATIVE R2
flag), gtex_tissue (matched GTEx v8 tissue or empty), col_order (display order: species → tissue → layer → kind →
group → sex → time).

## annotations.parquet

gene_symbol_human; `gtex_tpm|<GTEx tissue>` (median TPM, GTEx v8; matched on Ensembl gene id, then symbol);
mitocarta3 (MitoCarta3.0); go_mitochondrion; stable_complex (GO:CC respiratory chain, ribosome, proteasome,
spliceosome); go_any_complex (any GO:CC term containing "COMPLEX"); secreted_proxy (GO:CC extracellular proxy of
scripts/08: external encapsulating structure, basement membrane, collagen trimer, blood microparticle, platelet
alpha granule lumen, lipoprotein particles); go_contractile_fiber. GO:CC terms come from msigdbr
(model/ext/gocc_msigdbr.csv).

## genesets.parquet

data/raw/genesets_v1.csv (KEGG/Hallmark OXPHOS and glycolysis), fig/tables/genesets_extra.csv (KEGG TCA, all 50
Hallmark, GO:BP mitochondrial translation, Reactome mitochondrial biogenesis; msigdbr 26.1.1),
fig/tables/genesets_pkg.csv (MotrpacHumanPreSuspensionAnalysis MOLECULAR_SIGNATURES: GOBP OXPHOS, complement,
GO:CC contractile fibre, sarcoplasmic reticulum, MITOCARTA_ALL), model/ext/gocc_msigdbr.csv (all GO:CC), and
MYH_SET (MYH1/2/4/7).

## id_map.parquet

id_type `uniprot` (MoTrPAC human proteomics/Olink feature maps; accession with and without isoform suffix),
`ensembl` (ENSG without version from the human RNA feature map and RGD; ENSRNOG via RGD), `rat_symbol` (RGD
RAT_SYMBOL → HUMAN_ORTHOLOG_SYMBOL).

## provenance.parquet

built_at, build_seconds, repo_git_sha, python/pandas/numpy versions, MotrpacRatTraining6moData,
MotrpacHumanPreSuspensionAnalysis, MotrpacRatTraining6mo versions (read from the installed R packages),
join_table_version, row counts, crosscheck_join_vs_rebuilt (JSON), GTEx/MitoCarta sources, sha256 of the legacy
scripts whose code is executed, and `input:<path>` = {bytes, sha256} for every input file.

## METAB layer

Built by `store.build_metab()` from `store/raw_metab/` (provenance in `raw_metab/SOURCES.md`); read with
`store.load_metab()`, `store.load_metab_features()`, `store.load_pathway_map()`, and used through
`motrpac_probe.metab`.

| file | rows | what |
|---|---|---|
| `metab.parquet` | 168,401 | METAB rows, **same columns and dtypes as contrasts.parquet**, `layer = "METAB"`, `source = "metab_v1"` |
| `metab_columns.parquet` | 217 | column index, same fields as columns.parquet; `col_order` continues after the last columns.parquet column |
| `metab_features.parquet` | 3,511 | one row per RefMet name used in the rows: refmet_name, super/main/sub_class, kegg, hmdb, pubchem, class_source, refmet_id, inchi_key, species (rat / human / human/rat), name_key |
| `pathway_map.csv` | 647 | many-to-many pathway map: pathway, side (gene / metabolite), member, rule, source, measured_in (see `PATHWAY_MAP.md`) |

Differences from the gene layers:

* `gene_symbol_human` holds the **RefMet name** of the metabolite, not a gene. `feature_id` is
  `<platform>:<feature_ID>` (rat) or `<platform>:<RefMet name>` (human).
* `stat` is the published t (rat `tscore`, human mixed-model `t`); `stat_type = "t"`. `fdr` is the published
  adjusted p (rat BY across datasets; human package adj_p_value).
* `fdr_bh` = BH within dataset × tissue × contrast over all tested platform features (rat internal standards
  included), computed **before** collapsing. Rat: one row per RefMet name per column = the platform feature with
  max |t| (`legacy.collapse`), `n_collapsed` = number of platform features for that name; `[iSTD]` rows are not
  kept. Human: the package already reports one lowest-CV platform per RefMet name (`n_collapsed = 1`); KET and NEFA
  (no RefMet name) are dropped.
* `training_fdr` = rat package selection_fdr; `baseline_expr` = rat reference_average_intensity / human AveExpr
  (platform scale, not comparable across platforms); `prot_n_missing` is NA.
* Rat and human RefMet notations differ for some names (rat `CAR(3:0)`, human `CAR 3:0`); `name_key` (lower case,
  spaces, brackets and `;` removed) links them. Metabolite signatures map to every variant.

Rows per tissue (columns in brackets): human VL 18,711 (21), ADIPOSE 14,007 (21), BLOOD 37,653 (33); rat ADRNL
1,736, BAT 9,034, COLON 1,872, CORTEX 1,656, HEART 9,584, HIPPOC 8,584, HYPOTH 1,720, KIDNEY 10,072, LIVER 10,144,
LUNG 10,023, PLASMA 8,512, SKM-GN 8,288, SKM-VL 1,568, SMLINT 1,784, SPLEEN 1,800, WAT-SC 8,619 (8 each), OVARY 904
and TESTES 864 (4 each), VENACV 1,266 (6).

Rules for use: genes and metabolites are never joined one-to-one; a pathway is compared across layers only as
set-level scores (cameraPR of the gene set or metabolite set vs the rest of the same column). Metabolite pools are
concentrations, not flux.
