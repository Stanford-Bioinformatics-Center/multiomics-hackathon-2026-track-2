# raw_metab sources

All files were produced on 2026-09-26 by `tool/scripts/dump_metab.R` (run from `tool/`; Slurm job 1521 on danilogin), using R 4.4.3, MotrpacRatTraining6moData 2.0.0, MotrpacHumanPreSuspensionAnalysis 2.0.8, MotrpacBicQC 2.0.0 and data.table 1.18.6.1. The script is write-once: it skips any output file that already exists. Values are copied from the package objects without modification; only columns were selected and a few bookkeeping columns (`stat_type`, `source_object`, annotation columns) were added.

## Rat (MoTrPAC PASS1B, 6-month endurance training)

The rat metabolomics data are deposited at Metabolomics Workbench, project PR001020.

### rat_metab_DA_v1.csv.gz (113,048 rows)
Row-bind of the 19 `METAB_<TISSUE>_DA` objects from MotrpacRatTraining6moData (loaded with `data(list = name, package = "MotrpacRatTraining6moData")`). Package title: "Differential analysis of merged metabolomics datasets". One row per tissue x platform (`dataset`) x `feature_ID` x `sex` x `comparison_group` (1w, 2w, 4w, 8w trained vs sex-matched sedentary control).

- Tissues (19): ADRNL, BAT, COLON, CORTEX, HEART, HIPPOC, HYPOTH, KIDNEY, LIVER, LUNG, OVARY (female only), PLASMA, SKM-GN, SKM-VL, SMLINT, SPLEEN, TESTES (male only), VENACV, WAT-SC. ADRNL, COLON, CORTEX, HYPOTH, OVARY, SKM-VL, SMLINT, SPLEEN, TESTES and VENACV have only one platform. The other tissues have 10 to 13 platforms.
- Platforms (`dataset`): metab-t-acoa, metab-t-amines, metab-t-etamidpos, metab-t-ka, metab-t-nuc, metab-t-oxylipneg, metab-t-tca, metab-u-hilicpos, metab-u-ionpneg, metab-u-lrpneg, metab-u-lrppos, metab-u-rpneg, metab-u-rppos.
- 14,419 unique (tissue, dataset, feature_ID) features; 2,808 unique feature_ID strings; 2,489 unique RefMet names. `feature_ID` is not unique within a tissue, because the same metabolite can be measured on several platforms. Use (tissue, dataset, feature_ID) as the key.
- Columns: tissue, assay (always METAB), dataset, site, is_targeted, feature_ID, feature, metabolite_refmet, metabolite, sex, comparison_group, logFC, logFC_se, tscore, p_value, adj_p_value, selection_fdr, comparison_average_intensity, reference_average_intensity, stat_type, source_object.
- Stat type: `tscore` (t statistic) for every row, so `stat_type` is always "tscore". These objects have no `zscore` or `numNAs` column.
- `adj_p_value` is BY-adjusted across all datasets within the assay. `selection_fdr` is the IHW-adjusted training p-value used to select training-regulated analytes (tissue as covariate). `feature` (format METAB;TISSUE;feature_ID) is filled only for training-regulated features at 5% IHW FDR (3,735 unique); it is NA otherwise.
- `metabolite_refmet` (RefMet name) is non-missing in 100% of rows.

### rat_metab_DA_metareg_v1.csv.gz (100,012 rows)
Row-bind of the 19 `METAB_<TISSUE>_DA_METAREG` objects (package title: "Meta-regression of metabolomics differential analysis results"). For features that were measured on more than one platform, the per-platform rows are replaced by one meta-regression row set with `dataset == "meta-reg"`. The `site` field for these rows lists the combined site_platform pairs. Meta-regression rows carry `zscore` (tscore NA) plus `meta_reg_het_p` and `meta_reg_pvalue`. All other rows carry `tscore`. `stat_type` records which one applies: 91,908 tscore rows and 8,104 zscore rows. There are 12,783 unique (tissue, dataset, feature_ID) and the same 19 tissues. Columns are the same as the DA file plus zscore, meta_reg_het_p and meta_reg_pvalue. Meta-regression rows have NA intensities.

### rat_metab_feature_map_v1.csv.gz (14,420 rows)
`METAB_FEATURE_ID_MAP` ("Metabolite feature IDs and metadata"), all columns: tissue, dataset, metabolite_name, metabolite_refmet, feature_ID_sample_data, feature_ID_da, feature_ID_metareg, dataset_metareg, feature, rt, mz, neutral_mass, formula. Every (tissue, dataset, feature_ID) in rat_metab_DA_v1 is found in this map as (tissue, dataset, feature_ID_da). `metabolite_refmet` is non-missing for 100% of rows (2,489 unique names).

RefMet annotation columns were added by joining on `metabolite_refmet`: refmet_name_current, refmet_id, super_class, main_class, sub_class, refmet_formula, pubchem_cid, kegg_id, hmdb_id, chebi_id, lipidmaps_id, inchi_key and class_source. The column is named `refmet_formula` because `formula` already exists in the map. `class_source` records where the annotation came from:
- RefMet_bulk_download: 5,594 rows, exact or case-insensitive name match to the current RefMet table (see below).
- MotrpacBicQC_named_20210505: 8,725 rows. When the current RefMet table had no match, the row was filled from the MotrpacBicQC extdata file `motrpac-metabolomics-named-revised-20210505.csv`. That file is the MoTrPAC named-metabolite dictionary, with RefMet names, classes and IDs as of 2021-05-05. Most of these rows are lipids whose names use the older RefMet notation (for example `PC(36:5)_lp_b`), which the current RefMet release renamed. This source has no refmet_id or lipidmaps_id.
- NA: 101 rows. These rows have no class (18 names, for example `CAR(4:0)_rp_a`, `LPC(16:0)_rp_b`, `eucramide`).

Of the 2,489 RefMet names, 2,471 have super/main/sub class. KEGG IDs are present in 88% of rows, HMDB in 42% and PubChem CID in 44%.

## Human (MoTrPAC PreSuspension acute exercise)

### human_metab_DA_v1.csv.gz (70,437 rows)
Row-bind of the package data objects `MUSCLE_METAB_DA`, `BLOOD_METAB_DA` and `ADIPOSE_METAB_DA` from MotrpacHumanPreSuspensionAnalysis. These are the same objects that `load_differential_analysis(selected_omes = "metab")` returns. All metabolomics platforms are collapsed into the ome name "metab" (`assay` column), and the platform is in `platform`.

- Columns: tissue, assay, platform, contrast_short, contrast_category, contrast_type, randomGroupCode, Timepoint, feature_id, logFC, CI.L_calculated, CI.R_calculated, t, z.std, degrees_of_freedom, AveExpr, p_value, adj_p_value, source_object.
- Stat type: `t` (t statistic from the mixed model; `full_model` in the package object is `~ 0 + group_timepoint + BMI + calculatedAge + codedsiteid + raw_intensity_post + Sex + (1 | pid)` for 67,389 rows and the same model without raw_intensity_post for 3,048 rows; full_model is not carried into the CSV). `z.std` is also provided.
- muscle: 18,711 rows, 891 features, 21 contrasts, 10 platforms.
- blood: 37,719 rows, 1,143 features, 33 contrasts, 10 platforms. These include metab-t-conv, which has only KET and NEFA.
- adipose: 14,007 rows, 667 features, 21 contrasts, 12 platforms.
- Platforms across tissues: metab-t-acoa, metab-t-amines, metab-t-conv, metab-t-ka, metab-t-nuc, metab-t-oxylipneg, metab-t-tca, metab-u-hilicpos, metab-u-ionpneg, metab-u-lrpneg, metab-u-lrppos, metab-u-rpneg, metab-u-rppos.
- Contrasts: all three tissues share 21 contrasts: {Endur, Resist, Control} x {post_15_30_45_min, post_3.5_4_hr, post_24_hr} vs pre_exercise; Endur/Resist minus Control delta-deltas; Endur minus Resist delta-deltas; and the pre_exercise between-group contrasts. Blood has 12 more contrasts for the during_20_min, during_40_min and post_10_min timepoints. contrast_category is one of EE-CON, RE-CON, EE-RE, EE-EE, RE-RE, CON-CON.
- (feature_id, contrast_short) is unique within a tissue. `feature_id` is already a RefMet name: the package picks the lowest-CV platform measurement per metabolite and labels it with its RefMet name.

### human_metab_clinical_DA_v1.csv.gz (198 rows)
`BLOOD_METAB_T_CLINICAL_DA` (platform metab-t-clinical: Lactate, Glycerol, KET, NEFA, Cortisol, Glucose; 33 contrasts). This is a clinical-chemistry panel. `load_differential_analysis` only returns it with `load_clinical = TRUE`. It is kept separate from the main file, and its columns are the same as the main DA file.

### human_metab_feature_map_v1.csv.gz (2,701 rows)
One row per unique (tissue, platform, feature_id) in human_metab_DA_v1. There are 1,465 unique feature_id values.
- `refmet_name` = feature_id when feature_id appears in `METABOLOMICS_CVS$refmet_name`. This covers 1,463 of 1,465 feature_ids (99.9%). The unmapped ones are KET and NEFA (metab-t-conv).
- `platform_feature_id`, `site` and `feature_cv` come from the `METABOLOMICS_CVS` row with `lowest_CV == "yes"` for that tissue, platform and RefMet name. This is the original platform feature name.
- `f2g_refmet_id` and `f2g_kegg_id` come from `HUMAN_FEATURE_TO_GENE` rows with assay "metab", joined on refmet_name. `f2g_kegg_id` is present in 82% of rows.
- The RefMet annotation columns (refmet_name_current, refmet_id, super_class, main_class, sub_class, formula, pubchem_cid, kegg_id, hmdb_id, chebi_id, lipidmaps_id, inchi_key, class_source) are added the same way as for rat. Sources: 2,689 rows from RefMet_bulk_download and 9 from MotrpacBicQC_named_20210505. 3 rows have no class. 1,462 of 1,465 feature_ids have super/main/sub class. RefMet-download kegg_id is present in 33% of rows, hmdb_id in 39% and pubchem_cid in 41%.

## RefMet (web)

### refmet_classes_v1.csv.gz (3,563 rows)
Downloaded on 2026-09-26 from the Metabolomics Workbench RefMet bulk download, https://www.metabolomicsworkbench.org/databases/refmet/refmet_download.php (208,170 RefMet entries at retrieval). The download was then subset to the 3,563 distinct query names used by the files above. Query names come from rat metabolite_refmet, human DA feature_id, METABOLOMICS_CVS refmet_name and HUMAN_FEATURE_TO_GENE refmet_name. Names were matched in this order: exact refmet_name (1,580), case-insensitive name (2), then refmet_id from HUMAN_FEATURE_TO_GENE (0). 1,981 query names did not match (match_type NA). These are mostly lipid names in the older notation used by MoTrPAC in 2021, such as `PC(36:3)`, `SM(d42:1)` and `TG(51:2)_lp_b`. Columns: query_name, match_type, refmet_id, refmet_name, super_class, main_class, sub_class, formula, exactmass, pubchem_cid, chebi_id, hmdb_id, lipidmaps_id, kegg_id, inchi_key, retrieved, source_url. The full 208k-row download was not kept.
