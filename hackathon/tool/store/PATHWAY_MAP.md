# Pathway map (store/pathway_map.csv)

Built by `python -m motrpac_probe.metab` (`metab.build_pathway_map`) from store files only. Columns: `pathway`, `side` (gene | metabolite), `member` (human gene symbol or MoTrPAC RefMet name), `rule` (why the member is in the pathway), `source` (citation), `measured_in` (gene: layers of contrasts.parquet where the gene is measured; metabolite: species in store/metab_features.parquet).

The map is **many-to-many**. Each pathway has a gene side and a metabolite side, and each side is tested only inside its own layer (cameraPR of the set vs every other feature of the same column). Genes and metabolites are never matched one-to-one, and a metabolite can belong to more than one pathway (CAR(3:0) and CAR(5:0) are in fatty-acid beta-oxidation by RefMet class and in BCAA catabolism by curation). Metabolite pools are concentrations, not flux.

## Gene side

One msigdbr gene set per pathway, from `store/pathway_genesets_msigdbr.csv` (msigdbr 26.1.1, collection C2:CP:KEGG_LEGACY). KEGG legacy `KEGG_FATTY_ACID_METABOLISM` is KEGG hsa00071 (fatty acid degradation, i.e. beta-oxidation). `KEGG_PURINE_METABOLISM` and `KEGG_PYRIMIDINE_METABOLISM` include the RNA/DNA polymerase subunits that KEGG places in those maps; they are kept as published.

| pathway | gene set | KEGG | genes in set | measured in the store |
|---|---|---|---|---|
| TCA cycle | KEGG_CITRATE_CYCLE_TCA_CYCLE | hsa00020 | 32 | 31 |
| Glycolysis / gluconeogenesis | KEGG_GLYCOLYSIS_GLUCONEOGENESIS | hsa00010 | 62 | 60 |
| Fatty-acid beta-oxidation | KEGG_FATTY_ACID_METABOLISM | hsa00071 | 42 | 42 |
| BCAA catabolism | KEGG_VALINE_LEUCINE_AND_ISOLEUCINE_DEGRADATION | hsa00280 | 44 | 44 |
| Purine metabolism | KEGG_PURINE_METABOLISM | hsa00230 | 159 | 159 |
| Pyrimidine metabolism | KEGG_PYRIMIDINE_METABOLISM | hsa00240 | 98 | 98 |

## Metabolite side

Members are RefMet names that exist in the MoTrPAC metabolomics rows (store/metab_features.parquet). Two rules are used:

* **RefMet class**: every MoTrPAC RefMet name with that `sub_class` (RefMet bulk download 2026-09-26, or the MoTrPAC 2021-05-05 named-metabolite dictionary for older names; see raw_metab/SOURCES.md). TCA cycle uses sub_class `TCA acids`; fatty-acid beta-oxidation uses sub_class `Acyl carnitines` restricted to names beginning `CAR` (this drops 3-Dehydroxycarnitine, a carnitine precursor).
* **Curated**: a short list of intermediates of the corresponding KEGG reference map (map00020, map00010, map00280, map00230, map00240). `build_pathway_map` refuses a curated name that is not in the MoTrPAC vocabulary. Each curated name is expanded to its notation variants (rat RefMet-2021 `CAR(3:0)` and human current-RefMet `CAR 3:0` share `name_key`).

The `source` column records each member's RefMet sub_class and KEGG compound id as annotated in store/raw_metab. Rat and human use different notations for acylcarnitines and some TCA acids (rat `Oxoglutaric acid`, human `2-Oxoglutaric acid`), so a set can hold two names for one compound; only one of them is measured in any column.

| pathway | metabolites | rule |
|---|---|---|
| TCA cycle | 10 | RefMet sub_class = TCA acids |
| TCA cycle | 2 | curated: TCA-cycle acyl-CoA (KEGG map00020 compound) |
| Glycolysis / gluconeogenesis | 10 | curated: glycolysis/gluconeogenesis intermediate (KEGG map00010 compound) |
| Fatty-acid beta-oxidation | 128 | RefMet sub_class = Acyl carnitines (names CAR(...) / CAR ...) |
| Fatty-acid beta-oxidation | 1 | curated: free carnitine (carnitine shuttle) |
| BCAA catabolism | 17 | curated: BCAA, branched-chain keto acid or BCAA-derived acyl-CoA/acylcarnitine (KEGG map00280) |
| Purine metabolism | 24 | curated: purine base/nucleoside/nucleotide or degradation product (KEGG map00230) |
| Pyrimidine metabolism | 18 | curated: pyrimidine base/nucleoside/nucleotide or degradation product (KEGG map00240) |

### Members

* **TCA cycle** (12): `2-Oxoglutaric acid`, `Acetyl-CoA`, `Citric acid`, `Citric acid/Isocitric acid`, `Fumaric acid`, `Isocitric acid`, `Malic acid`, `Oxaloacetic acid`, `Oxoglutaric acid`, `Succinic acid`, `Succinyl-CoA`, `cis-Aconitic acid`
* **Glycolysis / gluconeogenesis** (10): `Dihydroxyacetone phosphate`, `Fructose 1,6-bisphosphate`, `Fructose 6-phosphate`, `Glucose`, `Glyceraldehyde 3-phosphate`, `Hexose 6-phosphate`, `Lactic acid`, `Phosphoenolpyruvic acid`, `Phosphoglyceric acid`, `Pyruvic acid`
* **Fatty-acid beta-oxidation** (129): `CAR 10:0`, `CAR 10:0;OH`, `CAR 10:1`, `CAR 10:2`, `CAR 11:0`, `CAR 11:1`, `CAR 12:0`, `CAR 12:0;OH`, `CAR 12:1`, `CAR 12:2`, `CAR 12:3`, `CAR 13:0`, `CAR 13:1`, `CAR 14:0`, `CAR 14:0;OH`, `CAR 14:1`, `CAR 14:1;OH`, `CAR 14:2`, `CAR 14:2;OH`, `CAR 16:0`, `CAR 16:0;OH`, `CAR 16:1`, `CAR 16:2`, `CAR 16:3`, `CAR 17:0`, `CAR 18:0`, `CAR 18:0;OH`, `CAR 18:1`, `CAR 18:1;OH`, `CAR 18:2`, `CAR 18:3`, `CAR 20:0`, `CAR 20:1`, `CAR 20:2`, `CAR 20:3`, `CAR 20:4`, `CAR 24:0`, `CAR 26:0`, `CAR 26:1`, `CAR 2:0`, `CAR 3:0`, `CAR 4:0`, `CAR 4:0;3Me`, `CAR 4:0;OH`, `CAR 5:0`, `CAR 5:0;OH`, `CAR 5:1`, `CAR 6:0`, `CAR 6:0;OH`, `CAR 7:0`, `CAR 8:0`, `CAR 8:0;OH`, `CAR 8:1`, `CAR 9:0`, `CAR DC3:0`, `CAR DC3:0;2Me`, `CAR DC5:0`, `CAR DC6:0`, `CAR(10:0(OH))`, `CAR(10:0)`, `CAR(10:1)`, `CAR(10:2)`, `CAR(11:0)`, `CAR(12:0(OH))`, `CAR(12:0)`, `CAR(12:0)_lp_a`, `CAR(12:0)_lp_b`, `CAR(12:1)`, `CAR(13:0)`, `CAR(13:1)`, `CAR(14:0(OH))`, `CAR(14:0)`, `CAR(14:1(OH))`, `CAR(14:1)`, `CAR(14:2(OH))`, `CAR(14:2)`, `CAR(15:0)_lp_a`, `CAR(15:0)_lp_b`, `CAR(16:0(OH))`, `CAR(16:0)`, `CAR(16:1)`, `CAR(16:1)_lp_a`, `CAR(16:1)_lp_b`, `CAR(16:2(OH))`, `CAR(16:2)`, `CAR(16:2)_lp_a`, `CAR(16:3)`, `CAR(17:0)_lp_a`, `CAR(17:0)_lp_b`, `CAR(17:1)`, `CAR(18:0(OH))`, `CAR(18:0)`, `CAR(18:1(OH))`, `CAR(18:1)`, `CAR(18:1)_lp_a`, `CAR(18:2)`, `CAR(18:2)_lp_a`, `CAR(18:3)_lp_a`, `CAR(19:0)`, `CAR(20:0)`, `CAR(20:1(OH))`, `CAR(20:1)`, `CAR(20:2(OH))`, `CAR(20:2)_lp_a`, `CAR(20:4)`, `CAR(21:2)`, `CAR(22:2)`, `CAR(22:4)`, `CAR(24:0)`, `CAR(24:1)`, `CAR(26:0)`, `CAR(26:1)`, `CAR(2:0)`, `CAR(3:0)`, `CAR(4:0(OH))`, `CAR(4:0)`, `CAR(5:0(OH))`, `CAR(5:0)`, `CAR(5:1)`, `CAR(6:0)`, `CAR(7:0)`, `CAR(8:0(OH))`, `CAR(8:0)`, `CAR(8:1)`, `CAR(9:0)`, `CAR(DC3:0(2Me))`, `CAR(DC3:0)`, `CAR(DC5:0)`, `Carnitine`
* **BCAA catabolism** (17): `3-Methyl-2-oxovaleric acid`, `CAR 3:0`, `CAR 5:0`, `CAR 5:0;OH`, `CAR 5:1`, `CAR(3:0)`, `CAR(5:0(OH))`, `CAR(5:0)`, `CAR(5:1)`, `Isoleucine`, `Isovaleryl-CoA`, `Ketoisovaleric acid`, `Ketoleucine`, `Leucine`, `Leucine/Isoleucine`, `Propionyl-CoA`, `Valine`
* **Purine metabolism** (24): `3',5' cyclic AMP`, `ADP`, `AMP`, `ATP`, `Adenine`, `Adenosine`, `Adenylsuccinic acid`, `Allantoin`, `Deoxyadenosine`, `Deoxyguanosine`, `GDP`, `GMP`, `GTP`, `Guanine`, `Guanosine`, `Hypoxanthine`, `IMP`, `Inosine`, `Phosphoribosyl pyrophosphate`, `Uric acid`, `XMP`, `Xanthine`, `Xanthosine`, `dATP`
* **Pyrimidine metabolism** (18): `CDP`, `CMP`, `CTP`, `Cytidine`, `Cytosine`, `Deoxycytidine`, `Deoxyuridine`, `Thymidine`, `UDP`, `UMP`, `UTP`, `Uracil`, `Ureidopropionic acid`, `Uridine`, `beta-Alanine`, `dCTP`, `dTTP`, `dUMP`

## Scores

`metab.pathway_panel`: per pathway and per comparison column that has RNA, protein and metabolite data at the same contrast (rat F/M 1, 2, 4, 8 wk trained vs sedentary; human VL EE/RE minus control at 15-45 min, 3.5-4 h, 24 h), `legacy.camera_pr` (port of limma::cameraPR, inter-gene correlation 0.01) of the set vs every other measured feature of that column. t > 0 = the set is up with exercise relative to the rest of the column. Cells with fewer than 3 measured members are NaN.
