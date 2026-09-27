# Ranked-feature pathway enrichment

This module accepts a **complete measured-universe ranking** for one disease contrast and a GMT collection. It produces an exploratory competitive pathway score and a canonical pathway signature that can be submitted to `query_core` when the MoTrPAC reference contains the same named collection. The input may be RNA (`HGNC`, `RatGeneSymbol`, or `RGD`) or metabolite (`RefMet` or `RefMetName`). Pathway matching requires the same collection and name; it does not establish identical membership unless versions are verified independently.

```powershell
python pathways/ranked_enrichment.py `
  --signature full_rank_disease.csv.gz `
  --gmt gene_sets.gmt.gz `
  --universe-provenance universe.json `
  --gmt-namespace HGNC `
  --pathway-namespace GO `
  --collection C5/GOBP `
  --score-field statistic `
  --out-dir results/pathways
```

The signature uses the shared canonical columns: `study_id,species,tissue,layer,feature_id,id_namespace,log2_fc,statistic,p_value,q_value,contrast`. It must contain one study, species, tissue, layer, namespace, and contrast, with one nonempty identifier and finite signed ranking score per **measured feature**. Selected hits or only significant rows are invalid. For multiple contrasts, run each full ranking separately. `--score-field log2_fc` is available when a signed test statistic is unavailable. The GMT member namespace must equal the signature namespace; map identifiers upstream and record the mapping method in `feature_id_policy`. Matching only trims/collapses whitespace and ignores case. It does **not** infer metabolite synonyms, isomers, adducts, or RefMet IDs from names.

`universe.json` must declare how completeness was established. Example:

```json
{
  "source": "study accession and analysis release",
  "complete_measured_universe": true,
  "expected_feature_count": 19593,
  "inclusion_rule": "all QC-passing measured genes, including nonsignificant genes",
  "feature_id_policy": "one declared gene symbol per gene after a documented probe rule",
  "gmt_source": "gene-set provider and download location",
  "gmt_version": "release identifier or access date"
}
```

The code checks shape, uniqueness, finite scores, namespace and an optional expected count. **It cannot prove from a CSV that no nonsignificant features were omitted**; the user declaration is preserved in `manifest.json` with source file hashes.

For each GMT set, the test compares its member ranks with ranks of all other measured features. Its signed statistic is an asymptotic Wilcoxon rank-sum Z with tie correction; positive means members tend to rank higher in the disease contrast. `rank_biserial` is a separate descriptive rank effect. A two-sided p value is calculated under exchangeability of feature labels within the measured universe. BH q values are adjusted over **all eligible sets in the supplied GMT for that contrast**, including nonsignificant ones. Defaults require 10–500 measured members; excluded sets and overlap counts remain in `pathway_tests.csv`. `pathway_signature.csv` contains only eligible tests and keeps `log2_fc` blank because a pathway rank shift is not a molecular fold change. Test type, score field, set collection, q scope and source hashes are recorded.

**Interpretation limit:** feature scores can be correlated, especially genes in the same pathway. This simple competitive null does not account for that correlation, so its p and q values are exploratory and may be anti-conservative. Wu and Smyth's [CAMERA paper](https://academic.oup.com/nar/article/40/17/e133/2411151) shows why even modest inter-gene correlation can inflate simpler competitive test results; a [Mann–Whitney–Wilcoxon gene-set methods paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9140214/) describes the rank-sum approach used here. This module is neither CAMERA nor the older PAH PAGE analysis, and their score magnitudes should not be subtracted or compared. For a defensible inferential disease pathway analysis, rerun from sample-level expression with a method that models inter-gene correlation, using an appropriate design and full GMT. The same limitation applies to correlated metabolites, and RefMet-name matches can be especially incomplete or ambiguous.

## Optional R adapter: `limma::cameraPR`

`camera_pr.R` runs `limma::cameraPR` on the same full-rank RNA signature, GMT and universe-provenance inputs (RNA only; requires R with `limma`, `jsonlite`, `digest`):

```powershell
Rscript pathways/camera_pr.R `
  --signature full_rank_disease.csv.gz --gmt gene_sets.gmt.gz `
  --universe-provenance universe.json --gmt-namespace HGNC `
  --pathway-namespace GO --collection C5/GOBP --out-dir results/camera_pr `
  [--score-field statistic] [--min-size 10] [--max-size 500] [--inter-gene-cor 0.01]
```

It writes the same `pathway_tests.csv`, `pathway_signature.csv` and hashed `manifest.json` layout, with the cameraPR two-sided P converted to a signed Z (sign from `Direction`) and BH q over eligible sets. cameraPR penalizes the test for a **preset** inter-gene correlation (default 0.01) because a pre-ranked table carries no sample-level covariance; it does not estimate set-specific correlation, so treat the preset as a sensitivity parameter and consider rerunning with larger values. It is therefore closer to CAMERA than the rank-sum test above, but it is **not** the legacy sample-level CAMERA analysis, and its Z is not comparable to PAGE or rank-sum scores.

**Validation status:** only a synthetic smoke test (in `tests/`, skipped when R or the packages are absent) has been run. No real disease ranking with an official GOBP GMT has been validated through this adapter yet.

Test with `python -m unittest discover -s pathways/tests -v` from the package root.
