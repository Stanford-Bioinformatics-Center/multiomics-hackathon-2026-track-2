# Reusable MoTrPAC query core

This Python package matches a normalized disease signature against published
MoTrPAC contrast summaries. It is independent of PAH study IDs, feature counts,
and file locations. The output is an auditable cross-study comparison, not an
exercise effect measured in patients.

## Build a reference from the bundled exports

From the containing `generalized/` folder, when it is placed inside the existing `MoTrPAC Hackathon` folder:

```powershell
python -m query_core.adapter_legacy `
  --legacy-root .. `
  --output query_core/motrpac_reference.csv.gz `
  --include-biocarta
```

This reads the existing human muscle RNA/protein, blood RNA/Olink,
metabolomics, and full blood GOBP tables. The optional BioCarta table has only
the **six PAH-selected sets**; its rows carry `reference_scope=partial...`.
Neither table requires R at query time. The legacy muscle RNA export has only
24-hour rows, and the bundled exports contain no PTM results. A separate full
muscle export can be merged when available:

```powershell
python -m query_core.merge_reference `
  --input query_core/motrpac_reference.csv.gz `
  --input path/to/full_muscle_reference.csv.gz `
  --output query_core/motrpac_reference_full.csv.gz `
  --audit query_core/replaced_reference_rows.csv `
  --expected-package-version 2.0.8 --expected-collection c2.0
```

The last input wins on a duplicate source feature/assay/tissue/time/contrast
key; **every replacement is logged**. Missing or mixed package/collection
versions fail. The full muscle export must use this package's normalized
reference schema and carry `source_feature_id`, `source_package_version`, and
`source_collection`.

An optional `--metabolite-feature-name-fallback` treats a MoTrPAC metabolite
`feature_id` as a RefMet name only when both `refmet_id` and `refmet_name` are
blank. Its rows carry `match_basis=feature_id_name_fallback`, and query matches
are marked `name_only`. It is off by default. The bundled c2.0 metabolite
export has 13,455 rows where both RefMet fields contain the literal missing
token `NA`. The adapter normalizes these to blanks. Only with the opt-in flag
are those rows offered as lower-confidence names.

## Input contract

CSV and CSV.gz are accepted. Both disease and reference rows require:

| Column | Meaning |
| --- | --- |
| `study_id` | Source study or release label |
| `species` | `human` or `rat` |
| `tissue` | Tissue label, normalized to lower case |
| `layer` | `rna`, `protein`, `phosphosite`, `metabolite`, or `pathway` |
| `feature_id` | Identifier in `id_namespace` |
| `id_namespace` | `HGNC` (human gene symbol), `RatGeneSymbol`/`RGD`, `UniProt`, `RefMet`, `RefMetName`, `GO`, or `BIOCARTA`; opaque IDs may be retained but do not match directly |
| `log2_fc` | Signed numerator minus denominator log2 fold change; blank is allowed for pathways without this estimate |
| `statistic` | Signed test statistic, if known; blank allowed |
| `p_value`, `q_value` | Original study values; blank allowed. Do not fill `q_value` with a value from a different testing family |
| `contrast` | Explicit numerator-minus-denominator text |

The reference also requires `timepoint` and `contrast_category` (for example
`EE-CON`). Optional typed mapping columns are `gene_symbol`, `uniprot`,
`refmet_id`, and `refmet_name`. Optional context columns include `biospecimen`,
`effect_scale`, `statistic_type`, `source_feature_id`, `source_q_value`,
`source_q_scope`, `assay`, `platform`, `source_package_version`, and
`source_collection`. Pathways require `pathway_collection` as well as a
collection-specific `feature_id`; for example, `C2/BIOCARTA` plus
`BIOCARTA_TOB1_PATHWAY`. A pathway name alone is not joined across collections.
MoTrPAC's GOBP `feature_id` is its set label, not a stable GO accession.
`pathway_membership_version` may be supplied on both sides to assert identical
gene-set membership; matching names alone do not assert this.

The interface does not interpret `HGNC:12345` accessions as gene symbols;
`HGNC` here means an approved human symbol. Map other identifier forms
upstream. It does not infer orthology from similar names.

## Query

```powershell
python -m query_core.cli `
  --disease examples/signatures/pah_blood_rna.csv.gz `
  --reference query_core/motrpac_reference.csv.gz `
  --out query_core/demo_pah_blood `
  --tissue blood --reference-contrast-category EE-CON
```

Optional filters: `--species`, `--reference-species`, `--tissue`,
`--disease-contrast`, `--reference-contrast-category`, `--reference-contrast`,
and repeatable `--timepoint`. By default all available disease contrasts and
reference groups are analyzed separately. `--rank-field auto` uses signed
statistics only when every unique pair in a group has one; otherwise it uses
signed log2 fold change. It reports Spearman rho only for at least three
nonduplicated, nonconstant pairs, with no inferential P value. Choose
`--rank-field statistic` or `--rank-field log2_fc` explicitly if desired.

For human-to-rat queries, provide an **explicit one-to-one** map and do not
relabel the disease rows:

```powershell
python -m query_core.cli --disease disease.csv --reference rat_reference.csv.gz `
  --out results/rat --species human --reference-species rat `
  --ortholog-map orthologs.csv
```

`orthologs.csv` needs `source_species,source_gene_symbol,target_species,target_gene_symbol`.
Conflicting one-to-many or many-to-one rows fail. Cross-species matches are
marked `ortholog_gene_link`; metabolite names are never used across species.

## Outputs

- `matches.csv.gz`: one row per candidate reference match, with matched typed
  identifiers, mapping status, biospecimen comparison, source effects and
  statistics, signed relationship, and source provenance.
- `coverage.csv`: one row for **every disease input row** after filtering.
  `unmatched.csv` gives explicit reasons; `ambiguous.csv.gz` retains every
  one-to-many candidate instead of choosing a convenient result.
- `summary_by_group.csv`: counts unique, ambiguous, and unmatched disease rows
  by reference layer, time, and contrast, including groups with zero matches.
- `rank_correlation.csv`: descriptive Spearman results and all exclusion counts.
- `summary.json`: source paths/hashes, parameters, per-study/species release
  audit, counts, and interpretation limits. Different known package or
  collection releases within one reference study/species fail; blank release
  metadata stays visible as unknown.

`mapping_status=unique` means one reference feature in that **group**, not a
unique biological mechanism. `match_confidence` describes identifier
resolution: `high`, `name_only`, or `ortholog`. `match_relation` distinguishes
`same_analyte`, `gene_linked_same_layer`, `gene_linked_cross_layer`,
`protein_site_link`, `protein_linked_site`, `metabolite_name_only`,
`ortholog_gene_link`, and `pathway_name_overlap`. A protein linked to RNA by a
gene symbol is a cross-layer gene link, not the same measured analyte. Two
different proteins or phosphosites linked only by a gene symbol are
`gene_linked_same_layer`, while sites sharing only a parent UniProt accession
are `protein_linked_site`. A
matching BioCarta name is `pathway_name_overlap` unless both sources explicitly
declare the same pathway membership version. `reference_timepoint_order` is
supplied for the bundled MoTrPAC times; generic references may leave it blank.
`match_confidence=high` describes identifier resolution only; it is not a
statement that two assays, cohorts, or pathway definitions are biologically
equivalent.

`sign_relation=same/opposite` is descriptive for signed estimates. The
`supported_direction` label requires both source q values below `--q-threshold`
(default 0.05) and both nonzero effects exceeding `--min-abs-log2-fc` (default
0). Missing q gives `evidence_status=missing_q`; q above threshold means
`not_jointly_supported`, **not** no change. For pathway-only signed statistics,
the source log2 fold change stays blank and `sign_basis=signed_statistic`.

Only same-species, same-tissue typed identities match by default. Different
`biospecimen` values are retained and flagged. Gene-to-metabolite pathway
bridges and PTM occupancy are outside this direct-identity query.

## Test

```powershell
python -m unittest discover -s query_core/tests -v
```
