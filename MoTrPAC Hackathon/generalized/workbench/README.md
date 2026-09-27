# Named Metabolomics Workbench study adapter

`workbench_adapter.py` converts **local named-metabolite** Workbench `/factors` and `/data` JSON files to a signed disease-feature table for the shared MoTrPAC comparison pipeline. It does not choose a disease group, a comparator, a tissue, a measurement scale, or a pairing rule from study names. A researcher records those choices in a reviewed design JSON.

## Fetch by accession

```sh
python fetch_workbench.py ST000763 --output-dir downloaded_workbench
```

This optional command accepts only a validated `ST` accession, retrieves the official Workbench named `/factors` and `/data` JSON endpoints, and writes `downloaded_workbench/ST000763/{factors,data,source_manifest}.json`. The manifest records exact URLs, retrieval time, byte counts, and SHA256 hashes. Existing files require `--overwrite`. The fetcher does not choose a contrast or calculate statistics. Use the offline adapter below after reviewing the study factors and design.

## CLI contract

```sh
python workbench_adapter.py --factors factors.json --data data.json \
  --design study_design.json --output disease_signature.csv
```

The two inputs are the JSON forms of the official `study/study_id/{ST_ID}/factors` and `study/study_id/{ST_ID}/data` endpoints, downloaded separately or obtained with `fetch_workbench.py`. No network request is made by the adapter. Paths may be absolute or relative to the command's working directory. The output directory is created if needed. The main CSV contains rows with a RefMet ID or standardized RefMet name. A sibling `disease_signature_unmatched.csv` contains tested features without either identifier, for coverage review. A `disease_signature_duplicates.csv` sidecar records all source features that map to repeated canonical IDs and which one was retained. All output rows retain `source_feature_id = analysis_id:metabolite_id`. The adapter does not guess a chemical identity for unmatched features.

Required output columns are `study_id,species,tissue,layer,feature_id,id_namespace,log2_fc,statistic,p_value,q_value,contrast`. `layer` is `metabolite`; `species` is `human` or `rat`. `tissue` uses the MoTrPAC anatomical value, such as `blood`; optional `biospecimen` records a more specific source such as `plasma`. A RefMet ID, when present in the source record, becomes `feature_id` with `id_namespace=RefMet`. Otherwise a standardized `refmet_name` becomes `feature_id` with `id_namespace=RefMetName`. Additional columns carry source IDs, names, units, arm sample counts, test status, and a duplicate canonical ID flag. **A `RefMetName` match is an exact standardized-name link, not proof that two platforms measured the same isomer.** Repeated canonical IDs require an explicit `duplicate_policy`: `error` (default), `keep_all` (requires downstream many-to-one handling), or `highest_coverage` (test-eligible feature first, then largest observed arm sample count, then source ID). The last rule does not use effect size or P value, and the duplicate sidecar preserves alternatives for review.

The signed `log2_fc` is numerator minus denominator. For `unpaired`, it is mean log2 case abundance minus mean log2 control abundance at the selected time. For `paired_change`, it is mean within-subject endpoint minus baseline in the selected case group. For `difference_in_changes`, it is the case mean within-subject change minus the control mean within-subject change. The `statistic` is the corresponding t statistic; `p_value` comes from Welch's two-group test or a one-sample t test of paired differences. `q_value` is Benjamini–Hochberg across all eligible source features in this run, **before** dropping unmatched IDs or duplicate names. Features below the configured sample-coverage gate remain in the CSV with a blank statistic/P/q and `status=insufficient_coverage`. For paired designs, the coverage denominator counts only subjects with **both** designated samples, regardless of whether a particular feature was measured in those samples; a lone rest or peak sample cannot enter a paired contrast. `n_design_case` and `n_design_control` expose those denominators alongside feature-specific measured-pair counts.

## Study design JSON

The supplied [ST000763 resting example](st000763_rest.example.json) exercises the unpaired route. Required top-level fields are `study_id`, `species`, `tissue`, `contrast`, `design`, and `normalization`. `species` currently accepts `human` and `rat`, the MoTrPAC reference species. Optional `biospecimen` records the actual source (`plasma` in the example). Optional `analysis_ids` selects specific Workbench assays. Mixed source `units` are rejected unless the reviewer explicitly sets `allow_mixed_units: true`; the safer choice is to analyze compatible assays separately.

`design.type` accepts:

- `unpaired`: `group_factor`, `case_values`, `control_values`, and optionally `time_factor` plus `time_value`.
- `paired_change`: `group_factor`, `case_values`, `time_factor`, `baseline_value`, `endpoint_value`, and a verified subject identifier.
- `difference_in_changes`: the paired fields plus `control_values`.

Group values are nonempty arrays, so a reviewer can define a compound comparator. Factor names are matched without regard to letter case; factor **values** are exact. For paired designs, provide exactly one of `subject_id_factor` (a Workbench factor name) or `subject_map_csv` (a CSV with `sample_id,subject_id`, relative to the design file or absolute). Sample order and adjacent sample numbers are never used to infer pairs. Missing subject IDs, a subject with multiple samples at the selected time, or a subject that changes group cause a validation error. Subjects without both selected times contribute no paired difference. Sample counts and coverage are exposed in the output.

`normalization.input_scale` must be `linear` or `log2`. Linear values are transformed with `log2(value + pseudocount)`; nonpositive transformed inputs are rejected. Already log2 values are used as supplied and cannot receive a pseudocount. Optional `sample_center` is `none` (default) or `median_within_analysis`; centering happens on the log2 scale within each Workbench analysis. The adapter does not claim to remove batch effects. `min_samples_per_arm` defaults to 3 and `min_fraction_per_arm` to 0.7. These are eligibility gates, not a substitute for a study-specific statistical plan.

## Local validation

Run the synthetic tests with:

```sh
python -m unittest discover -s tests -v
```

For the checked local ST000763 files, use `st000763_rest.example.json` with their `factors.json` and `data.json`. The example deliberately uses **resting samples only**: ST000763's local files do not provide verified subject IDs, so the adapter rejects a paired rest-to-peak analysis. In this study, the PAH samples were collected in the catheterization setting and healthy samples noninvasively. The contrast is therefore confounded by disease, systemic sclerosis, and collection setting. It must not be interpreted as an isolated PAH effect.

## Scope of “any Workbench study”

The [Workbench REST specification](https://www.metabolomicsworkbench.org/tools/MWRestAPIv1.2.pdf) provides stable endpoints for study summary, factors, named metabolites and data, plus separate untargeted-data endpoints. This adapter supports the **named `/data` JSON shape** only. The [mwTab format specification](https://www.metabolomicsworkbench.org/data/mwTab_specification.pdf) says subject IDs are optional and sample factors use study-defined `NAME:VALUE` entries. A generic parser can read the envelope; it cannot safely infer every study's biological contrast, paired design, normalization or metabolite identity. Some studies expose only untargeted m/z features or data files outside mwTab, which need separate ingestion and chemical-annotation steps.

The [CFDE documentation](https://cfde.cloud/data/documentation) describes C2M2 as a standard for **discovery metadata**. The [NIH CFDE FAQ](https://commonfund.nih.gov/dataecosystem/faqs) says CFDE does not replicate or directly serve the underlying program datasets. CFDE can help find a study; the analysis-ready named data should be obtained from Metabolomics Workbench and reviewed with a study-specific design file.
