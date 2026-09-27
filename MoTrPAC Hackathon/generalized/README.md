# Reusable MoTrPAC cross-omic analysis engine

This folder is the app-facing backend for the Track 2 project. It accepts processed disease-study effects, matches them to compatible MoTrPAC results, writes machine-readable audits and plots, and builds a separate within-MoTrPAC RNA–protein discordance catalog and predictive pilot. The earlier numbered PAH scripts remain reproducible source-specific examples in the parent folder.

**Current boundary:** A generic analysis engine can accept any correctly described processed signature; it cannot infer an arbitrary study's biological contrast, sample pairing, normalization, or metabolite identity. Unsupported comparisons appear in coverage and mapping outputs. The Workbench adapter supports named `/data` JSON with an explicit study-design file. No GUI is included here; an app can invoke these CLIs or import their Python functions.

## Setup

Use Python 3.12 and `python -m pip install -r requirements.txt`. The query core itself uses only the Python standard library; the Workbench adapter, plots, and predictor use the packages in `requirements.txt`. The human MoTrPAC full-layer exporter needs R 4.4+ and `MotrpacHumanPreSuspensionAnalysis` 2.0.8/c2.0. The bundled-reference adapter runs without R.

## Processed signature schema

One row is one measured feature in one signed comparison. Required CSV headers are:

`study_id,species,tissue,layer,feature_id,id_namespace,log2_fc,statistic,p_value,q_value,contrast`

`layer` is `rna`, `protein`, `phosphosite`, `metabolite`, or `pathway`. The positive direction is the named numerator minus denominator; `log2_fc` must not be fabricated for studies reporting only a signed statistic. Numeric cells may be blank when unavailable. Declare `id_namespace` (`HGNC`, `UniProt`, `RatGeneSymbol`, `RefMet`, `RefMetName`, or a pathway collection) and supply optional `gene_symbol`, `uniprot`, `refmet_id`, `refmet_name`, `pathway_collection`, `biospecimen`, `effect_scale`, `statistic_type`, and provenance/caveat columns as available. `pathway_collection` is required for pathway matching. MoTrPAC reference rows add `timepoint` and `contrast_category`.

Examples generated from the existing repo are in `examples/signatures/`. `examples/make_pah_signatures.py` reconstructs them and keeps caveats visible. The GSE33463 gene table retains its original **probe-level** adjusted P in `source_q_value`, while canonical `q_value` is blank because it is not formal gene-level FDR. The Malenfant nine-protein table keeps the paper-era symbols and two explicit old-to-current mappings in `examples/pah_gene_aliases.csv`. ST000763 examples include both PAH-versus-healthy and same-setting PAH-versus-normal-pressure comparisons; the former confounds disease, systemic sclerosis, and sampling setting.

## Quick start: query the bundled reference

Run from this folder (the parent `..` is `MoTrPAC Hackathon`):

```sh
python run_config.py --config examples/pah_jobs.example.json --out demo_all
```

The example config runs five independent PAH demonstration jobs using the included full human muscle and legacy blood/metabolite reference, plus the rat training protein reference. It writes `demo_all/manifest.json`, one query folder per job, and plots where requested. An app can generate the same JSON from uploaded processed signatures; its paths are resolved relative to the config file. Each job reports its own failure without hiding completed jobs. `tests/test_run_config.py` uses an unrelated synthetic study ID to check this interface.
Use a fresh `--out` directory for each run; the manifest identifies current results, but older files can remain if an existing directory is reused.

To rebuild the bundled portion of the human reference from the older scripts, run `python -m query_core.adapter_legacy --legacy-root .. --output query_core/motrpac_reference.csv.gz --include-biocarta`. The checked-in `query_core/motrpac_reference_full.csv.gz` then adds human muscle RNA, total protein, and phosphosites from the official c2.0 package, with duplicate replacements audited in `query_core/reference_duplicate_audit.csv.gz`. The bundled portion contains blood RNA/Olink, blood and muscle metabolomics, full blood GOBP, and a **partial six-set** BioCarta table. Neither file contains every human MoTrPAC tissue/assay/contrast. Run the package exporter below for additional compatible results. The reference rows retain original assay feature IDs and release metadata.

The query writes `matches.csv.gz`, `coverage.csv`, `unmatched.csv`, `ambiguous.csv.gz`, `summary_by_group.csv`, `rank_correlation.csv`, and `summary.json`. Candidate matches include typed relation (`same_analyte`, `gene_linked_cross_layer`, `ortholog_gene_link`, pathway-name overlap, or metabolite-name-only), mapping confidence, biospecimen compatibility, signed relation, and evidence status. An ambiguous feature is never silently collapsed into a unique hit. The plot renderer writes audited bubble, rank, and coverage figures plus underlying plot data.

## Full human muscle layers and the Track 2 pilot

The included full reference contains the official summary-level human muscle EE-CON RNA, total protein, and phosphosite effects. Rebuild the 24-hour pilot with:

```sh
python discordance/build_discordance.py --input query_core/motrpac_reference_full.csv.gz --out-dir discordance_ee_24h --tissue muscle --contrast-category EE-CON --target-time post_24_hr --earlier-times post_15_30_45_min post_3.5_4_hr
```

This reproduces the included demo metrics exactly. Use `discordance/export_motrpac_human.R` to export another compatible assay, tissue, or exercise mode. See [discordance/README.md](discordance/README.md) for the R command, model rules, and output definitions. The catalog calls a pair *supported concordant* or *supported opposite* only under declared effect-size, q, and CI rules. A large protein P value is never called no change; an RNA-supported/protein-equivalent class requires the protein CI wholly within a predeclared small-effect band. The predictor compares held-out protein log2 effects from zero, RNA-only, and RNA plus earlier RNA/PTM features, using gene-grouped folds and out-of-fold metrics.

In the current human-muscle EE-CON 24-hour pilot, 5,586 unambiguous RNA/protein pairs enter the model. The temporal model's out-of-fold R² is approximately 0.012 versus approximately 0.010 for RNA-only, with essentially identical MAE. These are **weak** predictions. The catalog's 285 RNA-supported/protein-equivalent events are descriptive candidates, not demonstrated synthesis or degradation mechanisms. The PAH examples must not be used to train or select this model. See [TRACK2_DESIGN.md](TRACK2_DESIGN.md) for the scientific pitch and validation criteria.

The companion phosphosite-to-parent audit finds 735 early EE-CON phosphosite feature rows with supported changes and an equivalently small measured parent-protein effect, spanning 350 genes under strict one-parent mapping. This is a stronger *descriptive* Track 2 demonstration of layer divergence; it is not measured phosphorylation occupancy or a validated predictor. See `discordance/demo_ptm_parent_ee/` and the audit CLI in `discordance/README.md`.

## Workbench study ingestion

```sh
python workbench/fetch_workbench.py ST000763 --output-dir workbench_downloads
python workbench/workbench_adapter.py --factors workbench_downloads/ST000763/factors.json --data workbench_downloads/ST000763/data.json --design workbench/st000763_rest.example.json --output st000763_signature.csv
```

The fetcher accepts only an ST accession, downloads the official named-data/factors endpoints, and saves URLs, retrieval time, and SHA256 hashes. The adapter requires a reviewed contrast/design JSON. Paired analyses require real subject IDs; the existing ST000763 public files lack them, so the example is an unpaired resting contrast. Untargeted m/z features without an audited chemical annotation stay unmatched. See [workbench/README.md](workbench/README.md).

## Full-ranking pathway analysis

`pathways/ranked_enrichment.py` accepts a **complete measured feature universe** from RNA or RefMet-identified metabolomics plus a GMT collection and a reviewed universe-provenance JSON. It returns an auditable pathway table and canonical pathway signature that can be queried against a compatible MoTrPAC pathway collection. Its competitive rank-sum P/q values are exploratory because correlated genes or metabolites violate the simple label-exchangeability assumption; this is not the legacy CAMERA result. For RNA, an optional R adapter `pathways/camera_pr.R` runs `limma::cameraPR` on the same inputs with a preset inter-gene correlation; it has passed only a synthetic smoke test and has not been validated on a real disease ranking. See [pathways/README.md](pathways/README.md) for required inputs and method choices.

## Rat training reference

`examples/export_rat_training_reference.R` exports the **full** bundled gastrocnemius protein table and a one-to-one human-to-rat ortholog map from the official rat objects. The cached example files are `examples/rat_training_reference.csv.gz` and `examples/human_to_rat_orthologs.csv`. Human disease signatures can query rat training only with `--reference-species rat --ortholog-map examples/human_to_rat_orthologs.csv`; the output retains both species and reports the ortholog relation. Rat training spans weeks and is context, not direct validation of the human acute model. The current bundled rat reference has protein only; other rat omes require separate official exports.

## Five former analyses, one reusable structure

| Former example | Shared capability | Study-specific adapter that remains |
| --- | --- | --- |
| Human muscle nine PAH proteins | Typed protein/RNA/time match and discordance context | Paper Table 2 transcription and old symbol crosswalk |
| Rat PAH protein comparison | Explicit ortholog mapping and training timecourse lookup | Rat package source files and intervention metadata |
| Cheadle BioCarta pathways | Exact database/collection pathway-name comparison and bubble plot | Published selected PAGE rows; no inferred fold change |
| GSE33463 blood RNA | Full-rank RNA/Olink matching, time-specific Spearman relation, and optional full-universe pathway testing | GEO platform probe mapping and IPAH-minus-healthy fit |
| ST000763 metabolomics | Reviewed Workbench design, RefMet/name match, drift/contrast context | Study factors, measurement scale, comparator and confounder choices |

The generic engine does **not** fit raw RNA/protein matrices. To add a new platform, implement a small source adapter that emits the schema above, rather than embedding disease names and fixed sample counts in analysis kernels. The old counts remain regression checks for the old examples, not validity rules for new studies.

## Verify

```sh
python -m unittest discover -s query_core/tests -v
python -m unittest discover -s workbench/tests -v
python -m unittest discover -s plots/tests -v
python -m unittest discover -s discordance -p 'test_*.py' -v
python -m unittest discover -s pathways/tests -v
python -m unittest discover -s tests -v
```

For a quick falsification check, replace `examples/signatures/pah_muscle_protein.csv.gz` with a different processed signature using the same schema. The engine should report its measured overlap and unsupported analyses without PAH-specific assertions. The PAH blood rank test reproduces the six legacy Spearman values to numerical precision; this verifies the adapter and ranking, not the biological interpretation.
