# Stanford Multi-omics Hackathon 2026 Track 2

## Omic Discordance Explained

This project answers the Track 2 **"Omic Discordance Explained"** challenge: omic layers within the MoTrPAC data are often discordant, so the goal is to compare two or more compatible layers — RNA, total protein, and post-translational modifications (PTM) — within a defined tissue and identify a model that explains when they agree and when they disagree.

The headline deliverable is the **discordance catalog and its predictive model**: a per-feature catalog that classifies each event as `supported_concordant`, `supported_opposite`, `rna_response_protein_equivalent`, or `indeterminate`, paired with an out-of-fold ridge regression that compares a zero-change baseline, an RNA-only model, and a temporal model. The out-of-fold R² is near zero, so the model is reported honestly as a **weak prediction — a result, not a success claim**. A committed muscle EE-CON demo is served read-only through the app so the catalog and model can be explored without recomputing any statistics.

Reported results are **cross-cohort association findings, not a claim of PAH treatment efficacy**. Cohorts are **analyzed separately and are not merged** (muscle protein, plasma metabolomics, and blood are kept distinct). The **PTM occupancy caveat applies to all reported findings**: PTM signal is not automatically a measure of modification occupancy, enzyme activity, or functional consequence.

## Deliverables

1. **Discordance catalog + predictive model (headline).** The concordant/discordant event catalog and the out-of-fold ridge model of discordance, surfaced as the app's default Discordance view (query builder → data visualization → discordance interpretation) over the committed muscle EE-CON demo.
2. **PTM-parent divergence audit.** A phosphosite-versus-parent-protein divergence summary per timepoint, carrying the PTM occupancy caveat.
3. **Directed and generalized signature analyses.** The directed male-rat SKM-GN muscle-protein signature and the generalized query_core path, each rendered with fully labeled, cohort-separated results.

## Data

Compatible MoTrPAC transcriptomic, proteomic, and PTM measurements with protein, site, pathway, and interaction annotations, analyzed within a defined tissue. The discordance demo currently served is a single committed muscle EE-CON `post_24_hr` run plus its PTM-parent audit; live recompute of the catalog from an arbitrary tissue/contrast/timepoint is a documented follow-on, not part of this MVP.

## Honest framing

The directed headline result — male rat SKM-GN protein at 8 weeks — is **q = 0.0584 over a frozen 52-column BH multiplicity family and is labeled not significant**. The 16-column `q = 0.0413` value appears only as an explicitly labeled sensitivity analysis, never as a headline. All statistics originate in the engines (mprobe, query_core, and the standalone discordance modules); the API and web layers only serve, shape, and render.

> [!IMPORTANT]
> PTM signal is not automatically a measure of modification occupancy, enzyme activity, or functional consequence.

## Documentation map

- [WORKFLOWS.md](WORKFLOWS.md) — every analysis workflow traced through seven stages, with the API request→response contract and the point→evidence→source chain.
- [DECISIONS.md](DECISIONS.md) — architecture decision records explaining the reasons behind each choice.
- [VERIFICATION.md](VERIFICATION.md) — the clean-clone reproduction log with exact commands, counts, and hashes per gate.
- [METABOLOMICS.md](METABOLOMICS.md) — the metabolomics case-study module, its relationship to the engine, and the convergence cross-check.
- [RUN_LOCAL.md](RUN_LOCAL.md) — how to run the API and web app locally, including the exact ports.

## Setup

For a fresh clone, [BOOTSTRAP.md](BOOTSTRAP.md) documents the one-command bootstrap (`bash scripts/bootstrap.sh`) that editable-installs the engine and API, fetches the mprobe store, and runs a smoke check within 600 seconds; it also documents the `parents[3]` path-coupling constraint. The `RUN_LOCAL.md` guide (linked above) then covers the two-process local run — the FastAPI service plus the Vite web app — and the exact port wiring. The engine and API install as editable Python packages, and the web app installs from the committed lockfile; the reproduction log linked in the documentation map records the counts reproduced at each step.
