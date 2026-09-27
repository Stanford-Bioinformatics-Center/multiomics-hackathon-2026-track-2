# Handoff — next feature development (priority order)

Hands off from the integration work (branch `integration/exercise-signature-explorer`, pushed to
origin) to the next round. Grounded in the actual code so the next planning session starts from facts.

---

## Current state (what exists and works)

**Branch:** `integration/exercise-signature-explorer` (pushed; upstream
`origin/integration/exercise-signature-explorer`; 29 commits ahead of `main`). Do NOT rebase onto or
push over `origin/t3code/build-motrpac-probe-tool` (teammate's branch).

**Three pieces, one repo:**
- `hackathon/tool/` — the mprobe engine (scientific source of truth for the directed-signature path).
  122 pass / 34 skip (allowlist). Store at `hackathon/tool/store` (hash d615c2f0b9aec182). Dev venv at
  `hackathon/tool/.venv`.
- `apps/api/` — FastAPI service wrapping the engine + two adapters. 44 tests pass.
- `apps/web/` — React (Vite) frontend, 8 nav views, 18 tests pass, builds clean.
- `MoTrPAC Hackathon/` — analysis modules (rat comparison, BioCarta transcriptomics, Metabolomics, and
  the generalized `query_core` backend), merged from main.

**Three analysis types exposed through API + UI (all data-backed, no mocks):**
1. directed_signature — mprobe run of the PAH muscle 19-protein signature (view 05 Live results).
2. metabolomics_case_study — read-only ST000763 plasma case study (view 06 Metabolomics).
3. generalized_query — query_core engine, run a bundled signature vs MoTrPAC (view 07 Generalized
   query); PAH examples reproduce the legacy hardcoded results.

**Key API facts the next work depends on:**
- AnalysisRequest (schema.py) ALREADY accepts signature_rows, signature_csv_text, OR example_name,
  plus target_species (rat|human), selected_omics, tissue, sex, timepoint, fdr_threshold. So
  user-chosen signatures/species/omics are a FRONTEND gap, not a backend one.
- ColumnResult ALREADY carries camera_p (raw/unadjusted p), camera_fdr (BH over the 52-col family),
  dataset, species, tissue, layer, sex, timepoint, plus camera_t_up/down and camera_fdr_up/down. So
  "show raw p / species / dataset in the table" is a DISPLAY fix, not a backend gap.
- /api/catalog returns tissues/layers/sexes/timepoints per species (store-derived) + the capability
  matrix — everything a signature/omics picker needs.
- Mapping-confirmation endpoints exist: /api/signatures/validate, /api/mappings/preview (surface
  ambiguous candidates; built-in fixture is unambiguous).
- Endpoints: /api/health, /catalog, /availability, /signatures/validate, /mappings/preview,
  /comparisons (+/{id}, /features, /provenance, /report, /export), /metabolomics/casestudy
  (+/convergence, /export), /generalized/signatures, /generalized/query/{id}. CORS enabled for
  localhost (added late; restart uvicorn to pick it up).

**Frontend gaps (the substance of the requested work):**
- LiveDashboard.tsx HARDCODES example_name="pah_muscle_malenfant2015", target_species="rat",
  selected_omics=["transcriptomics","proteomics"], fdr_threshold=0.05. No user input — this is the
  "bundled queries" to unbundle. (GeneralizedQuery.tsx has a bundled-signature picklist but no
  upload/paste.)
- LiveDashboard set-level table shows: Comparison, opposed/measured, cameraPR t, BH q (family),
  significant, verdict. It does NOT show the raw/nominal p (camera_p) nor a clear species/dataset
  column.
- Views 01-04 + 08 (Architecture, Technical flow, Researcher workflow, Results dashboard [mock], Judge
  slide) are static design/handoff artifacts with mock data (e.g. hardcoded genes PPARGC1A/SOD2/COL1A1
  in the mock Dashboard). These are the "old system-design stuff" to migrate off the frontend to docs.

**Docs in place (extend, don't duplicate):** DECISIONS.md (ADR-0001..0012), VERIFICATION.md,
REQUIREMENTS_TRACEABILITY.md, METABOLOMICS.md, RUN_LOCAL.md.

**Known open items (carry forward):**
- Some tools (fs_write / shell) were intermittently unavailable in prior sessions; several files were
  written via shell heredoc + verified by tests/read-back. Re-read via the file tool if in doubt.
  Note: ~/.kiro/settings/ is a protected path an agent cannot write (by design) — the user edits
  permissions.yaml.
- Gallery under hackathon/tool/site/ records a stale commit (5e2cde3); regen procedure in
  scripts_regenerate_gallery.md (own commit).
- pnpm vs npm not standardized in CI (app ships pnpm-lock.yaml; dev used npm).
- Browser point-to-source enforced at API level only; no Playwright e2e yet.

---

## Feature 1 (priority 1) — Trace the workflows

Goal: a clear, interpretable trace of each end-to-end workflow so the next developer (and judges) can
follow input -> mapping -> engine -> result -> provenance without reverse-engineering code.

Scope:
- One trace per analysis type: (a) directed_signature (mprobe), (b) metabolomics_case_study,
  (c) generalized_query. For each: entry point -> module/function -> data read -> transform -> output
  fields -> provenance/run_id -> which UI view renders it.
- Request->response contract for each API endpoint used, and the point->evidence->source chain
  (already testable: test_point_to_evidence_to_source_chain).
- Include the CLI path for query_core (from MoTrPAC Hackathon/generalized/) alongside the API path.

Deliverable: WORKFLOWS.md (or docs/workflows/) with a mermaid diagram per workflow + a numbered
narrative + exact files/functions touched. Documentation only; no code change. Feeds Features 2 and 3.

Why first: cheapest, de-risks the other two (you'll know exactly what the frontend calls and what the
mock views hide before you touch them), and a direct judging win (Methods & Approach, interpretability).

---

## Feature 2 (priority 2) — Clean up the frontend

Four sub-tasks. All additive to a backend that already supports them.

2a. Unbundle the queries (user chooses the analysis).
- Replace LiveDashboard's hardcoded example_name/species/omics with a real query builder: signature
  source (built-in example | paste rows | upload CSV -> signature_rows/signature_csv_text), target
  species (rat|human), omics, tissue/sex/timepoint (from /api/catalog), FDR threshold.
- Reuse QueryBuilder.tsx domain logic (already models species->study-design and availability); wire it
  to the live API instead of the in-memory catalog.
- Enforce the mandatory mapping-confirmation gate for uploads via /api/mappings/preview before
  /api/comparisons (surface ambiguous candidates; never silently pick).
- FDR threshold = analysis param (changes run_id); "include nonsignificant" = presentation-only (no
  new run). Viz mode from RETURNED layers (already implemented).

2b. Let disease signatures AND protein sets be user-chosen.
- Directed path: any gene/protein signature (API + mapper accept gene_symbol/uniprot/ensembl/
  rat_symbol + direction). Expose built-in examples (hackathon/tool/examples/*.csv) as a picklist plus
  paste/upload.
- Generalized path (GeneralizedQuery.tsx): add upload/paste beyond the bundled picklist (the CLI
  accepts arbitrary signatures matching the query_core schema; the API adapter currently only runs
  bundled ones — likely needs a new endpoint accepting an uploaded signature).

2c. Fix the tables — add raw/nominal p; make enrichment stringency explicit.
- Add a raw p (nominal, camera_p) column next to BH q, so stringency is visible (raw p describes the
  enrichment; BH q is the family-adjusted call). Both already in ColumnResult.
- Make clear BH q is over the frozen 52-col family and raw p is nominal (not adjusted). Consider
  showing up/down-half cameraPR (camera_t_up/down, camera_fdr_up/down) for enrichment detail.

2d. Tables must clearly show species compared and source datasets.
- Add explicit species and dataset columns (both in ColumnResult) to every results table (Live,
  Generalized, Metabolomics), plus the contrast (disease + exercise numerator/denominator are
  first-class in the response). Never render a bare number without species+dataset+contrast.
- For cross-species (human signature -> rat reference), show the ortholog-link relation.

Deliverable: updated apps/web components + client, new tests (component/domain), all builds/tests
green. Backend likely unchanged except possibly a "run uploaded generalized signature" endpoint (2b).
Update RUN_LOCAL.md if the UX changes.

---

## Feature 3 (priority 3) — Remove old system-design views; migrate to docs; document cleanly

Goal: the frontend shows only real, data-backed analyses; the design/mock material becomes docs.

Scope:
- Remove or gate the static design views from the nav: 01 Architecture, 02 Technical flow,
  03 Researcher workflow, 04 Results dashboard (mock), 08 Judge slide. The mock Dashboard (hardcoded
  PPARGC1A/SOD2/COL1A1) must not read as a real result.
- Migrate their content into documentation (reuse the mermaid diagrams from Feature 1). Preserve the
  Figma-Make design provenance (source repo/commit).
- Result: nav is just the live analyses (Live results, Metabolomics, Generalized query) plus maybe a
  single About/Methods page linking to docs.
- Document cleanly + interpretably: a top-level README that actually describes the project (verify the
  root README isn't still a stub), links to WORKFLOWS/DECISIONS/VERIFICATION/METABOLOMICS/RUN_LOCAL,
  and states the honest framing (cross-cohort association, not a PAH treatment effect; separate
  cohorts kept separate).

Deliverable: trimmed apps/web nav (design views removed/relocated), new/updated docs, root README
rewritten, all tests green. Largely deletion + documentation; low risk once Features 1-2 are done.

Order note: do Feature 3 AFTER 1 and 2 — the workflow traces (1) supply the doc content that replaces
the removed views (3), and don't delete the design views until the real query UI (2) fully covers what
a user needs.

---

## Non-negotiables to preserve across all three features

- mprobe / query_core / the standalone modules stay the source of truth. No statistics in the frontend
  or the FastAPI layer.
- Species-only study selection (rat->chronic, human->acute), derived in the backend. No acute/chronic
  input control.
- Both contrasts first-class; never "healthy gene set." Show species + dataset + contrast on every
  result.
- Frozen 52-column BH multiplicity family for the directed path; the male rat SKM-GN protein 8-wk
  headline is q=0.0584, NOT significant — never surface the 16-column q=0.0413 as a headline (only as
  a labeled sensitivity analysis).
- Honest framing + guardrails + PTM caveat rendered in the UI; cohorts kept separate (muscle protein
  vs plasma metabolomics vs blood).
- Traceability: each feature gets ADR entries in DECISIONS.md, a VERIFICATION.md checkpoint (commands,
  counts, hashes), and REQUIREMENTS_TRACEABILITY.md rows.
- Verify before claiming: engine (pytest -q -rs, expect 122/34), api (pytest apps/api/tests), web
  (pnpm test + pnpm run build); confirm python analysis.py stays byte-identical to the golden baseline
  whenever anything near the engine changes.

## Suggested first steps for the next session
1. Trace the three analysis types into WORKFLOWS.md (Feature 1).
2. Prototype the unbundled query builder in LiveDashboard/QueryBuilder against the live API (2a/2b),
   then the table columns (2c/2d).
3. Only then remove the design views and finish the docs (Feature 3).
