# VERIFICATION — clean-clone / CI reproduction log

Purpose: prove the integrated project reproduces from a clean state. Records exact
commands, versions, checksums, and passed/skipped tests. Append per gate.

## Environment baseline
- OS: macOS (darwin, arm64)
- git: 2.47.0
- git-lfs: 3.7.0 (present but NOT required — see below)

## Repo state at start of integration
- Integration branch base: `origin/t3code/build-motrpac-probe-tool` @ `9962a43`
- Merged: `origin/main` @ `57e01e7` (merge commit `09ae0d7`)
- Merge base of tool branch and main: `993f746`
- React app source repo: `github.com/StanchPillow55/ExerciseSignatureExplorerHandoff` @ `d9fd564`

## LFS finding (Gate 1)
- The React repo's `.gitattributes` declares LFS rules, but `git lfs ls-files` is EMPTY.
- The only tracked binary, `src/imports/Screenshot_2026-09-26_at_6.48.13_PM.png`, is a
  real 152,756-byte PNG git blob (`git cat-file -p HEAD:...` begins with PNG magic bytes,
  not a `version https://git-lfs...` pointer).
- Conclusion: there is NO LFS object to materialize; a subtree import carries the real
  blob. The parent repo has no LFS dependency.

---

## Gate checkpoints
(Appended as each gate completes.)

### Gate 1 checkpoint — Repository & integration (PASS)
Files changed / added:
- `DECISIONS.md`, `VERIFICATION.md`, `REQUIREMENTS_TRACEABILITY.md` (new; traceability trail)
- `apps/web/**` (React app imported via `git subtree add --prefix=apps/web`, no `--squash`)
- `apps/web/.gitattributes` (LFS rules replaced with plain `binary` attributes)

Why each change was necessary:
- Branch off tool branch + merge main → single integration line with engine + analysis
  modules (ADR-0002).
- Subtree import (not submodule, not rm-rf) → self-contained repo, preserved history
  (ADR-0003).
- `.gitattributes` neutralized → fresh clone works without Git LFS.

Commands run (key):
- `git switch -c integration/exercise-signature-explorer origin/t3code/build-motrpac-probe-tool`
- `git merge --no-ff origin/main`  → merge commit `09ae0d7`
- `git subtree add --prefix=apps/web ese-src main`  (source `d9fd564`)  → commit `7863879`
- `GIT_LFS_SKIP_SMUDGE=1 git clone --branch integration/... <repo>.git /tmp/ese-freshclone`

Results:
- Integration branch HEAD after Gate 1: `728f36f`
- `hackathon/tool/` (mprobe) and `MoTrPAC Hackathon/Rat Comparison PAH Proteins/` coexist.
- React app at `apps/web`, tracked as real files (not a gitlink).
- PNG `apps/web/src/imports/Screenshot_2026-09-26_at_6.48.13_PM.png`: real PNG blob
  (152,756 bytes), NOT an LFS pointer, verified in a fresh clone with LFS smudge skipped.
- In-repo build: `npm run build` OK (JS 269.10 kB / gzip 81.42 kB), `npm test` 16/16 pass.
- Fresh clone (no original repo, LFS skipped): install OK, build OK, tests 16/16 pass.

Skips / assumptions / limitations:
- The React app ships both `pnpm-lock.yaml` (original) and now installs cleanly with npm;
  the stray `package-lock.json` from an npm install was removed to avoid dual lockfiles.
  Toolchain (`.mise.toml`) prefers pnpm 10.34.3; CI should standardize on one. (Open item.)
- Source repo URL: `github.com/StanchPillow55/ExerciseSignatureExplorerHandoff` @ `d9fd564`.

### Gate 2 checkpoint — Fixtures & mapping (PASS)
Files added / changed:
- `hackathon/tool/examples/pah_muscle_malenfant2015.provenance.csv` (new; lossless superset
  schema: paper vs canonical ids, ratio, disease_log2fc, reported p-value text, correction notes)
- `hackathon/tool/examples/pah_muscle_malenfant2015.mapping_manifest.{csv,json}` (new; committed
  pre-approved mapping for the built-in fixture)
- `hackathon/tool/examples/pah_muscle_sensitivity_views.json` (new; full19/lower9/upper10/fibre)
- `hackathon/tool/src/motrpac_probe/signature.py` (added `candidates_for()`, `mapping_audit()`,
  `_all_maps()`; existing `map_genes`/`load` UNCHANGED so golden outputs are preserved)
- `hackathon/tool/tests/test_pah_fixture_table2.py` (new; 8 tests)
- `hackathon/tool/tests/test_mapping_candidates.py` (new; 4 tests)
- `hackathon/tool/.gitignore` (added `.venv/`)
- NOT changed: `pah_muscle_malenfant2015.csv`, `pah_muscle_lower9_malenfant2015.csv` (ADR-0005).

Why each change was necessary:
- Enrich via companion table, not by editing the fixture → preserves the live input + golden
  outputs while adding full paper↔canonical provenance (ADR-0005).
- `candidates_for`/`mapping_audit` → surface all candidates (176 Ensembl ids map to >1 human
  symbol) for the mandatory confirmation gate; never `distinct()` away duplication; conserve rows.
- Sensitivity views → planned subgroup analyses (incl. fibre/myosin competing explanation), not
  post-hoc gene removal.

Environment / setup:
- Created `hackathon/tool/.venv` (gitignored), `pip install -e ".[test]"`, `mprobe store fetch`
  (store bundle `mprobe-store-v1`, 10 files verified, store hash `d615c2f0b9aec182`).

Commands run (key):
- `python analysis.py` before edits → `/tmp/golden_full19_baseline.txt`
- `python -m pytest -q -rs` before edits → 110 passed, 34 skipped (allowlist baseline)
- `python -m pytest tests/test_pah_fixture_table2.py tests/test_mapping_candidates.py -q` → 12 passed
- `python analysis.py` after edits → byte-identical to baseline (diff empty) ⇒ NO scientific drift
- `python -m pytest -q -rs` after edits → 122 passed, 34 skipped

Golden numbers reproduced (full-19, unchanged):
- 19 genes counted, 25/25 rows mapped.
- Rat gastrocnemius male protein 8 wk: 13/19 opposed, cameraPR t +2.23 (matches QA).
- Rat gastrocnemius female protein 8 wk: 15/19 opposed, cameraPR t +4.62.

Expected-skip allowlist (34): test_golden.py (30), test_camera.py (3: 2 deck + 1 R env),
test_metab.py (1), test_store.py (1) — all "needs local hackathon/ data & deck numbers (not in
git)" or "R env not available". CI must use `pytest -q -rs` and allowlist exactly these.

Fixture provenance / integrity notes:
- Fixture `pah_muscle_malenfant2015.csv` sha256_16 = `0c80d82e5546ed9b`.
- All 25 rows map unambiguously to a single measured MoTrPAC symbol (n_ambiguous=0) under the
  released ID map; new uploads still require interactive confirmation.
- MYH1 correction recorded: paper prints Q9UKX2 (which is MYH2); canonical MYH1 = P12882. Both
  preserved. Isoform suffixes dropped for GLO1 (Q04760-2→Q04760) and ATP2A1 (O14983-2→O14983).

### Gate 3 checkpoint — Scientific service run_analysis() (PASS)
Files added:
- `apps/api/motrpac_probe_service/{__init__,schema,service}.py` (the service adapter)
- `apps/api/pyproject.toml`, `apps/api/README.md`, `apps/api/tests/test_service.py` (8 tests)

Why:
- One UI-independent entry point `run_analysis(request) -> AnalysisResponse` wrapping
  `motrpac_probe.run.compute()`; performs no statistics; returns only JSON-safe values.
- Content-addressed input (temp file named by sha256 of the signature bytes); the service never
  accepts a client-supplied server path. This also fixes the confirmed KeyError('_path') by setting
  `R["opts"]["_path"]` before calling `run.provenance()`.
- Deterministic `run_id = sha256(signature bytes + analysis params + ordered family column IDs +
  mapping version + store hash + code version + seed + schema version)[:24]`. Presentation-only
  options (include_nonsignificant, display tissue/sex/timepoint) are excluded; FDR threshold is
  included.
- Multiplicity family reported from the engine's core columns; BH already applied across the whole
  family inside compute(); the service never re-runs BH. Provenance records method, family_size,
  n_tests, threshold, and all ordered column IDs.

Verified numbers:
- Family size = 52 columns (engine `core_cols`); n_tests = non-missing cameraPR tests.
- Headline male rat SKM-GN protein 8 wk: camera_t = +2.2295, camera_p = 0.0258,
  camera_fdr (52-col BH) = 0.0584 -> NOT significant at q<0.05; 13/19 opposed. COMPUTED by the
  service and asserted to |.-0.0584|<5e-4 (not hardcoded).
- run_id identical across repeated identical requests and across presentation-only changes;
  differs when fdr_threshold changes 0.05 -> 0.10.
- Whole response serializes with `json.dumps(..., allow_nan=False)`.

Commands:
- `pip install -e apps/api[test]` (editable, alongside the editable engine)
- `python -m pytest apps/api/tests/test_service.py -v` -> 8 passed

### Gate 4 checkpoint — Catalog, lineage, availability from store.columns (PASS)
Files added / changed:
- `apps/api/motrpac_probe_service/catalog.py` (new; catalog + availability + capability matrix)
- `apps/api/motrpac_probe_service/__init__.py`, `schema.py`, `service.py` (lineage fields wired)
- `apps/api/tests/test_catalog.py` (new; 8 tests), `tests/test_service.py` (+2 lineage tests)

Why:
- Catalog and availability are derived from `store.load_columns()`, never hardcoded constants, so
  the store is the single source of truth for what exists.
- Fixes the React catalog bug: rat SKM-VL is RNA-only; SKM-GN has RNA+PROT+PHOSPHO.
- Capability matrix per layer (engine/api/react/demo_validated): transcriptomics & proteomics are
  fully exposed; phosphoproteomics is engine+api but not react; metabolomics is engine-only (not
  api/react); genomics/epigenomics unsupported. So an engine-supported/UI-not-exposed layer is not
  mislabeled globally "planned".
- Source-species gating: mapper supports human + rat only; mouse/other are explicitly unsupported.
- Lineage: every measured FeatureEvidence carries source_feature_id, n_collapsed,
  aggregation_method, mapping_decision_id -> point -> evidence -> source is reconstructable/testable.

Verified numbers (from the store, 422 columns total):
- rat SKM-GN exercise layers = {RNA, PROT, PHOSPHO}; SKM-VL = {RNA} only.
- Collapse rule = max |stat|: in rat SKM-GN protein, MYH7 collapses 2 features (feature_id
  NP_058936.1), PDLIM3 collapses 3 (XP_006253181.1), NDUFA9 = 1 (NP_001094222.1).
- Availability: (rat, [transcriptomics, proteomics], SKM-VL) -> partially_available (proteomics
  unavailable); SKM-GN -> available; non-existent tissue -> no_matching_context (not an error).

Commands:
- `python -m pytest apps/api/tests -v` -> 18 passed (8 service + 2 lineage + 8 catalog)

### Gate 5 checkpoint — End-to-end (PASS)
Files added / changed:
- `apps/api/motrpac_probe_service/app.py` (FastAPI; OpenAPI = transport truth), `export.py` (bundle
  + HTML report), `pyproject.toml` (test extra += fastapi, httpx)
- `apps/api/tests/test_api.py` (10 incl. point->source), `test_parity.py` (3)
- `apps/web/src/api/client.ts` (typed client), `src/components/LiveDashboard.tsx` (live view),
  `src/App.tsx` (+ "Live results (API)" view), `src/domain/analysis.ts`
  (getVisualizationModeFromReturnedLayers), `src/domain/analysis.test.ts` (+2 tests)
- `.github/workflows/ci.yml` (python + web jobs, skip allowlist)
- `hackathon/tool/scripts_regenerate_gallery.md` (gallery regen procedure)
- `hackathon/tool/src/motrpac_probe/signature.py` (mapping_audit now enforces the id-column
  requirement so validate/preview return 422 on malformed CSVs)

Why:
- OpenAPI/Pydantic app wraps the service; endpoints /health, /catalog, /availability,
  /signatures/validate, /mappings/preview, /comparisons(+/{id}/features,/provenance,/report,/export).
- Scientific empty states (no_matching_context, no_compatible_data, zero-mapped) return HTTP 200;
  malformed input returns 422.
- React: live client + LiveDashboard show real numbers (summary, set-level table with family
  q-values, headline, both contrasts, guardrails, PTM caveat, export/report links). Static design
  views retained. Visualization mode derives from RETURNED layers (new domain fn + tests) so a
  requested-but-unreturned layer cannot fabricate a matrix.

Verified end-to-end (live uvicorn on 127.0.0.1:8137):
- GET /api/health -> ok, store_hash d615c2f0b9aec182.
- GET /api/catalog -> rat SKM-VL layers ["RNA"] (store-derived).
- POST /api/comparisons {example pah_muscle_malenfant2015, rat} -> status ok, 19 counted genes,
  returned_layers [PROT, RNA], headline male rat SKM-GN protein 8wk q=0.0584 significant=False.
- GET .../export -> zip with exactly the 9 bundle files.
- Parity: service set-level results match a direct run.compute() cell-for-cell (n_opposed,
  n_measured, camera_t, camera_fdr to 6 dp) at q<0.05 and q<0.10; store+counted genes match.
- point -> evidence -> source: every measured feature's source_feature_id exists in its column's
  store rows (test_api.py::test_point_to_evidence_to_source_chain).

Test totals:
- apps/api: 31 passed (8 service + 2 lineage + 8 catalog + 10 api + 3 parity).
- apps/web: 18 passed (16 original + 2 viz-mode-from-returned-layers); production build OK
  (JS 276.58 kB / gzip 83.29 kB).
- engine suite unchanged at 122 passed / 34 skipped (allowlist).

Gallery: committed reports under hackathon/tool/site record stale repo_git_sha 5e2cde3 (predates
this branch). Regeneration on this branch records the current HEAD (verified a1dbfc4). Left the
committed teammate-authored reports in place; regen procedure documented in
scripts_regenerate_gallery.md (should be its own reviewed commit).

Headline policy: male rat SKM-GN protein 8wk reported as q=0.0584 (NOT significant at q<0.05,
directionally concordant only) under the 52-column family. The 16-column q=0.0413 is not used as a
headline anywhere; it would only appear as an explicitly labeled sensitivity analysis.

### Post-Gate-5 merge checkpoint — origin/main (metabolomics pipeline) (PASS, no conflicts)
- Fetched origin/main which advanced to `2140f9b` "Add PAH metabolomics pipeline mapped to MoTrPAC".
- Merged into the integration branch: merge commit `fb161d9` (parents f756b6c + 2140f9b). Working
  tree clean; NO merge conflicts. Branch is 22 ahead / 0 behind origin/main (fully contains 2140f9b).
- What 2140f9b added: a self-contained module under `MoTrPAC Hackathon/Metabolomics/` (ST000763 PAH
  metabolomics: scripts 00–07, data, outputs, README). Additive analysis module only.
- Impact assessment: 2140f9b does NOT touch `hackathon/tool/` — the mprobe engine, examples, tests,
  and store are untouched. Therefore this is the additive-module case, not a scientific-output STOP
  condition; no app change required and no golden re-verification forced.
- Post-merge verification (run anyway):
  - `python analysis.py` output byte-identical to /tmp/golden_full19_baseline.txt (no scientific drift).
  - engine: 122 passed / 34 skipped (allowlist unchanged).
  - api: 31 passed. web: 18 passed + production build OK.
- Note: the new Metabolomics module is the deferred "Track 2 blood/pathway"-adjacent work maturing on
  main; it is NOT yet wired into the app (metabolomics remains engine-supported / UI-not-exposed per
  the capability matrix, ADR-0009). Wiring it into the API/UI is a future, separate task.

### Metabolomics integration checkpoint — case-study adapter + convergence (PASS)
Files added / changed:
- `apps/api/motrpac_probe_service/metab_casestudy.py` (Option-1 adapter), `convergence.py`
  (cross-check), `app.py` (+3 endpoints), `catalog.py` (+ANALYSIS_TYPES), `__init__.py` (exports)
- `apps/api/tests/test_metabolomics.py` (7 tests)
- `apps/web/src/api/client.ts` (metab types+methods), `src/components/MetabolomicsCaseStudy.tsx`
  (separate top-level view), `src/App.tsx` (+ "Metabolomics (ST000763)" view)
- `METABOLOMICS.md` (module<->engine relationship, integration, cross-check, caveat)

Why: surface the merged standalone metabolomics module (2140f9b) as a SEPARATE read-only case-study
analysis type, without blending it with the muscle/rat PAH work and without recomputing it in the
service (module stays source of truth; ADR-0011).

Verified:
- Adapter payload matches `06_context_summary.json`: 41 hits, 27 matched blood, hit_labels
  {stable 24, no-match 14, drifts 3}, exercise-sensitive hits 0% (Fisher p 0.609), median drift 3.09.
  JSON-safe (json.dumps allow_nan=False).
- Convergence cross-check: module 04 MoTrPAC export vs engine METAB store, EE-CON, on the unambiguous
  RefMet-name subset -> IDENTICAL, 1156/1156 cells, max abs diff 0.0 (both 2.0.8 / c2.0). Ambiguous
  blank/duplicate-RefMet cells excluded and counted (key ambiguity, not disagreement).
- API endpoints (via FastAPI TestClient): /api/metabolomics/casestudy 200, /convergence identical,
  /casestudy/export zip {summary,hits,conclusion,provenance}; catalog analysis_types includes
  metabolomics_case_study (api+react true, "SEPARATE cohort").
- Test totals after this pass: api 38 passed; web 18 passed + build OK; engine 122/34 UNCHANGED;
  `python analysis.py` byte-identical to the pre-existing golden baseline (adapter-only).
- UI: separate top-level "Metabolomics (ST000763)" view; honest conclusion + both contrasts + EE-EE
  caveat + separate-cohort banner + convergence badge + export all rendered. Engine live metabolite
  scoring (Option 2) remains a documented follow-on.

### Merge + generalized-pipeline replication checkpoint (PASS)
- Merged origin/main 0882687 "Add generalized Track 2 backend" (merge e61c9a8). Additive module
  `MoTrPAC Hackathon/generalized/` (query_core engine, discordance, pathways, plots, signatures);
  does NOT touch hackathon/tool/ or apps/. Clean merge, README auto-merged.
- App suites after merge: engine 122 passed / 34 skipped; api 38 passed; web 18 + build OK;
  `python analysis.py` byte-identical to golden baseline.
- generalized query_core tests: 8 passed.
- REPLICATION (hardcoded -> generalized): ran the PAH blood signature through query_core against the
  bundled MoTrPAC reference:
    python -m query_core.cli --disease examples/signatures/pah_blood_rna.csv.gz \
      --reference query_core/motrpac_reference.csv.gz --out <out> --tissue blood \
      --reference-contrast-category EE-CON
  RNA-layer Spearman rho (13,055 shared genes) reproduces the legacy hardcoded pipeline exactly:
    during_20_min -0.0788, during_40_min -0.0286, post_10_min +0.0601,
    post_15_30_45_min +0.1247, post_3.5_4_hr +0.1161, post_24_hr +0.0899.
  Matches the legacy values (-0.079, -0.029, +0.060, +0.125, +0.116, +0.090) to the decimal:
  same input + same MoTrPAC reference + generalized engine -> identical result. This is the
  "generalize the hardcoded pipeline" confirmation.


### Generalized query wiring + RUN_LOCAL checkpoint (PASS)
Files added/changed:
- `apps/api/motrpac_probe_service/generalized_query.py` (adapter over query_core.engine.run_query),
  `app.py` (+2 endpoints), `catalog.py` (+generalized_query analysis type), `__init__.py` (export)
- `apps/api/tests/test_generalized.py` (5 tests)
- `apps/web/src/api/client.ts` (types+methods), `src/components/GeneralizedQuery.tsx` (view),
  `src/App.tsx` (+ "Generalized query" view, now 8 views)
- `RUN_LOCAL.md` (two-process run instructions + port guidance)

Why: expose the generalized query_core engine (CLI-only in the repo) through the API/UI as a third
analysis type, and document local run (a port clash with an unrelated ":8000 aegis" app was causing
confusion).

Verified:
- Adapter reproduces the legacy blood result THROUGH the service: rna-layer Spearman
  during_20_min -0.0788, during_40_min -0.0286, post_10_min +0.0601, post_15_30_45_min +0.1247,
  post_3.5_4_hr +0.1161, post_24_hr +0.0899 (matches the hardcoded pipeline). JSON-safe.
- Endpoints (TestClient): /api/generalized/signatures -> 5; /api/generalized/query/pah_blood_rna ->
  200 with the reproduced rho; unknown id -> 404; catalog analysis_types includes generalized_query.
- Test totals: api 43 passed (was 38 + 5 generalized); web 18 + build OK (8 views); engine
  122/34 UNCHANGED; golden byte-identical (adapter-only).

Local-run diagnosis: `:8000` is an unrelated "aegis" app; this app's web UI is Vite on 8443/5173 and
the API is uvicorn on a chosen port (docs use 8765). A stale Vite was found on :5173. RUN_LOCAL.md
records the exact two-terminal commands and VITE_API_BASE_URL wiring.


### CORS fix checkpoint (PASS)
- Symptom: browser showed "TypeError: Failed to fetch"; API logged OPTIONS /api/comparisons 405.
  Root cause: FastAPI app had no CORS middleware; the web origin (8443) differs from the API origin
  (8765), so preflight OPTIONS were rejected and every browser fetch failed. Server-side/TestClient
  calls worked, which is why tests were green but the UI was not.
- Fix: added CORSMiddleware to app.py, allowing the default Vite ports plus any localhost/127.0.0.1
  port (regex), overridable via MPROBE_CORS_ORIGINS.
- Verified (TestClient): preflight OPTIONS /api/comparisons -> 200 with access-control-allow-origin
  http://localhost:8443 (was 405); GET carries the header. Added regression test
  test_cors_preflight_allows_web_origin. api suite 44 passed.

### R1 checkpoint — WORKFLOWS.md authored (documentation-first) (PASS)
Files added / changed:
- `WORKFLOWS.md` (new; root) — one mermaid + numbered seven-stage narrative per five analysis types,
  API request→response contract table, point→evidence→source ordered sequence, CLI-only labels.
- `DECISIONS.md` (ADR-0013), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R1 rows) — traceability trail (R13).
- NO runtime source file added, removed, or altered (R1 AC7).

Why each change was necessary:
- Documentation-first (ADR-0013): trace every analysis path before building features so success
  criteria are analyzed against the real code, and so the mermaid diagrams exist for reuse in R6.

How the workflows were traced (read-only investigation):
- API endpoints + directed path: `apps/api/motrpac_probe_service/app.py`, `service.py`.
- Read-only adapters: `metab_casestudy.py`, `generalized_query.py`, `catalog.py`.
- Discordance/PTM CLI + demo outputs: `MoTrPAC Hackathon/generalized/discordance/`
  (`build_discordance.py`, `audit_ptm_parent.py`, `README.md`, `demo_muscle_ee/*`,
  `demo_ptm_parent_ee/*`).
- UI views: `apps/web/src/App.tsx`, `components/*.tsx`, `api/client.ts`.
- point→evidence→source test: `apps/api/tests/test_api.py::test_point_to_evidence_to_source_chain`.

Structure grep check (commands run, in repo root):
- ``grep -c '```mermaid' WORKFLOWS.md`` → 5 (one per analysis type; R1 AC2).
- `grep -E '^## [0-9]\.' WORKFLOWS.md` → the five analysis-type sections
  (directed_signature, metabolomics_case_study, generalized_query, discordance catalog + predictive
  model, PTM-parent audit; R1 AC2).
- Seven-stage labels (`Entry point.`/`Module / function.`/`Data read.`/`Transform.`/`Output fields.`/
  `Provenance / run_id.`/`UI view.`) present across all sections (35 = 7×5, +1 template mention;
  R1 AC3), with "not applicable" used for missing stages (`grep -c 'not applicable'` → 4; R1 AC3).
- `grep -c 'CLI-only' WORKFLOWS.md` → 8 (discordance catalog, predictive model, PTM-parent labeled
  CLI-only; R1 AC6).
- `grep -c 'test_point_to_evidence_to_source_chain' WORKFLOWS.md` → 1 (R1 AC5).
- `grep -c 'PLANNED' WORKFLOWS.md` → 12; the contract table covers existing endpoints plus the three
  planned discordance endpoints (`/api/discordance/{catalog,model,ptm-parent}`) and the planned
  generalized upload endpoint (`POST /api/generalized/query`) (R1 AC4).

Path-coupling confirmation (read-only):
- `python3 -c "from pathlib import Path; print(Path('apps/api/motrpac_probe_service/metab_casestudy.py').resolve().parents[3])"`
  → repo root `.../multiomics-hackathon-2026-track-2` (confirms the `parents[3]` constraint the
  read-only adapters rely on; documented in WORKFLOWS.md).

Result: R1 deliverable is `WORKFLOWS.md` only; no runtime source changed. Structure checks PASS.

### R2 checkpoint — discordance deliverables via read-only API (PASS, with one deferred 2.4 test)
Files added / changed:
- `apps/api/motrpac_probe_service/discordance_casestudy.py` (new; Option-1 read-only adapter:
  `build_catalog_summary()`, `build_model()`, `build_ptm_parent()`; `catalog_available()` /
  `model_available()` / `ptm_available()` guards; event-class domain gate over `KNOWN_EVENT_CLASSES`).
- `apps/api/motrpac_probe_service/app.py` (+3 GET-only endpoints: `/api/discordance/catalog|model|ptm-parent`,
  each 404 when unavailable, `FileNotFoundError` → 404).
- `apps/api/motrpac_probe_service/catalog.py` (additive `ANALYSIS_TYPES["discordance_catalog"]` entry).
- `apps/api/tests/test_discordance.py` (Properties 1, 3, 4 + parity/endpoint unit tests; the
  task-2.4 Property 2 test is present but 2.4 is still in progress).
- `DECISIONS.md` (ADR-0014), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R2 rows) — traceability trail (R13).

Why each change was necessary:
- Surface the three CLI-only Track 2 deliverables through the proven Option-1 adapter pattern
  (ADR-0014) without recomputing statistics in the transport layer (R11.1 source-of-truth boundary).
- GET-only endpoints + `available()` guards enforce the read-only, no-write contract (R2 AC1) and
  return 404 (never fabricated/partial data) when a committed output is missing/invalid (R2 AC3, AC8).

Commands run (venv: `hackathon/tool/.venv`, from repo root) — verify-before-claiming (R13.4):
- `hackathon/tool/.venv/bin/python -m pytest apps/api/tests/test_discordance.py -q`
  → `1 failed, 14 passed` (the 1 failure is `test_property2_read_only_idempotence_and_byte_stability`).
- `... -m pytest apps/api/tests/test_discordance.py --collect-only -q` → 15 tests collected.
- `... "apps/api/tests/test_discordance.py::test_property2_read_only_idempotence_and_byte_stability"`
  run IN ISOLATION → `1 passed` (Property 2 passes on its own; fails only under full-file interaction).

Committed counts asserted (read straight from the committed source files, then confirmed by the tests):
- `run_summary.json.catalog_classes` (`demo_muscle_ee/run_summary.json`):
  `supported_concordant`=1, `supported_opposite`=0 (zero-filled by the adapter; absent from the file),
  `rna_response_protein_equivalent`=285, `indeterminate`=5642, `no_protein_measurement`=9228,
  `no_rna_measurement`=255.
- Adapter `classification_counts` = {supported_concordant 1, supported_opposite 0,
  rna_response_protein_equivalent 285, indeterminate 5642}; `coverage_counts` =
  {no_protein_measurement 9228, no_rna_measurement 255} (kept separate; 4 classes + 2 coverage states).
- `model_metrics.csv` = exactly 6 rows = {all_mapped, rna_responsive} × {zero, rna_only, temporal}.
- `demo_ptm_parent_ee/ptm_parent_summary.csv` = 3 timepoint rows;
  `ptm_parent_candidates.csv` = 749 rows → `candidate_count`=749.
- `Path(discordance_casestudy.py).resolve().parents[3]` == repo root (path-coupling asserted by
  `test_parents3_path_coupling`).

Passing R2 tests (14): `test_parents3_path_coupling`, `test_committed_catalog_class_counts_exact`,
`test_run_summary_catalog_classes_exact`, `test_build_catalog_summary_classification_and_coverage_counts`,
`test_build_model_has_exactly_six_rows_over_the_fixed_grid`, `test_build_ptm_parent_summary_and_candidate_count`,
`test_property1_read_fidelity_and_json_safety`, `test_property3_count_parity_and_fixed_model_grid`,
`test_property4_event_class_domain_gate`, `test_endpoint_discordance_catalog`,
`test_endpoint_discordance_model`, `test_endpoint_discordance_ptm_parent`,
`test_endpoints_are_read_only`, `test_endpoint_catalog_unavailable_returns_404`.

Deferred / open item:
- `test_property2_read_only_idempotence_and_byte_stability` belongs to task 2.4 (Property 2), which is
  still in progress (`[-]`) and is NOT part of task 2.6. It PASSES in isolation (`1 passed`) and fails
  only under full-file test interaction — a test-isolation matter in the not-yet-finalized 2.4
  (Property 4 monkeypatches the adapter's module-level `DEMO_DIR`/`CATALOG_CSV`/`RUN_SUMMARY_JSON`
  constants, and the module-scoped TestClient fixture couples across Hypothesis examples). It is NOT
  a defect in the R2 adapter or the three endpoints, all of which the 14 functional/property tests
  above confirm. Finalizing 2.4's isolation is tracked with task 2.4.

Result: R2 adapter + endpoints + catalog entry verified against the committed demo outputs; all R2
functional and property assertions PASS. The single non-passing test is the in-progress task-2.4
Property 2 (deferred, not a task-2.6 concern).

### R3 checkpoint — Discordance three-stage main view + layer_discordance in Live (PASS)
Files added / changed:
- `apps/web/src/api/client.ts` (task 4.1) — `discordanceCatalog()`, `discordanceModel()`,
  `discordancePtmParent()` methods + `DiscordanceCatalogResponse`, `DiscordanceModelResponse`,
  `DiscordanceModelMetricRow`, `PtmParentResponse` interfaces (mirror the R2 adapter payload shapes).
- `apps/web/src/components/Discordance.tsx` (task 4.2/4.3) — the three-stage view: Stage 1 query
  builder over the committed demo (UI text states live recompute is out of scope); Stage 2
  visualization (class-count chart, model-metrics comparison chart, PTM-parent timing, optional
  RNA/protein scatter ONLY if per-feature rows are delivered); Stage 3 interpretation (all four
  classification classes incl. zero, two coverage states kept separate, R²/MAE/RMSE model-metrics
  table for zero/rna_only/temporal, PTM-parent divergence summary, PTM occupancy caveat, visible
  weak-prediction statement). Single exclusive empty/error branch
  (`data-testid="discordance-empty-or-error"`) rendering NONE of counts/metrics/PTM on failure or an
  empty catalog (R3 AC9).
- `apps/web/src/components/LiveDashboard.tsx` (task 4.4) — renders the `layer_discordance` rows
  already returned by `run_analysis` as a "Layer discordance" table (R3 AC8).
- `apps/web/src/components/Discordance.property.test.tsx` (task 4.5) — Property 14.
- `apps/web/src/components/Discordance.test.tsx` (task 4.6, authored in parallel; present at
  verification time) — three-stage order, PTM summary + occupancy caveat, weak-prediction text.
- `DECISIONS.md` (ADR-0015), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R3 rows) — traceability trail (R13).

Why each change was necessary:
- Surface the R2-adapter discordance deliverables in the UI as a guided query → visualize → interpret
  flow (R3 AC2), keeping the source-of-truth boundary (R11.1: the client computes no statistics; it
  renders the read-only adapter payloads verbatim). See ADR-0015.
- Render all four classification classes zero-filled (R3 AC3) and the two coverage states separately;
  render R²/MAE/RMSE per model variant (R3 AC4); render the PTM-parent divergence summary (R3 AC6),
  the PTM occupancy caveat as visible text (R3 AC7), and a visible weak-prediction statement alongside
  the metrics (R3 AC5, honest framing R12.4).
- Render `layer_discordance` in the Live view (R3 AC8) — previously computed but shown in zero views.
- A single exclusive empty/error return so failure and delivered-but-empty catalogs render only the
  empty-state and none of the results (R3 AC9).

Commands run (from `apps/web`) — verify-before-claiming (R13.4):
- `npm run test:domain` (`node --experimental-strip-types --test src/domain/analysis.test.ts`)
  → `# tests 18  # pass 18  # fail 0` (18 passed, 0 failed).
- `npm run test:web` (`vitest run`) → `Test Files  3 passed (3)  Tests  7 passed (7)`:
  - `src/components/Discordance.test.tsx` — 4 passed (task 4.6 component tests).
  - `src/components/Discordance.property.test.tsx` — 2 passed (Property 14; each ≥100 fast-check
    iterations: success branch = all four class counts incl. zero + each present variant's
    R²/MAE/RMSE render, empty/error absent; failure/empty branch = empty/error block present and
    NONE of counts/metrics/PTM render).
  - `src/components/ResultsTable.property.test.tsx` — 1 passed (Property 8, unrelated to R3).

Property 14 (`# Feature: finalize-discordance-mvp, Property 14`) — result completeness +
empty-state exclusivity — PASS (2 tests, ≥100 iterations each). fast-check drives the two mutually
exclusive branches; the API client module is mocked so no network is touched.

Component tests (`Discordance.test.tsx`, task 4.6) — PASS (4 tests): three-stage marker order is
`STAGE 1 · QUERY` → `STAGE 2 · VISUALIZE` → `STAGE 3 · INTERPRET` (asserted by both text order and
`compareDocumentPosition`) (R3 AC2); PTM-parent divergence summary present with candidate-gene count
and timepoint row (R3 AC6); PTM occupancy caveat present as visible text under the "Occupancy caveat"
heading (R3 AC7); weak-prediction note present alongside the metrics inside the
"Model metrics (R² / MAE / RMSE)" card, labeled by a visible "WEAK PREDICTION" marker (R3 AC5).

Scope note (honest): R3 AC1 (Discordance is the DEFAULT view on initial load) is NOT delivered by
this task — the default-view/nav wiring is task 8.1 (R6). The `Discordance.tsx` component is ready;
its default-nav wiring is tracked with 8.1. All other R3 ACs (2–9) are delivered and verified above.

Result: R3 Discordance three-stage view + `layer_discordance` Live rendering verified against the
mocked R2-adapter payload shapes; all present R3 tests PASS (domain 18/18, web 7/7 across 3 files).
Engine and `apps/api` behavior unchanged (web-only changes).

### R4 checkpoint — Live query builder rewired to the live API + mandatory mapping-preview gate (PASS)
Files added / changed:
- `apps/web/src/components/QueryBuilder.tsx` (task 5.1) — rewired to `/api/catalog` + live comparison,
  reusing `getStudyContext`, `reconcileQueryForSpecies`, `getVisualizationModeFromReturnedLayers`,
  `layerCapabilities` from `domain/analysis.ts`. Signature source is EXACTLY ONE of built-in example /
  pasted text / uploaded CSV (example/paste → `signature_rows`|`example_name`; upload →
  `signature_csv_text`); invalid upload (empty / over-max-size / unparseable) is rejected, the prior
  source is retained unchanged, and an error is shown. Catalog-driven selectors: total `/api/catalog`
  failure disables ALL selectors + catalog-unavailable error; partial degradation disables ONLY the
  empty dimension's selector. FDR presets over `0 < t <= 1` with the UI-inclusive-upper vs API
  `gt=0`/`lt=1` boundary note. Mandatory mapping-preview gate blocks `/api/comparisons` until confirmed
  (preview failure blocks + retains selections + errors; every multi-candidate term surfaced with all
  candidates, none auto-selected). Species-derived study shown read-only (no acute/chronic input); both
  contrasts first-class; no "healthy gene set". Also fixed the pre-existing tsc error
  (`replaceAll` → `split`/`join`).
- `apps/web/src/components/LiveDashboard.tsx` (task 5.2) — removed the hardcoded request
  (`example_name="pah_muscle_malenfant2015"`, `target_species="rat"`, `fdr_threshold=0.05`); results
  driven by the confirmed QueryBuilder run contract; honest q=0.0584 headline preserved.
- `apps/web/src/components/QueryBuilder.source.property.test.tsx` (Property 7),
  `QueryBuilder.selectors.property.test.tsx` (Property 9),
  `QueryBuilder.mappinggate.property.test.tsx` (Property 10) — new.
- `apps/web/src/domain/analysis.test.ts` — extended with FDR-boundary + reconciliation cases.
- `apps/api/tests/test_service.py` — Property 5 (run_id) + Property 11 (species-derived study) present.
- `DECISIONS.md` (ADR-0016), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R4 rows) — traceability trail (R13).
- `service._run_id` NOT re-implemented — already records `fdr_threshold` and excludes
  `include_nonsignificant` (ADR-0008), verified by Property 5.

Why each change was necessary:
- Drive the UI from `/api/catalog` so the store is the single source of truth for available targets
  (ADR-0009, ADR-0016) — removes the last in-memory-catalog drift surface (R4 AC5, AC9).
- Exactly-one signature source with `signature_rows` vs `signature_csv_text` mapping (R4 AC1–AC3) and
  invalid-upload rejection that retains the prior source (R4 AC4).
- Total vs partial catalog degradation (R4 AC6, AC7): total failure disables all + errors; a single
  missing dimension disables only its own selector.
- Mandatory mapping-preview confirmation gate so ambiguous mappings are a researcher decision, never a
  silent first-wins (R4 AC10, AC11, AC12; ADR-0007 boundary).
- FDR constrained to `0 < t <= 1` with the boundary note (R4 AC8); FDR in run_id / include-nonsig out
  of run_id kept a backend concern (R4 AC13, AC14; ADR-0008 / Property 5).

Commands run — verify-before-claiming (R13.4):
- From `apps/web`: `npm run test:web` (`vitest run`)
  → `Test Files  6 passed (6)  Tests  14 passed (14)`:
  - `QueryBuilder.source.property.test.tsx` — 2 passed (Property 7: exactly one source active;
    example/paste → `signature_rows`|`example_name`, upload → `signature_csv_text`; invalid upload
    rejected, prior source retained, error surfaced, no source switch).
  - `QueryBuilder.selectors.property.test.tsx` — 2 passed (Property 9: each per-dimension selector
    enabled iff its dimension has ≥1 value — empty dimension disables only its own selector; total
    catalog failure disables ALL selectors + catalog-unavailable error).
  - `QueryBuilder.mappinggate.property.test.tsx` — 3 passed (Property 10: comparison gated —
    `onRun` fires exactly once, only after confirm, with the preview attached; every
    ambiguous/multi-candidate row surfaced with all candidates and none pre-selected; preview failure
    blocks the comparison, shows an error, and retains selections).
  - `Discordance.test.tsx` (4), `Discordance.property.test.tsx` (2), `ResultsTable.property.test.tsx`
    (1) — unrelated to R4, still passing.
- From `apps/web`: `npm run test:domain`
  (`node --experimental-strip-types --test src/domain/analysis.test.ts`)
  → `# tests 26  # pass 26  # fail 0` (extended from 18 → 26 with the FDR-boundary `0 < t <= 1` and
  catalog-driven reconciliation cases; R4 AC8).
- From `hackathon/tool` (venv `hackathon/tool/.venv`):
  `.venv/bin/python -m pytest ../../apps/api/tests/test_service.py -q` → `13 passed`, and
  `--collect-only -q` → 13 tests. Includes:
  - `test_property5_run_id_determinism` (Property 5, `# Feature: finalize-discordance-mvp, Property 5`,
    Validates Requirements 4.13, 4.14): identical inputs → identical run_id; changing `fdr_threshold`
    changes the run_id; toggling `include_nonsignificant` (presentation-only) does NOT — confirms the
    FDR-in / include-nonsig-out boundary without re-implementing `_run_id`.
  - `test_property11_species_derived_study` + `test_property11_unsupported_species_has_no_design`
    (Property 11): study design derived from species alone (rat → chronic, human → acute); unsupported
    source species never fabricate a design.

Boundary note recorded (R4 AC8): the UI presets treat the FDR upper bound as INCLUSIVE (`0 < t <= 1`,
so `1.0` is selectable) while the API request validation expresses it as `gt=0` / `lt=1` (exclusive at
`1.0`). This UI-inclusive-upper vs API-exclusive-upper difference is documented, not silently
reconciled.

Scope note (honest): R4 AC13/AC14 (FDR-in-run_id, include-nonsignificant-out) are a `service._run_id`
behavior already delivered under ADR-0008 and re-verified here by Property 5; no new run_id code was
written for R4. The QueryBuilder simply treats FDR as an analysis parameter and include-nonsignificant
as presentation-only in the request it sends.

Result: `QueryBuilder.tsx` catalog-driven rewire + exactly-one signature source + total/partial
selector degradation + FDR presets/boundary note + mandatory mapping-preview confirmation gate, and
the `LiveDashboard.tsx` de-hardcoding, verified against mocked API-client payload shapes. All R4 tests
PASS: web 14/14 across 6 files, domain 26/26, service 13/13. `service._run_id` unchanged; `apps/api`
engine behavior unchanged (web-only changes + the pre-existing tsc fix).

### R5 checkpoint — Generalized query upload + shared ResultsTable labeling invariant (PASS)
Files added / changed:
- `apps/api/motrpac_probe_service/generalized_query.py` (task 7.2) — added `run_uploaded(signature_csv_text,
  *, tissue, reference_contrast_category, reference="default")` mirroring `run_bundled`: guards on
  `available()`, validates against the query_core schema via `_validate_query_core_schema` (REUSES
  query_core's own loader — no reimplemented schema) returning `(ok, failing_requirement)`, raises
  `SchemaError(which)` on non-conforming input WITHOUT running (prior signature unchanged), writes a
  CONTENT-ADDRESSED sha256 temp file on success, runs `_run_query`, normalizes identically to
  `run_bundled`, and cleans up the temp file. Computes NO statistics (query_core is source of truth).
- `apps/api/motrpac_probe_service/app.py` (task 7.3) — `POST /api/generalized/query` accepting
  `{ signature_csv_text, tissue, reference_contrast_category, reference }`; returns the normalized
  `run_bundled`-shaped payload on success and 422 NAMING the failing schema requirement on
  non-conforming input (translated from `SchemaError`).
- `apps/web/src/api/client.ts` (task 7.4) — typed `generalizedQueryUpload(body)` + `GeneralizedQueryUploadInput`
  → `POST /generalized/query`.
- `apps/web/src/components/GeneralizedQuery.tsx` (task 7.5) — upload + paste controls alongside the
  bundled picklist; a pasted signature that cannot be parsed is rejected with a parse-failure
  indication WITHOUT clearing the pasted text.
- `apps/web/src/components/ResultsTable.tsx` (task 7.6) — shared component: `camera_p` (raw nominal p)
  immediately adjacent to `camera_fdr` (BH q), each with a visible labeled header (q = adjusted over
  the frozen 52-column family; p = nominal/unadjusted); distinctly labeled `species` and `dataset`
  columns; the row-model TYPE REQUIRES species + dataset + contrast so no bare value can render;
  cross-species rows show the ortholog-link relation; where a headline is surfaced it is the
  q = 0.0584 not-significant result, with the 16-column q = 0.0413 only ever a labeled sensitivity
  analysis (R12.2, R12.3). Reads existing `ColumnResult` fields — no backend change.
- Adoption (task 7.7): `ResultsTable` adopted in the Live set-level `ColumnResult` table (where a
  genuine results-row model exists). Generalized/Discordance/Metabolomics KEPT their bespoke tables
  because their row shapes (rank correlations / class counts / metabolite hits) lack the
  `camera_p`/`camera_fdr` + species/dataset/contrast row model; forcing the shared component onto them
  would fabricate fields. Cohorts kept separate in presentation (R12.5).
- `apps/api/tests/test_generalized.py` (task 7.1) — Property 6 as both a unit test
  (`test_property6_upload_schema_enforcement`) and an endpoint test
  (`test_property6_endpoint_rejects_nonconforming`).
- `apps/web/src/components/ResultsTable.property.test.tsx` (task 7.8) — Property 8.
- `DECISIONS.md` (ADR-0017), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R5 rows) — traceability trail (R13).

Why each change was necessary:
- Put the generalized pipeline behind an upload/paste path in-process (ADR-0017), reusing the
  ADR-0012 adapter pattern so query_core stays the sole source of truth and the upload path performs
  no statistics (R11.1). Validate-then-run so non-conforming input is rejected WITHOUT running and the
  response names the unmet schema requirement (R5 AC1, AC2); the content-addressed temp file never
  trusts a client path (ADR-0008).
- Parse failure on a pasted signature retains the analyst's text (R5 AC4).
- A shared `ResultsTable` whose row-model TYPE requires species + dataset + contrast makes the
  never-a-bare-value invariant hold by construction (R5 AC7, AC8), with `camera_p` adjacent to
  `camera_fdr` under labeled headers (R5 AC5), distinct species + dataset columns (R5 AC6), and the
  ortholog-link relation on cross-species rows (R5 AC9); the honest q=0.0584 headline is preserved and
  the 16-column q=0.0413 is only a labeled sensitivity analysis (R12.2, R12.3).

Commands run — verify-before-claiming (R13.4):
- From `hackathon/tool` (venv `hackathon/tool/.venv`):
  `.venv/bin/python -m pytest ../../apps/api/tests/test_generalized.py -q` → `7 passed` (was 5 + 2
  Property-6 tests). `--collect-only -q` → 7 tests, incl.:
  - `test_property6_upload_schema_enforcement` (Property 6, `# Feature: finalize-discordance-mvp,
    Property 6`, Validates Requirements 5.1, 5.2): the endpoint runs the query iff the signature
    conforms to the query_core schema; non-conforming input is rejected WITHOUT running, returns an
    error naming the unmet schema requirement, prior signature unchanged (`query_core.engine.run_query`
    mocked to keep iterations cheap).
  - `test_property6_endpoint_rejects_nonconforming` (Property 6 at the endpoint: non-conforming → 422
    naming the failing schema requirement, no run).
- From `hackathon/tool`: `.venv/bin/python -m pytest ../../apps/api/tests -q` → `64 passed`
  (full api suite; the generalized additions do not regress the other tests).
- From `apps/web`: `npm run test:web` (`vitest run`) → `Test Files  7 passed (7)  Tests  19 passed (19)`:
  - `ResultsTable.property.test.tsx` — 1 passed (Property 8, `# Feature: finalize-discordance-mvp,
    Property 8`: every rendered row shows species / dataset / contrast; no row dropped; cross-species
    rows show the ortholog-link relation — ≥100 fast-check iterations).
  - `GeneralizedQuery.test.tsx` — 5 passed (upload + paste controls; parse-failure preserves pasted
    text).
  - `Discordance.test.tsx` (4), `Discordance.property.test.tsx` (2 — Property 14),
    `QueryBuilder.selectors.property.test.tsx` (2 — Property 9),
    `QueryBuilder.source.property.test.tsx` (2 — Property 7),
    `QueryBuilder.mappinggate.property.test.tsx` (3 — Property 10) — unrelated to R5, still passing.
- From `apps/web`: `npm run build` → OK (JS 279.69 kB / gzip 83.07 kB; CSS 42.34 kB / gzip 9.32 kB;
  25 modules transformed).

Scope note (honest): R5 AC5/AC6/AC9 labeling invariants are enforced by the shared `ResultsTable`
in the Live set-level `ColumnResult` table, where a genuine `camera_p`/`camera_fdr` +
species/dataset/contrast results-row model exists. The Generalized (rank correlations), Discordance
(class/coverage counts), and Metabolomics (metabolite hits) views KEEP their bespoke tables because
their row shapes lack that row model; forcing the shared component onto them would fabricate
species/dataset/contrast/`camera_*` fields those results do not have (an integrity violation). This
is the honest adoption scope, not a partial delivery — the invariant holds by construction wherever a
real results-row model is rendered, and cohorts remain separate in presentation (R12.5).

Result: `run_uploaded()` + `SchemaError` + `_validate_query_core_schema()` + `POST
/api/generalized/query` (422 naming the failing schema requirement), `generalizedQueryUpload()`,
`GeneralizedQuery.tsx` upload/paste (parse-failure preserves pasted text), and the shared
`ResultsTable` (adopted in Live; the other three views correctly kept bespoke) verified. All R5 tests
PASS: api generalized 7/7, full api 64/64, web 7 files / 19 tests, web build OK. Engine golden
unchanged (additive `run_uploaded` + endpoint; the shared table is web-only).

### R6 checkpoint — Retire mock design views, consolidate navigation, migrate content with Figma-Make attribution (PASS)
Files added / changed:
- `apps/web/src/App.tsx` (task 8.1) — removed the FIVE mock design views (architecture/01,
  technical/02, workflow/03, dashboard/04, judge-slide/08) and all their helpers; the `View` type now
  names ONLY the five live views. Discordance is the DEFAULT view on initial load; the nav presents
  EXACTLY five entries — Discordance (default, eyebrow "01"), Live results (API) ("02"), Generalized
  query ("03"), Metabolomics (ST000763) ("04"), and a single combined About / Methods ("05"). Added a
  PURE redirect resolver: `resolveView(requested, lastLiveView)` (a live identifier → itself; any
  retired/unknown identifier → the last-used live view when that is itself live, else the default
  Discordance) and `initialView()` (runs before first render, reads the requested identifier and the
  last-used live view from sessionStorage key `esx.lastView`). A retired view is NEVER returned, so a
  direct navigation to a retired identifier cannot mount retired content. `AboutMethods` links to the
  migrated docs (`DESIGN_PROVENANCE.md`, `WORKFLOWS.md`). Removed the banned mock identifiers PPARGC1A
  / SOD2 / COL1A1 from navigable views (the `QueryBuilder` placeholder → neutral GENE1/GENE2).
- `DESIGN_PROVENANCE.md` (task 8.2) — migrated the content of the removed views 01–04+08, each with a
  per-view Figma-Make design-origin attribution line preserved (plus a global Figma-Make provenance
  section), referencing the R1 `WORKFLOWS.md` mermaid diagrams by name with NO diagram-source change.
- `apps/web/src/App.property.test.tsx` (task 8.3) — Property 15 (retired-view redirect resolution).
- `apps/web/src/App.test.tsx` (task 8.4) — nav = exactly five live entries + Discordance default +
  aggregate banned-identifier grep over navigable component source.
- `DECISIONS.md` (ADR-0018), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R6 rows) — traceability trail (R13).

Why each change was necessary:
- Removing the five mock views from navigation is the direct satisfaction of R6 AC1/AC6 — a user is
  never shown static mock content presented as real results.
- A PURE `resolveView`/`initialView` resolver makes R6 AC2 hold by construction: because it can never
  return a retired identifier, a direct navigation to a retired route cannot mount retired content
  (stronger than a runtime guard that could be bypassed or forgotten). It degrades gracefully to the
  last-used live view (sessionStorage `esx.lastView`) or the default Discordance.
- Migrating the content into `DESIGN_PROVENANCE.md` with per-view Figma-Make attribution preserves the
  design story and its authorship (R6 AC4/AC5); referencing the R1 diagrams by name avoids diagram
  drift (ADR-0013). `AboutMethods` links to the migrated docs (R6 AC7).
- Removing PPARGC1A/SOD2/COL1A1 from navigable views satisfies R6 AC3 — the fabricated identifiers
  never appear where a user could read them as data.

Commands run — verify-before-claiming (R13.4):
- From `apps/web`: `npm run test:web` (`vitest run`) → `Test Files  9 passed (9)  Tests  34 passed (34)`:
  - `App.property.test.tsx` — 4 passed, including Property 15 (`# Feature: finalize-discordance-mvp,
    Property 15`, retired-view redirect resolution, Requirement 6.2): `resolveView` never returns a
    retired view and always returns a live view; a retired/unknown identifier resolves to the last-used
    live view or, absent one, to Discordance.
  - `App.test.tsx` — 11 passed: the nav presents exactly the five live entries (Discordance, Live,
    Generalized, Metabolomics, About/Methods), Discordance is the default view on load, and the
    aggregate banned-identifier test reads each navigable component source from disk and asserts none
    contains PPARGC1A/SOD2/COL1A1.
  - `GeneralizedQuery.test.tsx` (5), `ResultsTable.property.test.tsx` (1 — Property 8),
    `Discordance.test.tsx` (4), `Discordance.property.test.tsx` (2 — Property 14),
    `QueryBuilder.selectors.property.test.tsx` (2 — Property 9),
    `QueryBuilder.source.property.test.tsx` (2 — Property 7),
    `QueryBuilder.mappinggate.property.test.tsx` (3 — Property 10) — unrelated to R6, still passing.
  - Note: one full-suite run showed a single flaky `waitForWrapper` async-timeout in `App.test.tsx`
    under concurrent load; a re-run of the full suite passed 34/34, and the two App test files run in
    isolation (`vitest run src/App.test.tsx src/App.property.test.tsx`) pass 15/15 deterministically
    (11 + 4). The failure was a test-runner timing artifact, not an App or resolver defect.
- From `apps/web`: `npm run build` → OK (JS 279.69 kB / gzip 83.07 kB; CSS 42.40 kB / gzip 9.33 kB;
  25 modules transformed).
- Banned-identifier grep over navigable source (R6 AC3), verify-before-claiming:
  `grep -rnE 'PPARGC1A|SOD2|COL1A1' apps/web/src --include=*.tsx --include=*.ts` returns matches ONLY
  in test files — `App.test.tsx` (the comment, the `BANNED_IDENTIFIERS` array, and the assertion) and
  `QueryBuilder.source.property.test.tsx` (`VALID_CSV`/`VALID_PASTE` fixtures). Excluding test files
  (`… | grep -vE '\.test\.'`) yields ZERO matches: no navigable component source contains any banned
  identifier. The remaining test-file matches assert/exercise their absence.

Result: `App.tsx` renders only the five live views (Discordance default; nav = exactly five entries);
the pure `resolveView`/`initialView` resolver (sessionStorage `esx.lastView`) maps any retired/unknown
identifier to the last-used live view or Discordance and never returns a retired view; the mock
content of views 01–04+08 lives in `DESIGN_PROVENANCE.md` with per-view Figma-Make attribution and
R1-diagram references (no diagram-source change); `AboutMethods` links to the migrated docs;
PPARGC1A/SOD2/COL1A1 removed from all navigable views. All R6 tests PASS: web 9 files / 34 tests, web
build OK, banned-identifier grep clean over navigable source. Engine and `apps/api` behavior unchanged
(web-only changes + docs migration).

### R7 checkpoint — root README rewritten, anchored to Track 2 (DOCS ONLY) (PASS)
Files added / changed:
- `README.md` (rewritten) — Track 2 anchor in the first section; discordance catalog + predictive
  model as the headline deliverable (also Deliverables item 1); honest framing (cross-cohort, not PAH
  efficacy; cohorts separate/not merged; PTM occupancy caveat); one navigable link each to
  `WORKFLOWS.md`, `DECISIONS.md`, `VERIFICATION.md`, `METABOLOMICS.md`, `RUN_LOCAL.md`; no empty
  heading. Replaced the prior empty section skeleton.
- `DECISIONS.md` (ADR-0019), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R7 rows) — traceability trail (R13).
- NO runtime source file added, removed, or altered (R7 is DOCS ONLY).

Why each change was necessary:
- Make the README the honest front door (ADR-0019): lead with the Track 2 challenge and the headline
  discordance deliverable, state the framing/caveats up front, and give a reader exactly one link to
  each key doc — closing the prior empty-section skeleton (which left headings with no body).

Structure + link-resolvability grep checks (commands run in repo root) — verify-before-claiming (R13.4):
- Heading order — `grep -nE '^#{1,6} ' README.md | head -6`:
  1 `# Stanford Multi-omics Hackathon 2026 Track 2`; 3 `## Omic Discordance Explained`;
  11 `## Deliverables`; 17 `## Data`; 21 `## Honest framing`; 28 `## Documentation map`.
  ⇒ `## Omic Discordance Explained` is the first section, before the second heading `## Deliverables`
  (R7 AC1). The Track 2 challenge is named in that first section.
- Headline (R7 AC2) — the discordance catalog + predictive model is named as the headline deliverable
  in the first section (line 7, "The headline deliverable is the discordance catalog and its
  predictive model") AND as Deliverables item 1 (line 13), before the PTM-parent audit and the
  directed/generalized analyses.
- Framing (R7 AC3, AC4) — `grep -nE 'cross-cohort association|not a claim of PAH treatment efficacy' README.md`,
  `grep -nE 'analyzed separately and are not merged' README.md`,
  `grep -niE 'PTM occupancy caveat applies' README.md` → all three present on line 9
  ("cross-cohort association findings, not a claim of PAH treatment efficacy"; "analyzed separately and
  are not merged"; "PTM occupancy caveat applies to all reported findings").
- Exactly one link per target (R7 AC5) — for each of `WORKFLOWS.md`, `DECISIONS.md`, `VERIFICATION.md`,
  `METABOLOMICS.md`, `RUN_LOCAL.md`: `grep -oE "\]\(<target>\)" README.md | wc -l` → 1, 1, 1, 1, 1.
  (First pass showed VERIFICATION.md=2 and RUN_LOCAL.md=2 because the Setup section re-linked them;
  corrected so Setup references those docs in prose, restoring exactly-one each.)
- Link targets resolve (R7 AC5) —
  `grep -oE '\]\([A-Za-z_]+\.md\)' README.md | sed -E 's/\]\(([^)]+)\)/\1/' | sort -u`
  then `test -f` each → `DECISIONS.md`, `METABOLOMICS.md`, `RUN_LOCAL.md`, `VERIFICATION.md`,
  `WORKFLOWS.md` all RESOLVE (files exist at repo root).
- No empty heading (R7 AC6) — an awk scan that treats a following DEEPER heading as valid subordinate
  content (so the H1 title followed by the H2 section is NOT empty) and otherwise requires non-blank
  body text before the next heading/EOF → reports NO empty headings. Every heading is followed by body
  text or a subordinate heading.

Result: `README.md` is documentation only; no runtime source changed. All six R7 acceptance criteria
(AC1–AC6) PASS by the grep/structure checks above. Engine, `apps/api`, and `apps/web` behavior
unchanged (docs only).

### R8 checkpoint — fresh-clone bootstrap + pnpm standardization + `parents[3]` docs (PASS)
Files added / changed (R8 work recorded here; this checkpoint is task 11.4, DOCS ONLY):
- `scripts/bootstrap.sh` (task 11.1) — one-command fresh-clone bootstrap: engine-first editable
  installs (`hackathon/tool[test]` then `apps/api[api,test]`), `mprobe store fetch`, a named smoke
  check (engine import → store present → service import → `catalog.build_catalog()` non-empty),
  `set -euo pipefail` + per-step guards naming the failed target/check, and a 600 s deadline guard.
- `.github/workflows/ci.yml` + `apps/web/package.json` (task 11.2) — standardized on pnpm: removed the
  stray `apps/web/package-lock.json`, added `"packageManager": "pnpm@10.34.3"`, every web CI step uses
  `pnpm/action-setup@v4` + `pnpm install --frozen-lockfile` + `pnpm test` / `pnpm run build`
  (honoring the committed `apps/web/pnpm-lock.yaml`); no standalone `npm` invocation remains.
- `BOOTSTRAP.md` (task 11.3) — documents the one-command invocation (linked once each from `README.md`
  and `RUN_LOCAL.md`) and the `parents[3]` path-coupling constraint (depth rule + a one-line verifier).
- `DECISIONS.md` (ADR-0020), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R8 rows) — traceability trail (R13).
- NO runtime source file or the bootstrap script was added, removed, or altered by task 11.4 (DOCS ONLY).

Why each change was necessary:
- One command with a named-target fail-fast contract is the honest form of R8 (ADR-0020): every
  failure mode names the specific target/check and exits non-zero rather than leaving a half-installed
  environment (R8 AC2–AC4); the engine-first ordering matches the real dependency direction (the
  service depends on the engine editable package); the end-to-end smoke check proves the app can build
  its catalog, not merely that pip exited zero (R8 AC1).
- A single declared package manager (pnpm) removes the dual-lockfile drift the Gate 1 checkpoint
  flagged as an open item (R8 AC6).
- Documenting `parents[3]` turns an invisible, silently-breakable coupling into something a reader can
  verify (R8 AC5, AC7).

Commands run (repo root) — verify-before-claiming (R13.4):
- Bootstrap syntax (R8 AC1–AC4 script well-formedness):
  - `bash -n scripts/bootstrap.sh && echo "SYNTAX_OK"` → `SYNTAX_OK` (exit 0).
- Embedded smoke-check heredoc compiles (R8 AC4 — the named-check Python block is valid):
  - extract: `awk "/<<'PYEOF'/{f=1;next} /^PYEOF$/{f=0} f" scripts/bootstrap.sh > /tmp/bootstrap_smoke.py`
    → 34 lines.
  - `python3 -m py_compile /tmp/bootstrap_smoke.py && echo "PY_COMPILE_OK"` → `PY_COMPILE_OK` (exit 0).
- Single package manager — no standalone `npm` invocation remains (R8 AC6):
  - word-boundary grep (excludes `pnpm`) over CI:
    `grep -rnE '(^|[^p[:alnum:]])npm (install|test|run|ci|exec)' .github/`
    → NO match (grep exit 1) → `NO_STANDALONE_NPM_IN_.github`.
  - over the setup docs:
    `grep -nE '(^|[^p[:alnum:]])npm (install|test|run|ci|exec)' BOOTSTRAP.md RUN_LOCAL.md README.md`
    → NO match (grep exit 1) → `NO_STANDALONE_NPM_IN_SETUP_DOCS`.
  - (A naive `grep -rnE 'npm (install|…)'` DOES match the `pnpm install`/`pnpm run …` lines because
    "pnpm install" contains "npm install" as a substring — those hits are pnpm, not npm; the
    word-boundary form above disambiguates and confirms zero genuine npm invocations.)
  - pnpm is present in CI: `grep -nE 'pnpm' .github/workflows/ci.yml` → `pnpm/action-setup@v4`,
    `cache: "pnpm"`, `pnpm install --frozen-lockfile`, `pnpm test`, `pnpm run build` (both web and e2e
    jobs).
  - `apps/web/package.json` line 6: `"packageManager": "pnpm@10.34.3"`.
  - `test -f apps/web/package-lock.json` → ABSENT (removed); `test -f apps/web/pnpm-lock.yaml` → PRESENT.
- `BOOTSTRAP.md` link resolvability + exactly-one-link invariant (R8 AC5):
  - link TO `BOOTSTRAP.md` from README: `grep -oE '\]\(BOOTSTRAP\.md\)' README.md | wc -l` → `1`
    (README.md:38 — "For a fresh clone, [BOOTSTRAP.md](BOOTSTRAP.md) documents the one-command bootstrap…").
  - link TO `BOOTSTRAP.md` from RUN_LOCAL: `grep -oE '\]\(BOOTSTRAP\.md\)' RUN_LOCAL.md | wc -l` → `1`
    (RUN_LOCAL.md:12 — "the one-command bootstrap documented in [BOOTSTRAP.md](BOOTSTRAP.md)").
  - target resolves: `test -f BOOTSTRAP.md` → EXISTS.
  - links FROM `BOOTSTRAP.md` resolve: its only `*.md` link target is `RUN_LOCAL.md` → RESOLVES
    (`test -f` passes).

Bootstrap contract recorded (from `scripts/bootstrap.sh`, R8 AC1–AC4):
- Order: create/reuse `hackathon/tool/.venv` → editable-install engine (`pip install -e
  "hackathon/tool[test]"`) → editable-install service (`pip install -e "apps/api[api,test]"`) →
  `mprobe store fetch` → named smoke check → 600 s deadline assertion → `bootstrap succeeded` (exit 0).
- Engine-first: the engine install precedes the service install because the service depends on the
  engine as an editable package (a failure in the engine install stops before the service is attempted).
- Named failure modes: engine-install failure names `hackathon/tool`; service-install failure names
  `apps/api`; store-fetch failure emits a store-fetch error; a failed smoke check names the specific
  check (`import motrpac_probe (engine)` / `engine store present` / `import motrpac_probe_service
  (service)` / `service catalog builds`); overrunning `BOOTSTRAP_DEADLINE_S=600` reports elapsed vs
  deadline — every mode exits non-zero.
- Smoke check (the py_compile'd heredoc) asserts, in order: `import motrpac_probe`;
  `store.CONTRASTS.exists()`; `import motrpac_probe_service`; `catalog.build_catalog()` returns a
  non-empty result.

`parents[3]` path-coupling recorded (R8 AC7, from `BOOTSTRAP.md`): `Path(__file__).resolve().parents[3]`
resolves to the repo root ONLY when an adapter sits exactly three directories below the root
(`apps/api/motrpac_probe_service/<name>.py`); the doc gives the by-hand depth rule (`depth == 3`) and a
one-line `python3` verifier. The `apps/api` suite pins this with
`test_parents3_path_coupling` in `apps/api/tests/test_discordance.py` (see the R2 checkpoint).

Scope note (honest): task 11.4 is DOCS ONLY — it records the already-completed R8 work (11.1 bootstrap
script, 11.2 pnpm standardization, 11.3 `BOOTSTRAP.md`) and did NOT modify the bootstrap script or any
runtime source. The checks above are static (`bash -n`, `py_compile`, grep, `test -f`); a full runtime
`bash scripts/bootstrap.sh` execution (actual editable installs + `mprobe store fetch` + smoke check
within 600 s) is a live-environment run and is not part of this docs-only traceability task.

Result: R8 bootstrap contract, single-package-manager (pnpm) standardization, and `parents[3]`
documentation verified by the static checks above — `bash -n` SYNTAX_OK, smoke-check heredoc
PY_COMPILE_OK, zero standalone `npm` invocations in CI or setup docs, and `BOOTSTRAP.md` linked exactly
once each from `README.md` and `RUN_LOCAL.md` with all links resolving. Engine, `apps/api`, and
`apps/web` behavior unchanged (docs only).

### R9 checkpoint — automated verification (httpx smoke + Playwright e2e) wired into CI, full-battery gate + golden-hash guard (PASS)
Files added / changed (R9 work recorded here; this checkpoint is task 12.5, DOCS ONLY):
- `apps/api/tests/test_smoke_httpx.py` (task 12.1) — a SINGLE-execution httpx smoke check that boots
  the FastAPI app IN-PROCESS via `httpx.ASGITransport`, requests every documented endpoint (health,
  catalog, availability, mappings/preview, signatures/validate, directed `/api/comparisons` + run
  subresources, the three metabolomics endpoints, the generalized list + bundled + UPLOAD endpoints,
  and the three discordance endpoints), and asserts exact values — headline `q == 0.0584 ± 0.0001`
  AND not-significant; discordance per-class integer counts (supported_concordant=1,
  supported_opposite=0, rna_response_protein_equivalent=285, indeterminate=5642, each an `int`); the
  metabolomics null result (`hits_fraction == 0.0`, `fisher_p == 0.609`). It NAMES the offending
  endpoint/value on any non-success or mismatch (R9 AC3) and SKIPS (never fabricates a pass) when the
  store or a committed demo is absent.
- `apps/web/e2e/live-views.spec.ts` + `apps/web/playwright.config.ts` (task 12.2) — the config boots
  BOTH real processes (API via `uvicorn motrpac_probe_service.app:app` from `apps/api`, then Vite via
  `pnpm run dev` with `VITE_API_BASE_URL` pointed at the API) per the RUN_LOCAL.md port model
  (8765/8443, env-overridable), pnpm per R8 AC6. The spec asserts, per live view (Discordance, Live,
  Generalized, Metabolomics), that REAL data renders and the error/empty markers are absent (R9 AC4),
  and it drives the browser-level point→evidence→source chain from a rendered Live point through the
  page's run-scoped export/report links, `/features` evidence (`evidence_id` run-scoped, `column_id`
  shown), the store-backed `source_feature_id`, and the report's verbatim honest-framing source text
  (R9 AC5) — the browser counterpart to `test_point_to_evidence_to_source_chain`.
- `.github/workflows/ci.yml` (tasks 12.3 + 12.4) — three jobs: `python` (golden-hash guard → engine
  EXACTLY-122/34 gate → `apps/api` tests → httpx must-not-skip gate), `e2e` (installs engine+service
  with `[api]` uvicorn, fetches the store, installs web deps with pnpm `--frozen-lockfile` + chromium,
  runs `pnpm run test:e2e`), and `web` (`pnpm test` + `pnpm run build`), all on push to
  `integration/**`/`main` and on every pull request.
- `DECISIONS.md` (ADR-0021), `VERIFICATION.md` (this checkpoint),
  `REQUIREMENTS_TRACEABILITY.md` (R9 rows) — traceability trail (R13).
- NO runtime source, CI, test, or the golden-hash file was added, removed, or altered by task 12.5
  (DOCS ONLY).

Why each change was necessary:
- An in-process httpx single pass proves the whole HTTP surface is wired AND the exact committed
  headline/discordance/metabolomics values still flow through it, naming any drift (R9 AC1–AC3);
  a real-boot Playwright battery proves the four live views render real data and the
  point→evidence→source chain holds in the browser (R9 AC4/AC5); wiring both into CI as PR-triggered
  jobs makes the guarantee automatic (R9 AC6); the countable full-battery gate (exactly 122/34,
  `apps/api` + web tests/build zero failures) makes "green" an enforceable merge condition (R9
  AC7/AC8); the golden-hash guard fails closed on any change to the frozen `analysis.py` (R9 AC9).
  See ADR-0021.

Commands run for THIS docs-only checkpoint — verify-before-claiming (R13.4). These are cheap STATIC
checks re-run here plus results RECORDED by the prior tasks (12.1–12.4); the full LIVE CI battery is
NOT re-executed by this docs-only task (see the scope note below):

STATIC checks re-run here (real outputs):
- Golden-hash guard (R9 AC9), from `hackathon/tool`:
  `shasum -a 256 -c analysis.py.sha256` → `analysis.py: OK` (exit 0). The recorded reference in
  `hackathon/tool/analysis.py.sha256` is
  `1c1089723ec7c9d67c4c0f7ab51262717f80d95322c7da8df60187a06cab7372  analysis.py`, and the working
  `analysis.py` hashes byte-identical to it ⇒ the frozen scientific core is unchanged.
- CI gate literals (R9 AC7/AC8), `grep -nE 'EXPECTED_PASSED|EXPECTED_SKIPPED' .github/workflows/ci.yml`
  → line 76 `EXPECTED_PASSED, EXPECTED_SKIPPED = 122, 34` and line 77 the guard
  `if passed != EXPECTED_PASSED or skipped != EXPECTED_SKIPPED:` — the engine step fails the battery
  (and blocks the merge) on any deviation from EXACTLY 122 passed / 34 skipped.
- CI smoke must-not-skip gate (R9 AC1–AC3 wiring), in the `python` job: after `mprobe store fetch`, the
  step runs `pytest … test_smoke_httpx.py -rs` and fails if the output contains `no tests ran` /
  ` skipped` / `SKIPPED` or does not report `1 passed` (ci.yml lines 95–98) — so the smoke must
  actually execute in CI, not silently skip.
- ci.yml YAML validity + job set (R9 AC6): `ruby -ryaml -e 'd=YAML.load_file(".github/workflows/ci.yml")
  …'` → `YAML_OK` / `jobs: python, e2e, web` — the dedicated `e2e` job and the smoke-carrying `python`
  job are both present; `apps/web/package.json` defines `test:e2e` (`playwright test`) and
  `test:e2e:install` (`playwright install --with-deps chromium`), the scripts the `e2e` job invokes.
- Referenced files exist (R13.4): `apps/api/tests/test_smoke_httpx.py`, `apps/web/e2e/live-views.spec.ts`,
  `apps/web/playwright.config.ts`, `.github/workflows/ci.yml`, `hackathon/tool/analysis.py.sha256`, and
  `hackathon/tool/analysis.py` all EXIST (`test -f` each).

RECORDED from the prior tasks' runs (NOT re-executed here):
- Task 12.1 httpx smoke: `python -m pytest apps/api/tests/test_smoke_httpx.py` → `1 passed` (the single
  smoke pass across every documented endpoint with all exact values matching; this is the result the
  CI must-not-skip gate asserts).
- Task 12.2 Playwright e2e: `pnpm run test:e2e` locally → `5 passed` (the four per-view real-data specs
  under "R9 AC4 — every live view renders real data" + the one "R9 AC5 — browser-level
  point → evidence → source chain" spec).
- Engine full-battery gate (R9 AC7): the engine `pytest -q -rs` count matches the counts recorded
  across prior VERIFICATION checkpoints (Gate 2 and later: `122 passed, 34 skipped` over the
  store/deck/R skip allowlist) — the exact pair the ci.yml gate enforces. The `apps/api` and web
  `pnpm test` + `pnpm run build` legs are recorded green in the R2–R8 checkpoints above.

Exact expected values asserted by the smoke check (source-of-truth: engines + committed demos):
- headline `q = 0.0584 ± 0.0001`, labeled NOT significant (R12.2).
- discordance `classification_counts`: supported_concordant=1, supported_opposite=0,
  rna_response_protein_equivalent=285, indeterminate=5642 (each an integer).
- metabolomics null result: blood exercise-sensitive `hits_fraction = 0.0`, Fisher `p = 0.609`.

Honest static-vs-recorded split (R13.4): STATIC-VERIFIED here = the golden hash (`shasum -c` → OK),
the ci.yml `122/34` gate literals + smoke must-not-skip gate, ci.yml YAML validity / job set, and the
existence of all referenced files. RECORDED-FROM-PRIOR-RUNS (not re-run in this docs-only task) = the
httpx smoke `1 passed`, the e2e `5 passed`, and the engine `122 passed / 34 skipped` battery. This
checkpoint does NOT claim the full live CI battery (engine + apps/api + web + e2e all executing in CI)
ran as part of task 12.5 — it did not; that is a live-environment run, and only the cheap static checks
above were re-executed here.

Result: R9's automated verification is in place and wired into CI — the in-process httpx smoke (exact
headline/discordance/metabolomics values, names offender on failure, skips-not-fakes when store/demo
absent), the real-boot Playwright e2e (per-view real data + browser point→evidence→source chain), the
CI jobs (`python` with the golden-hash guard + engine 122/34 gate + `apps/api` + httpx must-not-skip
gate; dedicated `e2e`; `web` tests+build), the full-battery merge gate (EXACTLY 122/34, zero
`apps/api`/web failures), and the golden-hash guard on `analysis.py`. Static checks re-run here PASS
(golden hash OK at `1c1089723ec7c9d67c4c0f7ab51262717f80d95322c7da8df60187a06cab7372`; ci.yml gate =
`122, 34`; YAML valid, jobs `python, e2e, web`; all referenced files exist); the smoke `1 passed`,
e2e `5 passed`, and engine `122/34` are recorded from the prior tasks. Engine, `apps/api`, and
`apps/web` behavior unchanged (docs only).
### R10 checkpoint — merge-to-main plan (`MERGE_PLAN.md`) structure + no-execution confirmation (PASS)
Files added / changed (R10 work recorded here; this checkpoint is task 13.1, DOCS ONLY):
- `MERGE_PLAN.md` (task 13) — a DOCUMENTATION-ONLY merge-to-main plan at the repo root that executes
  nothing (R10 AC1). Section 1 documents the three pre-merge preconditions — all verification gates
  passing; `git status` reporting zero modified/staged/untracked/deleted files; the local branch zero
  commits ahead of upstream (R10 AC2), each with a copy-pasteable verifier command shown as text.
  Section 2 documents the block rules (any non-passing gate, OR any uncommitted/untracked file, OR
  branch ≥1 commit ahead → BLOCKED, and identify WHICH precondition failed — R10 AC3). Section 3
  documents the no-rebase-over-`origin/t3code/build-motrpac-probe-tool` constraint (R10 AC4). Section 4
  gives the EXACT `gh pr create` with base `main` / head `integration/exercise-signature-explorer`
  (R10 AC5). Section 5 gives the PR template with EXACTLY three labeled sections — `## Summary`,
  `## What was tested`, `## Deferred items` (R10 AC6). Every `git`/`gh` command lives inside a fenced
  code block (documented, not executed).
- `DECISIONS.md` (ADR-0022), `VERIFICATION.md` (this checkpoint), `REQUIREMENTS_TRACEABILITY.md`
  (R10 rows) — traceability trail (R13).
- NO runtime source, CI, or test was added, removed, or altered by task 13/13.1 (DOCS ONLY), and NO
  `git`, `gh`, or shell VCS command was run to author or verify the plan.

Why each change was necessary:
- A written, documentation-only merge plan lets the team see the exact merge procedure, its gating
  preconditions, its block semantics, the no-rebase constraint, the exact PR-create command, and the
  fixed PR template WITHOUT the feature itself taking any irreversible or shared-history-affecting
  action (R10 AC1–AC6). See ADR-0022.

Commands run for THIS docs-only checkpoint — verify-before-claiming (R13.4). These are cheap HARMLESS
STATIC checks only (grep / awk / test over the documentation); NO `git`, `gh`, or any VCS/shell command
was executed, and the plan itself checks no git state (it DOCUMENTS the preconditions as text):

STATIC checks re-run here (real outputs):
- Three PR-template sections present (R10 AC6): `grep -cE '^## (Summary|What was tested|Deferred
  items)$' MERGE_PLAN.md` → `3`; `grep -nE '^## (Summary|What was tested|Deferred items)$'
  MERGE_PLAN.md` → `132:## Summary`, `142:## What was tested`, `150:## Deferred items`.
- Exact `gh pr create` base/head (R10 AC5): `grep -n 'gh pr create --base main --head
  integration/exercise-signature-explorer' MERGE_PLAN.md` → line 115
  (`gh pr create --base main --head integration/exercise-signature-explorer --title "Finalize
  Discordance MVP" --body-file PR_BODY.md`).
- No-rebase clause + `t3code` reference (R10 AC4): `grep -niE 'no rebase' MERGE_PLAN.md` → line 95
  ("**No rebase may be performed that replays commits over"); `grep -n
  'origin/t3code/build-motrpac-probe-tool' MERGE_PLAN.md` → lines 96, 99, 103, 104 (the constraint +
  its rationale + the explicit "do not run any `git rebase` whose base or `--onto` target is …").
- Three documented preconditions (R10 AC2): `grep -niE 'verification gates report a passing|git
  status.*zero|zero commits ahead' MERGE_PLAN.md` → line 26 ("All verification gates report a passing
  status."), line 41 ("`git status` reports **zero modified, zero staged, zero untracked, and zero
  deleted**"), line 52 ("the local `integration/exercise-signature-explorer` branch is **zero commits
  ahead**").
- Docs-only / executes-nothing statement (R10 AC1): `grep -niE 'executes nothing|DOCUMENTATION ONLY'
  MERGE_PLAN.md` → line 1 (title), line 4 ("**This document is DOCUMENTATION ONLY. It executes
  nothing.**"), line 160 (closing scope note "This plan is documentation only …").
- NO-EXECUTION confirmation (R10 AC1) — every `git`/`gh` command is quarantined inside a fenced code
  block (documented, not executed): an indentation-aware awk scan that toggles on `^[[:space:]]*```
  fences and strips inline-backtick spans reports `bare git/gh command lines inside code fences: 3`
  and `bare git/gh command lines OUTSIDE code fences: 2`, where the 3 in-fence lines are the real
  commands (`git status --porcelain`, `git rev-list --count @{upstream}..HEAD`, `gh pr create …`) and
  the 2 "outside" hits are PROSE sentences that merely contain the word "git" ("…run the verification
  battery, or query any remote." and "…does not inspect the git working tree, does not…"), NOT
  executable command lines. ⇒ the document shows commands only inside code blocks and executes none.

Scope note (honest): task 13.1 is DOCS ONLY — it records the R10 merge plan authored in task 13 and
altered NO runtime source, CI, or test. The checks above are HARMLESS STATIC documentation checks
(`grep`, `awk`, no VCS access). Per the R10 constraint and R13.4, NO `git`, `gh`, or shell VCS command
was executed to author or verify this plan; the plan DOCUMENTS the pre-merge preconditions as text and
does not itself check git state, run the verification battery, or query any remote.

Result: R10's merge-to-main plan is in place as documentation only — `MERGE_PLAN.md` states it executes
nothing (AC1), documents the three preconditions (passing gates / clean `git status` / zero commits
ahead — AC2) and the block rules that identify the failed precondition (AC3), forbids any rebase over
`origin/t3code/build-motrpac-probe-tool` (AC4), gives the exact `gh pr create --base main --head
integration/exercise-signature-explorer` (AC5), and supplies the three-section PR template (AC6).
Static checks re-run here PASS (3 PR sections at lines 132/142/150; exact base/head at line 115;
no-rebase clause + 4 `t3code` mentions; 3 preconditions at lines 26/41/52; executes-nothing at lines
1/4/160; all 3 real `git`/`gh` commands inside code fences, zero bare commands outside). Engine,
`apps/api`, and `apps/web` behavior unchanged (docs only).

