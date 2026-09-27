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
