# Architecture Decision Records — Exercise Signature Explorer

Each entry: Context / Decision / Reasons / Alternatives rejected / Consequences.
These records exist to make the project interpretable and traceable ("reasons for
actions"), a stated judging criterion. Newest at the bottom.

---

## ADR-0001 — mprobe is the sole scientific engine; wrap, do not rebuild
- **Context.** A teammate branch (`t3code/build-motrpac-probe-tool`, tip `9962a43`)
  contains `mprobe`, a Python package with a parquet store, typed analysis functions
  (`core`, `layers`, `nulls`, `everywhere`), a report generator with `provenance.json`,
  a Marimo app, a CLI, and a test suite (110 passed / 34 skipped on a clean run). A
  separate React "Exercise Signature Explorer" prototype exists with mock data.
- **Decision.** Treat `mprobe` as the only scientific engine. The FastAPI layer and the
  React app perform no statistics; they call/render mprobe outputs.
- **Reasons.** Avoids a second source of scientific truth; preserves the reproducible,
  provenance-tracked pipeline; matches the hackathon's Documentation & Integrity and
  Methods & Approach criteria; fastest path to a truthful vertical slice.
- **Alternatives rejected.** (a) Rebuild a backend that re-runs the R/Python scripts and
  re-normalizes to Parquet — duplicates what the store already is; creates drift.
  (b) Keep the React app's in-memory logic as the engine — it has no real statistics.
- **Consequences.** All numeric results trace to mprobe. The API is thin. Parity between
  React, the service, and Marimo becomes a testable requirement.

## ADR-0002 — Integration branch off the tool branch, merge main into it
- **Context.** The React app lives on `main`; mprobe on the tool branch. `main`
  (`57e01e7`) adds additive analysis modules under `MoTrPAC Hackathon/` (rat comparison,
  BioCarta transcriptomics) with no overlap with `hackathon/tool/`.
- **Decision.** Base `integration/exercise-signature-explorer` on the tool branch and
  merge `origin/main` into it.
- **Reasons.** The tool is the center of gravity (ADR-0001); main only contributes
  additive modules, so the merge is low-risk and additive.
- **Alternatives rejected.** Base off `main` and merge the tool branch — would make the
  larger, more complex codebase the "incoming" side and invert the mental model.
- **Consequences.** Clean merge (verified). Both the engine and the analysis modules
  coexist on the integration branch.

## ADR-0003 — Import the React app via `git subtree` (not submodule), history preserved
- **Context.** The React app is a nested standalone git repo
  (`github.com/StanchPillow55/ExerciseSignatureExplorerHandoff`, commit `d9fd564`) sitting
  untracked inside the working tree. It declares LFS rules in `.gitattributes`, but
  `git lfs ls-files` is empty and its one PNG is a normal 152 KB git blob (not an LFS
  pointer) — so there is no LFS object to materialize.
- **Decision.** Import into `apps/web` with `git subtree add` (no `--squash`), preserving
  the source history. Neutralize the LFS rules in the imported `.gitattributes` so the
  parent repo (which has no LFS) never tries to smudge these paths. Record the source URL
  and imported SHA here and in VERIFICATION.md.
- **Reasons.** Subtree gives judges a self-contained repo with no recursive clone and no
  LFS dependency; `--no-squash` preserves the teammate's authorship (Collaboration
  criterion).
- **Alternatives rejected.** Submodule (requires recursive clone + intact upstream);
  `rm -rf .git` then copy (destroys history).
- **Consequences.** A fresh clone builds the frontend without access to the original repo
  or Git LFS.

## ADR-0004 — Multiplicity family frozen at mprobe's existing 53-column BH set (MVP)
- **Context.** Three entry points adjust p-values over different comparison families
  (`compute()` ~53 columns; `analysis.py` 8; Marimo dataset-scoped). This changes
  significance: full-19 male rat SKM-GN protein at 8 wk is q=0.0584 over 53 columns
  (not significant) vs q=0.0413 over 16 (significant).
- **Decision.** Freeze the 53-column BH family as the MVP analysis contract. BH is applied
  across the full 53 before any tissue/sex/timepoint/display filtering. Provenance records
  all 53 ordered column IDs and `n_tests`. Display filtering never re-runs BH.
- **Reasons.** One honest, reproducible family; makes React/service/Marimo parity
  achievable; avoids p-hacking by post-hoc family narrowing.
- **Alternatives rejected.** Context-scoped BH (16 columns) as default — cherry-picks the
  significant answer.
- **Consequences.** Headline: male rat SKM-GN protein 8 wk is q=0.0584, NOT significant at
  q<0.05, directionally concordant only. The 16-column q=0.0413 may appear only as an
  explicitly labeled sensitivity analysis. A parity test computes (not hardcodes) 0.0584.

## ADR-0005 — Enrich the existing full-19 fixture; never recreate a 9+10 source of truth
- **Context.** The full-19 directed signature already exists at
  `hackathon/tool/examples/pah_muscle_malenfant2015.csv` (25 rows = 19 counted + 6
  `group=phenotype` context rows). It stores the PAH/control ratio in a `weight` column,
  though mprobe analysis is unweighted. Symbols are canonicalized but paper-original
  identifiers are not recorded. The lower-9 file
  (`examples/pah_muscle_lower9_malenfant2015.csv`) is the committed regression fixture.
- **Decision.** Enrich the existing full-19 fixture (or a companion provenance table) with
  a lossless superset schema that preserves BOTH paper and canonical identifiers, the
  ratio, `disease_log2fc`, reported p-value text, and correction notes. Leave the lower-9
  file unchanged. Keep `weight` only for backward compatibility, documented as not a
  statistical weight.
- **Reasons.** "The repo is law": do not fork a second source of truth. Preserve paper
  provenance (integrity) and the regression baseline.
- **Alternatives rejected.** Build a fresh 9+10 CSV — creates a duplicate/competing fixture.
- **Consequences.** A Table-2 validation test checks all 19 records, ratios, directions,
  p-values, and identifier mappings.

## ADR-0006 — OpenAPI/Pydantic is the transport source of truth; generate the TS client
- **Context.** The React app already declares a typed contract in `domain/analysis.ts`
  that closely mirrors mprobe outputs.
- **Decision.** Define request/response schemas once as Pydantic models; generate the
  TypeScript client from the OpenAPI spec. Converge toward the existing `domain/analysis.ts`
  names so the UI change is minimal. mprobe remains the scientific source of truth.
- **Reasons.** One generator, no hand-mirrored drift between Python and TypeScript.
- **Alternatives rejected.** Hand-mirror Pydantic ↔ TS types — drifts silently.
- **Consequences.** Type changes flow from one place; the client is regenerated in CI.

## ADR-0007 — Enrich the fixture via a companion table; expose mapping candidates additively
- **Context.** `signature.load()` consumes the fixture's `weight` column and has committed golden
  outputs; editing the fixture risks scientific drift. Separately, `_maps()` does
  `drop_duplicates("id")` (first-wins) and per-row mapping picks the first candidate — 176 Ensembl
  ids in the released ID map resolve to >1 human symbol, so ambiguity is silently resolved.
- **Decision.** (a) Add provenance as a *companion* CSV (`*.provenance.csv`) and leave
  `pah_muscle_malenfant2015.csv` byte-identical; verified `python analysis.py` output is unchanged.
  (b) Add `candidates_for()` and `mapping_audit()` that expose ALL candidates and conserve input-row
  counts, WITHOUT modifying `map_genes()`/`load()` (so the scoring path and golden outputs are
  untouched). The audit is what the mapping-confirmation gate will consume.
- **Reasons.** Preserve the regression baseline and paper provenance simultaneously; make ambiguity
  visible for researcher confirmation without changing the science.
- **Alternatives rejected.** Editing the fixture in place (drift risk); changing `map_genes` to
  return candidates (would alter the golden path).
- **Consequences.** 12 new tests (Table-2 validation + mapping candidates + sensitivity partition);
  full suite 122 passed / 34 skipped; no drift. New uploads still require interactive confirmation;
  the built-in fixture ships a committed mapping manifest (all 25 rows unambiguous).

## ADR-0008 — Service wraps compute(); content-addressed input; run_id excludes presentation options
- **Context.** `run.compute()` returns a dict with DataFrames, a Store, a Timer, and NaN/NumPy
  values; `run.provenance()` reads `opts['_path']` which only `main_run` set (KeyError otherwise).
- **Decision.** `apps/api` exposes `run_analysis(request) -> AnalysisResponse` that: (a) materializes
  the signature to a temp file whose name embeds its sha256 (content-addressed; no client file
  paths); (b) sets `R["opts"]["_path"]` to that file before calling `provenance()`; (c) extracts only
  JSON-safe scalars (NaN/inf -> None, numpy -> python); (d) computes a deterministic `run_id` over
  signature bytes + analysis params + ordered family column IDs + mapping/store/code versions + seed
  + schema version, EXCLUDING presentation-only fields (include_nonsignificant, display filters).
- **Reasons.** A stable, safe, serializable contract with real provenance; the run_id changes only
  when the science changes; the path fix also closes a path-traversal hole.
- **Alternatives rejected.** Shelling out to the `mprobe` CLI (opaque, unserializable, slow);
  returning the raw compute() dict (non-serializable, leaks internals).
- **Consequences.** 8 service tests; family size recorded as 52 (engine `core_cols`), not a hardcoded
  53; headline q=0.0584 computed. The service is the single seam React (Gate 5) and Marimo call.

## ADR-0009 — Catalog & availability derived from store.columns; explicit capability matrix
- **Context.** The React app shipped an in-memory catalog that (a) listed proteomics for rat SKM-VL
  (the store has SKM-VL RNA only; SKM-GN has RNA+PROT+PHOSPHO) and (b) offered mouse/other source
  species the mapper cannot resolve. Layers were labeled a flat "implemented/planned".
- **Decision.** Generate the catalog and availability from `store.load_columns()`. Publish a
  capability matrix with four independent flags per layer (engine/api/react/demo_validated) so an
  engine-supported but UI-unexposed layer (phospho, metab) is not mislabeled globally "planned".
  Gate the source species to human+rat; mouse/other return an explicit unsupported status. Add
  lineage fields (source_feature_id, n_collapsed, aggregation_method="max_abs_stat" when >1,
  mapping_decision_id) to every feature so point -> evidence -> source is testable.
- **Reasons.** The store is the single source of truth for what exists; hardcoded catalogs drift and
  misreport availability (an integrity risk). Lineage makes the traceability claim verifiable.
- **Alternatives rejected.** Keeping the React constants (already wrong); a single "supported"
  boolean per layer (hides the engine-vs-UI distinction).
- **Consequences.** 8 catalog tests + 2 lineage tests. "No compatible data" is a valid HTTP-200
  outcome, not an error. React (Gate 5) consumes this catalog instead of its constants.

## ADR-0010 — Gate 5: FastAPI transport, live React view added (not a rewrite), gallery left as-is
- **Context.** The React app is a 5-view design handoff with mock data. The requirement is to connect
  it to real results without destroying the design artifact, and to expose the engine over HTTP.
- **Decision.** (a) Add a FastAPI app whose OpenAPI/Pydantic schema is the transport source of truth;
  endpoints wrap the service. Scientific empty states are HTTP 200; malformed input is 422. (b) Add a
  typed `src/api/client.ts` and a new "Live results (API)" view (`LiveDashboard.tsx`) that renders
  real numbers from one run object, deriving visualization mode from RETURNED layers via a new
  `getVisualizationModeFromReturnedLayers`. Keep the existing design views. (c) Report the male rat
  SKM-GN protein 8wk result as q=0.0584 (not significant); never surface the 16-column q=0.0413 as a
  headline. (d) Do NOT regenerate the committed teammate gallery in this work; document the stale
  `5e2cde3` provenance and a regeneration procedure instead (regen belongs in its own reviewed commit).
- **Reasons.** Preserves the teammate's design work and gallery (Collaboration/integrity), adds a real
  end-to-end path, and keeps the honest headline. Deriving viz mode from returned layers fixes the
  fake-multi-omic bug.
- **Alternatives rejected.** Rewriting the mock dashboard in place (destroys the design view and risks
  large churn); regenerating the gallery now (large teammate-owned diff mixed into feature work).
- **Consequences.** 31 api tests + 18 web tests; live round-trip verified over HTTP; CI runs both.
  Point -> evidence -> source enforced at the API level (browser e2e is a future add).

## ADR-0011 — Metabolomics wired in as a separate read-only case-study analysis type (Option-1 adapter)
- **Context.** origin/main merged a standalone `MoTrPAC Hackathon/Metabolomics/` module (ST000763
  SSc-PAH plasma, mapped to MoTrPAC). The engine also has a METAB layer (`metab.py`). These are
  complementary — case study vs reusable scorer — not competing sources of truth (they compute
  different comparisons on different data). The task was to surface metabolomics in the app.
- **Decision.** (a) Expose the module as a SEPARATE analysis type `metabolomics_case_study` via an
  Option-1 adapter (`metab_casestudy.py`): read the module's committed outputs, normalize to the
  service evidence/provenance schema, no stats and no R at request time; the module stays the source
  of truth. (b) Give it its own top-level UI view, never a tab in the muscle-protein dashboard, so
  the SSc-PAH plasma cohort is never blended with the Malenfant muscle signature or rat exercise data.
  (c) Add a convergence cross-check (`convergence.py`) reconciling the module's MoTrPAC export against
  the engine METAB store on shared EE-CON cells. (d) Keep the engine's live directed-signature
  metabolite scoring (Option 2) as a documented follow-on (capability matrix api/react False).
- **Reasons.** Fastest honest path to the metabolomics story without R-at-request-time; preserves the
  module's authorship and its cautionary conclusion; the cross-check is an integrity win (two
  independent pipelines agree on the MoTrPAC side).
- **Alternatives rejected.** Route metabolite signatures through the engine live (Option 2) first —
  answers a different question and needs more work; reimplement the module in mprobe — forks a second
  source of truth.
- **Consequences.** 7 metabolomics tests (adapter counts match 06_context_summary.json; convergence
  identical 1156/1156 cells, max diff 0.0; endpoints; export bundle; catalog analysis type). api suite
  38 passed; web 18 passed + build OK; engine golden byte-identical (adapter-only, no engine change).
  Honest conclusion (setting/scleroderma, not PAH; null acute response) carried verbatim to UI+export.


## ADR-0012 — Generalized query_core wired in as a third analysis type; RUN_LOCAL added
- **Context.** The merged `generalized/` module (`query_core`) is CLI-only; the frontend did not
  expose it. A user opening the app could not reach the "run any signature vs MoTrPAC" flow, and a
  port clash (an unrelated "aegis" app on :8000) made the app hard to find locally.
- **Decision.** (a) Add `generalized_query.py`, an adapter that imports `query_core.engine.run_query`
  in-process (adds `generalized/` to sys.path), runs a bundled signature into a temp dir, and returns
  a JSON-safe summary + rank correlations + provenance. Endpoints `/api/generalized/signatures` and
  `/api/generalized/query/{id}`; catalog gains the `generalized_query` analysis type. (b) Add a
  `GeneralizedQuery.tsx` frontend view (8th nav item). (c) Add `RUN_LOCAL.md` documenting the two
  processes and exact ports, and that :8000 is an unrelated app.
- **Reasons.** Puts the generalized pipeline behind the UI (no longer CLI-only) using the same
  adapter pattern as mprobe and the metabolomics case study; `query_core` stays the source of truth.
  RUN_LOCAL removes the port confusion.
- **Alternatives rejected.** Reimplement query_core inside the service (forks a source of truth);
  leave it CLI-only (the user explicitly wanted it in the frontend).
- **Consequences.** 5 generalized tests (adapter reproduces the legacy six blood Spearman values
  through the service; endpoints; catalog). api suite 43 passed; web 18 + build OK (8 views). Engine
  golden unchanged (adapter-only).

## ADR-0013 (R1) — Documentation-first: author WORKFLOWS.md before any feature source change
- **Context.** The Track 2 "Omic Discordance Explained" deliverables (discordance catalog +
  predictive model, PTM-parent audit) already exist as CLI-only pilots with committed demo outputs,
  and the FastAPI layer already computes `layer_discordance` that no UI renders. Before exposing any
  of this, the team needs one traced map of every analysis path so success criteria are analyzed
  against reality, not assumptions (Requirement 1, documentation-first).
- **Decision.** Author a single `WORKFLOWS.md` at the repository root **before** any
  feature-implementation source file is touched. Trace all five analysis types
  (`directed_signature`, `metabolomics_case_study`, `generalized_query`, discordance catalog +
  predictive model, PTM-parent audit) through the same seven stages in order (entry point →
  module/function → data read → transform → output fields → provenance/run_id → UI view), writing
  "not applicable" for any missing stage rather than omitting it. Include exactly one mermaid diagram
  + one numbered narrative per type; a request→response contract table covering existing endpoints
  plus the three planned discordance endpoints and the planned generalized upload endpoint; the
  point→evidence→source chain as an ordered sequence referencing `test_point_to_evidence_to_source_chain`
  by name; and CLI-only labels on the discordance catalog, predictive model, and PTM-parent paths.
  The five mermaid diagrams are authored for reuse in R6 (mock-view content migration) without
  modification to their diagram source.
- **Reasons.** Documentation-first makes the build phase deliberate and reviewable; a single traced
  map is the shared reference every subsequent feature (R2–R7) points back to; authoring the mermaid
  diagrams now lets R6 migrate mock-view content by reusing them verbatim (no diagram drift).
- **Alternatives rejected.** (a) Trace workflows inline as each feature is built — fragments the map
  and lets the build outrun understanding. (b) Author diagrams per-feature and reconcile later —
  guarantees diagram drift and duplicate sources; R6 could not reuse them verbatim.
- **Consequences.** `WORKFLOWS.md` is documentation only and altered **no runtime source file**
  (R1 AC7). Structure verified by grep (5 mermaid diagrams, five analysis-type sections, the
  seven-stage labels, "not applicable" present, CLI-only labels, the test referenced by name, all
  four planned endpoints in the contract table). The `parents[3]` path-coupling used by the read-only
  adapters was confirmed to resolve to the repo root.

## ADR-0014 (R2) — Discordance deliverables exposed via a read-only Option-1 adapter (GET-only, event-class domain gate)
- **Context.** The Track 2 headline deliverables — the concordant/discordant event catalog, the
  out-of-fold ridge predictive model of discordance, and the PTM phosphosite-vs-parent divergence
  audit — already exist as CLI-only pilots under `MoTrPAC Hackathon/generalized/discordance/`
  (`build_discordance.py`, `audit_ptm_parent.py`) with committed demo outputs
  (`demo_muscle_ee/{catalog.csv,model_metrics.csv,model_predictions.csv,run_summary.json}`,
  `demo_ptm_parent_ee/{ptm_parent_summary.csv,ptm_parent_candidates.csv}`). No API surfaced them and
  no UI rendered them. A wrinkle: `run_summary.json.catalog_classes` carries SIX labels — four true
  classification classes plus two coverage states — while R2 AC6 fixes exactly FOUR supported
  classification classes (`supported_concordant`, `supported_opposite`, `rna_response_protein_equivalent`,
  `indeterminate`). `supported_opposite` has zero rows in this committed demo, so a Counter over the
  file never emits it.
- **Decision.** Expose the three deliverables through a SEPARATE read-only adapter
  `apps/api/motrpac_probe_service/discordance_casestudy.py`, mirroring the metabolomics Option-1
  pattern (ADR-0011): read the committed demo outputs, normalize to JSON-safe dataclasses, compute NO
  statistics at request time; the standalone discordance modules remain the sole source of truth.
  Three GET-ONLY endpoints in `app.py` — `/api/discordance/catalog`, `/api/discordance/model`,
  `/api/discordance/ptm-parent` — each 404 when its `available()` guard is False, with no `_RUNS`
  cache entry and no write/create/update/delete verb honored (write verbs return 405). Reconcile the
  six-vs-four wrinkle by surfacing `classification_counts` as the four supported classes
  ZERO-FILLED (`supported_opposite` = 0 always present) and `coverage_counts` as the two coverage
  states (`no_protein_measurement`, `no_rna_measurement`) kept SEPARATE. Gate the catalog on the full
  event-class domain (`KNOWN_EVENT_CLASSES` = 4 classes ∪ 2 coverage states): `catalog_available()`
  is True exactly when every row's `event_class` is in that set; one missing/malformed/out-of-set
  label makes the catalog unavailable and `build_catalog_summary()` raises `FileNotFoundError`
  (R2 AC8, Property 4). Carry the PTM occupancy caveat and the honest weak-prediction note verbatim
  in the payloads. Add an additive `ANALYSIS_TYPES["discordance_catalog"]` entry in `catalog.py`
  (label/cohort/engine/api/react flags + "serves committed demo; live recompute is a follow-on").
- **Reasons.** The Option-1 adapter is the proven, fastest honest path to surface committed pilot
  outputs without recomputing statistics in the transport layer (source-of-truth boundary, R11.1);
  GET-only + `available()` guards make the read-only contract enforceable and testable; zero-filling
  the four classes keeps `supported_opposite` first-class even at zero rows (R2 AC6) while the
  separate `coverage_counts` avoids inflating the classification family with coverage bookkeeping;
  the event-class domain gate turns a malformed catalog into an honest "unavailable" (404) rather
  than a silently truncated or fabricated result (R2 AC3, AC8).
- **Alternatives rejected.** (a) Recompute the catalog/model live in the service — forks a second
  source of truth and violates R11.1. (b) Fold the two coverage states into the four classification
  classes — overstates the classification family and breaks the AC6 four-class contract. (c) Drop or
  silently skip malformed catalog rows — hides integrity problems; the gate must fail closed. (d) A
  single combined `/api/discordance` endpoint — couples three independently-available resources into
  one availability decision.
- **Consequences.** Committed counts served and asserted by tests (verified by running them, R13.4):
  `classification_counts` = {supported_concordant 1, supported_opposite 0,
  rna_response_protein_equivalent 285, indeterminate 5642}; `coverage_counts` =
  {no_protein_measurement 9228, no_rna_measurement 255}; `model_metrics` = exactly 6 rows over
  {all_mapped, rna_responsive} × {zero, rna_only, temporal}; PTM-parent summary = 3 timepoint rows;
  `candidate_count` = 749. `Path(discordance_casestudy.py).resolve().parents[3]` resolves to the repo
  root (the `parents[3]` path-coupling the adapter relies on). `apps/api/tests/test_discordance.py`:
  14 passed of 15 (see VERIFICATION R2 checkpoint) — the one non-pass is the task-2.4 Property-2
  idempotence/byte-stability test (task 2.4 still in progress), which passes in isolation and fails
  only under full-file test interaction (a 2.4 test-isolation matter), not a defect in the R2 adapter
  or endpoints. Engine and existing api behavior unchanged (adapter + additive endpoints only).

## ADR-0015 (R3) — Discordance is the three-stage main view; four-class + two-coverage reconciliation rendered honestly, predictions framed weak
- **Context.** The Track 2 headline deliverables are now served read-only by the R2 adapter
  (ADR-0014) but no UI rendered them, and `run_analysis` already returns `layer_discordance` that
  zero views showed. R3 makes discordance the app's main view, organized as a three-stage flow
  (query builder → data visualization → discordance interpretation). Two wrinkles carry over from
  the adapter: (a) the catalog reconciles into FOUR classification classes
  (`supported_concordant`, `supported_opposite`, `rna_response_protein_equivalent`, `indeterminate`)
  plus TWO separate coverage states (`no_protein_measurement`, `no_rna_measurement`), and
  `supported_opposite` is zero in the committed demo; (b) the out-of-fold ridge model is a weak
  prediction (near-zero R²), a result and not a success claim, so the UI must say so in visible text.
  A per-feature RNA/protein scatter is desirable but the adapter does not deliver per-feature catalog
  rows by default, so the scatter must be omitted (never faked) when rows are absent.
- **Decision.** Build `apps/web/src/components/Discordance.tsx` as a single three-stage view rendered
  in order — Stage 1 (query builder over the committed demo option(s), with UI text stating live
  recompute is out of scope), Stage 2 (visualization: class-count chart from `classification_counts`,
  model-metrics comparison chart from `model.metrics`, PTM-parent timing from `ptm.summary_rows`, and
  an optional same-time RNA/protein scatter rendered ONLY if per-feature catalog rows are delivered),
  Stage 3 (interpretation: all four classification classes INCLUDING any zero class, the two coverage
  states kept SEPARATE and labeled "coverage state (not a class)", a model-metrics table of R²/MAE/RMSE
  for the zero/rna_only/temporal variants, the PTM-parent divergence summary, the PTM occupancy caveat
  as visible text, and a visible weak-prediction statement alongside the metrics). Render a SINGLE
  exclusive empty/error branch (`data-testid="discordance-empty-or-error"`): the catalog fetch failing
  OR a delivered-but-empty catalog (`total_events === 0`, or all four classes summing to zero) renders
  ONLY the empty-state and NONE of the class counts, model metrics, or PTM-parent summary (R3 AC9).
  Model and PTM are best-effort resources fetched after the anchor catalog; each may be independently
  unavailable and the view still renders the stages it can. Separately, render `layer_discordance` in
  `LiveDashboard.tsx` as a "Layer discordance" table (R3 AC8). The client computes NO statistics; it
  renders the read-only adapter payloads verbatim.
- **Reasons.** The three-stage flow gives the analyst a guided path from query to charts to classified
  results (R3 AC2) while keeping the source-of-truth boundary intact (R11.1 — no client stats).
  Rendering all four classes zero-filled keeps `supported_opposite` first-class at zero rows (R3 AC3)
  and the separate coverage cards avoid inflating the classification family with coverage bookkeeping,
  matching the adapter's reconciliation (ADR-0014). The visible weak-prediction statement is the honest
  framing R12.4 requires — the near-zero R² is reported as a result, not upgraded to a success claim.
  Omitting the scatter when per-feature rows are absent avoids fabricating data points (integrity). A
  single exclusive empty/error return lets Property 14 (task 4.5) and the component tests (task 4.6)
  assert crisp mutual exclusivity: on success everything renders; on failure/empty nothing does.
- **Alternatives rejected.** (a) Fold the two coverage states into the four classification classes —
  overstates the classification family and breaks the AC3 four-class contract. (b) Compute a
  live catalog/model in the client — forks a source of truth and violates R11.1. (c) Fabricate or
  interpolate scatter points when per-feature rows are absent — an integrity violation; the scatter is
  omitted instead. (d) Render metrics without a weak-prediction statement — would overstate a near-zero
  R² model (violates R12.4). (e) Separate empty and error into two independent branches that could both
  partially render — would let a stray count/metric leak onto the empty path (breaks R3 AC9).
- **Consequences.** `Discordance.tsx` (three-stage view + exclusive empty/error branch) and the
  `LiveDashboard.tsx` `layer_discordance` table are in place. Verified by running the web tests
  (R13.4): `npm run test:domain` → 18 passed / 0 failed; `npm run test:web` (vitest) → 3 files, 7
  tests passed — `Discordance.property.test.tsx` (Property 14: result completeness + empty-state
  exclusivity, 2 tests, ≥100 fast-check iterations each) and `Discordance.test.tsx` (task 4.6
  component tests, 4 tests: three-stage QUERY→VISUALIZE→INTERPRET order, PTM-parent divergence summary,
  PTM occupancy caveat, weak-prediction statement) are BOTH present and passing;
  `ResultsTable.property.test.tsx` (Property 8, 1 test) also passes. R3 AC1 (Discordance as the
  DEFAULT view) is NOT delivered here — the default-view/nav wiring is task 8.1 (R6); the Discordance
  component is ready and its default-nav wiring is tracked with 8.1. Engine and `apps/api` behavior
  unchanged (web-only changes).

## ADR-0016 (R4) — Live query builder rewired to `/api/catalog`; mandatory mapping-preview gate; single-source signature; species-derived study
- **Context.** `QueryBuilder.tsx` drove the analyst flow from an in-memory catalog and
  `LiveDashboard.tsx` issued a single hardcoded request (`example_name="pah_muscle_malenfant2015"`,
  `target_species="rat"`, `fdr_threshold=0.05`). This meant: (a) the selectors could offer target
  combinations the store does not actually have (the same class of drift ADR-0009 fixed on the
  backend), (b) an uploaded/pasted signature had no first-class path, and (c) mappings were resolved
  silently — the 176 Ensembl ids that map to >1 human symbol (ADR-0007) were never surfaced for
  confirmation before a comparison ran. Three boundaries carry over from the backend: the run_id
  already includes `fdr_threshold` and excludes `include_nonsignificant` (ADR-0008, verified by
  Property 5, so it is NOT re-implemented here); the FDR domain is `0 < t <= 1`, but the API request
  validation expresses the bound as `gt=0` / `lt=1` — so a UI preset of exactly `1.0` is inclusive on
  the UI side while the API's `lt=1` is exclusive at the endpoint; and study design is derived from
  species alone (rat → chronic, human → acute; Property 11), so there is no acute/chronic input.
- **Decision.** Rewire `QueryBuilder.tsx` to source catalog and comparison data from the live API
  (`/api/catalog`), reusing the already-valid domain logic (`getStudyContext`,
  `reconcileQueryForSpecies`, `getVisualizationModeFromReturnedLayers`, `layerCapabilities`) from
  `domain/analysis.ts` rather than re-deriving it. (a) The signature source is EXACTLY ONE of built-in
  example, pasted text, or uploaded CSV: example/paste → `signature_rows` (or `example_name`), upload →
  `signature_csv_text`. (b) An invalid upload (empty, over max size, or unparseable CSV) is rejected,
  the previously selected source is retained unchanged, and an error is shown — the source never
  silently switches. (c) Target selectors (species, omics, tissue, sex, timepoint) are catalog-driven:
  a TOTAL `/api/catalog` failure disables ALL selectors and shows a catalog-unavailable error, while a
  PARTIAL response (some dimensions present, one dimension missing/empty) disables ONLY the empty
  dimension's selector and leaves the populated ones enabled. (d) FDR is chosen from presets over
  `0 < t <= 1`, with the UI-inclusive-upper vs API `gt=0`/`lt=1` boundary recorded as a note (not a
  silent mismatch). (e) A MANDATORY mapping-preview gate blocks the `/api/comparisons` call until the
  previewed mappings are explicitly confirmed: a preview failure blocks the comparison, retains the
  current selections, and shows an error; every multi-candidate term is surfaced with all its
  candidates and NONE is auto-selected. (f) Species-derived study is shown read-only (no acute/chronic
  control); both contrasts remain first-class and the phrase "healthy gene set" never appears. Rewire
  `LiveDashboard.tsx` to remove its hardcoded request and drive results from the confirmed
  QueryBuilder run contract, preserving the honest q=0.0584 headline. `service._run_id` is left as-is
  (already correct per ADR-0008 / Property 5). Also fixed the pre-existing `QueryBuilder.tsx` tsc error
  (`replaceAll` → `split`/`join`) so the file type-checks.
- **Reasons.** The store/catalog is the single source of truth for what exists (ADR-0009); driving the
  UI from `/api/catalog` removes the last in-memory-catalog drift surface. A single, explicit signature
  source prevents ambiguous "which input won" states. The mandatory preview gate makes mapping
  ambiguity a researcher decision instead of a silent first-wins pick (ADR-0007) — the integrity crux
  of R4. Reusing the existing domain functions keeps the client statistics-free (R11.1) and avoids a
  second copy of the reconciliation logic. Recording the FDR boundary note keeps the UI-inclusive vs
  API-exclusive difference honest rather than hidden.
- **Alternatives rejected.** (a) Keep the in-memory catalog in the UI — reintroduces the ADR-0009
  drift the backend already fixed. (b) Allow more than one signature source to be simultaneously
  "active" and let the API disambiguate — ambiguous and untestable; exactly-one is the AC1 contract.
  (c) Auto-select the first mapping candidate to skip the gate — the silent first-wins behavior R4
  exists to remove; violates AC10/AC12. (d) Fail the comparison silently (or clear selections) on a
  preview failure — loses the analyst's work and hides the error; AC11 requires retain + error. (e)
  Disable ALL selectors on any partial catalog gap — needlessly blocks usable dimensions; AC7 requires
  disabling only the empty one. (f) Add an acute/chronic input — study design is species-derived
  (Property 11); an input would let the UI contradict the backend. (g) Re-implement run_id handling in
  the client — it is a backend concern already verified by Property 5.
- **Consequences.** `QueryBuilder.tsx` is catalog-driven with the exactly-one-source rule, the
  partial/total selector-degradation behavior, the FDR presets + boundary note, and the mandatory
  mapping-preview confirmation gate; `LiveDashboard.tsx` is driven by the confirmed run contract with
  the honest headline preserved. Verified by running the tests (R13.4): `npm run test:web` → 6 files /
  14 tests passed (incl. Property 7 source exclusivity + field mapping, Property 9 catalog-driven
  selector enablement, Property 10 mapping-preview confirmation gate); `npm run test:domain` → 26
  tests passed / 0 failed (extended with the FDR-boundary + reconciliation paths); `test_service.py`
  → 13 passed (incl. Property 5 run_id determinism confirming FDR-in / include_nonsignificant-out, and
  Property 11 species-derived study). `service._run_id` unchanged (already correct per ADR-0008);
  `apps/api` engine behavior unchanged (web-only changes plus the pre-existing tsc fix).

## ADR-0017 (R5) — Generalized upload wired in-process; content-addressed temp file; a shared ResultsTable enforces the never-a-bare-value labeling invariant where a genuine results-row model exists
- **Context.** The generalized `query_core` pipeline was reachable only as bundled built-in
  signatures (ADR-0012); an analyst could not upload or paste their own signature and run it against
  the MoTrPAC reference. Separately, each view rendered its own bespoke results table, so the
  labeling invariants R5 requires — raw nominal p (`camera_p`) shown immediately adjacent to the BH q
  (`camera_fdr`) with each header stating that the q is adjusted over the frozen 52-column
  multiplicity family and the p is nominal/unadjusted; a distinctly labeled `species` column and a
  distinctly labeled `dataset` column; species + dataset + contrast on EVERY rendered value; the
  ortholog-link relation on cross-species results; and the honest headline being the `q = 0.0584`
  not-significant result with the 16-column `q = 0.0413` only ever a labeled sensitivity analysis
  (R12.2, R12.3) — were not consistently enforceable. `ColumnResult` already carries `camera_p`,
  `camera_fdr`, `species`, and `dataset`, so no backend change is needed to render the set-level
  table; the gap was a shared component and an upload path, not new statistics. A wrinkle: not every
  view produces the same row shape — rank correlations (generalized), class/coverage counts
  (discordance), and metabolite hits (metabolomics) are NOT `camera_p`/`camera_fdr` +
  species/dataset/contrast rows, so a single table cannot be force-fit onto all of them without
  fabricating fields.
- **Decision.** (a) Add `run_uploaded(signature_csv_text, *, tissue, reference_contrast_category,
  reference="default")` to `generalized_query.py`, mirroring `run_bundled`: guard on `available()`;
  validate the uploaded signature against the query_core schema via `_validate_query_core_schema`
  (which REUSES query_core's own loader — the engine remains the sole source of truth, no
  reimplemented schema), returning `(ok, failing_requirement)`; on failure raise `SchemaError(which)`
  WITHOUT running the query and leaving any prior signature unchanged; on success write a
  CONTENT-ADDRESSED temp file whose name embeds the sha256 of the signature bytes (no client-supplied
  server path, matching the ADR-0008 pattern), run `_run_query`, normalize identically to
  `run_bundled`, and clean up the temp file. Compute NO statistics. (b) Add `POST
  /api/generalized/query` in `app.py` accepting `{ signature_csv_text, tissue,
  reference_contrast_category, reference }`, returning the normalized `run_bundled`-shaped payload on
  success and 422 NAMING the unmet schema requirement on non-conforming input (translated from
  `SchemaError`). (c) Add a typed `generalizedQueryUpload()` + `GeneralizedQueryUploadInput` to
  `client.ts`; add upload + paste controls to `GeneralizedQuery.tsx` where a pasted signature that
  cannot be parsed is rejected with a parse-failure indication WITHOUT clearing the pasted text.
  (d) Create a shared `ResultsTable` component whose row-model TYPE requires `species`, `dataset`,
  and `contrast` on every row (so the component can never render a bare value), places `camera_p`
  immediately adjacent to `camera_fdr` with the two labeled headers, renders the ortholog-link
  relation for cross-species rows, and — where a headline is surfaced — surfaces the `q = 0.0584`
  not-significant result, using the 16-column `q = 0.0413` only as an explicitly labeled sensitivity
  analysis. (e) ADOPTION SCOPE (recorded honestly): `ResultsTable` is adopted in the Live set-level
  `ColumnResult` table, where a genuine `camera_p`/`camera_fdr` + species/dataset/contrast results-row
  model exists. The Generalized (rank correlations), Discordance (class/coverage counts), and
  Metabolomics (metabolite hits) views KEEP their bespoke tables because their row shapes are
  different and lack that row model; forcing the shared component onto them would require fabricating
  species/dataset/contrast/`camera_*` fields those results do not have. Cohorts remain separate in
  presentation (muscle protein / plasma metabolomics / blood distinct; R12.5).
- **Reasons.** In-process `run_uploaded` reuses the proven adapter pattern (ADR-0012) so query_core
  stays the single source of truth and the upload path performs no statistics (R11.1); the
  content-addressed temp file keeps inputs safe/serializable and never trusts a client path (ADR-0008)
  and also closes a path-traversal surface. Type-enforcing species/dataset/contrast on the row model
  is the strongest possible form of R5 AC7/AC8 — the component literally cannot compile a bare value —
  rather than relying on a runtime check that could be forgotten. Adopting the shared table where a
  real results-row model exists (the Live set-level table) makes the labeling invariant hold there by
  construction, while keeping the other three views' bespoke tables is the HONEST choice: their rows
  are genuinely a different shape, and inventing `camera_p`/`camera_fdr`/species/dataset/contrast for
  a rank correlation or a metabolite-hit count would fabricate scientific fields (an integrity
  violation) purely to satisfy a uniformity that the data does not support.
- **Alternatives rejected.** (a) Reimplement the query_core schema check in the adapter — forks a
  second source of truth; `_validate_query_core_schema` reuses the engine's own loader instead. (b)
  Accept a client-supplied server file path for the uploaded signature — a path-traversal risk; the
  content-addressed temp file is written from the posted bytes. (c) Run the query first and validate
  after — R5 AC1/AC2 require rejecting non-conforming input WITHOUT running and returning the unmet
  requirement; validate-then-run is the contract. (d) Clear the analyst's pasted text on a parse
  failure — loses their work; R5 AC4 requires retaining it. (e) A runtime-only "every row has
  species/dataset/contrast" assertion instead of a type-enforced row model — weaker; the type makes a
  bare value unrepresentable (R5 AC8). (f) Force the single `ResultsTable` onto ALL four views
  (Generalized/Discordance/Metabolomics included) — would fabricate `camera_p`/`camera_fdr` +
  species/dataset/contrast fields for rank correlations, class counts, and metabolite hits that do not
  have them; a fake-uniformity integrity violation. (g) Merge cohorts into one table for visual
  consistency — violates R12.5 cohort separation.
- **Consequences.** `generalized_query.run_uploaded()` + `SchemaError` +
  `_validate_query_core_schema()` and the `POST /api/generalized/query` endpoint (422 naming the
  failing schema requirement) are in place; `client.ts` has `generalizedQueryUpload()`;
  `GeneralizedQuery.tsx` has upload + paste controls with parse-failure text preserved; the shared
  `ResultsTable` (camera_p adjacent to camera_fdr with labeled headers, type-required
  species/dataset/contrast, ortholog link for cross-species, q=0.0584 headline / 16-col only as
  sensitivity) is adopted in the Live set-level table and the other three views correctly kept their
  bespoke tables (rank correlations / class counts / metabolite hits are different shapes, not forced;
  cohorts kept separate per R12.5). Verified by running the tests (R13.4): from `hackathon/tool`
  (venv `hackathon/tool/.venv`), `.venv/bin/python -m pytest ../../apps/api/tests/test_generalized.py
  -q` → `7 passed` (including Property 6 as both a unit test `test_property6_upload_schema_enforcement`
  and an endpoint test `test_property6_endpoint_rejects_nonconforming`), and the full
  `.venv/bin/python -m pytest ../../apps/api/tests -q` → `64 passed`. From `apps/web`, `npm run
  test:web` (vitest) → `7 files / 19 tests passed`, including `ResultsTable.property.test.tsx`
  (Property 8: every row shows species/dataset/contrast, no row dropped, cross-species rows show the
  ortholog link — ≥100 fast-check iterations) and `GeneralizedQuery.test.tsx` (5 tests); `npm run
  build` OK (JS 279.69 kB / gzip 83.07 kB). Engine golden unchanged (backend change is the additive
  `run_uploaded` + endpoint only; the shared table is web-only). Scope note (honest): R5 AC5/AC6/AC9
  labeling invariants are enforced by the shared `ResultsTable` in the Live set-level `ColumnResult`
  table where a genuine results-row model exists; the Generalized/Discordance/Metabolomics views keep
  their bespoke tables because their row shapes (rank correlations, class counts, metabolite hits)
  lack the `camera_p`/`camera_fdr` + species/dataset/contrast row model, and forcing the shared
  component onto them would fabricate fields.

## ADR-0018 (R6) — Retire the five mock design views; pure redirect resolver maps any retired/unknown identifier to a live view; mock content migrated to DESIGN_PROVENANCE.md with per-view Figma-Make attribution
- **Context.** The React app shipped as a Figma-Make design handoff whose navigation included FIVE
  mock design views — architecture/01, technical/02, workflow/03, dashboard/04, and judge-slide/08 —
  rendering static hardcoded content (including the fabricated gene identifiers PPARGC1A, SOD2, and
  COL1A1) presented alongside the real analysis views. This is the integrity risk R6 exists to close:
  a user could be shown illustrative mock content indistinguishable from computed results. Three
  wrinkles carry into the change: (a) the mock views hold genuine design/narrative content authored in
  Figma Make, so deleting them outright would lose the design story and its provenance; (b) a user (or
  a stale bookmark/session) can navigate DIRECTLY to a retired identifier, so simply removing the nav
  entries is not enough — a direct hit must not render a retired view; and (c) R1 already authored the
  five mermaid diagrams for reuse (ADR-0013), so the migrated documentation should reference them by
  name without touching their diagram source.
- **Decision.** (a) Remove the five mock views and all their helpers from `apps/web/src/App.tsx`; the
  `View` type now names ONLY the five live views. Discordance is the default view on initial load, and
  the nav presents EXACTLY five entries — Discordance (default), Live, Generalized, Metabolomics, and a
  single combined About/Methods page (`AboutMethods`). (b) Add a PURE redirect resolver — `resolveView`
  (a live identifier resolves to itself; any retired or otherwise unknown identifier resolves to the
  last-used live view when that is itself live, otherwise to the default Discordance view) and
  `initialView` (runs before first render, reading the requested identifier and the last-used live view
  from sessionStorage key `esx.lastView`). A retired view is NEVER returned, so callers render the
  result directly and retired content is never mounted. The last-used live view is persisted to
  `esx.lastView` on navigation. (c) Migrate the content of the removed views 01–04+08 into
  `DESIGN_PROVENANCE.md` (task 8.2), preserving a per-view Figma-Make design-origin attribution line on
  each migrated view (plus a global Figma-Make provenance section) and referencing the R1
  `WORKFLOWS.md` diagrams by name — NO diagram-source change. (d) `AboutMethods` links to the migrated
  docs (`DESIGN_PROVENANCE.md` and `WORKFLOWS.md`). (e) Remove the banned mock identifiers PPARGC1A /
  SOD2 / COL1A1 from every navigable view (the `QueryBuilder` placeholder now uses neutral GENE1/GENE2).
- **Reasons.** Removing the mock views from navigation is the direct satisfaction of R6 AC1/AC6 — a
  user is never shown static mock content presented as real results (the integrity crux of R6). A PURE
  `resolveView`/`initialView` resolver is the strongest form of R6 AC2: because it can never return a
  retired identifier, a direct navigation to a retired route cannot mount retired content by
  construction (not merely by a runtime guard that could be bypassed), and it degrades gracefully to
  the last-used live view or Discordance. Migrating the content to `DESIGN_PROVENANCE.md` with per-view
  Figma-Make attribution preserves the design story and its authorship (R6 AC4/AC5, Collaboration/
  integrity) instead of destroying it, and reusing the R1 diagrams by name avoids diagram drift
  (ADR-0013). Dropping PPARGC1A/SOD2/COL1A1 from navigable views satisfies R6 AC3 — the fabricated
  identifiers never appear where a user could read them as data.
- **Alternatives rejected.** (a) Delete the mock views outright with no migration — loses the
  Figma-Make design story and its provenance (violates R6 AC4/AC5). (b) Only hide the nav entries but
  leave the views routable — a direct navigation to a retired identifier would still render mock
  content (violates R6 AC2). (c) A runtime guard inside each view component that redirects on mount —
  weaker than a pure resolver: the retired component would still be imported/mounted before redirect,
  and the guard could be forgotten on a new view; the pure resolver makes a retired render
  unrepresentable. (d) Rewrite the mock diagrams fresh in the migrated doc — guarantees drift from the
  R1 diagram source; referencing them by name keeps a single source. (e) Keep PPARGC1A/SOD2/COL1A1 as
  harmless placeholders — they are fabricated gene identifiers that R6 AC3 bans from navigable views.
- **Consequences.** `App.tsx` renders only the five live views with Discordance as the default; the
  nav has exactly five entries; `resolveView`/`initialView` (pure, sessionStorage `esx.lastView`) map
  any retired/unknown identifier to the last-used live view or Discordance and never return a retired
  view; `AboutMethods` links to `DESIGN_PROVENANCE.md` and `WORKFLOWS.md`; the mock content of views
  01–04+08 lives in `DESIGN_PROVENANCE.md` with per-view Figma-Make attribution and R1-diagram
  references (no diagram-source change); PPARGC1A/SOD2/COL1A1 no longer appear in any navigable view
  (QueryBuilder placeholder → GENE1/GENE2). Verified by running the tests (R13.4): from `apps/web`,
  `npm run test:web` (`vitest run`) → 9 files / 34 tests passed, including `App.property.test.tsx`
  (Property 15: `resolveView` never returns a retired view and always returns a live view — retired/
  unknown → last-used live view or Discordance) and `App.test.tsx` (11 tests: nav presents exactly the
  five live entries, Discordance is the default view, and the aggregate banned-identifier grep asserts
  no navigable component source contains PPARGC1A/SOD2/COL1A1); `npm run build` → OK (JS 279.69 kB /
  gzip 83.07 kB; CSS 42.40 kB / gzip 9.33 kB; 25 modules transformed). The banned-identifier grep over
  navigable source (`grep -rnE 'PPARGC1A|SOD2|COL1A1' apps/web/src --include=*.tsx --include=*.ts`
  excluding test files) is CLEAN — the only remaining matches are in test files (`App.test.tsx`
  assertions + `QueryBuilder.source.property.test.tsx` CSV fixtures), which assert/exercise their
  absence. Engine and `apps/api` behavior unchanged (web-only changes + docs migration).

## ADR-0019 (R7) — Root README rewritten to anchor on Track 2; discordance catalog + model as the headline deliverable; one navigable link per key doc
- **Context.** The root `README.md` opened with the raw challenge prompt but then trailed off into an
  empty section skeleton (`Project Snapshot`, `Research Question`, `Workflow`, `Setup`, `Inputs &
  Outputs`, `Methods`, `Validation`, `Reuse`) — headings with no body content — and it did not lead
  with the honest scope or point a new reader to the traceability docs. R7 makes the README the honest
  front door: name the Track 2 "Omic Discordance Explained" challenge in the first section, present the
  discordance catalog + predictive model as the headline deliverable before any other, state the
  cross-cohort (not PAH-efficacy) framing plus the cohorts-separate and PTM-occupancy caveats, link
  once to each of the five key docs, and leave no empty heading. This is documentation only (R7 is a
  DOCS-ONLY task) and touches no runtime source.
- **Decision.** Rewrite `README.md` with: (a) an H1 title immediately followed by the H2
  `## Omic Discordance Explained` section whose body (before the second H2 `## Deliverables`) names the
  Track 2 challenge as the anchoring topic and states the honest framing (cross-cohort association, not
  a PAH treatment-efficacy claim; cohorts analyzed separately and not merged; PTM occupancy caveat
  applies) (R7 AC1, AC3, AC4). (b) The discordance catalog + predictive model presented as the headline
  deliverable in that first section AND as item 1 of a `## Deliverables` list, before the PTM-parent
  audit and the directed/generalized analyses (R7 AC2). (c) The honest directed headline recorded as
  `q = 0.0584` over the frozen 52-column BH family, labeled not significant, with the 16-column
  `q = 0.0413` only ever a labeled sensitivity analysis (R12 discipline). (d) A `## Documentation map`
  with EXACTLY one navigable Markdown link to each of `WORKFLOWS.md`, `DECISIONS.md`, `VERIFICATION.md`,
  `METABOLOMICS.md`, and `RUN_LOCAL.md` (R7 AC5); the Setup section references those docs in prose
  rather than re-linking them, to preserve the exactly-one-link invariant. (e) No heading left without
  body text or subordinate content (R7 AC6).
- **Reasons.** Anchoring the README on Track 2 and leading with the discordance catalog + model makes a
  new reader immediately understand the honest scope and the headline deliverable (the R7 user story).
  Exactly-one-link-per-doc keeps the documentation map unambiguous and satisfies AC5 as written (a
  second link to the same target would violate the "one navigable link to each" count). Replacing the
  empty section skeleton with real body content closes AC6 (no heading followed by empty content) and
  removes the misleading impression that those sections were authored. Carrying the honest framing and
  the `q = 0.0584` headline into the README keeps the front door consistent with R12 and with the rest
  of the docs.
- **Alternatives rejected.** (a) Keep the empty section skeleton and only prepend framing — leaves
  headings with empty content (violates AC6). (b) Link the key docs from multiple sections (e.g. both a
  documentation map and a setup section) — would push `VERIFICATION.md`/`RUN_LOCAL.md` to two links
  each and violate AC5's exactly-one count (this was caught and corrected during verification). (c) Lead
  with the directed muscle-protein result — buries the Track 2 headline deliverable and contradicts AC2.
- **Consequences.** `README.md` is documentation only; no runtime source changed. Structure and
  link-resolvability verified by grep BEFORE claiming (R13.4; see the VERIFICATION R7 checkpoint):
  heading order confirms `## Omic Discordance Explained` is the first section before the second heading
  `## Deliverables` (AC1); the headline discordance catalog + model appears in that first section and as
  Deliverables item 1 (AC2); the three framing statements are present (AC3, AC4); each of the five doc
  links appears exactly once and every `*.md` link target resolves to an existing file (AC5); and a
  subordinate-heading-aware scan finds no heading followed by empty content (AC6). Engine, `apps/api`,
  and `apps/web` behavior unchanged (docs only).

## ADR-0020 (R8) — One-command bootstrap with a named-target fail-fast contract; pnpm as the single package manager; the `parents[3]` path-coupling documented
- **Context.** A fresh clone needed a reliable, single-command path to a working state, but the setup
  was implicit: a reader had to know to create a venv, editable-install the engine BEFORE the service
  (the service depends on the engine as an editable package, not a pinned dependency), fetch the
  mprobe store, and somehow confirm the wiring was intact — with no timing expectation and no honest
  failure signal if a step silently half-completed. Two related integrity hazards sat alongside this.
  (a) The imported React app (ADR-0003) carried a `pnpm-lock.yaml`, but a stray
  `apps/web/package-lock.json` from an incidental `npm install` had crept in (flagged as an open item
  back in the Gate 1 checkpoint), so the repo had TWO competing lockfiles and CI/setup instructions
  could drift between package managers — the exact ambiguity R8 AC6 exists to close. (b) The
  read-only API adapters (ADR-0011, ADR-0014) locate their committed demo files by climbing
  `Path(__file__).resolve().parents[3]` to the repo root; this coupling is invisible and silently
  breaks if an adapter is moved shallower or deeper, yet nothing documented the required
  directory-depth relationship so a reader could verify a given file location.
- **Decision.** (a) Author `scripts/bootstrap.sh` as the single fresh-clone command
  (`bash scripts/bootstrap.sh`, runnable from anywhere — it resolves the repo root from its own
  location via `SCRIPT_DIR`/`REPO_ROOT`). It runs, in order: create/reuse `hackathon/tool/.venv` →
  editable-install the engine (`pip install -e "hackathon/tool[test]"`) → editable-install the service
  (`pip install -e "apps/api[api,test]"`) → `mprobe store fetch` → a named smoke check that imports
  `motrpac_probe`, confirms the store is present (`store.CONTRASTS.exists()`), imports
  `motrpac_probe_service`, and confirms `catalog.build_catalog()` builds a non-empty catalog. It is
  ENGINE-FIRST by design: the engine editable install precedes the service install because the service
  depends on the engine editable package, so an engine failure must stop before the service is
  attempted. Fail-fast is `set -euo pipefail` plus explicit per-step guards that emit a message NAMING
  the failed target/check: a failed editable install names `hackathon/tool` or `apps/api` and runs no
  later step (R8 AC2); a failed/unavailable `mprobe store fetch` emits a store-fetch error (R8 AC3);
  the smoke check exits zero iff every named check passes and otherwise exits non-zero naming the
  failed check (R8 AC4). A `BOOTSTRAP_DEADLINE_S=600` timing guard records elapsed time and fails if
  the run exceeds 600 s (R8 AC1). (b) Standardize on pnpm as the SINGLE declared package manager:
  remove the stray `apps/web/package-lock.json`, add `"packageManager": "pnpm@10.34.3"` to
  `apps/web/package.json`, and use `pnpm/action-setup@v4` + `pnpm install --frozen-lockfile` + `pnpm
  test` / `pnpm run build` (honoring the committed `apps/web/pnpm-lock.yaml`) in every CI job, leaving
  NO standalone `npm` invocation in either CI or the setup docs (R8 AC6). (c) Document the bootstrap in
  `BOOTSTRAP.md`, linked by exactly one resolvable link from BOTH `README.md` and `RUN_LOCAL.md`
  (R8 AC5), and document the `parents[3]` path-coupling constraint there: `parents[3]` means "go up
  four directory levels from the adapter file", which resolves to the repo root ONLY when an adapter
  sits exactly three directories below the root (`apps/api/motrpac_probe_service/<name>.py`), with a
  by-hand depth rule (`depth == 3`) and a one-line `python3` check a reader can run to verify a given
  file location satisfies the constraint (R8 AC7).
- **Reasons.** A single command with a NAMED-target fail-fast contract is the honest form of R8:
  every failure mode produces a specific, actionable message and a non-zero exit rather than a
  half-installed environment that looks fine (R8 AC2–AC4); engine-first ordering matches the real
  dependency direction so the error points at the true first failure. The smoke check verifies the
  wiring END-TO-END (engine import → store present → service import → catalog builds) so "bootstrap
  succeeded" means the app can actually build its catalog, not merely that pip exited zero. Collapsing
  to a single package manager removes the dual-lockfile drift the Gate 1 checkpoint flagged as an open
  item: pnpm was already the toolchain's preference (`.mise.toml` pinned pnpm 10.34.3) and the imported
  app shipped `pnpm-lock.yaml`, so pnpm is the lockfile that reflects the real dependency tree — the
  `packageManager` field + `--frozen-lockfile` make that choice enforceable and reproducible.
  Documenting `parents[3]` turns an invisible, silently-breakable coupling into something a reader can
  verify by counting directory segments or running one command, so a future file move that would break
  the adapters is caught by understanding rather than by a mysterious FileNotFoundError.
- **Alternatives rejected.** (a) A README prose list of manual steps instead of a script — every reader
  re-derives the ordering and there is no timing/fail-fast contract or end-to-end smoke check; drift is
  guaranteed. (b) Install the service before the engine (or in parallel) — inverts the real dependency
  and makes an engine failure surface as a confusing service-install error. (c) Keep both lockfiles and
  let CI pick one — the dual-lockfile ambiguity R8 AC6 exists to remove; installs could resolve
  differently under npm vs pnpm. (d) Standardize on npm instead of pnpm — would discard the committed
  `pnpm-lock.yaml` and the `.mise.toml` pnpm pin, forcing a lockfile regeneration and losing the exact
  resolved tree the app was authored against. (e) A smoke check that only imports the packages without
  building the catalog — would pass on a broken store/wiring; building the catalog is what proves the
  app is actually usable. (f) Leave `parents[3]` undocumented and rely on the test alone — the test
  (`test_parents3_path_coupling`) catches a regression but does not TEACH a reader the depth rule they
  need to place a new adapter correctly.
- **Consequences.** `scripts/bootstrap.sh` is the one-command fresh-clone entry (engine-first editable
  installs → `mprobe store fetch` → named end-to-end smoke check, `set -euo pipefail` + per-step guards,
  600 s deadline guard); pnpm is the single declared package manager (`"packageManager":
  "pnpm@10.34.3"` in `apps/web/package.json`; the stray `apps/web/package-lock.json` is removed; every
  CI job uses `pnpm/action-setup@v4` + `pnpm install --frozen-lockfile`); `BOOTSTRAP.md` documents the
  invocation and the `parents[3]` depth rule and is linked once each from `README.md` and `RUN_LOCAL.md`.
  Verified by running the checks BEFORE claiming (R13.4; see the VERIFICATION R8 checkpoint):
  `bash -n scripts/bootstrap.sh` → SYNTAX_OK; the embedded smoke-check heredoc extracted and
  `python3 -m py_compile`'d → PY_COMPILE_OK (34 lines); a word-boundary grep for standalone `npm
  (install|test|run|ci|exec)` over `.github/` and the setup docs (`BOOTSTRAP.md`, `RUN_LOCAL.md`,
  `README.md`) → ZERO matches (all web steps are `pnpm …`), and `apps/web/package-lock.json` is ABSENT
  while `apps/web/pnpm-lock.yaml` is present; `BOOTSTRAP.md` is linked exactly once from each of
  `README.md` and `RUN_LOCAL.md` (both links resolve to the existing file) and its own single `*.md`
  link (`RUN_LOCAL.md`) resolves. This ADR/checkpoint is documentation only (task 11.4, DOCS ONLY): it
  records the R8 work completed in tasks 11.1 (bootstrap script), 11.2 (pnpm standardization), and 11.3
  (`BOOTSTRAP.md` docs) and altered NO runtime source or the bootstrap script. Engine, `apps/api`, and
  `apps/web` behavior unchanged.

## ADR-0021 (R9) — Automated verification wired into CI: in-process httpx smoke with exact-value assertions, real-boot Playwright e2e with the browser point→evidence→source chain, a full-battery merge gate (exactly 122/34), and a golden-hash guard on `analysis.py`
- **Context.** Through R8 the app was verified piecemeal — the engine `pytest` battery, the
  `apps/api` tests, and the web `vitest`/build each ran, and the point→evidence→source chain was
  asserted server-side by `apps/api/tests/test_api.py::test_point_to_evidence_to_source_chain` — but
  nothing (a) exercised every documented HTTP endpoint end-to-end and asserted the exact committed
  headline/discordance/metabolomics values in one pass, (b) proved the live views render real data in
  a real browser against a real API, or (c) turned "the battery is green" into an enforceable,
  countable merge gate. R9 closes those gaps and wires the result into CI so a regression is caught on
  every pull request. Several boundaries carry in: the honest headline is `q = 0.0584 ± 0.0001`
  labeled not-significant (R12.2, ADR-0004/ADR-0010); the discordance per-class integer counts are the
  committed `supported_concordant`=1 / `supported_opposite`=0 / `rna_response_protein_equivalent`=285
  / `indeterminate`=5642 (ADR-0014); the metabolomics story is a NULL result (0% blood
  exercise-sensitive hits, Fisher p 0.609; ADR-0011); the engine battery is EXACTLY 122 passed / 34
  skipped over an allowlist of known store/deck/R skips (recorded since Gate 2); and `analysis.py` is
  the frozen scientific core whose golden reference is `hackathon/tool/analysis.py.sha256`.
- **Decision.** (a) httpx SMOKE (task 12.1, `apps/api/tests/test_smoke_httpx.py`): a SINGLE-execution
  (not 100-iteration) check that boots the FastAPI app IN-PROCESS via `httpx.ASGITransport` (no live
  server needed, matching the rest of `apps/api/tests`), requests every documented endpoint —
  health, catalog, availability, mappings/preview, signatures/validate, the directed
  `/api/comparisons` + its run subresources, the three metabolomics endpoints, the generalized list +
  bundled + UPLOAD endpoints, and the three discordance endpoints — and asserts each returns success.
  It asserts the EXACT values: headline `q == 0.0584 ± 0.0001` AND `significant_at_threshold is False`;
  the four discordance per-class integer counts (with an `isinstance(int)` type check); and the
  metabolomics null result (`hits_fraction == 0.0`, `fisher_p == 0.609`). On any non-success OR value
  mismatch it raises `SmokeFailure` NAMING the offending endpoint/value (R9 AC3); when the store or a
  committed demo is absent it SKIPS with an explicit reason rather than fabricating a pass. (b)
  Playwright E2E (task 12.2, `apps/web/e2e/live-views.spec.ts` + `apps/web/playwright.config.ts`):
  `playwright.config.ts` boots BOTH real processes via its `webServer` array — the API (`uvicorn
  motrpac_probe_service.app:app` from `apps/api`) then Vite (`pnpm run dev`, `VITE_API_BASE_URL`
  pointed at the API) — per the RUN_LOCAL.md port model (8765/8443, env-overridable), pnpm per R8 AC6.
  The spec asserts, per live view (Discordance, Live, Generalized, Metabolomics), that REAL data
  renders (the frozen discordance class + coverage counts; a real Live comparison + the
  `layer_discordance` table; generalized rank-correlation rows; the metabolomics hits table +
  conclusion) and that the known error/empty markers are ABSENT; and it drives the browser-level
  point→evidence→source chain — select a rendered Live point, follow the page's own run-scoped
  export/report links and `/features` evidence (run-scoped `evidence_id`, shown `column_id`), and the
  store-backed `source_feature_id` + the report's verbatim honest-framing source text — as the
  browser counterpart to `test_point_to_evidence_to_source_chain`. (c) CI WIRING (task 12.3,
  `.github/workflows/ci.yml`): the httpx smoke runs inside the `python` job with a MUST-NOT-SKIP gate
  (after `mprobe store fetch`, assert the smoke reported `1 passed` and did not skip / "no tests ran"),
  and a dedicated `e2e` job installs the engine+service (with `[api]` for uvicorn), fetches the store,
  installs web deps (pnpm `--frozen-lockfile`) + the chromium browser, and runs `pnpm run test:e2e`,
  on every push to `integration/**`/`main` and every pull request. (d) FULL-BATTERY GATE + GOLDEN-HASH
  GUARD (task 12.4, in `ci.yml`): the engine step parses the `pytest` summary and FAILS unless
  `passed == 122 AND skipped == 34` (with an allowlist that rejects any non-allowlisted skip), so any
  deviation blocks the battery and the merge (R9 AC7/AC8); the `apps/api` step and the web `pnpm test`
  + `pnpm run build` steps must report zero failures; and a `Golden-hash guard` step runs `shasum -a
  256 -c analysis.py.sha256` on every push/PR, failing the moment `analysis.py` is not byte-identical
  to the golden reference (R9 AC9).
- **Reasons.** An in-process httpx pass is the fastest honest way to prove the whole HTTP surface is
  wired AND that the exact committed values still flow through it — one run, no server, naming the
  offending endpoint/value on any drift (R9 AC1–AC3) — while skipping (never faking) when the store or
  a demo is genuinely absent. A real-boot Playwright battery is the only way to prove the four live
  views render real data as a USER sees it and that the point→evidence→source chain holds in the
  browser, not merely at the API layer (R9 AC4/AC5). Wiring both into CI as PR-triggered jobs (the
  smoke inside `python` with a must-not-skip gate, e2e as its own job) makes the guarantee automatic on
  every change (R9 AC6). Encoding the battery as a COUNTABLE gate — exactly 122/34, `apps/api` zero
  failures, web tests+build zero failures — turns "green" into an enforceable, non-negotiable merge
  condition (R9 AC7/AC8) instead of a human eyeball. The golden-hash guard makes any change to the
  frozen scientific core `analysis.py` fail closed unless the reference is deliberately re-blessed
  (R9 AC9), protecting the honest headline from silent drift.
- **Alternatives rejected.** (a) Run the smoke against a live uvicorn server in CI — slower and
  flakier than the in-process ASGITransport; the in-process client hits the same app object. (b) Make
  the smoke a 100-iteration property test — R9 AC1/AC2 want wiring + exact values in ONE pass, not a
  distribution; the property tests already live elsewhere. (c) Let the smoke silently pass when the
  store/demo is absent — would hide a broken CI environment; it SKIPS with a reason and the CI gate
  additionally asserts it did NOT skip. (d) Assert the live views with component tests only (no
  browser) — cannot prove the real API↔web boot or the browser-level chain (R9 AC4/AC5). (e) A gate
  that accepts "≥122 passed" or ignores the skip count — would let a newly-skipped (silently broken)
  test or an unexpected extra pass slip through; the battery is complete ONLY at EXACTLY 122/34
  (R9 AC7). (f) Guard `analysis.py` with a test that recomputes numbers instead of a content hash —
  weaker; the byte-identical hash check is the strongest statement that the frozen core is unchanged
  (R9 AC9). (g) Fold e2e into the `python`/`web` jobs — couples two independently-installable
  toolchains (engine+store+uvicorn AND node+pnpm+chromium) into one job; a dedicated `e2e` job is
  cleaner and fails in isolation.
- **Consequences.** `apps/api/tests/test_smoke_httpx.py` (in-process httpx, single pass, exact-value
  asserts, names offender on failure, skips-not-fakes when store/demo absent);
  `apps/web/e2e/live-views.spec.ts` + `apps/web/playwright.config.ts` (real API+web boot, per-view
  real-data assertions, browser point→evidence→source chain); and `.github/workflows/ci.yml` with
  three jobs (`python`, `e2e`, `web`) — the `python` job carrying the golden-hash guard, the engine
  EXACTLY-122/34 gate, the `apps/api` step, and the httpx must-not-skip gate; the `e2e` job booting
  the real pair and running `pnpm run test:e2e`; the `web` job running `pnpm test` + `pnpm run build`.
  Verified BEFORE claiming (R13.4; see the VERIFICATION R9 checkpoint): the smoke test file, the e2e
  spec, the Playwright config, and `analysis.py.sha256` all exist; `shasum -a 256 -c
  hackathon/tool/analysis.py.sha256` → `analysis.py: OK` (hash
  `1c1089723ec7c9d67c4c0f7ab51262717f80d95322c7da8df60187a06cab7372`); the ci.yml gate literally
  asserts `EXPECTED_PASSED, EXPECTED_SKIPPED = 122, 34`; ci.yml parses as valid YAML with jobs
  `python, e2e, web`; and the smoke `1 passed` / e2e `5 passed` results are those recorded by tasks
  12.1/12.2. This ADR/checkpoint is documentation only (task 12.5, DOCS ONLY): it RECORDS the R9 work
  completed in tasks 12.1 (httpx smoke), 12.2 (Playwright e2e + config), 12.3 (CI wiring), and 12.4
  (full-battery gate + golden-hash guard) and altered NO runtime source, CI, test, or the golden-hash
  file. The full LIVE CI battery (engine 122/34 + apps/api + web + e2e all executing in CI) is a
  live-environment run and is NOT re-executed by this docs-only task; the checks above are static
  (`shasum -c`, YAML parse, ci.yml gate/job inspection, file-existence) or cited from the prior tasks'
  recorded runs. Engine, `apps/api`, and `apps/web` behavior unchanged (docs only).
## ADR-0022 (R10) — Merge-to-main plan is DOCUMENTATION ONLY: `MERGE_PLAN.md` records the pre-merge preconditions, the block rules, the no-rebase-over-`t3code` constraint, the exact `gh pr create` (base `main` / head `integration/exercise-signature-explorer`), and a three-section PR template — and executes nothing
- **Context.** With R1–R9 complete (docs, discordance API, the Discordance default view + three-stage
  flow, the live query builder, generalized upload + shared results table, mock-view retirement, the
  Track-2 README, fresh-clone bootstrap, and automated verification wired into CI), the only remaining
  step is the merge of `integration/exercise-signature-explorer` into `main`. R10 requires that this
  merge be **planned as a document, not performed** — nothing in this feature may trigger, invoke, or
  execute any `git`, `gh`, or shell command (R10 AC1). Several boundaries carry in: the verification
  battery is the R9 CI gate (engine EXACTLY 122 passed / 34 skipped, `apps/api` zero failures, web
  `pnpm test` + `pnpm run build` zero failures, httpx smoke exact values, Playwright e2e, golden-hash
  guard on `analysis.py`; ADR-0021); the integration branch was cut off the tool branch and
  `origin/t3code/build-motrpac-probe-tool` is shared upstream history whose commits must not be
  replayed over; and the honest headline `q = 0.0584` not-significant + cross-cohort-not-PAH framing
  (ADR-0004/ADR-0010/ADR-0011) is what the PR summary must carry.
- **Decision.** (a) MERGE PLAN AS A DOCUMENT (task 13, `MERGE_PLAN.md` at the repo root, alongside the
  other planning docs — `WORKFLOWS.md`, `DECISIONS.md`, `VERIFICATION.md`, `REQUIREMENTS_TRACEABILITY.md`,
  `BOOTSTRAP.md`, `METABOLOMICS.md`): a documentation-only plan that states at the very top, in an
  `> [!IMPORTANT]` callout, that it EXECUTES NOTHING and that authoring/reading it performs no
  version-control action (R10 AC1). (b) PRECONDITIONS (Section 1): the merge is permitted only when ALL
  three hold — every verification gate reports passing; `git status` reports zero modified/staged/
  untracked/deleted files (documented verifier: `git status --porcelain` must produce NO output); and
  the local branch is zero commits ahead of upstream (documented verifier: `git rev-list --count
  @{upstream}..HEAD` must be `0`) — R10 AC2. (c) BLOCK RULES (Section 2): if any gate is non-passing,
  OR the tree has ≥1 uncommitted/untracked file, OR the branch is ≥1 commit ahead of upstream, the
  merge is BLOCKED and the plan requires identifying WHICH precondition failed (a block-condition table
  records what to capture + how to resolve each) — R10 AC3. (d) NO-REBASE CONSTRAINT (Section 3): no
  rebase may replay commits over `origin/t3code/build-motrpac-probe-tool`; the merge is a merge/PR only,
  and no `git rebase` whose base/`--onto` is that branch may be run — R10 AC4. (e) EXACT PR-CREATE
  COMMAND (Section 4), shown as literal text in a fenced code block: `gh pr create --base main --head
  integration/exercise-signature-explorer --title "Finalize Discordance MVP" --body-file PR_BODY.md`
  (R10 AC5). (f) PR TEMPLATE (Section 5) with EXACTLY three labeled sections — `## Summary`, `## What
  was tested`, `## Deferred items` — R10 AC6. Every `git`/`gh` command in the doc lives INSIDE a fenced
  code block for a human to run deliberately later; the doc itself runs none of them.
- **Reasons.** A written plan lets the team see the exact merge procedure, its gating preconditions,
  and the failure/block semantics before anyone touches the branch, without the feature itself taking
  any irreversible or shared-history-affecting action (R10 AC1). Encoding the three preconditions as
  copy-pasteable verifier commands (`git status --porcelain` → empty; `git rev-list --count
  @{upstream}..HEAD` → 0) makes "clean tree" and "not ahead of upstream" objectively checkable by an
  operator rather than a vibe (R10 AC2). Spelling out the block rules with a "record which precondition
  failed" requirement guarantees a blocked merge is diagnosable, not silently retried (R10 AC3). The
  no-rebase-over-`t3code` clause protects shared upstream provenance from being rewritten (R10 AC4).
  Pinning the EXACT base/head on the `gh pr create` command removes any ambiguity about what merges
  into what (R10 AC5), and the fixed three-section PR template makes every PR self-describe its scope,
  its verification evidence, and its deferred work (R10 AC6).
- **Alternatives rejected.** (a) Actually run the merge / open the PR as part of this feature — directly
  violates R10 AC1 (documentation only) and is an irreversible shared-history action; rejected. (b)
  Have this plan *check* the preconditions (run `git status` / `git rev-list` / the battery to "verify"
  readiness) — still executes VCS/shell commands, violating R10 AC1; the plan DOCUMENTS the
  preconditions as text and leaves execution to the operator. (c) Allow a rebase onto
  `origin/t3code/build-motrpac-probe-tool` to linearize history — rewrites shared commits and breaks
  provenance; explicitly forbidden by R10 AC4. (d) Put the merge steps in a runnable script
  (`merge.sh`) — invites accidental execution and blurs the documentation-only boundary; a read-only
  `.md` with commands quarantined in fenced blocks is safer. (e) A free-form PR description instead of
  the fixed template — R10 AC6 requires EXACTLY the three labeled sections; a template guarantees it.
- **Consequences.** `MERGE_PLAN.md` (repo root) with the docs-only `> [!IMPORTANT]` header, Section 1
  (three preconditions + verifier commands), Section 2 (block rules + block-condition table), Section 3
  (no-rebase-over-`t3code`), Section 4 (exact `gh pr create --base main --head
  integration/exercise-signature-explorer`), Section 5 (the three-section PR template), and a closing
  honest scope note. Verified BEFORE claiming (R13.4; see the VERIFICATION R10 checkpoint) using only
  harmless STATIC checks (no `git`/`gh`/VCS command run): `grep -cE '^## (Summary|What was tested|
  Deferred items)$' MERGE_PLAN.md` → `3` (the three PR sections, at lines 132/142/150); `grep -n 'gh pr
  create --base main --head integration/exercise-signature-explorer' MERGE_PLAN.md` → line 115 (exact
  base/head); the no-rebase clause + all four `origin/t3code/build-motrpac-probe-tool` mentions present
  (lines 95–104); the three documented preconditions present (passing gates / `git status` zero files /
  zero commits ahead — lines 26, 41, 52); and a `git`/`gh` fence-containment scan (awk, indentation-
  aware, inline-backtick spans stripped) → all 3 real commands (`git status --porcelain`, `git rev-list
  --count @{upstream}..HEAD`, `gh pr create …`) sit INSIDE fenced code blocks, with ZERO bare
  executable command lines outside fences (the only two non-fence `git`-word hits are prose sentences —
  "run the verification battery, or query any remote" and "does not inspect the git working tree" — not
  commands). This ADR/checkpoint is documentation only (tasks 13 + 13.1, DOCS ONLY): it RECORDS the R10
  merge plan and altered NO runtime source, CI, or test; NO `git`, `gh`, or shell VCS command was run
  in authoring the plan or verifying it. Engine, `apps/api`, and `apps/web` behavior unchanged (docs
  only).


## ADR-0023 (merge) — Integrate `origin/motrpac-explorer`: Explorer as the default live view, LiveDashboard kept as backup, mocks stay retired; saved-vs-live cross-check made tolerant + artifact-aware; engine battery gate 122 → 123
- **Context.** After the finalize-discordance-mvp work (ADR-0013…ADR-0022) was committed and pushed,
  a parallel branch `origin/motrpac-explorer` added a generalized MoTrPAC Explorer — any molecule
  list (genes/proteins/metabolites/pathways) in, every matching MoTrPAC result by omic layer out,
  with the set-level cameraPR / rank-correlation computed on demand against the pre-processed
  `motrpac_probe` store (a later commit, `85f897a`, made bundled examples load committed saved
  results from `apps/web/public/examples/*.json` while custom uploads/pastes still run live). The
  two branches diverged from the same base and overlapped heavily in the web app. Integrating them
  raised a genuine product decision (not a mechanical merge): the explorer branch's `App.tsx`
  restored the mock design views (architecture/technical/workflow/dashboard/slide) that R6
  deliberately RETIRED, rebuilt section 05 as the Explorer, and DELETED `LiveDashboard.tsx`, whereas
  HEAD had a five-live-view model with Discordance as the default and `LiveDashboard.tsx` carrying
  the `layer_discordance` table. The team also pivoted away from "no live demo": the interactive
  live surface is now the intended primary, so Discordance no longer needs to be the default.
- **Decision.** Merge (not rebase — consistent with `MERGE_PLAN.md` R10 AC4; both branches are
  published). Resolutions: (a) **App.tsx** — keep HEAD's R6 navigation model (mock design views stay
  RETIRED and documented in `DESIGN_PROVENANCE.md`; the `resolveView`/`initialView` redirect resolver
  and sessionStorage persistence are preserved) and ADD the Explorer as a live view. The merged app
  has SIX live views; **Explorer is the default landing view**, **Live results (LiveDashboard) is the
  backup live view**, and **Generalized query (bundled example signatures) is a further backup**;
  Discordance is kept as the Track 2 "Omic Discordance Explained" headline showcase over committed
  demo outputs. Nav order: Explorer → Live results (API) → Discordance → Generalized → Metabolomics →
  About/Methods. (b) **LiveDashboard.tsx** — KEPT against the explorer branch's delete (modify/delete
  resolved in favor of HEAD); it is the directed live view plus the `layer_discordance` table.
  (c) **pyproject.toml** — union of both optional-dependency sets (`test` += `hypothesis>=6.100`;
  `api` += `openpyxl>=3.1` for the explorer's `.xlsx` demo-input upload). `openpyxl` was ALSO added
  to the `test` extra so `apps/api[test]` (what the CI `python` job installs) can run the explorer
  tests that read `.xlsx` inputs. (d) **MotrpacExplorer.tsx** — the input-format placeholder gene
  `PPARGC1A` was changed to `ATP5F1A` so the R6 AC3 banned-mock-identifier invariant
  (PPARGC1A/SOD2/COL1A1) holds across the merged navigable surface; `App.test.tsx` extends the
  banned-identifier grep to cover MotrpacExplorer/liveCharts/ui and was updated to the six-view
  Explorer-default model. (e) **`test_saved_example_results_match_the_live_analysis`** (added by the
  explorer branch) was changed from an EXACT byte-equality check to a **tolerance-based, artifact-aware**
  comparison: numeric fields agree within `rtol=1e-9`/`atol=1e-12` (last-ULP float differences arise
  from platform BLAS/summation order — the saved JSONs are generated on one machine, CI/dev on another),
  and an example whose live rerun cannot reproduce a saved-available layer because a required BUILT
  artifact is absent (the gitignored R-built pathway table `apps/api/data/motrpac_camera_pathways.csv.gz`,
  produced by `Rscript apps/api/scripts_build_pathways.R`) is skipped-in-loop with an explicit reason
  rather than failing. (f) **CI full-battery gate** raised from EXACTLY 122 → **123** passed / 34
  skipped: the explorer branch added `hackathon/tool/tests/test_multiple_testing.py` (covering the
  additive `camera_bonferroni` column in `core.py`), a legitimate +1 to the engine suite. The R9
  golden-hash guard on `analysis.py` is UNAFFECTED — the merge changed `core.py`, not `analysis.py`
  (hash still `analysis.py: OK`).
- **Reasons.** Keeping the R6 retirement while adding Explorer honors the finalize-discordance-mvp
  deliverable AND the new generalized tool with nothing lost (Option "integrate both"). Explorer as
  default matches the team's pivot to a live-first primary surface; LiveDashboard and Generalized
  remain as backups so a demo degrades gracefully (live tool → directed live view → bundled examples).
  A merge (not rebase) preserves both authors' history and matches the documented merge discipline.
  The exact-equality saved-vs-live test was not portable (bit-identical floats across machines, plus
  a gitignored R artifact); the tolerant, artifact-aware form preserves its intent (saved results
  reflect the real analysis) while passing deterministically anywhere — mirroring the engine suite's
  existing expected-skip-when-local-data-absent posture. Raising the gate to 123 records the real
  new count rather than masking it.
- **Alternatives rejected.** (a) Take the explorer branch's `App.tsx` wholesale — un-retires the mock
  views and contradicts the R6 deliverable. (b) Let Explorer REPLACE LiveDashboard (accept the delete)
  — drops the `layer_discordance` rendering (part of R3). (c) Rebase either branch — rewrites published
  history; rejected per R10 AC4 / git safety. (d) Regenerate the saved example JSONs locally to fix the
  float mismatch — only moves the last-ULP mismatch to the next machine and still can't reproduce the
  pathways example without the R table. (e) Leave the saved-vs-live test failing / build the R table in
  CI — leaves CI red or adds a heavy, non-portable R dependency to the gate. (f) Add `openpyxl` only to
  the `api` extra — the CI `python` job installs `[test]`, so the explorer xlsx tests would fail there.
- **Consequences.** Two merge commits (`08dec87` bringing in `c164747`, then `4babcbd` bringing in
  `85f897a`). Discordance and Explorer endpoints coexist on the backend (app/catalog/service/client
  auto-merged). Verified by running the battery BEFORE claiming (R13.4): engine **123 passed / 34
  skipped**; `apps/api` **81 passed** (includes the new explorer/story/metab_upload tests, with the
  saved-vs-live cross-check now tolerant/artifact-aware); web **38 passed** + `pnpm run build` OK;
  golden hash `analysis.py: OK`; `ci.yml` gate updated to 123/34 and YAML re-validated
  (`jobs: python, e2e, web`). Note: `apps/api[api]`/`[test]` now require `openpyxl` (a fresh clone /
  bootstrap must install it — already covered by the extras). The explorer's saved example JSONs are
  regenerated by `apps/api/scripts_save_examples.py` (a provenance artifact, analogous to the gallery
  regen procedure).

## ADR-0024 (multiomic-comparison-view) — Generalized multi-omic Compare-All view as a UI-only delta over the merged RNA↔protein pair co-view; props-based shared-filter lift; B2 store-sourced study/dataset labels
- **Context.** The `origin/motrpac-explorer` merge (ADR-0023) shipped, inside
  `apps/web/src/components/MotrpacExplorer.tsx`, a focused two-layer RNA↔protein co-view: a combined
  `pair` HeatGrid and a `pairTime` SmallMultiples that match molecules across transcriptomics and
  proteomics by the five-part key `species|tissue|category|time|sex` (the `pKey`/`pColIndex`/
  `pMolIndex`/`pValueAt`/`paired`/`pairGroups`/`pairCell`/`pairPanels` machinery). The Explorer fetches
  EVERY layer in one `ExplorerResponse` from the single existing `POST /explorer/analyse` endpoint and
  today shows one layer at a time via a tab index. Two gaps remained: (a) there was no way to see every
  available omic's existing view at once under one shared biological context — each `LayerView`/
  `PathwayView` owned its OWN filter state, so panels could drift apart; and (b) `ExplorerColumn`
  carried no explicit study/dataset provenance, so any per-panel study label would have to be invented
  in React. This feature is a deliberately **delta-only** addition on top of the merged pair co-view: it
  builds ON that reality rather than re-specifying it. `explorer.py` already reads `c.dataset` off store
  columns (values `rat_train` / `human_acute`) and `analyse()` already returns every layer in one
  payload — so neither a store-schema change nor an API request-shape change is needed.
- **Decision.** Ship two capabilities, both additive and read-only w.r.t. the statistical engine.
  - **(A′) UI-only Compare-All comparison grid + shared filter bar.** Add a mode toggle
    (`type Mode = "tabs" | "compare"`, default `tabs`) in the results header next to the existing
    `RuleSwitch`. Tabs Mode is unchanged (exactly one layer visible). Compare-All Mode renders every
    available layer's EXISTING view simultaneously in a responsive Comparison Grid — a sibling modifier
    class `.lr-cards.xp-grid` in `LiveDashboard.css` that switches from single-column to equal `1fr`
    tracks and reflows to fewer columns at the existing `@media (max-width: 860px)` breakpoint (even
    split = pure responsive CSS; no splitter/drag-to-resize). The grid renders a list of **units** where
    a unit is one grid cell: when BOTH transcriptomics and proteomics are available they collapse into
    ONE fused **RNA↔Protein Pair Unit** — a single `LayerView` given its `partner`, reusing the merged
    pair views — occupying one cell (Requirement 3, Reading A); phosphoproteomics is its own unit
    (`LayerView`, no partner); pathways is its own unit (`PathwayView`); a single available RNA-only /
    protein-only layer is its own unit; and an **unavailable layer contributes NO unit and reserves no
    grid space**. The even split is `units.length` equal tracks (the fused pair counting as one). No new
    chart type, no new statistic, no new classification is introduced.
  - **Props-based (NOT React context) shared-filter lift.** A single `Filters`
    (`{ species, tissue, contrast, sex, time, layer }`) is lifted OUT of each `LayerView`/`PathwayView`
    (which each held their own `useState<Filters>`) into a single `useState<Filters | null>` in the
    parent `MotrpacExplorer`, passed down as **props** (`filters` + `onFiltersChange`) — not via a
    React context/provider. The seeding (`defaults(...)`, kept a pure module function) and the change
    resolution are lifted into the parent's `onFiltersChange(next)` handler and applied **once,
    centrally**: the `pick("species", v)` rebuild (`defaults(cols.filter(c => c.species === v), …)`),
    the `tissue`-change `time: ALL` reset, and the human ⇒ `sex: ALL` rule are copied **verbatim** from
    the pre-refactor `LayerView` so behavior is preserved. Shared-bar option values are derived from the
    columns present in the current `ExplorerResponse` (reusing the existing `uniq(...)`/`timesOf(...)`
    derivations), never a static list. The pair machinery is kept exactly as shipped — the same
    `partner` expression (`l.input === layer.input` sibling) feeds both the Tabs-mode pair view and the
    Compare-All pair unit; there is **no second molecule-pairing implementation**.
  - **(B, variant B2) Backend store-sourced study/dataset labels.** `explorer.py` adds two fields to
    each `ExplorerColumn` built in `_layer()` and `_pathway_layer()`: `dataset` (from the store column's
    `c.dataset`; the human-only pathway path via a tiny `_dataset_for_species("human") → "human_acute"`
    helper matching the store's own convention) and `study_label` (from a new `_study_label(dataset)`
    helper). `_study_label` maps the store dataset code to a display phrase (`rat_train` → "endurance
    training", `human_acute` → "acute exercise") and appends the authoritative version read from
    `store.load_provenance()` under keys `MotrpacRatTraining6moData` (→ 2.0.0) and
    `MotrpacHumanPreSuspensionAnalysis` (→ 2.0.8) — deliberately the **PreSuspensionAnalysis** key, NOT
    the sibling `MotrpacRatTraining6mo` key. The provenance read is wrapped in an `@lru_cache`
    `_provenance()` that swallows load failures and returns `{}`; when the version is missing or the
    dataset is unknown, `study_label` falls back to the `dataset` code itself — never omitted, never a
    fabricated study name (R4.7). The `ExplorerColumn` TypeScript type in `client.ts` declares
    `dataset: string; study_label: string;`. The label is surfaced via a `provenanceLabel(c)` helper
    (dot-joined `species · study_label · tissue_label · contrast · sex · time`, human sex → `all`) in
    every panel's `Card` chips (each Compare-All unit shows its OWN provenance) and added to EVERY CSV
    row-builder (the `pair` card CSV, the `molecules` card CSV, the set-level `comparisons` CSV, and the
    `PathwayView` CSV, which also gains `species`). The saved example JSONs under
    `apps/web/public/examples/*.json` are regenerated by `apps/api/scripts_save_examples.py` (no script
    change) so every `ExplorerColumn` carries the two new label fields; the existing page footnote
    (published summary statistics; association not treatment effects; hardcoded citation versions) is
    preserved as-is.
- **Reasons.** Building A′ on the existing pair co-view means side-by-side display adds NO computation —
  `ExplorerResponse` already contains every layer — so the feature stays honestly UI-only (no API
  request-shape change, no new chart type, no new statistic/classification). Props (not context) suit a
  shallow, self-contained tree (`MotrpacExplorer → LayerView/PathwayView`, two levels) with few, known
  consumers: props keep data flow explicit and trivially testable (render `LayerView` with a fixed
  `filters` prop and assert its selection), whereas context would add a provider and indirection with no
  reuse benefit at this size. Copying the species/tissue/human-sex resolution verbatim into one central
  resolver preserves the shipped behavior while making one bar drive every panel. B2 keeps the Store the
  single source of truth for the study identity: the authoritative version comes from
  `store.load_provenance()`, never a React or Python literal, which is the whole reason B2 was chosen
  over deriving labels in the frontend. Fusing RNA+protein into one unit (Reading A) matches the biology
  the pair view already encodes and reuses that machinery rather than duplicating it.
- **Alternatives rejected.** (a) Take the upstream/other-branch approach wholesale — rejected; it would
  not honor the delta-only scope and risked re-introducing retired behavior. (b) React **context** for
  the shared filters — rejected; unnecessary indirection/provider for a two-level tree with few known
  consumers, and less directly testable than props. (c) Derive the study/dataset labels in **React** —
  rejected; the Store is the source of truth (R9.2), so labels are read from `store.load_provenance()`
  in the backend and never hardcoded in the client. (d) A **second** molecule-pairing implementation for
  the Compare-All pair unit — rejected; the fused unit REUSES the existing `partner` pair machinery
  (`pKey`/`pColIndex`/`pMolIndex`/`pValueAt`/`paired`/`pairGroups`/`pairCell`/`pairPanels`, the
  `pair`/`pairTime` cards) unchanged, so there is one pairing path, not two. (e) A general **N-way,
  user-selectable arbitrary-pair fusion** (pick any subset of layers and which two are jointly matched) —
  **deferred**, cited in the requirements Future Scope; this feature fixes the fused unit to RNA↔protein
  and shows other omics as independent panels, leaving the generalized selection UI + matching path to a
  future feature. (f) Sourcing the rat version from the sibling `MotrpacRatTraining6mo` key — rejected in
  favor of `MotrpacRatTraining6moData` (→ 2.0.0).
- **Consequences.** The RNA↔protein pair view is **preserved (generalized, not duplicated)**: Tabs Mode
  on the transcriptomics/proteomics tab renders it exactly as before, and Compare-All Mode renders it as
  the single fused pair unit — same five-part key, same values, same statistics. The feature **computes
  no statistics** and **does not change the API request shape** (`ExplorerResponse`/`ExplorerLayer`/
  `api.explorerAnalyse` unchanged); the Discordance, Live, Generalized, and Metabolomics views are
  untouched. Verified by running the battery BEFORE claiming (R11.4; see the VERIFICATION
  "multiomic-comparison-view" checkpoint): `apps/api` **81 passed** (2 deprecation warnings) — includes
  the explorer suite, with the tolerant, artifact-aware saved-vs-live cross-check
  `test_saved_example_results_match_the_live_analysis` passing after the saved JSONs were regenerated
  (task 7.1); web `pnpm run build` OK (dist emitted) and `pnpm test` **38 passed** (9 files). The
  `motrpac_probe` engine suite is UNTOUCHED by this feature (no engine files changed) and was therefore
  not re-run as part of this feature's battery — it remains at its last recorded 123 passed / 34 skipped
  from the ADR-0023 explorer-merge checkpoint. Environment: `openpyxl` (required for the explorer xlsx
  tests) is installed; the R-built pathway table is absent in this checkout, so `blood_pathway_6`'s
  pathway layer is skip-on-missing-artifact in the cross-check. The optional (`*`) property/component
  tests (tasks 4.4, 5.3, 6.4, 8.1) are not implemented (MVP); they are tracked as PLANNED.

## ADR-0025 (submission readiness) — README rewritten Explorer-first to the hackathon's eight-section scaffold; WORKFLOWS traces the Explorer; a JUDGING_CRITERIA map is added
- **Context.** After the `origin/motrpac-explorer` merge (ADR-0023/ADR-0024) the shipped app
  became **Explorer-default with six navigable views**, but the submission-facing docs still
  described the deprecated **discordance-default, five-view** world: the root `README.md` led with
  the discordance catalog as the headline and named the Discordance view as the default; `WORKFLOWS.md`
  traced only the five discordance-era analysis paths and did not mention the Explorer (the new
  headline); there was no README-level architecture/workflow diagram, no AI-usage disclosure, and no
  document mapping the project to the hackathon's judging criteria. The hackathon slides specify an
  **eight-section README scaffold** (Project snapshot, Research question, Workflow, Setup, Inputs &
  outputs, Methods & provenance, Validation & limitations, Reuse & continuation) and **seven judging
  criteria** (scientific novelty, impact, content, collaboration, methods & approach, documentation &
  integrity, technical complexity), with the explicit caveat that technical complexity alone is not
  rewarded. The user confirmed the pivot: the **Explorer is the single headline** (the paper→demo
  trace), and this work is **documentation only** with a merge to `main` as the next task.
- **Decision.** (a) Rewrite the root `README.md` to the eight-section scaffold, Explorer-first: lead
  with the Explorer (any molecule list in → every matching MoTrPAC result out), state it is the default
  view among six, keep the honest-framing block verbatim, and make every section **point to the
  authoritative doc** (BOOTSTRAP/RUN_LOCAL/SCHEMA/WORKFLOWS/DECISIONS/VERIFICATION/METABOLOMICS/HANDOFF/
  LICENSE) rather than duplicating it. Add two inline mermaid diagrams (a Data→Method→Result→User
  workflow and a system-design/architecture diagram showing the source-of-truth boundary), an AI-usage
  disclosure (Codex/Kiro/Claude/ChatGPT contributed across prompting + implementation; all numbers
  originate in the engines, enforced by the golden-hash guard and the test battery), contributor names
  (roles marked to-be-confirmed, since neither `main`'s README nor git history records explicit roles),
  citations, and a screenshot placeholder (figures are user-provided, out of scope here). (b) Add a
  new `## 0. explorer` section to `WORKFLOWS.md` tracing the Explorer through the same seven stages as
  the other paths with one mermaid diagram, grounded in `explorer.py` / `app.py` / `client.ts`; update
  the intro list to "six analysis types" (Explorer as default/headline) and add the three
  `/api/explorer/*` endpoints to the contract table; the existing five diagrams are unchanged. (c) Add
  `JUDGING_CRITERIA.md` mapping the seven criteria to concrete repo/presentation evidence plus the
  minimum-defensible-submission checklist with honest status (presentation-only items and the
  screenshot/roles gaps marked as gaps, not asserted). (d) Reconcile stale phrasing: fix `MERGE_PLAN.md`
  (engine gate 122→123; PR summary now Explorer-first) and confirm `README.md` carries no
  discordance-default/five-view/122 phrasing.
- **Reasons.** The front door and the interpretability reference must match the shipped reality or a
  judge lands on the Explorer while the README promises Discordance. Pointing to existing docs (rather
  than duplicating) keeps a single source of truth and avoids drift, matching the repo's
  documentation-integrity posture. Mapping to the exact scaffold + criteria makes the submission legible
  to judges and frames technical complexity per the slide caveat (coherence + reproducibility over raw
  complexity). Marking roles/screenshots as explicit gaps rather than inventing them preserves the
  verify-before-claiming discipline (R13.4). Historical VERIFICATION/ADR-0021 checkpoints that record
  `122/34` are **left intact** as accurate point-in-time logs (the 122→123 change is already documented
  in ADR-0023); only the forward-looking gate in `MERGE_PLAN.md` was updated to 123.
- **Alternatives rejected.** (a) Keep discordance-first framing — contradicts the shipped default and
  the user's confirmed pivot. (b) Restore Discordance as the default in `App.tsx` — a runtime change,
  out of the docs-only scope and against the user's deprecation decision. (c) Duplicate setup/method
  content into the README — creates drift; rejected in favor of links. (d) Invent contributor roles or
  fabricate an Explorer screenshot — violates integrity; left as marked gaps for the team. (e) Rewrite
  the historical `122/34` checkpoints to `123/34` — would falsify the reproduction log; rejected.
- **Consequences.** Docs-only change set: `README.md` (full rewrite), `WORKFLOWS.md` (new Explorer
  section + contract rows + intro), `JUDGING_CRITERIA.md` (new), `MERGE_PLAN.md` (123 + Explorer-first
  PR summary), `DECISIONS.md` (this ADR), `VERIFICATION.md` (checkpoint below). No runtime source, CI,
  test, or engine file changed; the golden-hash guard on `analysis.py` is unaffected. Two team actions
  remain before submission: insert Explorer/Discordance screenshots (README §7) and confirm per-person
  contributor roles (README §8). The next task is the merge to `main` per `MERGE_PLAN.md`.

## ADR-0026 (merge prep) — Relocate the stale, regenerable `tool/site` gallery to `legacy/site` before merging to main
- **Context.** Preparing the `integration/exercise-signature-explorer` → `main` merge (44 commits,
  ~1003 files), the largest non-functional payload was `hackathon/tool/site/` — **193 MB across 572
  files** (pre-rendered report HTML + multi-MB CSV tables, incl. a 60 MB `pathway_concordance_long.csv`).
  It is **not referenced by the app runtime** (no `apps/` usage), and HANDOFF flagged it as recording a
  **stale commit**. It is a **generated** artifact: the engine's `motrpac_probe site` CLI re-renders it
  into `tool/site/` from the example run folders (`scripts_regenerate_gallery.md`); the CLI does not
  read the committed copy.
- **Decision.** `git mv hackathon/tool/site legacy/site` (history preserved) and add `legacy/README.md`
  explaining that the folder holds set-aside, regenerable artifacts not used by the live app. Keep it
  tracked (per the user's instruction to relocate, not untrack), so the live `hackathon/tool/` tree
  carries source rather than a large regenerable output, while the artifact remains available for
  reference.
- **Reasons.** Relocating a stale, regenerable, unreferenced 193 MB artifact out of the tool source
  tree keeps `main`'s working tree honest about what is source vs. output, without losing the reference
  copy or rewriting history. A `git mv` preserves provenance (renames, not delete+add). Regeneration is
  documented, so the canonical way to get a fresh gallery is the CLI, not these files.
- **Alternatives rejected.** (a) Untrack via `git rm --cached` + `.gitignore` — the user chose to
  relocate to `legacy/`, not untrack; and it would drop the reference copy. (b) Leave it at
  `hackathon/tool/site/` — keeps a large stale output intermixed with source, the problem being solved.
  (c) History rewrite to shrink the repo — out of scope and destructive; NOTE that `git mv` does NOT
  shrink the repo (the 193 MB stays in history and in the tree, just relocated).
- **Consequences.** 572 files relocated `hackathon/tool/site/** → legacy/site/**`; new `legacy/README.md`
  and this ADR. No runtime source, engine, API, web, CI, or test file changed; the `motrpac_probe site`
  CLI still regenerates into `hackathon/tool/site/` (now untracked until regenerated). The repository
  size is unchanged (relocation, not removal). This is a merge-prep hygiene move; the merge itself
  follows MERGE_PLAN.md (merge commit, no squash, no rebase over the shared tool branch).
