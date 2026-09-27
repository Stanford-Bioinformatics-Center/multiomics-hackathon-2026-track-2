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
