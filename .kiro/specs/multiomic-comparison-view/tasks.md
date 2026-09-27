# Implementation Plan: Multi-Omic Comparison View

## Overview

This plan implements the delta-only Multi-Omic Comparison View on top of the already-merged RNA↔protein co-view. Work is sequenced so shared foundations land before the UI that consumes them, per the confirmed order:

1. **Backend labeling first (R4):** add `dataset` + `study_label` to every `ExplorerColumn` in `explorer.py`, sourced from the store, so the client has real provenance to render.
2. **Client type next (R4.6):** declare the two new fields on the `ExplorerColumn` TypeScript type so the UI can consume them.
3. **State-lift before Compare-All (R2 → R1/R3):** lift `Filters` out of `LayerView`/`PathwayView` into the parent `MotrpacExplorer` (shared filter bar), THEN add Compare-All Mode + the responsive Comparison Grid whose even-split units depend on the shared filters.
4. **Provenance surfacing (R5, R6, R10):** panel chips/headers and CSV row-builders gain `study_label`; human sex renders `all`; the footnote is preserved.
5. **Saved examples (R7):** regenerate the example JSONs so they carry the new fields and the tolerant cross-check stays green/skip-on-missing-artifact.
6. **Tests (R11):** new Explorer component tests + extended `apps/api/tests/test_explorer.py`; optional (`*`) property tests.
7. **Traceability + checkpoints last (R11):** per-feature ADR + `VERIFICATION.md` checkpoint + `REQUIREMENTS_TRACEABILITY.md` rows, with a final full-battery checkpoint.

The engine/store remains the single source of truth; **nothing new computes statistics** (R9.1). Feature A′ is UI-only — no new API request shape, no new chart types, no new statistics, no new classification (R8.3, R8.4). Feature B only *labels* columns from store-held provenance (R9.2).

**Grounding facts (from design.md — WIRE/RENDER/verify, do not re-implement):**
- `explorer.py` already reads `c.dataset` off store columns (`_category()`/`_columns()` branch on it); the store's `column_index()` groups a `dataset` column with values `human_acute` / `rat_train` — so `dataset` is a **confirmed store field** needing no discovery.
- Store columns carry **no per-column version string**; the study version lives in `store.load_provenance()` under keys `MotrpacRatTraining6moData` / `MotrpacHumanPreSuspensionAnalysis`. Confirm these keys exist in the running store build BEFORE wiring the version in (verify-before-claiming); if absent, the deterministic fallback (R4.7) applies.
- The existing `MotrpacExplorer.tsx` footnote hardcodes citation versions and is a page-level citation — **preserve it as-is** (R10.3). The per-panel `study_label` comes from the store (R5.5).
- The pair machinery (`pKey`, `pColIndex`, `pMolIndex`, `pValueAt`, `paired`, `pairGroups`, `pairCell`, `pairPanels`, the `pair`/`pairTime` cards) is kept **exactly as shipped** — no second pairing implementation (R3.1, R3.4).
- `.lr-cards` in `LiveDashboard.css` is today a single-column grid; the Comparison Grid adds a sibling `.lr-cards.xp-grid` modifier reusing the existing `@media (max-width: 860px)` breakpoint idiom.

**Legend:** Sub-tasks marked with `*` are optional (property/unit tests) and may be skipped for a faster MVP; the coding agent MUST NOT auto-implement `*` sub-tasks and MUST implement non-`*` sub-tasks. Documentation-only tasks are labeled **[DOCS ONLY]** and touch no runtime source.

---

## Tasks

- [x] 1. Add `dataset` + `study_label` to backend columns from the store (R4, R9.2)
  - [x] 1.1 Confirm store provenance keys before wiring the version (verify-before-claiming)
    - Confirm `motrpac_probe.store.load_provenance()` exposes the keys `MotrpacRatTraining6moData` and `MotrpacHumanPreSuspensionAnalysis` in the running store build before relying on them for `study_label` versions.
    - Record the confirmation (observed keys/values) so the version wiring is grounded, not assumed; if a key is absent, note that the R4.7 fallback path will apply for that dataset.
    - Do not fabricate a version; this task only observes the store.
    - _Requirements: 4.4, 11.4_
  - [x] 1.2 Add the `_study_label()` and `_dataset_for_species()` helpers in `explorer.py`
    - Add a module-level `_study_label(dataset)` that maps the store dataset code to a display phrase (`rat_train` → "endurance training", `human_acute` → "acute exercise") and appends the store provenance version from `store.load_provenance()` keyed by `MotrpacRatTraining6moData` / `MotrpacHumanPreSuspensionAnalysis` (e.g. "endurance training (2.0.0)").
    - Wrap the provenance read in an `@lru_cache` `_provenance()` that swallows load failures and returns `{}`; when the version is missing or the dataset is unknown, return the `dataset` code itself — never omit, never fabricate a study name (R4.7).
    - Add a tiny `_dataset_for_species(species)` returning `"human_acute"` for human and `"rat_train"` for rat, matching the store's own convention.
    - Compute NO statistics here (R9.1).
    - _Requirements: 4.3, 4.4, 4.5, 4.7, 9.1, 9.2_
  - [x] 1.3 Emit `dataset` + `study_label` from `_layer()`
    - In `explorer.py` `_layer()`, add `dataset=c.dataset` and `study_label=_study_label(c.dataset)` to the per-column `row = dict(...)`; leave all other fields and statistic values unchanged.
    - _Requirements: 4.1, 4.3, 4.5, 4.7_
  - [x] 1.4 Emit `dataset` + `study_label` from `_pathway_layer()`
    - In `explorer.py` `_pathway_layer()` (human-only path), add `dataset=_dataset_for_species("human")` (→ `"human_acute"`) and `study_label=_study_label(_dataset_for_species("human"))` to each column dict; leave all other fields unchanged.
    - _Requirements: 4.2, 4.3, 4.5, 4.7_

- [x] 2. Declare the new fields on the client `ExplorerColumn` type (R4.6)
  - [x] 2.1 Add `dataset` + `study_label` to `ExplorerColumn` in `client.ts`
    - In `apps/web/src/api/client.ts` (~line 412), add `dataset: string;` and `study_label: string;` to the `ExplorerColumn` interface. Leave `ExplorerValue`, `ExplorerLayer`, `ExplorerResponse`, and `api.explorerAnalyse` unchanged (no request-shape change, R8.3).
    - _Requirements: 4.6, 8.3_

- [~] 3. Checkpoint — backend fields + client type
  - Ensure `apps/api` tests pass and the web build/typecheck compiles with the new fields declared; ask the user if questions arise.

- [x] 4. Lift `Filters` into the parent as the Shared Filter Bar (R2, R3.5)
  - [x] 4.1 Hold `Filters` state in `MotrpacExplorer` and pass it down as props
    - In `apps/web/src/components/MotrpacExplorer.tsx`, add `const [filters, setFilters] = useState<Filters | null>(null)` in the parent; seed it from the active input's primary layer via the existing module-scope `defaults(...)` (hoisted, not re-derived) when a result loads or the active input changes.
    - Change `LayerView` and `PathwayView` to receive `filters` + `onFiltersChange` as props instead of owning `useState<Filters>`; remove their internal `const [f, setF]` and the `useEffect(() => setF(defaults(...)), [layer])`, mapping `f → filters` and `setF → onFiltersChange`.
    - Keep the shared-bar option derivations (`tissues`, `categories`, `timesOf(...)`, species scoping) sourced from the columns present in the current `ExplorerResponse`, reusing the existing `uniq(...)` helpers (R2.7).
    - _Requirements: 2.1, 2.2, 2.7_
  - [x] 4.2 Move the species/tissue/human-sex resolution into the parent resolver
    - Move the `pick("species", v)` rebuild (`defaults(cols.filter(c => c.species === v), …)`) and the `tissue`-change `time: ALL` reset out of `LayerView` into the parent's `onFiltersChange(next)` handler, copied **verbatim** so behavior is preserved (R3.4); apply the resolution once, centrally, then hand the resolved `Filters` to every panel.
    - Force `sex: ALL` whenever the resolved `species === "human"` (R2.5), matching today's `defaults()` and the human columns' `sex_label === "all"`.
    - _Requirements: 2.3, 2.4, 2.5, 3.4_
  - [x] 4.3 Preserve single-tab behavior including the RNA↔Protein Pair View
    - Keep Tabs Mode identical: the existing tab bar selects exactly one layer; the active `LayerView`/`PathwayView` renders, now driven by the shared `filters` prop. On the transcriptomics/proteomics tab, keep rendering the RNA↔Protein Pair View exactly as before, with the unchanged `partner` computation (`l.input === layer.input` sibling) and the unchanged `species|tissue|category|time|sex` key (R3.2, R3.3, R3.5).
    - Keep the existing empty-selection node (`lr-empty` "No MoTrPAC comparison matches these filters.") — never an error (R2.6, R10.2).
    - _Requirements: 3.2, 3.3, 3.5, 2.6, 10.2_
  - [ ]* 4.4 Write property tests for shared selection + option derivation + human sex + pairing
    - **Property 1: Shared filters select exactly the matching columns** — selection equals columns matching species/tissue/contrast/sex(or ALL)/time(or ALL); no match → empty state, no throw. Tag `Feature: multiomic-comparison-view, Property 1: shared filters select exactly the matching columns`. min 100 iterations, fast-check.
    - **Property 2: Shared-bar options derived from response columns** — offered option values equal the distinct column-field values in the response (respecting species scoping), never a value absent from columns. Tag `Feature: multiomic-comparison-view, Property 2: shared-bar options derived from response columns`. min 100 iterations, fast-check.
    - **Property 3: Human comparisons carry sex `all`** — effective sex for matching and label is `all` for any human column; no per-sex contrast on human results. Tag `Feature: multiomic-comparison-view, Property 3: human comparisons carry sex all`. min 100 iterations, fast-check.
    - **Property 4: RNA↔Protein pairing uses the unchanged five-part key** — join pairs by `species|tissue|category|time|sex`, producing the same paired set as the pre-refactor `LayerView` logic. Tag `Feature: multiomic-comparison-view, Property 4: RNA-protein pairing uses the unchanged five-part key`. min 100 iterations, fast-check.
    - New file: `apps/web/src/components/MotrpacExplorer.test.tsx` (or a sibling `*.property.test.tsx`). Use fast-check; never hand-roll PBT.
    - _Requirements: 2.3, 2.5, 2.6, 2.7, 3.1, 3.2, 3.3, 3.4, 5.3, 10.2, 10.4; Properties 1, 2, 3, 4_

- [x] 5. Add Compare-All Mode + the responsive Comparison Grid (R1, R3.6, R3.7)
  - [x] 5.1 Add the mode toggle and unit model in `MotrpacExplorer.tsx`
    - Add `type Mode = "tabs" | "compare"` with `const [mode, setMode] = useState<Mode>("tabs")` (default Tabs Mode, R1.2) and a toggle control in the results header next to the existing `RuleSwitch` (R1.1).
    - Compute the per-input `Unit` list: RNA+protein both available → one fused `{ kind: "pair" }` unit (a `LayerView` with `partner`); RNA-only/protein-only/phospho → `{ kind: "layer" }`; pathways → `{ kind: "pathways" }` (`PathwayView`); unavailable layers contribute NO unit (R1.7).
    - In Compare-All Mode, render each unit in the grid using the existing view for that layer — no new chart type (R1.6, R8.4); RNA and protein render as the single fused pair unit occupying one cell, never two independent panels (R3.6); phospho and pathways are their own cells (R3.7). Returning to Tabs Mode shows exactly one layer again (R1.8).
    - _Requirements: 1.1, 1.2, 1.3, 1.6, 1.7, 1.8, 3.6, 3.7, 8.4_
  - [x] 5.2 Add the responsive even-split grid CSS
    - In `apps/web/src/components/LiveDashboard.css`, add a `.lr-cards.xp-grid` modifier that switches from single-column to equal `1fr` tracks across the displayed units and reflows to fewer columns at the existing `@media (max-width: 860px)` breakpoint (e.g. 2-up → single column). Pure responsive CSS — no splitter/drag-to-resize (R1.4, R1.5, R8.2).
    - _Requirements: 1.4, 1.5, 8.2_
  - [ ]* 5.3 Write Compare-All example tests
    - Assert: default is Tabs Mode with one layer visible and switching renders all available units (R1.1, R1.2, R1.3, R1.8); an unavailable layer's cell is omitted with no reserved space (R1.7); with RNA+protein available exactly one fused pair-unit cell renders (not two panels) and phospho/pathways are their own cells (R3.6, R3.7).
    - Add to `apps/web/src/components/MotrpacExplorer.test.tsx`.
    - _Requirements: 1.1, 1.2, 1.3, 1.7, 1.8, 3.6, 3.7_

- [x] 6. Surface full provenance in headers/chips, tables, and CSV (R5, R6, R10)
  - [x] 6.1 Add the `provenanceLabel(column)` helper and panel chips/headers
    - Add a `provenanceLabel(c: ExplorerColumn)` helper returning the dot-joined `species · study_label · tissue_label · contrast · sex · time` (human sex → `all`), e.g. "Rat · endurance training (2.0.0) · gastrocnemius · endurance vs control · male · 8 wk"; the study/dataset token comes from `c.study_label`, never a client literal (R5.5).
    - Feed these tokens into each panel's `Card` `chips`; in Compare-All Mode each unit shows its own provenance (R5.4); the fused pair unit shows the shared species/study/tissue/contrast/sex/time both layers are matched on.
    - _Requirements: 5.1, 5.2, 5.3, 5.4, 5.5, 10.4_
  - [x] 6.2 Add `study_label` to every CSV row-builder
    - Add `study_label` (from `c.study_label`) to each existing Explorer CSV row object — the `pair` card CSV, the `molecules` card CSV, the set-level `comparisons` CSV, and the `PathwayView` CSV (which also gains `species`) — so every exported row carries species · study/dataset · tissue · contrast · sex · time (R6.1, R6.2, R6.3, R6.4). Render human sex as `all`. Do not add a new statistic; only the label field.
    - _Requirements: 6.1, 6.2, 6.3, 6.4, 10.1_
  - [x] 6.3 Preserve the page footnote
    - Leave the existing `MotrpacExplorer.tsx` footnote (published summary statistics; association not treatment effects; hardcoded citation versions) unchanged (R10.3).
    - _Requirements: 10.3_
  - [ ]* 6.4 Write the provenance CSV property test + chip/footnote examples
    - **Property 5: Every emitted result row carries full provenance from the column** — every CSV row has non-empty `species`, `study_label`, `tissue`, `contrast`, `sex`, `time`, and `study_label` equals the source `ExplorerColumn.study_label`. Tag `Feature: multiomic-comparison-view, Property 5: every emitted result row carries full provenance from the column`. min 100 iterations, fast-check.
    - Example smokes: panel header/chips contain species · study_label · tissue · contrast · sex · time and the Provenance Label matches the R5.2 shape (R5.1, R5.4, R5.5); the published-summary-statistics footnote still renders (R10.3); guard that `api.explorerAnalyse` request body is unchanged and no new chart component is imported (R8.3, R8.4, R8.6).
    - Add to `apps/web/src/components/MotrpacExplorer.test.tsx`.
    - _Requirements: 5.1, 5.2, 5.4, 5.5, 6.2, 6.3, 6.4, 8.3, 8.4, 8.6, 10.3; Property 5_

- [x] 7. Regenerate saved examples with the new fields (R7)
  - [x] 7.1 Regenerate the example JSONs
    - Run `apps/api/scripts_save_examples.py` (no script code change needed) so each `apps/web/public/examples/*.json` `ExplorerColumn` carries `dataset` + `study_label`; the regeneration adds only the two label fields and changes no statistic value (R7.1, R7.4).
    - _Requirements: 7.1, 7.4_
  - [ ]* 7.2 Confirm the tolerant saved-vs-live cross-check stays green/skip-on-missing-artifact
    - Run the existing `test_saved_example_results_match_the_live_analysis` in `apps/api/tests/test_explorer.py` unchanged; confirm it reports no mismatch within tolerance and skips a layer with an explicit reason when a required built artifact is absent (R7.2, R7.3). Do not modify the cross-check logic.
    - _Requirements: 7.2, 7.3_

- [ ] 8. Extend backend column-labeling tests (R4, R6.4, R7.1)
  - [ ]* 8.1 Extend `test_explorer.py` for column labeling + pathway path + saved-example fields
    - **Property 6: Every backend-built column is labeled** — for any `ExplorerResponse` from `analyse()`, every `ExplorerColumn` across every available layer (including the pathway layer) has `dataset` equal to its store code and a non-empty `study_label`; a column with no store version has `study_label == dataset`. Tag `Feature: multiomic-comparison-view, Property 6: every backend-built column is labeled with dataset and study_label`. min 100 iterations, Hypothesis.
    - Examples: a rat column's `study_label` matches `endurance training (<version>)`, a human column's matches `acute exercise (<version>)` (version from `store.load_provenance()`); with provenance monkeypatched empty, `study_label == dataset` (R4.7); `_pathway_layer()` columns carry `dataset == "human_acute"` and a non-empty `study_label` (R4.2); after running `scripts_save_examples.py`, every `ExplorerColumn` in each regenerated JSON has `dataset` + `study_label` (R7.1).
    - Extend `apps/api/tests/test_explorer.py`. Use Hypothesis (already present per `.hypothesis/`); never hand-roll PBT.
    - _Requirements: 4.1, 4.2, 4.3, 4.5, 4.7, 6.4, 7.1; Property 6_

- [x] 9. Record traceability artifacts **[DOCS ONLY]** (R11)
  - Add the feature ADR to `DECISIONS.md` recording the UI-only comparison grid, the props-based (not context) state lift, and the B2 store-sourced labeling approach (R11.1).
  - Add a `VERIFICATION.md` checkpoint recording the exact commands and resulting counts: web build + tests (`pnpm --dir apps/web build`; `pnpm --dir apps/web test`), `apps/api` pytest including the Saved-vs-Live Cross-Check (`pytest apps/api/tests`), and a note that the `motrpac_probe` engine suite is untouched by this feature (R11.2).
  - Add `REQUIREMENTS_TRACEABILITY.md` rows mapping each requirement (R1–R10) to its code path (`explorer.py`, `client.ts`, `MotrpacExplorer.tsx`, `LiveDashboard.css`, `scripts_save_examples.py`), its test, and status (R11.3).
  - Verify every count/pass by running the relevant command BEFORE claiming it (R11.4). Do NOT record counts that were not observed.
  - _Requirements: 11.1, 11.2, 11.3, 11.4_

- [x] 10. Final checkpoint — full verification battery
  - Run the web build + vitest suite and the `apps/api` pytest suite (including the Saved-vs-Live Cross-Check); mark the feature complete only when BOTH report zero failures. Do NOT claim counts before running the commands. The engine suite is not run/changed by this feature. Ask the user if questions arise.
  - _Requirements: 11.5_

## Notes

- Tasks marked with `*` are optional (property/unit tests) and can be skipped for a faster MVP; the coding agent MUST NOT auto-implement `*` sub-tasks and MUST implement non-`*` sub-tasks.
- Task 9 is labeled **[DOCS ONLY]** and touches no runtime source.
- Each task references specific requirement sub-clauses and, where relevant, the design property numbers for traceability.
- **Ordering intent:** backend field + client type (tasks 1–2) land before the UI that consumes them; the state-lift (task 4) precedes Compare-All (task 5) because the grid's even-split units depend on the shared filters; provenance surfacing (task 6) and saved examples (task 7) follow; traceability + the final checkpoint (tasks 9–10) are last.
- **Do NOT re-implement** the RNA↔Protein pair machinery (`pKey`/`pColIndex`/`pMolIndex`/`pValueAt`/`paired`/`pairGroups`/`pairCell`/`pairPanels`, the `pair`/`pairTime` cards) — tasks WIRE/preserve it; `dataset` is a confirmed store field (no discovery); `analyse()` already returns every layer in one payload (no API-shape change).
- **Out of scope (R8) — no task may implement any of these:** cross-omic concordance or classification summary (R8.1); drag-to-resize / splitter panes (R8.2); any Explorer API request-shape change (R8.3); any new or altered per-omic statistic/classification (R8.4); changes to the Discordance, Live, Generalized, or Metabolomics views (R8.5); re-specifying or reimplementing the RNA↔Protein Pair View beyond not regressing it (R8.6).
- Cross-cutting R9 (source-of-truth boundary; store-derived labels; contrast vs sex are distinct axes) and R10 (honest labeling; honest empty states; preserved footnote; no per-sex human contrast) are woven into tasks 1, 4, 5, and 6 rather than deferred.
- Each property test is a single test at min 100 iterations, tagged `Feature: multiomic-comparison-view, Property N: <text>`, using fast-check (web) or Hypothesis (`apps/api`); never hand-rolled.

## Task Dependency Graph

Ordering intent: backend labeling (R4) → client type (R4.6) → state-lift (R2) → Compare-All (R1/R3) → provenance surfacing (R5/R6) → saved examples (R7) → backend tests → traceability (R11). Checkpoints (3, 10) and top-level parent tasks are not scheduled; only leaf sub-tasks and standalone leaf tasks appear below. Optional (`*`) test tasks are included.

```json
{
  "waves": [
    { "id": 0, "tasks": ["1.1"] },
    { "id": 1, "tasks": ["1.2"] },
    { "id": 2, "tasks": ["1.3", "1.4"] },
    { "id": 3, "tasks": ["2.1"] },
    { "id": 4, "tasks": ["4.1"] },
    { "id": 5, "tasks": ["4.2", "4.3"] },
    { "id": 6, "tasks": ["4.4", "5.1"] },
    { "id": 7, "tasks": ["5.2", "6.1", "6.2", "6.3"] },
    { "id": 8, "tasks": ["5.3", "6.4", "7.1"] },
    { "id": 9, "tasks": ["7.2", "8.1"] },
    { "id": 10, "tasks": ["9"] }
  ]
}
```
