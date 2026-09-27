# Requirements Document

## Introduction

This feature extends the MoTrPAC Explorer with a generalized multi-omic comparison view built on top of the RNA↔protein side-by-side view that the upstream `motrpac-explorer` branch just shipped and that is now merged into this branch. It is a deliberately re-scoped, delta-only feature: the merge already delivered a focused two-layer co-view (a combined "pair" HeatGrid and a "pairTime" SmallMultiples that match molecules across transcriptomics and proteomics by `species|tissue|category|time|sex`), and this spec builds ON that reality rather than re-specifying it.

Two capabilities are added. Feature A' generalizes the existing RNA↔protein pairing into a "Compare all layers" grid that shows every available layer's existing view simultaneously, driven by ONE shared filter bar lifted out of each per-panel `LayerView` into the parent `MotrpacExplorer`, so all panels always show the same biological context. Feature B (variant B2) adds explicit study/dataset labeling: `dataset` and `study_label` fields are added to each `ExplorerColumn` in `explorer.py`, sourced from the motrpac_probe store per species/assay, and surfaced in panel headers/chips, on-screen tables, and CSV exports so every displayed result makes explicit WHAT is compared and IN WHOM.

Feature A' is UI-only: it introduces no new API request shape, no new chart types, no new statistics, and no new classification. The existing `ExplorerResponse` already contains every layer in one payload, so showing layers side by side adds no computation. Feature B is a small backend addition that only labels columns from store-held provenance; it computes no statistics.

The engine/store remains the single source of truth. No statistics are computed in the web client or the FastAPI/adapter layer; study/dataset labels are read from the store, never hardcoded in React (this is the whole reason B2 is chosen over deriving labels in the frontend). Honest labeling is a first-class constraint: every displayed result carries species, dataset/study, tissue, contrast (category), sex, and time. Note that `category` is the CONTRAST (EE-CON endurance-vs-control, EE-EE before→after, TRAIN-SED trained-vs-sedentary, CON-CON control drift) while `sex` is a SEPARATE axis (male/female for rat, "all" for human); the contrast runs WITHIN each sex, so "male vs female" is not a contrast.

## Glossary

- **Explorer**: The `MotrpacExplorer` React component (`apps/web/src/components/MotrpacExplorer.tsx`) that fetches all layers in one `ExplorerResponse` and today shows one layer at a time via a tab index (`tab`/`setTab`).
- **LayerView**: The Explorer sub-component that renders a single layer's charts (overview, molecules, over-time, set-level bars) and, when given a `partner` layer on the transcriptomics/proteomics tab, also renders the RNA↔protein pair views. It currently owns its own per-panel filter state.
- **PathwayView**: The Explorer sub-component that renders the pathways layer (`layer.id === "pathways"`) as a DotPlot; it owns its own filter state distinct from `LayerView`.
- **Supported Layer**: One of the four layers the Explorer can display for a given result: genes/transcriptomics, proteins/proteomics, phospho/phosphoproteomics, and pathways.
- **Tabs Mode**: The existing single-omic default in which exactly one supported layer's view is visible at a time, selected by the tab bar.
- **Compare-All Mode**: The new mode added by this feature in which every available supported layer's existing view is displayed simultaneously in a responsive grid.
- **Comparison Grid**: The responsive layout that renders the available layers' views simultaneously in Compare-All mode, splitting screen width evenly across the displayed layers and reflowing on narrow screens.
- **Shared Filter Bar**: The single filter/flag control set — species, tissue, contrast, sex, time — lifted out of each `LayerView` into the parent Explorer, driving all displayed panels together.
- **Filters**: The filter state object `{ species, tissue, contrast, sex, time, layer }` used to select which MoTrPAC comparisons a panel displays.
- **RNA↔Protein Pair View**: The already-shipped two-layer co-view rendered by `LayerView` when a `partner` layer is present — a combined "pair" HeatGrid and a "pairTime" SmallMultiples matching molecules across transcriptomics and proteomics by `species|tissue|category|time|sex`.
- **ExplorerColumn**: The per-comparison metadata record on the API response (`apps/web/src/api/client.ts`, ~line 412) carrying `id, species, tissue, tissue_label, layer, category, exercise, time, time_label, time_rank, sex, n_tested?, n_measured?, set_*, rho_*`. It currently carries no explicit dataset/study name.
- **category (Contrast)**: The contrast axis carried on `ExplorerColumn.category`: EE-CON (endurance vs control), EE-EE (endurance before→after), TRAIN-SED (trained vs sedentary), CON-CON (control drift), and related codes.
- **sex axis**: The separate `ExplorerColumn.sex` axis — `male`/`female` for rat, `all` for human — orthogonal to the contrast.
- **dataset**: The store-held dataset code for a comparison column (`rat_train` for rat endurance training, `human_acute` for the human acute bout), already present on store columns and to be carried onto each `ExplorerColumn`.
- **study_label**: A human-readable study/dataset label to be added to each `ExplorerColumn`, sourced from the store's provenance (e.g. "endurance training (2.0.0)" for `MotrpacRatTraining6moData 2.0.0`, "acute exercise (2.0.8)" for `MotrpacHumanPreSuspensionAnalysis 2.0.8`).
- **Provenance Label**: The full displayed description of a result: species · study/dataset · tissue · contrast (category) · sex · time — e.g. "Rat · endurance training (2.0.0) · gastrocnemius · endurance vs control · male · 8 wk".
- **Store**: The motrpac_probe store accessed by `explorer.py`, holding the published MoTrPAC summary statistics and provenance (dataset codes and package versions). The single source of truth for all statistics and study/dataset labels.
- **Saved Example Results**: The static JSON files in `apps/web/public/examples/*.json`, generated by `apps/api/scripts_save_examples.py`, that the website loads for built-in examples and that `apps/api/tests/test_explorer.py` cross-checks against a live run.
- **Saved-vs-Live Cross-Check**: The tolerant, artifact-aware test in `apps/api/tests/test_explorer.py` (`test_saved_example_results_match_the_live_analysis`) that deep-compares saved example JSON against a fresh live analysis within a numeric tolerance and skips a layer when a required built artifact is absent.

## Requirements

### Requirement 1: Compare-All Mode Toggle and Comparison Grid

**User Story:** As an analyst, I want to switch the Explorer from single-tab viewing to a "Compare all layers" grid, so that I can see every available omic's existing view at once for scroll-free visual comparison.

#### Acceptance Criteria

1. THE Explorer SHALL provide a control that toggles between Tabs Mode and Compare-All Mode.
2. WHEN the Explorer loads a result, THE Explorer SHALL default to Tabs Mode with a single supported layer's view visible.
3. WHEN the user selects Compare-All Mode, THE Explorer SHALL render, in the Comparison Grid, the existing view of every available supported layer simultaneously.
4. WHEN Compare-All Mode is active, THE Comparison Grid SHALL split the available width evenly across the displayed layers.
5. WHERE the viewport is narrow, THE Comparison Grid SHALL reflow the panels into fewer columns (for example from a two-by-two arrangement to a stacked single column).
6. WHEN Compare-All Mode renders a supported layer's panel, THE Explorer SHALL use that layer's existing view (transcriptomics, proteomics, and phosphoproteomics via LayerView; pathways via PathwayView) without introducing a new chart type.
7. IF a supported layer is unavailable for the current result, THEN THE Comparison Grid SHALL omit that layer's panel and SHALL NOT reserve grid space for it.
8. WHEN the user returns to Tabs Mode from Compare-All Mode, THE Explorer SHALL again display exactly one supported layer's view at a time.

### Requirement 2: Shared Filter Bar Across All Panels

**User Story:** As an analyst, I want one filter bar to drive every panel, so that all omics always show the same biological context instead of each panel tracking its own filters.

#### Acceptance Criteria

1. THE Explorer SHALL present a single Shared Filter Bar exposing species, tissue, contrast, sex, and time controls.
2. THE Explorer SHALL hold the Filters state (species, tissue, contrast, sex, time) in the parent Explorer component rather than within each LayerView.
3. WHEN the user changes a value in the Shared Filter Bar, THE Explorer SHALL apply the same species, tissue, contrast, sex, and time selection to every displayed panel.
4. WHILE Compare-All Mode is active, THE Explorer SHALL render every displayed layer's panel using the biological context defined by the Shared Filter Bar.
5. WHERE the selected species is human, THE Shared Filter Bar SHALL treat sex as `all` for every panel, consistent with the human dataset carrying no per-sex contrast.
6. IF the current filter selection matches no MoTrPAC comparison for a given layer, THEN that layer's panel SHALL display its existing empty or no-match state and SHALL NOT display an error.
7. WHEN the Shared Filter Bar offers option values for a control, THE Explorer SHALL derive those values from the columns present in the current `ExplorerResponse`.

### Requirement 3: Generalization Subsumes and Preserves the RNA↔Protein Pair View

**User Story:** As a maintainer, I want the comparison grid to generalize the existing RNA↔protein pairing rather than duplicate or break it, so that the already-shipped co-view keeps working — and so that in Compare-All Mode RNA and protein are always presented as one jointly-matched unit while other omics appear as their own panels.

#### Acceptance Criteria

1. THE Comparison Grid SHALL treat the RNA↔Protein Pair View as the two-layer case of the general multi-layer comparison and SHALL NOT introduce a second, parallel implementation of molecule pairing.
2. WHEN both transcriptomics and proteomics are available for a result, THE Explorer SHALL continue to render the RNA↔Protein Pair View (the combined pair HeatGrid and pairTime SmallMultiples) matching molecules by `species|tissue|category|time|sex`.
3. WHEN the Shared Filter Bar drives the RNA↔Protein Pair View, THE Explorer SHALL match molecules across the two layers using the same `species|tissue|category|time|sex` key used before this feature.
4. THE Explorer SHALL NOT change the statistics, values, or matching keys of the RNA↔Protein Pair View as part of adding Compare-All Mode or the Shared Filter Bar.
5. WHEN Tabs Mode is active on the transcriptomics or proteomics tab, THE Explorer SHALL render the RNA↔Protein Pair View as it did before this feature.
6. WHEN Compare-All Mode is active and both transcriptomics and proteomics are available, THE Comparison Grid SHALL render RNA and protein together as the SINGLE fused RNA↔Protein Pair Unit (the molecule-matched pair views) occupying one cell of the grid, and SHALL NOT render RNA and protein as two independent, unjoined panels.
7. WHEN Compare-All Mode is active, THE Comparison Grid SHALL render phosphoproteomics and pathways as their own separate panels (not molecule-fused into the RNA↔Protein Pair Unit), and SHALL split the available width evenly across the displayed units (the fused RNA↔Protein Pair Unit counting as one unit).

### Requirement 4: Add dataset and study_label to ExplorerColumn from the Store

**User Story:** As a frontend developer, I want each ExplorerColumn to carry an explicit dataset and study_label sourced from the store, so that the UI can display provenance without hardcoding study names in React.

#### Acceptance Criteria

1. WHEN `explorer.py` builds an `ExplorerColumn` in `_layer()`, THE Explorer_Backend SHALL include a `dataset` field and a `study_label` field on that column.
2. WHEN `explorer.py` builds an `ExplorerColumn` in `_pathway_layer()`, THE Explorer_Backend SHALL include a `dataset` field and a `study_label` field on that column.
3. THE Explorer_Backend SHALL source the `dataset` value from the store column's dataset code (for example `rat_train` for rat and `human_acute` for human).
4. THE Explorer_Backend SHALL source the `study_label` value from the store's provenance for the column's species and assay and SHALL NOT hardcode study labels in the React client.
5. THE study_label for a rat endurance-training column SHALL identify the rat endurance-training study, and the study_label for a human acute-bout column SHALL identify the human acute-exercise study.
6. THE `ExplorerColumn` TypeScript type in `apps/web/src/api/client.ts` SHALL declare the `dataset` and `study_label` fields.
7. IF the store does not provide a study label for a column's species and assay, THEN THE Explorer_Backend SHALL set `study_label` to the column's `dataset` value rather than omit the field or fabricate a study name.

### Requirement 5: Surface Full Provenance in Panel Headers and Chips

**User Story:** As an analyst, I want each panel to state what is compared and in whom, so that I never read a chart without knowing its species, study, tissue, contrast, sex, and time.

#### Acceptance Criteria

1. WHEN a panel renders a result, THE Explorer SHALL display, in that panel's header or chips, the species, study/dataset (study_label), tissue, contrast (category), sex, and time for the displayed comparison.
2. THE Explorer SHALL render the Provenance Label in a form that makes explicit what is compared and in whom (for example "Rat · endurance training (2.0.0) · gastrocnemius · endurance vs control · male · 8 wk").
3. WHERE the displayed species is human, THE Explorer SHALL render the sex component of the Provenance Label as `all`.
4. WHEN Compare-All Mode displays multiple panels, THE Explorer SHALL show each panel's own species, study/dataset, tissue, contrast, sex, and time.
5. THE Explorer SHALL derive the study/dataset component of the Provenance Label from the `study_label` field on `ExplorerColumn` rather than from a value hardcoded in the client.

### Requirement 6: Provenance Columns in On-Screen Tables and CSV Exports

**User Story:** As an analyst, I want provenance visible in the tables and exports, so that every row states its species, study/dataset, tissue, contrast, sex, and time.

#### Acceptance Criteria

1. WHERE the Explorer renders an on-screen table of results, THE Explorer SHALL include visible columns for species, study/dataset (study_label), tissue, contrast (category), sex, and time.
2. WHEN the Explorer exports a CSV of results, THE Explorer SHALL include the species, study/dataset (study_label), tissue, contrast (category), sex, and time fields for each exported row.
3. THE Explorer SHALL NOT render a result row in a table without an accompanying species, study/dataset, tissue, contrast, sex, and time for that row.
4. THE Explorer SHALL populate the study/dataset column and CSV field from the `study_label` field on `ExplorerColumn`.

### Requirement 7: Regenerate Saved Example Results with the New Fields

**User Story:** As a maintainer, I want the saved example JSONs regenerated to carry the new fields, so that examples work offline and the tolerant saved-vs-live cross-check stays green.

#### Acceptance Criteria

1. WHEN `apps/api/scripts_save_examples.py` is run after the `dataset` and `study_label` fields are added, THE Saved Example Results SHALL include `dataset` and `study_label` on every `ExplorerColumn` in each regenerated JSON.
2. WHEN the Saved-vs-Live Cross-Check runs after regeneration, THE Saved-vs-Live Cross-Check SHALL report no mismatch between each saved example JSON and a fresh live analysis, within the test's existing numeric tolerance.
3. IF a saved example layer cannot be reproduced live because a required built artifact is absent, THEN THE Saved-vs-Live Cross-Check SHALL skip that layer with an explicit reason rather than fail.
4. THE regeneration SHALL NOT change any statistic value in the Saved Example Results beyond adding the `dataset` and `study_label` label fields.

### Requirement 8: Scope Boundaries (Explicitly Out of Scope)

**User Story:** As a maintainer, I want the boundaries of this feature stated explicitly, so that reviewers can confirm the delta added nothing beyond the two intended features.

#### Acceptance Criteria

1. THE Explorer SHALL NOT add any cross-omic concordance or classification summary as part of this feature.
2. THE Explorer SHALL NOT add drag-to-resize or splitter panes for the Comparison Grid.
3. THE Explorer_Backend SHALL NOT change the Explorer API request shape as part of this feature.
4. THE Explorer_Backend SHALL NOT add or alter any per-omic statistic or classification as part of this feature; the existing `ExplorerResponse` already contains every layer, so side-by-side display adds no computation.
5. THE feature SHALL NOT modify the Discordance, Live, Generalized, or Metabolomics views.
6. THE feature SHALL NOT re-specify or reimplement the already-shipped RNA↔Protein Pair View beyond requiring that it not be regressed.

### Requirement 9: Source-of-Truth and Backend-Derivation Boundary (Cross-Cutting)

**User Story:** As a maintainer, I want the statistical source-of-truth boundary and store-sourced labeling enforced, so that the presentation layers never compute or invent analysis.

#### Acceptance Criteria

1. THE Web_Client and Explorer_Backend SHALL NOT compute statistics; the Store (motrpac_probe) SHALL remain the source of truth for all statistics.
2. THE Explorer_Backend SHALL derive the `dataset` and `study_label` values from the Store and SHALL NOT hardcode study/dataset labels in the React client.
3. THE Explorer SHALL show species, study/dataset, contrast, sex, and time on every displayed result.
4. THE feature SHALL treat the contrast (category) and the sex axis as distinct dimensions and SHALL NOT present a "male vs female" contrast.

### Requirement 10: Honest Labeling and No-Data Handling (Cross-Cutting)

**User Story:** As a reader, I want honest labeling and honest empty states, so that no panel implies data that does not exist and no missing data reads as an error.

#### Acceptance Criteria

1. THE Explorer SHALL render, alongside every displayed result, its species, study/dataset, contrast, sex, and time.
2. IF a filter selection has no data for a given layer, THEN that layer's panel SHALL show its existing empty or no-match state as an honest outcome and SHALL NOT present it as a failure.
3. THE Explorer SHALL preserve the existing page footnote stating that MoTrPAC values are published summary statistics and that comparisons describe association, not treatment effects.
4. THE Explorer SHALL NOT label any human result with a per-sex contrast, consistent with human sex being `all`.

### Requirement 11: Traceability Discipline and Verify-Before-Claiming (Cross-Cutting)

**User Story:** As a maintainer, I want this feature traced through an ADR, a verification checkpoint, and a traceability matrix, so that every claim is reproducible and every decision is recorded.

#### Acceptance Criteria

1. WHERE this feature is implemented, THE Team SHALL add a corresponding ADR entry to `DECISIONS.md` recording the choice of a UI-only comparison grid plus the B2 backend labeling approach.
2. WHERE this feature is implemented, THE Team SHALL record a `VERIFICATION.md` checkpoint containing the exact commands run and the resulting test counts for the web and `apps/api` suites.
3. WHERE this feature is implemented, THE Team SHALL add rows to `REQUIREMENTS_TRACEABILITY.md` mapping each requirement to its code path, test, and status.
4. THE Team SHALL NOT claim a count, a passing test, or a behavior without first running the relevant test or command.
5. WHEN the verification battery runs, THE Team SHALL run the web tests and build and the `apps/api` tests, including the Saved-vs-Live Cross-Check, and SHALL mark the feature complete only when those report zero failures.


## Future Scope (Not This Feature)

The following is explicitly deferred and is recorded here so the intent is not lost:

- **User-selectable N-way comparison with arbitrary joint pairing.** Ideally the user could choose how
  many omics to compare (any subset of the supported layers) and select ANY TWO of them to be
  molecule-matched/fused jointly (as RNA↔protein is fused today), rather than the fusion being fixed to
  the RNA↔protein pair. In this feature the fused unit is fixed to RNA↔protein (Requirement 3, Reading A)
  and the other omics appear as independent panels; a future feature would generalize the fusion so the
  user picks which two layers are jointly matched and which others accompany them as panels. This would
  require new selection UI and a generalized molecule-matching/charting path beyond the current
  `species|tissue|category|time|sex` RNA↔protein join, and is therefore out of scope here.
