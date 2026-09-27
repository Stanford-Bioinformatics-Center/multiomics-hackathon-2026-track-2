# Requirements Document

## Introduction

This feature finalizes the Track 2 "Omic Discordance Explained" MVP. The Track 2 challenge is: omic layers within MoTrPAC data are often discordant, so the goal is to compare two or more compatible layers (RNA, total protein, PTM) within a defined tissue and identify a model that explains when they agree and when they disagree. The intended outputs are (a) a catalog of concordant and discordant events and (b) a predictive model of discordance. A standing caveat applies: PTM signal is not automatically a measure of modification occupancy, enzyme activity, or functional consequence.

The core of these Track 2 deliverables already exists in the repository but is hidden as CLI-only pilots with committed demo outputs (the discordance event catalog, the out-of-fold ridge predictive model, and the PTM phosphosite-vs-parent divergence audit). The FastAPI layer already computes `layer_discordance` but renders it in zero UI views. This feature surfaces those deliverables, makes discordance the application's main/default view organized as a three-stage flow (query builder → data visualization → discordance interpretation), hardens fresh-clone reproducibility, adds automated end-to-end and API verification wired into CI, and documents the merge-to-main plan without executing it.

The scope is deliberately sequenced. Workflow tracing (documentation) comes first so the team understands what each analysis path does before analyzing success criteria. Feature work follows. A final consolidated verification and hardening pass closes the effort. Throughout, the statistical engines (mprobe, query_core, and the standalone discordance modules) remain the single source of truth; no statistics are computed in the FastAPI layer or the React client.

The honest-framing posture is a first-class constraint: results represent cross-cohort association, not a PAH treatment claim; cohorts are kept separate; and the headline directed result for male rat SKM-GN protein at 8 weeks is q = 0.0584 (not significant) under a frozen 52-column BH multiplicity family.

## Glossary

- **Discordance**: The condition in which two or more compatible omic layers (RNA, total protein, PTM) within a defined tissue disagree in direction or magnitude of change for the same gene or feature.
- **Concordant/Discordant Event Catalog**: The per-feature classification table (`catalog.csv`) produced by `build_discordance.py`, assigning each event one of the classes `supported_concordant`, `supported_opposite`, `rna_response_protein_equivalent`, or `indeterminate`.
- **Predictive Model of Discordance**: The out-of-fold ridge regression comparison (`model_metrics.csv`, `model_predictions.csv`) across three model variants — zero-change baseline, RNA-only, and temporal — reported with out-of-fold R2, MAE, and RMSE, framed honestly as a weak prediction (a result, not a success claim).
- **cameraPR (camera_p)**: The nominal (raw) p-value produced by the CAMERA pre-ranked gene-set test, carried on `ColumnResult` as `camera_p`, distinct from the family-adjusted BH q-value.
- **BH Family**: The Benjamini-Hochberg multiplicity-adjustment family. The frozen 52-column family is the canonical adjusted set for the directed path; the 16-column family is a labeled sensitivity analysis only.
- **EE-CON**: The endurance-exercise-versus-control contrast context used in the demo datasets (e.g., the c2.0 EE-CON demo).
- **PTM / Phosphosite-Parent**: A post-translational modification measurement (e.g., phosphosite) compared against its parent protein; divergence between the two is audited by `audit_ptm_parent.py`.
- **Ortholog Link**: The cross-species mapping relation connecting a rat feature to its human counterpart (or vice versa), surfaced as an explicit relation in cross-species results.
- **run_id**: The provenance identifier that captures the analysis parameters (including the FDR threshold) for a given computation run.
- **Point → Evidence → Source Chain**: The traceable path from a rendered data point in the UI to its underlying evidence record to its source-of-truth origin (module/CLI output and provenance), asserted by `test_point_to_evidence_to_source_chain`.
- **Source-of-Truth Boundary**: The architectural rule that all statistics originate in the engines (mprobe, query_core, standalone discordance modules); the FastAPI layer and React client only serve, shape (JSON-safe), and render, never compute statistics.
- **Layer Discordance (LayerDiscordanceRow)**: The per-row discordance record already computed and exported by `apps/api` `service.py`, not yet rendered in any UI view.
- **Option-1 Adapter Pattern**: The metabolomics-established pattern of a read-only API adapter that serves committed demo module outputs, guarded by an `available()` check, returning JSON-safe payloads, with no statistics computed in the adapter.

## Requirements

### Requirement 1: Trace All Workflows into WORKFLOWS.md (Documentation First)

**User Story:** As a contributor, I want a single traced map of every analysis workflow, so that I understand what each path does before analyzing success criteria or building features.

#### Acceptance Criteria

1. WHEN work on Requirement 1 begins, THE Finalization_Team SHALL create a single `WORKFLOWS.md` document at the repository root before any commit that modifies feature-implementation source files.
2. THE `WORKFLOWS.md` document SHALL include exactly one mermaid diagram and one numbered narrative for each of the following five analysis types: directed_signature, metabolomics_case_study, generalized_query, the discordance catalog and predictive model, and the PTM-parent audit path.
3. WHERE an analysis type is described, THE `WORKFLOWS.md` document SHALL trace, in order, all seven of these stages: entry point, module or function name, data read (input file or dataset), transform, output fields, provenance identifier or run_id, and the UI view that renders the output; and IF any of the seven stages does not exist for that analysis type, THEN THE document SHALL state "not applicable" for that stage rather than omit it.
4. WHERE an API endpoint exists, THE `WORKFLOWS.md` document SHALL document its request-to-response contract, including request inputs, response outputs, and the triggering condition for each endpoint.
5. THE `WORKFLOWS.md` document SHALL document the point-to-evidence-to-source chain as an ordered sequence from point to evidence to source, and SHALL include a reference to the `test_point_to_evidence_to_source_chain` test by name.
6. WHERE an analysis path is currently CLI-only, THE `WORKFLOWS.md` document SHALL describe the discordance catalog, predictive model, and PTM-parent audit paths and SHALL label each as CLI-only.
7. THE Requirement 1 deliverable SHALL consist only of the `WORKFLOWS.md` document and SHALL NOT add, remove, or alter any runtime source file.

### Requirement 2: Expose Discordance Deliverables via Read-Only API

**User Story:** As a frontend developer, I want read-only API endpoints for the discordance catalog, predictive model, and PTM-parent audit, so that the UI can render Track 2 deliverables without recomputing statistics.

#### Acceptance Criteria

1. WHEN a client issues an HTTP GET request to a discordance endpoint, THE Discordance_API SHALL return the committed demo output for that resource among the event catalog (catalog.csv), the model metrics (model_metrics.csv), the model predictions (model_predictions.csv), the run summary (run_summary.json), and the PTM-parent divergence audit (audit_ptm_parent outputs), without accepting or honoring any write, create, update, or delete operation.
2. WHEN the Discordance_API serves any committed demo output, THE Discordance_API SHALL return the values exactly as present in the committed source files without computing, aggregating, or deriving any statistic within the adapter.
3. IF a required committed demo output file is not present when a client requests it, THEN THE Discordance_API SHALL report the resource as unavailable via the `available()` guard, SHALL NOT return fabricated, placeholder, or partial data, and SHALL return a response indicating the resource is unavailable.
4. WHEN the Discordance_API returns any payload, THE Discordance_API SHALL serialize it as JSON in which every value is a JSON-representable type (string, number, boolean, null, object, or array) with no NaN, Infinity, or non-serializable objects.
5. THE Discordance_API SHALL NOT modify, regenerate, or replace the standalone discordance modules and CLI, which remain the sole source of truth for all statistics.
6. WHEN the catalog endpoint returns the event catalog, THE Discordance_API SHALL return each event classified as exactly one of the supported classes: supported_concordant, supported_opposite, rna_response_protein_equivalent, or indeterminate.
7. THE Finalization_Team SHALL add automated tests that assert, for each served resource, that the adapter-returned record counts and per-class counts are equal to the counts in the corresponding committed demo output files.
8. IF the committed catalog output contains an event whose classification label is missing, malformed, or not one of the four supported classes (supported_concordant, supported_opposite, rna_response_protein_equivalent, indeterminate), THEN THE Discordance_API SHALL report the catalog resource as unavailable via the `available()` guard and SHALL NOT return the catalog, rather than excluding the event or failing silently.

### Requirement 3: Build the Discordance Main View with Three-Stage Flow

**User Story:** As an analyst, I want discordance to be the app's main view organized as query builder → visualization → interpretation, so that I can move from a query to charts to classified discordance results in one guided flow.

#### Acceptance Criteria

1. WHEN the web application loads without a specific view requested, THE Web_Application SHALL display the Discordance view as the default view.
2. THE Discordance_View SHALL organize interaction as a three-stage flow presented in order: query builder, then data visualization, then discordance interpretation.
3. WHEN a discordance query completes successfully, THE Discordance_View SHALL render a count for each of the four event classes supported_concordant, supported_opposite, rna_response_protein_equivalent, and indeterminate, including classes whose count is zero.
4. WHEN a discordance query completes successfully, THE Discordance_View SHALL render, for each of the zero-change, RNA-only, and temporal model variants, its out-of-fold R2, MAE, and RMSE values.
5. WHEN the predictive-model metrics are rendered, THE Discordance_View SHALL display a visible text statement identifying the predictions as weak alongside the metrics.
6. WHEN a discordance query completes successfully, THE Discordance_View SHALL render the PTM-parent divergence summary.
7. THE Discordance_View SHALL display the PTM occupancy caveat as visible text.
8. WHEN the Live view displays results, THE Live_View SHALL render the `layer_discordance` data returned by the API.
9. IF a discordance query fails or returns no results, THEN THE Discordance_View SHALL display an error or empty-state indication and SHALL NOT render event-class counts, model metrics, or the PTM-parent summary.

### Requirement 4: Unbundle the Live Query Builder Against the Live API

**User Story:** As an analyst, I want the Live query builder wired to the live API with an explicit mapping-confirmation gate, so that I control the signature, target, and thresholds and never get silently mismatched mappings.

#### Acceptance Criteria

1. THE Live_Query_Builder SHALL allow the signature source to be selected as exactly one of: a built-in example, pasted text content, or an uploaded CSV file.
2. WHEN the signature source is a built-in example or pasted text content, THE Live_Query_Builder SHALL provide the signature to the API as signature_rows.
3. WHEN the signature source is an uploaded CSV file, THE Live_Query_Builder SHALL provide the signature to the API as signature_csv_text.
4. IF an uploaded CSV file is empty, exceeds the maximum accepted size, or cannot be parsed as CSV, THEN THE Live_Query_Builder SHALL reject the upload, retain any previously selected signature source unchanged, and display an error message indicating the file is invalid.
5. WHEN the target selectors are displayed, THE Live_Query_Builder SHALL allow selection of target species, omics, tissue, sex, and timepoint from values sourced from `/api/catalog`.
6. IF the request to `/api/catalog` fails entirely, THEN THE Live_Query_Builder SHALL disable all target selectors and display an error message indicating the catalog is unavailable.
7. IF the `/api/catalog` response provides values for some target dimensions but is missing or empty for a specific dimension (species, omics, tissue, sex, or timepoint), THEN THE Live_Query_Builder SHALL disable only the selector(s) lacking values and SHALL leave the selectors that have values enabled.
8. THE Live_Query_Builder SHALL allow selection of an FDR threshold constrained to a value between 0.0 (exclusive) and 1.0 (inclusive).
9. THE Finalization_Team SHALL rewire `QueryBuilder.tsx` to source catalog and comparison data from the live API instead of the in-memory catalog.
10. WHEN a comparison is requested, THE Live_Query_Builder SHALL call `/api/mappings/preview` and SHALL block the call to `/api/comparisons` until the previewed mappings are explicitly confirmed by the user.
11. IF the request to `/api/mappings/preview` fails, THEN THE Live_Query_Builder SHALL not call `/api/comparisons`, SHALL retain the current selections unchanged, and SHALL display an error message indicating the mapping preview is unavailable.
12. IF a signature term maps to more than one candidate in the `/api/mappings/preview` response, THEN THE Live_Query_Builder SHALL present all such candidates to the user for confirmation and SHALL NOT select any candidate automatically.
13. THE Live_Query_Builder SHALL treat the FDR threshold as an analysis parameter recorded in the run_id.
14. THE Live_Query_Builder SHALL treat the include-nonsignificant option as presentation-only and SHALL NOT record it in the run_id.

### Requirement 5: Generalized Query Upload and Results-Table Consistency

**User Story:** As an analyst, I want to upload or paste a generalized signature and see consistent, fully labeled results tables everywhere, so that no result is ambiguous about its statistics, species, dataset, or contrast.

#### Acceptance Criteria

1. WHEN an analyst submits an uploaded generalized signature that conforms to the query_core schema, THE Generalized_API SHALL accept it through the new endpoint and run the query.
2. IF an uploaded or pasted generalized signature does not conform to the query_core schema, THEN THE Generalized_API SHALL reject the submission without running the query and return an error response indicating which schema requirement was not met, while retaining any previously loaded signature unchanged.
3. THE Generalized_View SHALL provide both an upload control and a paste control for supplying a generalized signature, in addition to the bundled-signature picklist.
4. WHEN an analyst submits a pasted signature that cannot be parsed as a query_core signature, THE Generalized_View SHALL reject the submission and display an error indication identifying the parse failure, without clearing the analyst's pasted text.
5. THE Results_Tables SHALL, across the Live, Generalized, Discordance, and Metabolomics views, display the raw nominal p-value (camera_p) in a column immediately adjacent to the BH q-value (camera_fdr) column, with each column carrying a visible header label stating that the BH q-value is adjusted over the frozen 52-column multiplicity family and that the raw p-value is nominal (unadjusted).
6. THE Results_Tables SHALL include a distinctly labeled species column and a distinctly labeled dataset column.
7. THE Results_Tables SHALL render, alongside every displayed result value, the species, dataset, and contrast associated with that value.
8. THE Results_Tables SHALL NOT render any result value without an accompanying species, dataset, and contrast for that value.
9. WHERE a result maps a query feature to a target of a different species, THE Results_Tables SHALL display the ortholog-link relation for that result.

### Requirement 6: Retire Mock Design Views and Consolidate Navigation

**User Story:** As a user, I want the navigation to show only live analyses plus a single methods page, so that I am never shown static mock content presented as real results.

#### Acceptance Criteria

1. THE Web_Application SHALL exclude the mock design views 01, 02, 03, 04, and 08 from the `App.tsx` navigation such that no navigation entry, route, or link renders any of these five views.
2. IF a user navigates directly to a route corresponding to any of the mock design views 01, 02, 03, 04, or 08, THEN THE Web_Application SHALL NOT render the mock view content and SHALL redirect to the last-used live analysis view, or to the default Discordance view when no live analysis view has been used in the session.
3. THE Web_Application SHALL NOT display the hardcoded mock gene identifiers PPARGC1A, SOD2, and COL1A1 in any navigable view of the application.
4. THE Finalization_Team SHALL migrate the content of the removed mock views 01, 02, 03, 04, and 08 into documentation, reusing the mermaid diagrams defined in Requirement 1 without modification to their diagram source.
5. THE Finalization_Team SHALL preserve the Figma-Make design provenance in the migrated documentation such that each migrated mock view retains an attribution identifying it as Figma-Make design origin.
6. THE Web_Application navigation SHALL present exactly five entries: Discordance as the default view displayed on initial load, Live, Generalized, Metabolomics, and a single combined About or Methods page.
7. WHEN a user selects the About or Methods navigation entry, THE Web_Application SHALL display a link that navigates to the migrated documentation.

### Requirement 7: Rewrite the Root README Anchored to Track 2

**User Story:** As a new reader, I want the root README to lead with the Track 2 challenge and the discordance deliverables, so that I immediately understand the honest scope and where to find details.

#### Acceptance Criteria

1. THE Root_README SHALL name the Track 2 "Omic Discordance Explained" challenge within the first section (the content preceding the second heading) as its anchoring topic.
2. THE Root_README SHALL identify the discordance catalog and the predictive model as the headline deliverable, presenting them before any other deliverable is described.
3. THE Root_README SHALL state the framing that reported results are cross-cohort association findings and SHALL state that they are not a claim of PAH treatment efficacy.
4. THE Root_README SHALL state that cohorts are analyzed separately and are not merged, and SHALL state that the PTM occupancy caveat applies to the reported findings.
5. THE Root_README SHALL contain one navigable Markdown link to each of the following five targets: `WORKFLOWS.md`, `DECISIONS.md`, `VERIFICATION.md`, `METABOLOMICS.md`, and `RUN_LOCAL.md`.
6. THE Root_README SHALL contain no heading that is followed by empty content (no body text and no subordinate content before the next heading or end of file).

### Requirement 8: Fresh-Clone Hardening and Bootstrap

**User Story:** As a new contributor, I want a one-command bootstrap and clear setup docs, so that a fresh clone runs reliably without hidden manual steps.

#### Acceptance Criteria

1. WHEN the bootstrap script is invoked with a single command on a fresh clone, THE Finalization_Team SHALL perform editable installs of `hackathon/tool` and `apps/api`, fetch the mprobe store, and complete a smoke check within 600 seconds.
2. IF any editable install (`hackathon/tool` or `apps/api`) fails, THEN THE bootstrap script SHALL halt, exit with a non-zero status code, and emit an error message identifying the failed install target without performing subsequent steps.
3. IF fetching the mprobe store fails or the store is unavailable, THEN THE bootstrap script SHALL halt, exit with a non-zero status code, and emit an error message indicating the mprobe store could not be fetched.
4. WHEN the smoke check completes, THE bootstrap script SHALL exit with a zero status code if all checks pass and a non-zero status code if any check fails, emitting an error message identifying the failed check.
5. THE Bootstrap_Documentation SHALL be referenced by a resolvable link from both the README and `RUN_LOCAL.md`.
6. THE Finalization_Team SHALL configure exactly one package manager (either pnpm or npm) as the single declared package manager used in all project setup instructions and all CI workflow definitions, with no invocation of the non-selected package manager in either location.
7. THE Bootstrap_Documentation SHALL document the `parents[3]` path-coupling constraint, stating the required directory depth relationship between the adapter files and the referenced paths such that a reader can verify whether a given file location satisfies the constraint.

### Requirement 9: Automated Verification Wired into CI

**User Story:** As a maintainer, I want scripted API smoke checks and browser end-to-end tests wired into CI, so that live views and key values are verified automatically before merge.

#### Acceptance Criteria

1. THE Verification_Suite SHALL include a scripted httpx API smoke check that issues a request to every documented API endpoint and asserts each returns a success response.
2. THE httpx API smoke check SHALL assert the headline q value equals 0.0584 within a tolerance of plus or minus 0.0001, the discordance class counts match their expected per-class integer values, and the metabolomics result is the expected null result.
3. IF any endpoint in the httpx API smoke check returns a non-success response or any asserted value does not match its expected value, THEN THE Verification_Suite SHALL fail and indicate which endpoint or value failed.
4. THE Verification_Suite SHALL include a Playwright end-to-end test that boots the API and web application and asserts that each live view renders its expected real data rather than empty, placeholder, or error content.
5. THE Playwright end-to-end test SHALL assert the browser-level point-to-evidence-to-source chain by selecting a data point, verifying the linked evidence is displayed, and verifying the linked source is displayed.
6. THE Finalization_Team SHALL wire both the API smoke check and the Playwright end-to-end test into CI so that both run automatically on each merge request.
7. WHEN the full verification battery is run, THE Finalization_Team SHALL run the engine pytest suite, the `apps/api` tests, and the web `npm test` and build; and THE verification battery SHALL be considered successfully completed only when the engine suite reports exactly 122 passed and 34 skipped, the `apps/api` tests report zero failures, and the web tests and build report zero failures.
8. IF the engine suite counts differ from 122 passed and 34 skipped, OR the `apps/api` tests report any failure, OR the web tests or build report any failure, THEN THE verification battery SHALL NOT be marked complete and the merge SHALL be blocked.
9. IF the engine is modified, THEN THE Finalization_Team SHALL confirm `analysis.py` is byte-identical to the golden reference by comparing their content hashes.
10. THE Finalization_Team SHALL record a `VERIFICATION.md` checkpoint containing the exact commands run, the resulting test counts, and the content hashes for the verification battery.

### Requirement 10: Merge-to-Main Plan and Pre-Merge Checklist (Documentation Only)

**User Story:** As a maintainer, I want a documented merge-to-main plan and pre-merge checklist, so that the merge can be executed safely later without surprises.

#### Acceptance Criteria

1. THE Merge_Plan SHALL be documentation only and SHALL NOT trigger, invoke, or execute any git, gh, or shell command as part of this feature.
2. THE Merge_Plan SHALL require that all verification gates report a passing status, that `git status` reports zero modified, staged, untracked, or deleted files, and that the local branch has zero commits ahead of its upstream remote before merge is permitted.
3. IF any verification gate reports a non-passing status, OR the working tree contains one or more uncommitted or untracked files, OR the local branch is one or more commits ahead of its upstream remote, THEN THE Merge_Plan SHALL document that the merge is blocked and SHALL identify which precondition failed.
4. THE Merge_Plan SHALL require that no rebase is performed that replays commits over `origin/t3code/build-motrpac-probe-tool`.
5. THE Merge_Plan SHALL include the exact `gh pr create` command that sets the base branch to `main` and the head branch to `integration/exercise-signature-explorer`.
6. THE Merge_Plan SHALL include a pull request description template containing exactly three labeled sections: a summary section, a what-was-tested section, and a deferred-items section.

### Requirement 11: Source-of-Truth and Backend-Derivation Boundary (Cross-Cutting)

**User Story:** As a maintainer, I want the statistical source-of-truth boundary and backend-derived study selection enforced everywhere, so that the presentation layers never compute or misroute analysis.

#### Acceptance Criteria

1. THE FastAPI_Layer and React_Client SHALL NOT compute statistics; mprobe, query_core, and the standalone discordance modules SHALL remain the source of truth.
2. THE Backend SHALL derive study selection from species, mapping rat to the chronic study and human to the acute study.
3. THE Web_Application SHALL NOT present an acute-or-chronic input control.
4. THE Web_Application SHALL treat both contrasts as first-class and SHALL NOT present a "healthy gene set."
5. THE Web_Application SHALL show species, dataset, and contrast on every result.

### Requirement 12: Honest Framing, Multiplicity Discipline, and Guardrails (Cross-Cutting)

**User Story:** As a reader, I want honest framing and multiplicity discipline enforced in the UI and headlines, so that no result is overstated.

#### Acceptance Criteria

1. THE Directed_Path SHALL use the frozen 52-column BH multiplicity family.
2. THE Web_Application SHALL report the headline male rat SKM-GN protein 8-week result as q = 0.0584 and SHALL label it as not significant.
3. THE Web_Application SHALL NOT surface the 16-column q = 0.0413 result as a headline and SHALL present it only as a labeled sensitivity analysis.
4. THE Web_Application SHALL render the honest framing, guardrails, and PTM occupancy caveat in the user interface.
5. THE Web_Application SHALL keep cohorts separate, treating muscle protein, plasma metabolomics, and blood as distinct.

### Requirement 13: Traceability Discipline and Verify-Before-Claiming (Cross-Cutting)

**User Story:** As a maintainer, I want each feature traced through ADRs, verification checkpoints, and a traceability matrix, so that every claim is reproducible and every action has a recorded reason.

#### Acceptance Criteria

1. WHERE a feature is implemented, THE Finalization_Team SHALL add corresponding ADR entries to `DECISIONS.md`.
2. WHERE a feature is implemented, THE Finalization_Team SHALL record a `VERIFICATION.md` checkpoint with commands, counts, and hashes.
3. WHERE a feature is implemented, THE Finalization_Team SHALL add corresponding rows to `REQUIREMENTS_TRACEABILITY.md` mapping requirement to code path to test to status.
4. THE Finalization_Team SHALL NOT claim a count or behavior without first running the relevant test or command.
