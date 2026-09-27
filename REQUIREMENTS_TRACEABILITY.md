# Requirements Traceability Matrix

Each requirement → code path → test → status. Keeps the implementation interpretable and
traceable. Status: TODO / IN PROGRESS / DONE / BLOCKED.

| ID | Requirement | Code path | Test | Status |
|----|-------------|-----------|------|--------|
| R-ENG-1 | mprobe is the only scientific engine (no stats in API/React) | `hackathon/tool/src/motrpac_probe/*`; `apps/api` calls only | parity test (R-PAR-1) | IN PROGRESS |
| R-BRANCH-1 | Integration branch off tool branch, merge main | git history (`09ae0d7`) | VERIFICATION.md | DONE |
| R-IMP-1 | React imported via subtree, history preserved, LFS-independent | `apps/web/` | fresh-clone build (R-FC-1) | DONE |
| R-FIX-1 | Full-19 fixture enriched (paper+canonical ids, ratio, log2fc, p-text, corrections) | `examples/pah_muscle_malenfant2015.provenance.csv` (companion) | `test_pah_fixture_table2.py` | DONE |
| R-FIX-2 | Table-2 validation test for all 19 records | `tests/test_pah_fixture_table2.py` | pytest (8 tests) | DONE |
| R-FIX-3 | Lower-9 regression fixture unchanged | `examples/pah_muscle_lower9_malenfant2015.csv` | `test_lower9_regression_fixture_unchanged` | DONE |
| R-MAP-1 | Mapper exposes all candidates before selection (no silent first-wins) | `signature.candidates_for` / `mapping_audit` | `test_mapping_candidates.py` | DONE |
| R-MAP-2 | Mapping audit reconciles row counts before/after every join/collapse | `signature.mapping_audit` | `test_audit_conserves_rows_with_duplicates` | DONE |
| R-FIX-4 | Enrichment causes NO scientific drift (analysis.py byte-identical) | golden baseline diff | manual (VERIFICATION Gate 2) | DONE |
| R-SENS-1 | Predefined sensitivity views (full19/lower9/upper10/fibre); planned not post-hoc | `examples/pah_muscle_sensitivity_views.json` | `test_sensitivity_views_partition_and_membership` | DONE |
| R-MAN-1 | Committed mapping manifest for built-in fixture; uploads still need confirmation | `examples/*.mapping_manifest.{csv,json}` | manifest generation | DONE |
| R-SVC-1 | `run_analysis(request)` returns JSON-safe typed response (no DataFrames/Store/Timer) | `apps/api/.../service.py` | `test_runs_and_is_json_safe` | DONE |
| R-SVC-2 | `_path` KeyError fixed via content-addressed input; no client file paths | `service._materialize_signature` | `test_path_keyerror_is_fixed_provenance_present`, `test_signature_rows_input_path_is_content_addressed` | DONE |
| R-SVC-3 | Deterministic run_id (content+params+family+cols+mapping+store+code+seed+schema) | `service._run_id` | `test_run_id_deterministic_and_excludes_presentation_options` | DONE |
| R-PROV-1 | Provenance records family column IDs, family size, n_tests, method, threshold, per-feature p/q/effect | `run.provenance` + `service` | `test_multiplicity_family_recorded` | DONE |
| R-HEAD-1 (part) | Headline male rat SKM-GN prot 8wk = q=0.0584 (computed), not significant | `service._headline` | `test_headline_male_rat_skmgn_protein_8wk_not_significant` | DONE (service; UI in Gate 5) |
| R-RUN-1 (part) | FDR threshold in run_id; include-nonsignificant NOT in run_id | `service._run_id` | `test_run_id_deterministic_...` | DONE (service; UI in Gate 5) |
| R-LIN-1 (part) | evidence_id on every feature row; traces to a column | `service` FeatureEvidence | `test_features_have_lineage_and_trace_to_columns` | DONE (service; browser test Gate 5) |
| R-UI-2 (part) | Both contrasts first-class; never "healthy gene set" | `service` ContrastPair | `test_guardrails_and_contrasts_present` | DONE (service; UI in Gate 5) |
| R-CAT-1 | Catalog/availability from store.columns, not React constants; SKM-VL RNA-only fixed | `catalog.build_catalog` / `resolve_availability` | `test_catalog.py` (SKM-VL/SKM-GN) | DONE |
| R-CAT-2 | Unsupported source species (mouse/other) disabled or explicit unsupported | `catalog.SUPPORTED_SOURCE_SPECIES` | `test_source_species_gating`, `test_unknown_species_no_matching_context` | DONE |
| R-CAT-3 | Capability matrix (engine/api/react/demo) per layer | `catalog.CAPABILITY_MATRIX` | `test_capability_matrix_layers` | DONE |
| R-LIN-1 | evidence_id, mapping_decision_id, source feature, n_collapsed, aggregation method on every evidence row | `service` FeatureEvidence | `test_feature_lineage_fields_present_and_collapse_recorded`, `test_mapping_decision_id_ties_to_audit` | DONE (service; browser test Gate 5) |
| R-MULT-1 | 53-column BH family frozen; BH before filtering; display never re-runs BH | service | parity + q=0.0584 compute test | TODO |
| R-HEAD-1 | Headline male rat SKM-GN prot 8wk = q=0.0584 (computed), not significant; 16-col only as labeled sensitivity | `service._headline` + LiveDashboard | `test_service`/`test_api` headline; live verified | DONE |
| R-VIZ-1 | Visualization mode from RETURNED layers, not requested | `analysis.getVisualizationModeFromReturnedLayers` + LiveDashboard | `analysis.test.ts` (+2) | DONE |
| R-RUN-1 | FDR threshold in run_id; include-nonsignificant is presentation-only | `service._run_id` | `test_run_id_deterministic_...` | DONE |
| R-UI-1 | Guardrails Safe/Unsafe table + PTM warning render in UI | `LiveDashboard.tsx` | live view + `test_guardrails_and_contrasts_present` | DONE |
| R-UI-2 | Both contrasts first-class in req/resp/charts/export; never "healthy gene set" | `schema.ContrastPair` + LiveDashboard + export | `test_guardrails_and_contrasts_present` | DONE |
| R-EXP-1 | Export bundle = 9 files | `export.bundle_files` | `test_export_bundle_contents` | DONE |
| R-EMPTY-1 | No-compatible-data / zero-mapped → HTTP 200; malformed → 422 | `app.py` | `test_availability_empty_state_is_200`, `test_validate_malformed_is_422` | DONE |
| R-PAR-1 | Service matches engine cell-for-cell (feature counts, classes, adjusted p, store+git) | `test_parity.py` | 3 parity tests | DONE |
| R-PTS-1 | Every plotted point resolves to an exported evidence row → source | `test_api.py` | `test_point_to_evidence_to_source_chain` | DONE (API-level; browser e2e future) |
| R-FC-1 | Fresh-clone: install, build, tests pass; image + routes load without original repo/LFS | `.github/workflows/ci.yml` web job | CI | DONE |
| R-CI-1 | `pytest -q -rs` with expected-skip allowlist | `.github/workflows/ci.yml` python job | CI | DONE |
| R-GAL-1 | Gallery regen from clean commit records current SHA (stale 5e2cde3 documented) | `scripts_regenerate_gallery.md` | manual verify (a1dbfc4) | DONE (procedure; regen is its own commit) |
