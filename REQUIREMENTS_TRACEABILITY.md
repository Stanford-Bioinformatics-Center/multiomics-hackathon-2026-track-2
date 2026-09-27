# Requirements Traceability Matrix

Each requirement → code path → test → status. Keeps the implementation interpretable and
traceable. Status: TODO / IN PROGRESS / DONE / BLOCKED.

| ID | Requirement | Code path | Test | Status |
|----|-------------|-----------|------|--------|
| R-ENG-1 | mprobe is the only scientific engine (no stats in API/React) | `hackathon/tool/src/motrpac_probe/*`; `apps/api` calls only | parity test (R-PAR-1) | IN PROGRESS |
| R-BRANCH-1 | Integration branch off tool branch, merge main | git history (`09ae0d7`) | VERIFICATION.md | DONE |
| R-IMP-1 | React imported via subtree, history preserved, LFS-independent | `apps/web/` | fresh-clone build (R-FC-1) | IN PROGRESS |
| R-FIX-1 | Full-19 fixture enriched (paper+canonical ids, ratio, log2fc, p-text, corrections) | `hackathon/tool/examples/pah_muscle_malenfant2015.csv` (+companion) | Table-2 validation (R-FIX-2) | TODO |
| R-FIX-2 | Table-2 validation test for all 19 records | `hackathon/tool/tests/` | pytest | TODO |
| R-FIX-3 | Lower-9 regression fixture unchanged | `examples/pah_muscle_lower9_malenfant2015.csv` | golden test | TODO |
| R-MAP-1 | Mapper exposes all candidates before selection (no silent first-wins) | `signature.py` mapping | mapping-audit test | TODO |
| R-MAP-2 | Mapping audit reconciles row counts before/after every join/collapse | audit module | reconciliation test | TODO |
| R-SVC-1 | `run_analysis(request)` returns JSON-safe typed response (no DataFrames/Store/Timer) | `apps/api` service | schema test | TODO |
| R-SVC-2 | `_path` KeyError fixed via content-addressed input; no client file paths | service adapter | unit test | TODO |
| R-SVC-3 | Deterministic run_id (content+params+family+cols+mapping+store+code+seed+schema) | service | unit test | TODO |
| R-PROV-1 | Provenance records all 53 column IDs, family size, n_tests, BH n, method, threshold, per-feature p/q/effect | `run.provenance` + service | provenance-completeness test | TODO |
| R-CAT-1 | Catalog/availability from store.columns, not React constants; SKM-VL RNA-only fixed | `apps/api` catalog | availability test | TODO |
| R-CAT-2 | Unsupported source species (mouse/other) disabled or explicit unsupported | catalog + UI | unit test | TODO |
| R-LIN-1 | evidence_id, mapping_decision_id, selected+candidate features, aggregation method on every evidence row | service | point→source test (R-PTS-1) | TODO |
| R-MULT-1 | 53-column BH family frozen; BH before filtering; display never re-runs BH | service | parity + q=0.0584 compute test | TODO |
| R-HEAD-1 | Headline male rat SKM-GN prot 8wk = q=0.0584 (computed), not significant; 16-col only as labeled sensitivity | service + UI | q-compute test (R-MULT-1) | TODO |
| R-VIZ-1 | Visualization mode from RETURNED layers, not requested | React | component test | TODO |
| R-RUN-1 | FDR threshold in run_id; include-nonsignificant is presentation-only | service + React | unit test | TODO |
| R-UI-1 | Guardrails Safe/Unsafe table + PTM warning render in UI | React | component test | TODO |
| R-UI-2 | Both contrasts (IPAH/control, trained/sedentary) first-class in req/resp/charts/export; never "healthy gene set" | schema + React | test | TODO |
| R-EXP-1 | Export bundle = input_signature, mapping_audit, feature_evidence, layer_comparison, summary, analysis_parameters, provenance, methods_and_limitations, report.html | service | export test | TODO |
| R-EMPTY-1 | No-compatible-data / zero-mapped → HTTP 200 scientific-empty; malformed → 422 | API | contract test | TODO |
| R-PAR-1 | React/service/Marimo parity (feature count, classes, RNA/PROT states, adjusted p, store+git) | tests | parity test | TODO |
| R-PTS-1 | Every plotted point resolves to an exported evidence row | browser test | e2e | TODO |
| R-FC-1 | Fresh-clone: install, build, tests pass; image + routes load without original repo/LFS | CI | CI job | TODO |
| R-CI-1 | `pytest -q -rs` with expected-skip allowlist | CI | CI job | TODO |
| R-GAL-1 | Gallery regenerated from clean commit; provenance records dirty state | `hackathon/tool/site` | build check | TODO |
