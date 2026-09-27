# apps/api — Exercise Signature Explorer backend adapter

A thin, UI-independent layer over the `motrpac_probe` scientific engine.

- `motrpac_probe_service/service.py` — `run_analysis(request) -> AnalysisResponse`. Wraps
  `motrpac_probe.run.compute()`, returns JSON-safe typed data (no DataFrames, Store, Timer,
  NaN, or NumPy scalars), and builds full provenance including the frozen multiplicity family.
- `motrpac_probe_service/schema.py` — dataclasses for the request/response contract (the Pydantic
  models and OpenAPI spec in Gate 5 are generated from / converge to these).

Design rules (DECISIONS.md):
- No statistics here. All numbers come from `motrpac_probe`.
- Study design is derived from the target species in the engine's catalog (rat→chronic,
  human→acute); there is no acute/chronic input.
- Multiplicity family = the engine's full core-column BH set (52 columns for the MVP). BH is
  applied across the whole family by the engine *before* any display filtering; the service never
  re-runs BH. All family column IDs and `n_tests` are recorded in provenance.
- Input is content-addressed (a sha256 of the signature bytes); the service never accepts an
  arbitrary server file path from a client (fixes the `_path` KeyError and closes a path-traversal
  hole).

Run the tests: `python -m pytest apps/api/tests -q` (needs the mprobe store; `mprobe store fetch`).
