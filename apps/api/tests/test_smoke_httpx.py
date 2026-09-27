"""R9 httpx API smoke check (task 12.1) — single-execution wiring + exact-value verification.

This is NOT a property test and NOT 100 iterations. It performs ONE pass that:

  * issues a request to every documented API endpoint (including the three discordance
    endpoints and the generalized upload endpoint) and asserts each returns a success
    response (R9 AC1);
  * asserts the headline `q == 0.0584 ± 0.0001` and that it is labeled not-significant
    (R12.2 / R9 AC2);
  * asserts the discordance per-class integer counts:
    supported_concordant=1, supported_opposite=0, rna_response_protein_equivalent=285,
    indeterminate=5642 (R9 AC2);
  * asserts the metabolomics null result (0% exercise-sensitive hits, Fisher p 0.609)
    exactly as committed in `metab_casestudy` / `test_metabolomics.py` (R9 AC2);
  * on any non-success response OR any value mismatch, FAILS and NAMES the offending
    endpoint/value (R9 AC3).

Transport: the check runs the FastAPI app IN-PROCESS via `httpx.ASGITransport`, so it needs
no live server (matching the in-process pattern the rest of `apps/api/tests` uses). It requires
the mprobe store (for the directed `/api/comparisons` headline) and the committed demo outputs
for the discordance/metabolomics/generalized adapters; when those are absent the check skips
with an explicit reason rather than fabricating a pass.

Run it directly (standalone, single execution):
    python -m apps.api.tests.test_smoke_httpx        # from repo root, or
    python apps/api/tests/test_smoke_httpx.py
or as part of the pytest battery wired into CI:
    python -m pytest apps/api/tests/test_smoke_httpx.py
"""
from __future__ import annotations

import asyncio
import csv
import io
import sys

# --------------------------------------------------------------------------- #
# Expected committed values (source-of-truth: the engines + committed demos,
# cross-checked against apps/api/tests/{test_api,test_discordance,test_metabolomics}.py).
# --------------------------------------------------------------------------- #
_HEADLINE_Q = 0.0584
_HEADLINE_TOL = 1e-4  # "0.0584 ± 0.0001"

_EXPECTED_DISCORDANCE_CLASS_COUNTS = {
    "supported_concordant": 1,
    "supported_opposite": 0,
    "rna_response_protein_equivalent": 285,
    "indeterminate": 5642,
}

# The metabolomics null result: 0% blood exercise-sensitive hits, Fisher p 0.609
# (metab_casestudy.build_case_study().null_result["blood_exercise_sensitive"]).
_EXPECTED_METAB_HITS_FRACTION = 0.0
_EXPECTED_METAB_FISHER_P = 0.609
_EXPECTED_METAB_FISHER_TOL = 1e-6

# A minimal disease signature that conforms to the query_core schema (11 required
# DISEASE_COLUMNS), so POST /api/generalized/query runs the query and returns 200.
_QUERY_CORE_DISEASE_COLUMNS = (
    "study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
    "log2_fc", "statistic", "p_value", "q_value", "contrast",
)


class SmokeFailure(AssertionError):
    """Raised when an endpoint returns non-success or an asserted value mismatches.

    The message always NAMES the offending endpoint and/or value (R9 AC3).
    """


def _conforming_signature_csv_text() -> str:
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow(_QUERY_CORE_DISEASE_COLUMNS)
    w.writerow([
        "StudyX", "human", "muscle", "protein", "Q16795", "UniProt",
        "-0.49", "", "0.004", "", "resting PAH minus control",
    ])
    return buf.getvalue()


def _require_success(resp, endpoint: str) -> dict | None:
    """Assert a 2xx response, NAMING the endpoint on failure. Returns parsed JSON (or None)."""
    if resp.status_code >= 400:
        raise SmokeFailure(
            f"endpoint {endpoint!r} returned non-success status "
            f"{resp.status_code}: {resp.text[:300]!r}"
        )
    try:
        return resp.json()
    except ValueError:
        return None


def _check_value(condition: bool, endpoint: str, value_name: str, detail: str) -> None:
    """Assert an exact-value condition, NAMING the endpoint + value on mismatch (R9 AC3)."""
    if not condition:
        raise SmokeFailure(
            f"value mismatch at endpoint {endpoint!r}, value {value_name!r}: {detail}"
        )


async def _run_smoke() -> list[str]:
    """Execute the single-pass smoke check. Returns a list of human-readable check lines.

    Raises SmokeFailure (naming endpoint/value) on the first non-success or mismatch.
    """
    import httpx

    from motrpac_probe_service.app import app

    log: list[str] = []
    transport = httpx.ASGITransport(app=app)
    async with httpx.AsyncClient(transport=transport, base_url="http://smoke") as client:

        # ---- 1. health -------------------------------------------------------------
        ep = "GET /api/health"
        body = _require_success(await client.get("/api/health"), ep)
        _check_value(body.get("status") == "ok", ep, "status",
                     f"expected 'ok', got {body.get('status')!r}")
        log.append(f"OK  {ep} -> status=ok")

        # ---- 2. catalog ------------------------------------------------------------
        ep = "GET /api/catalog"
        body = _require_success(await client.get("/api/catalog"), ep)
        _check_value("analysis_types" in body, ep, "analysis_types",
                     "missing analysis_types in catalog payload")
        log.append(f"OK  {ep} -> analysis_types present")

        # ---- 3. availability (scientific empty state is a 200) ---------------------
        ep = "POST /api/availability"
        _require_success(
            await client.post("/api/availability", json={
                "target_species": "rat", "selected_omics": ["transcriptomics"]}), ep)
        log.append(f"OK  {ep}")

        # ---- 4. mappings preview ---------------------------------------------------
        ep = "POST /api/mappings/preview"
        _require_success(
            await client.post("/api/mappings/preview",
                              json={"example_name": "pah_muscle_malenfant2015"}), ep)
        log.append(f"OK  {ep}")

        # ---- 5. signatures validate ------------------------------------------------
        ep = "POST /api/signatures/validate"
        _require_success(
            await client.post("/api/signatures/validate",
                              json={"example_name": "pah_muscle_malenfant2015"}), ep)
        log.append(f"OK  {ep}")

        # ---- 6. directed comparison + headline q (R12.2 / R9 AC2) ------------------
        ep = "POST /api/comparisons"
        body = _require_success(
            await client.post("/api/comparisons",
                              json={"example_name": "pah_muscle_malenfant2015",
                                    "target_species": "rat"}), ep)
        _check_value(body.get("status") == "ok", ep, "status",
                     f"expected 'ok', got {body.get('status')!r}")
        run_id = body["run_id"]
        headline = body["headline"]["male_rat_skm_gn_protein_8wk"]
        q = headline["camera_fdr_family"]
        _check_value(
            abs(q - _HEADLINE_Q) <= _HEADLINE_TOL, ep, "headline q (camera_fdr_family)",
            f"expected {_HEADLINE_Q} ± {_HEADLINE_TOL}, got {q}")
        _check_value(
            headline["significant_at_threshold"] is False, ep, "headline significance",
            "headline must be labeled NOT significant (R12.2)")
        log.append(f"OK  {ep} -> headline q={q:.4f} (±{_HEADLINE_TOL}), not significant; "
                   f"run_id={run_id}")

        # ---- 7. run subresources ---------------------------------------------------
        for sub in ("", "/features", "/provenance", "/report", "/export"):
            ep = f"GET /api/comparisons/{{run_id}}{sub}"
            _require_success(await client.get(f"/api/comparisons/{run_id}{sub}"), ep)
            log.append(f"OK  {ep}")

        # ---- 8. metabolomics case study + NULL RESULT (R9 AC2) ---------------------
        ep = "GET /api/metabolomics/casestudy"
        body = _require_success(await client.get("/api/metabolomics/casestudy"), ep)
        _check_value(body.get("analysis_type") == "metabolomics_case_study", ep,
                     "analysis_type",
                     f"expected 'metabolomics_case_study', got {body.get('analysis_type')!r}")
        es = body["null_result"]["blood_exercise_sensitive"]
        _check_value(
            es["hits_fraction"] == _EXPECTED_METAB_HITS_FRACTION, ep,
            "null_result.blood_exercise_sensitive.hits_fraction",
            f"expected {_EXPECTED_METAB_HITS_FRACTION} (0% exercise-sensitive), "
            f"got {es['hits_fraction']}")
        _check_value(
            abs(es["fisher_p"] - _EXPECTED_METAB_FISHER_P) < _EXPECTED_METAB_FISHER_TOL, ep,
            "null_result.blood_exercise_sensitive.fisher_p",
            f"expected {_EXPECTED_METAB_FISHER_P}, got {es['fisher_p']}")
        log.append(f"OK  {ep} -> null result: hits_fraction={es['hits_fraction']}, "
                   f"fisher_p={es['fisher_p']}")

        ep = "GET /api/metabolomics/casestudy/export"
        _require_success(await client.get("/api/metabolomics/casestudy/export"), ep)
        log.append(f"OK  {ep}")

        ep = "GET /api/metabolomics/convergence"
        _require_success(await client.get("/api/metabolomics/convergence"), ep)
        log.append(f"OK  {ep}")

        # ---- 9. generalized: bundled list, bundled query, and UPLOAD (R9 AC1) ------
        ep = "GET /api/generalized/signatures"
        _require_success(await client.get("/api/generalized/signatures"), ep)
        log.append(f"OK  {ep}")

        ep = "GET /api/generalized/query/pah_blood_rna"
        body = _require_success(await client.get("/api/generalized/query/pah_blood_rna"), ep)
        _check_value(body.get("analysis_type") == "generalized_query", ep, "analysis_type",
                     f"expected 'generalized_query', got {body.get('analysis_type')!r}")
        log.append(f"OK  {ep} -> analysis_type=generalized_query")

        ep = "POST /api/generalized/query"
        body = _require_success(
            await client.post("/api/generalized/query", json={
                "signature_csv_text": _conforming_signature_csv_text(),
                "tissue": "muscle",
                "reference_contrast_category": "EE-CON",
                "reference": "default",
            }), ep)
        _check_value(body.get("analysis_type") == "generalized_query", ep, "analysis_type",
                     f"expected 'generalized_query', got {body.get('analysis_type')!r}")
        log.append(f"OK  {ep} (uploaded conforming signature) -> analysis_type=generalized_query")

        # ---- 10. discordance: catalog (per-class counts), model, ptm-parent --------
        ep = "GET /api/discordance/catalog"
        body = _require_success(await client.get("/api/discordance/catalog"), ep)
        counts = body["classification_counts"]
        for cls, expected in _EXPECTED_DISCORDANCE_CLASS_COUNTS.items():
            actual = counts.get(cls)
            _check_value(
                actual == expected, ep, f"classification_counts.{cls}",
                f"expected {expected}, got {actual}")
            _check_value(
                isinstance(actual, int), ep, f"classification_counts.{cls} type",
                f"expected an integer count, got {type(actual).__name__} ({actual!r})")
        log.append(f"OK  {ep} -> classification_counts={counts}")

        ep = "GET /api/discordance/model"
        body = _require_success(await client.get("/api/discordance/model"), ep)
        _check_value(len(body.get("metrics", [])) == 6, ep, "metrics",
                     f"expected 6 model-metric rows, got {len(body.get('metrics', []))}")
        log.append(f"OK  {ep} -> 6 metric rows")

        ep = "GET /api/discordance/ptm-parent"
        body = _require_success(await client.get("/api/discordance/ptm-parent"), ep)
        _check_value(body.get("analysis_type") == "ptm_parent_audit", ep, "analysis_type",
                     f"expected 'ptm_parent_audit', got {body.get('analysis_type')!r}")
        log.append(f"OK  {ep} -> analysis_type=ptm_parent_audit")

    return log


def run_smoke() -> list[str]:
    """Synchronous entry point: run the single-pass smoke check and return the check log."""
    return asyncio.run(_run_smoke())


# --------------------------------------------------------------------------- #
# pytest wrapper — runs the single smoke pass as one test. Skips (never fabricates
# a pass) when the mprobe store or a required committed demo is absent.
# --------------------------------------------------------------------------- #
def _skip_reason() -> str | None:
    try:
        from motrpac_probe.paths import CONTRASTS
    except Exception as exc:  # engine not installed
        return f"mprobe engine not importable: {exc}"
    if not CONTRASTS.exists():
        return "no mprobe store: run `mprobe store fetch`"
    from motrpac_probe_service import discordance_casestudy as _dc
    from motrpac_probe_service import generalized_query as _gq
    from motrpac_probe_service import metab_casestudy as _mcs
    if not _dc.catalog_available():
        return "discordance demo outputs not present"
    if not _mcs.available():
        return "metabolomics module outputs not present"
    if not _gq.available():
        return "generalized query_core module not present"
    return None


def test_api_smoke_single_pass():
    import pytest

    reason = _skip_reason()
    if reason:
        pytest.skip(reason)
    log = run_smoke()
    # Every documented endpoint + every exact value passed; surface the trail for the record.
    print("\n".join(log))
    assert log, "smoke check produced no results"


# --------------------------------------------------------------------------- #
# Standalone execution (single pass, no pytest). Prints each check and exits
# non-zero with the offending endpoint/value named on any failure.
# --------------------------------------------------------------------------- #
def _main() -> int:
    reason = _skip_reason()
    if reason:
        print(f"SKIP: {reason}")
        return 0
    try:
        log = run_smoke()
    except SmokeFailure as exc:
        print(f"SMOKE FAILED: {exc}")
        return 1
    print("\n".join(log))
    print(f"\nSMOKE PASSED: {len(log)} checks across all documented endpoints.")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
