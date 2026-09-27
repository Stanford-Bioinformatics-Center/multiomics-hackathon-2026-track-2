"""Generalized query_core adapter + endpoints: the third analysis type."""
import json

import pytest

from motrpac_probe_service import generalized_query as gq

pytestmark = pytest.mark.skipif(not gq.available(),
                                reason="generalized query_core module not present")


@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient
    from motrpac_probe_service.app import app
    return TestClient(app)


def test_lists_bundled_signatures():
    sigs = {s["id"] for s in gq.list_signatures()["signatures"]}
    assert "pah_blood_rna" in sigs and "pah_muscle_protein" in sigs
    assert len(sigs) == 5


def test_blood_query_reproduces_legacy_spearman():
    """The generalized engine must reproduce the hardcoded pipeline's six blood Spearman values."""
    r = gq.run_bundled("pah_blood_rna")
    json.dumps(r, allow_nan=False)  # JSON-safe
    assert r["analysis_type"] == "generalized_query"
    rna = {x["timepoint"]: x["spearman_rho"] for x in r["rank_correlation"] if x["layer"] == "rna"}
    # legacy values from the hardcoded pipeline
    legacy = {"during_20_min": -0.0788, "during_40_min": -0.0286, "post_10_min": 0.0601,
              "post_15_30_45_min": 0.1247, "post_3.5_4_hr": 0.1161, "post_24_hr": 0.0899}
    for t, v in legacy.items():
        assert abs(rna[t] - v) < 1e-3, f"{t}: {rna[t]} != legacy {v}"


def test_query_provenance_present():
    r = gq.run_bundled("pah_muscle_protein")
    assert r["provenance"]["disease_sha256"] and r["provenance"]["reference_sha256"]
    assert "source of truth" in r["provenance"]["engine"].lower()
    assert r["counts"]["disease_rows"] == 9  # 9-protein muscle signature


def test_api_generalized_endpoints(client):
    sigs = client.get("/api/generalized/signatures")
    assert sigs.status_code == 200 and len(sigs.json()["signatures"]) == 5

    q = client.get("/api/generalized/query/pah_blood_rna")
    assert q.status_code == 200
    body = q.json()
    assert body["analysis_type"] == "generalized_query"
    rna = [x for x in body["rank_correlation"] if x["layer"] == "rna" and x["timepoint"] == "during_20_min"]
    assert rna and abs(rna[0]["spearman_rho"] - (-0.0788)) < 1e-3

    assert client.get("/api/generalized/query/does_not_exist").status_code == 404


def test_catalog_exposes_generalized_analysis_type(client):
    at = client.get("/api/catalog").json()["analysis_types"]
    assert "generalized_query" in at
    assert at["generalized_query"]["api"] and at["generalized_query"]["react"]

# --------------------------------------------------------------------------------------
# Property 6 — generalized upload schema enforcement (test-first for tasks 7.2 / 7.3)
# Feature: finalize-discordance-mvp, Property 6
#
# The endpoint runs the query IFF the signature conforms to the query_core schema:
#   * conforming input  -> query_core.engine.run_query is called exactly once, returns a stub
#   * non-conforming     -> query is NEVER run, an error names the unmet schema requirement,
#                           and any previously loaded signature is left unchanged.
#
# This targets the intended `run_uploaded(signature_csv_text, *, tissue,
# reference_contrast_category, reference="default")` API plus `_validate_query_core_schema`
# returning (ok, failing_requirement) and raising `SchemaError(which)` on failure — none of
# which exist yet. It is EXPECTED to fail/error until tasks 7.2/7.3 implement them.
#
# query_core is source of truth; we mock `query_core.engine.run_query` so iterations stay cheap
# and so a "run" is observable without touching real reference data.
# --------------------------------------------------------------------------------------
import io
import csv as _csv

from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

# The 11 columns query_core requires of a *disease* signature (schema.DISEASE_COLUMNS).
_REQUIRED_DISEASE_COLUMNS = (
    "study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
    "log2_fc", "statistic", "p_value", "q_value", "contrast",
)


def _conforming_csv_text() -> str:
    """A minimal disease signature that satisfies the query_core schema."""
    buf = io.StringIO()
    w = _csv.writer(buf)
    w.writerow(_REQUIRED_DISEASE_COLUMNS)
    # one valid row: human/muscle/protein, UniProt id, finite numbers, p/q in [0,1]
    w.writerow([
        "StudyX", "human", "muscle", "protein", "Q16795", "UniProt",
        "-0.49", "", "0.004", "", "resting PAH minus control",
    ])
    return buf.getvalue()


# Strategies -----------------------------------------------------------------------------

_conforming = st.just(_conforming_csv_text())


@st.composite
def _nonconforming(draw):
    """Signature CSV text that violates the query_core schema in a structural way.

    Every variant here makes `query_core.schema.load_csv(..., reference=False)` raise, so the
    validator must reject it and name the unmet requirement.
    """
    kind = draw(st.sampled_from(["missing_column", "empty", "not_csv", "bad_species", "bad_layer"]))

    if kind == "empty":
        return ""
    if kind == "not_csv":
        # a non-tabular blob with no recognizable header row
        return draw(st.text(alphabet="abc123 ", min_size=1, max_size=20))

    if kind == "missing_column":
        drop = draw(st.sampled_from(_REQUIRED_DISEASE_COLUMNS))
        cols = [c for c in _REQUIRED_DISEASE_COLUMNS if c != drop]
        buf = io.StringIO()
        w = _csv.writer(buf)
        w.writerow(cols)
        w.writerow(["v"] * len(cols))
        return buf.getvalue()

    # header present + one row, but a required biological axis is invalid
    buf = io.StringIO()
    w = _csv.writer(buf)
    w.writerow(_REQUIRED_DISEASE_COLUMNS)
    species = "martian" if kind == "bad_species" else "human"
    layer = "protein" if kind == "bad_species" else "quux"
    w.writerow([
        "StudyX", species, "muscle", layer, "Q16795", "UniProt",
        "-0.49", "", "0.004", "", "resting PAH minus control",
    ])
    return buf.getvalue()


def _install_run_query_spy(monkeypatch):
    """Patch query_core.engine.run_query with a counting stub; return the call counter list."""
    import importlib
    import sys

    p = str(gq.GENERALIZED_DIR)
    if p not in sys.path:
        sys.path.insert(0, p)
    engine = importlib.import_module("query_core.engine")

    calls = []

    def _stub(**kwargs):
        calls.append(kwargs)
        # a minimal summary shaped like the real engine's return so normalization can proceed
        return {"counts": {}, "filters": {}, "thresholds": {}, "interpretation": "",
                "disease_input": "stub", "disease_sha256": "0" * 64,
                "reference_input": "stub", "reference_sha256": "0" * 64,
                "reference_releases": []}

    monkeypatch.setattr(engine, "run_query", _stub, raising=True)
    return calls


@settings(max_examples=100, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(conforming=_conforming, nonconforming=_nonconforming())
def test_property6_upload_schema_enforcement(monkeypatch, conforming, nonconforming):
    """Property 6: run_uploaded runs query_core IFF the signature conforms; else rejects
    it, names the unmet schema requirement, runs nothing, and leaves prior state unchanged.

    Validates: Requirements 5.1, 5.2
    """
    SchemaError = getattr(gq, "SchemaError")  # fails until task 7.2 defines it

    # --- validator agrees with the engine's own schema expectations ---
    ok_c, which_c = gq._validate_query_core_schema(conforming)
    assert ok_c is True and which_c in (None, "")
    ok_n, which_n = gq._validate_query_core_schema(nonconforming)
    assert ok_n is False
    assert isinstance(which_n, str) and which_n.strip(), "must name the unmet schema requirement"

    # --- conforming input: query runs exactly once and a normalized payload comes back ---
    calls = _install_run_query_spy(monkeypatch)
    result = gq.run_uploaded(conforming, tissue="muscle",
                             reference_contrast_category="EE-CON")
    assert len(calls) == 1, "conforming input must run the query exactly once"
    assert result["analysis_type"] == "generalized_query"

    # --- non-conforming input: NEVER runs, raises naming the requirement, prior state intact ---
    calls_before = list(calls)
    with pytest.raises(SchemaError) as exc:
        gq.run_uploaded(nonconforming, tissue="muscle",
                        reference_contrast_category="EE-CON")
    assert calls == calls_before, "non-conforming input must not run the query"
    assert which_n in str(exc.value) or which_n in repr(exc.value), \
        "the raised error must identify the unmet schema requirement"


@settings(max_examples=100, deadline=None,
          suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(nonconforming=_nonconforming())
def test_property6_endpoint_rejects_nonconforming(monkeypatch, client, nonconforming):
    """Property 6 (endpoint): POST /api/generalized/query returns 422 naming the failing
    schema requirement for non-conforming input, and never runs the query.

    Validates: Requirements 5.1, 5.2
    """
    calls = _install_run_query_spy(monkeypatch)
    resp = client.post("/api/generalized/query", json={
        "signature_csv_text": nonconforming,
        "tissue": "muscle",
        "reference_contrast_category": "EE-CON",
        "reference": "default",
    })
    assert resp.status_code == 422, "non-conforming upload must be rejected with 422"
    assert calls == [], "non-conforming upload must not run the query"
    # the response body must identify which schema requirement was unmet
    assert resp.text.strip(), "422 response must name the unmet schema requirement"
