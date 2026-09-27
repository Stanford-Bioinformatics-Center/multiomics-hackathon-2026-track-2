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
