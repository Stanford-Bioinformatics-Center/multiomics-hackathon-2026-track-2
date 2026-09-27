"""API contract tests: endpoints, HTTP-200 scientific-empty vs 422 malformed, export bundle."""
import io
import zipfile

import pytest

pytestmark = pytest.mark.skipif(
    __import__("motrpac_probe.paths", fromlist=["CONTRASTS"]).CONTRASTS.exists() is False,
    reason="no mprobe store: run `mprobe store fetch`")


@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient
    from motrpac_probe_service.app import app
    return TestClient(app)


def test_health(client):
    r = client.get("/api/health")
    assert r.status_code == 200 and r.json()["status"] == "ok"
    assert r.json()["store_hash"]


def test_catalog_endpoint(client):
    cat = client.get("/api/catalog").json()
    rat = next(c for c in cat["contexts"] if c["species"] == "rat")
    skmvl = next(t for t in rat["tissues"] if t["tissue"] == "SKM-VL")
    assert skmvl["layers"] == ["RNA"]           # store-derived, not React constants
    assert "capability_matrix" in cat


def test_availability_empty_state_is_200(client):
    r = client.post("/api/availability", json={"target_species": "rat",
                                                "selected_omics": ["transcriptomics"], "tissue": "NOPE"})
    assert r.status_code == 200                  # scientific empty state, not an error
    assert r.json()["status"] == "no_matching_context"


def test_validate_malformed_is_422(client):
    # a CSV with no id column and no direction column is malformed input
    r = client.post("/api/signatures/validate", json={"signature_csv_text": "foo,bar\n1,2\n"})
    assert r.status_code == 422


def test_mappings_preview_flags_confirmation(client):
    r = client.post("/api/mappings/preview", json={"example_name": "pah_muscle_malenfant2015"})
    assert r.status_code == 200
    body = r.json()
    assert body["counts"]["n_input_rows"] == 25
    assert body["requires_confirmation"] is False   # built-in fixture is unambiguous


def test_comparisons_and_subresources(client):
    r = client.post("/api/comparisons", json={"example_name": "pah_muscle_malenfant2015",
                                              "target_species": "rat"})
    assert r.status_code == 200
    body = r.json()
    assert body["status"] == "ok"
    run_id = body["run_id"]
    # headline computed, not significant
    h = body["headline"]["male_rat_skm_gn_protein_8wk"]
    assert abs(h["camera_fdr_family"] - 0.0584) < 5e-4 and h["significant_at_threshold"] is False
    # subresources
    assert client.get(f"/api/comparisons/{run_id}").status_code == 200
    feats = client.get(f"/api/comparisons/{run_id}/features").json()["features"]
    assert feats and all(f["evidence_id"].startswith(run_id) for f in feats[:20])
    prov = client.get(f"/api/comparisons/{run_id}/provenance").json()
    assert prov["run_id"] == run_id and prov["multiplicity_family"]["family_size"] >= 40
    rep = client.get(f"/api/comparisons/{run_id}/report")
    assert rep.status_code == 200 and "not evidence that exercise treats PAH" in rep.text


def test_export_bundle_contents(client):
    run_id = client.post("/api/comparisons", json={"example_name": "pah_muscle_malenfant2015",
                                                   "target_species": "rat"}).json()["run_id"]
    r = client.get(f"/api/comparisons/{run_id}/export")
    assert r.status_code == 200 and r.headers["content-type"] == "application/zip"
    z = zipfile.ZipFile(io.BytesIO(r.content))
    names = set(z.namelist())
    expected = {"input_signature.csv", "mapping_audit.csv", "feature_evidence.csv", "layer_comparison.csv",
                "summary.json", "analysis_parameters.json", "provenance.json",
                "methods_and_limitations.md", "report.html"}
    assert names == expected, f"bundle missing/extra: {names ^ expected}"
    # feature_evidence.csv carries lineage columns
    fe = z.read("feature_evidence.csv").decode()
    assert "source_feature_id" in fe and "aggregation_method" in fe and "evidence_id" in fe


def test_comparisons_unknown_example_422(client):
    r = client.post("/api/comparisons", json={"example_name": "does_not_exist", "target_species": "rat"})
    assert r.status_code == 422


def test_request_requires_exactly_one_signature_source(client):
    r = client.post("/api/comparisons", json={"target_species": "rat"})  # no signature at all
    assert r.status_code == 422


def test_point_to_evidence_to_source_chain(client):
    """Every plotted point (a column in the dashboard table, and each feature within it) must resolve
    to an exported evidence row and then to a real store feature_id. This is the API-level guarantee
    behind the browser 'every point traces to an exported row' requirement."""
    from motrpac_probe import core
    S = core.get_store()

    run_id = client.post("/api/comparisons", json={"example_name": "pah_muscle_malenfant2015",
                                                   "target_species": "rat"}).json()["run_id"]
    body = client.get(f"/api/comparisons/{run_id}").json()
    feats = client.get(f"/api/comparisons/{run_id}/features").json()["features"]

    col_ids = {c["column_id"] for c in body["columns"]}
    # every feature traces to a shown column, carries a run-scoped evidence_id, and (if measured) a
    # source_feature_id that actually exists in that column's store rows.
    store_features_by_col = {cid: set(S.frame(cid).feature_id.astype(str)) for cid in col_ids}
    measured = [f for f in feats if f["measured"]]
    assert measured
    for f in measured:
        assert f["column_id"] in col_ids                       # point -> column
        assert f["evidence_id"].startswith(run_id)             # -> exported evidence row
        assert f["source_feature_id"] in store_features_by_col[f["column_id"]]  # -> source record
