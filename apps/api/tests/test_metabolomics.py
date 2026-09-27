"""Metabolomics case-study adapter + convergence cross-check tests."""
import io
import json
import zipfile

import pytest

from motrpac_probe_service import metab_casestudy as mcs
from motrpac_probe_service import convergence

pytestmark = pytest.mark.skipif(
    not mcs.available(),
    reason="metabolomics module outputs not present")


@pytest.fixture(scope="module")
def client():
    from fastapi.testclient import TestClient
    from motrpac_probe_service.app import app
    return TestClient(app)


def test_case_study_matches_committed_summary():
    r = mcs.build_case_study()
    assert r.analysis_type == "metabolomics_case_study"
    assert r.case_study_id == "st000763_rest_pah_vs_healthy"
    assert "SSc-PAH plasma" in r.cohort_label and "ST000763" in r.cohort_label
    # numbers must match 06_context_summary.json
    assert r.n_hits == 41
    assert r.n_hits_matched_blood == 27
    assert r.hit_labels == {"stable in MoTrPAC blood": 24, "no MoTrPAC blood match": 14,
                            "drifts without exercise": 3}
    # the null result (0% exercise-sensitive hits, Fisher p 0.609)
    es = r.null_result["blood_exercise_sensitive"]
    assert es["hits_fraction"] == 0.0 and abs(es["fisher_p"] - 0.609) < 1e-6
    assert abs(r.median_pah_effect_over_control_drift - 3.0915882707822293) < 1e-6
    assert len(r.hits) == 41


def test_case_study_is_json_safe():
    from dataclasses import asdict
    r = mcs.build_case_study()
    json.dumps(asdict(r), allow_nan=False)  # must not raise


def test_case_study_preserves_conclusion_and_contrasts():
    r = mcs.build_case_study()
    # honest conclusion carried verbatim: setting/scleroderma, not PAH; not a treatment claim
    assert "sampling setting and/or scleroderma" in r.conclusion.lower()
    assert "not evidence exercise treats pah" in r.conclusion.lower()
    # both contrasts first-class; the EE-EE like-for-like caveat present; never "healthy gene set"
    assert "EE-EE" in r.contrasts["motrpac_exercise"]
    assert "no more MoTrPAC-like" in r.contrasts["motrpac_exercise"]["ee_ee_caveat"]
    assert r.contrasts["disease"]["denominator"].startswith("healthy plasma at rest")
    # separate-cohort caveat is explicit
    assert any("do not blend" in c.lower() for c in r.caveats)


def test_case_study_provenance_records_source():
    r = mcs.build_case_study()
    p = r.provenance
    assert p["motrpac_package_version"] == "2.0.8"
    assert "st000763_data" in p["source_datasets"]
    assert p["output_files"]["06_context_table_hits.csv"]  # sha recorded
    assert "read-only" in p["source_of_truth"].lower()


def test_convergence_identical_on_unambiguous_subset():
    """Two independently built pipelines (module 04 export vs engine METAB store) must produce
    IDENTICAL EE-CON logFC where the RefMet name is unambiguous."""
    res = convergence.cross_check()
    assert res["contrast"] == "EE-CON"
    assert res["motrpac_package_version"] == "2.0.8" and res["source_collection"] == "c2.0"
    assert res["unambiguous_matched_cells"] > 500      # a substantial shared subset
    assert res["identical"] is True
    assert res["max_abs_diff"] == 0.0
    assert res["n_within_tol"] == res["unambiguous_matched_cells"]


def test_api_metab_endpoints(client):
    cs = client.get("/api/metabolomics/casestudy")
    assert cs.status_code == 200
    body = cs.json()
    assert body["n_hits"] == 41 and body["analysis_type"] == "metabolomics_case_study"

    conv = client.get("/api/metabolomics/convergence")
    assert conv.status_code == 200 and conv.json()["identical"] is True

    exp = client.get("/api/metabolomics/casestudy/export")
    assert exp.status_code == 200 and exp.headers["content-type"] == "application/zip"
    z = zipfile.ZipFile(io.BytesIO(exp.content))
    assert {"summary.json", "hits.json", "conclusion.txt", "provenance.json"} == set(z.namelist())
    assert "scleroderma" in z.read("conclusion.txt").decode().lower()


def test_catalog_exposes_metabolomics_analysis_type(client):
    cat = client.get("/api/catalog").json()
    at = cat["analysis_types"]
    assert "metabolomics_case_study" in at
    assert at["metabolomics_case_study"]["api"] and at["metabolomics_case_study"]["react"]
    assert "SEPARATE cohort" in at["metabolomics_case_study"]["cohort"]
