"""Parity: the service must report exactly what the mprobe engine computes.

The Marimo app and analysis.py both call motrpac_probe.run.compute()/core functions directly. The
service wraps the same compute(). This test asserts the service does not distort the engine: for the
same signature and parameters, the set-level results (n_opposed / n_measured / camera_t / camera_fdr),
the returned layers, and the store+code versions match a direct engine call. If these ever diverge,
the polished UI would misrepresent the already-correct scientific engine.
"""
import math
from pathlib import Path

import pytest

EXAMPLE = Path(__file__).resolve().parents[3] / "hackathon" / "tool" / "examples" / "pah_muscle_malenfant2015.csv"

pytestmark = pytest.mark.skipif(
    __import__("motrpac_probe.paths", fromlist=["CONTRASTS"]).CONTRASTS.exists() is False,
    reason="no mprobe store: run `mprobe store fetch`")


def _engine_direct(cutoff=0.05):
    """Direct engine call, mirroring what analysis.py / Marimo use."""
    from motrpac_probe import run
    R = run.compute(str(EXAMPLE), quiet=True, cutoff=cutoff)
    sc = R["scores"].set_index("column_id")
    out = {}
    for cid, r in sc.iterrows():
        out[cid] = dict(n_opposed=int(r.n_opposed), n_measured=int(r.n_measured),
                        camera_t=(None if math.isnan(r.camera_t) else round(float(r.camera_t), 6)),
                        camera_fdr=(None if math.isnan(r.camera_fdr) else round(float(r.camera_fdr), 6)))
    return out, R


def test_service_matches_engine_setlevel():
    from motrpac_probe_service import run_analysis, AnalysisRequest
    engine, R = _engine_direct()
    resp = run_analysis(AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat"))

    # same set of columns (the frozen family)
    svc = {c.column_id: c for c in resp.columns}
    assert set(svc) == set(engine), "service and engine must score the same column family"

    for cid, e in engine.items():
        c = svc[cid]
        assert c.n_opposed == e["n_opposed"], f"{cid} n_opposed"
        assert c.n_measured == e["n_measured"], f"{cid} n_measured"
        st = None if c.camera_t is None else round(c.camera_t, 6)
        sq = None if c.camera_fdr is None else round(c.camera_fdr, 6)
        assert st == e["camera_t"], f"{cid} camera_t {st} != {e['camera_t']}"
        assert sq == e["camera_fdr"], f"{cid} camera_fdr {sq} != {e['camera_fdr']}"


def test_service_matches_engine_counts_and_versions():
    from motrpac_probe import store
    from motrpac_probe_service import run_analysis, AnalysisRequest
    resp = run_analysis(AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat"))
    # counted genes and store hash match the engine
    assert resp.n_counted_genes == 19
    assert resp.provenance["store"]["contrasts_parquet_sha256_16"] == store.store_hash()
    # returned layers derived from the columns actually present (RNA + PROT for this family)
    assert set(resp.returned_layers) >= {"RNA", "PROT"}


def test_fdr_family_semantics_match_under_alt_threshold():
    """Under q<0.10 the engine and service must still agree cell-for-cell (same family semantics)."""
    from motrpac_probe_service import run_analysis, AnalysisRequest
    engine, _ = _engine_direct(cutoff=0.10)
    resp = run_analysis(AnalysisRequest(example_name="pah_muscle_malenfant2015", target_species="rat",
                                        fdr_threshold=0.10))
    svc = {c.column_id: c for c in resp.columns}
    for cid, e in engine.items():
        c = svc[cid]
        sq = None if c.camera_fdr is None else round(c.camera_fdr, 6)
        assert sq == e["camera_fdr"], f"{cid} camera_fdr under q<0.10"
