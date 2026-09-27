"""CLI smoke test: `mprobe run` on the five signature examples, plus the positive / negative / pathway-class controls.

The five runs are started in parallel once per module (about 30 s wall time). tca_intermediates_demo.csv is not
included: it has only a refmet_name column, which signature.load rejects (no metabolite path yet).
"""
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from motrpac_probe import paths, render

MPROBE = str(Path(sys.executable).parent / "mprobe")
EXAMPLES = {
    "pah_muscle_malenfant2015": ["--pool-sets", "KEGG_OXIDATIVE_PHOSPHORYLATION,KEGG_CITRATE_CYCLE_TCA_CYCLE",
                                 "--pool-groups", "OXPHOS-down,TCA/transport-down"],
    "pah_blood_cheadle2012_eds": [],
    "hostrup2022_hiit_proteome": [],
    "random_matched": [],
    "random_mito9": [],
}
NBOOT = 200


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    base = tmp_path_factory.mktemp("mprobe_cli")
    env = {**os.environ, "MPLBACKEND": "Agg", "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"}
    procs = {}
    for name, extra in EXAMPLES.items():
        out = base / name
        cmd = [MPROBE, "run", "--signature", str(paths.EXAMPLES / f"{name}.csv"), "--quiet", "--nboot", str(NBOOT),
               "--out", str(out)] + extra
        procs[name] = (out, subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env))
    res = {}
    for name, (out, p) in procs.items():
        so, se = p.communicate(timeout=600)
        res[name] = dict(out=out, rc=p.returncode, stdout=so, stderr=se)
    return res


def _scores(runs, name):
    return pd.read_csv(runs[name]["out"] / "tables" / "column_scores.csv").set_index("column_id")


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_run_outputs(runs, name):
    r = runs[name]
    assert r["rc"] == 0, r["stderr"][-3000:]
    out = r["out"]
    report = out / "report.html"
    assert report.is_file() and report.stat().st_size > 10_000
    assert r["stdout"].strip().endswith("report.html")
    pngs = sorted((out / "figures").glob("*.png"))
    assert {p.stem for p in pngs} >= {"coverage", "grid", "camera", "nulls", "everywhere", "layers"}
    for p in pngs:
        assert p.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", p
    assert render.validate_report(report) == []
    sc = pd.read_csv(out / "tables" / "column_scores.csv")
    for c in ["comparison", "column_id", "camera_t", "camera_fdr", "pct_abund", "pct_class", "verdict"]:
        assert c in sc.columns
    assert len(sc) == 52  # VL + SKM-GN + SKM-VL + HEART, RNA + PROT, exercise comparisons


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_provenance_json(runs, name):
    out = runs[name]["out"]
    p = json.loads((out / "provenance.json").read_text(encoding="utf-8"))
    for k in ["repo_git_sha", "packages", "store", "signature", "timestamp", "command"]:
        assert k in p, k
    sig_file = paths.EXAMPLES / f"{name}.csv"
    assert p["signature"]["sha256_16"] == hashlib.sha256(sig_file.read_bytes()).hexdigest()[:16]
    assert p["command"].startswith("mprobe run")
    assert p["parameters"]["nboot"] == NBOOT
    assert {"numpy", "pandas", "scipy"} <= set(p["packages"])
    assert p["store"]["contrasts_parquet_sha256_16"]
    assert len(p["repo_git_sha"]) >= 7


def test_pah_pool_table(runs):
    t = pd.read_csv(runs["pah_muscle_malenfant2015"]["out"] / "tables" / "pool_null.csv")
    assert set(t.pool_size) == {124, 120}


# ------------------------------------------------------------------------------------------------ controls
def test_positive_control_hostrup(runs):
    """Human HIIT proteome moves the same way as rat 8-wk gastrocnemius protein training (same direction = t < 0)."""
    r = _scores(runs, "hostrup2022_hiit_proteome").loc["rat_train|SKM-GN|PROT|F_8w"]
    assert r.camera_t < -3
    assert r.verdict == "same direction"


def test_negative_control_random_matched_no_verdicts(runs):
    sc = _scores(runs, "random_matched")
    assert not (sc.camera_fdr < 0.05).any()
    assert set(sc.verdict) == {"no set-level shift"}
    assert 25 <= sc.pct_class.median() <= 75


def test_negative_control_random_matched_class_percentiles(runs):
    """Negative control lands at the null: median class-null percentile within [5, 95] and at most 10% of the core
    columns outside [2, 98]. (Columns are correlated, e.g. the heart protein columns, so a few extreme
    percentiles are expected; none reaches a set-level verdict, see above.)"""
    sc = _scores(runs, "random_matched")
    assert 5 <= sc.pct_class.median() <= 95
    assert (~sc.pct_class.between(2, 98)).mean() <= 0.10, sc.loc[~sc.pct_class.between(2, 98), "pct_class"].to_dict()


def test_pathway_class_control_random_mito9(runs):
    """Any mito set marked 'down' looks opposed against abundance-matched genes but not against its own class."""
    r = _scores(runs, "random_mito9").loc["human_acute|VL|RNA|EE_vs_CON_24h"]
    assert r.pct_abund >= 90
    assert r.pct_class < 90
