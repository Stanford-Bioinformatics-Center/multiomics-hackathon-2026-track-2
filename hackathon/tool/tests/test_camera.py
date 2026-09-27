"""cameraPR: the vectorised port equals the scripts' port, and the scripts' port equals limma.

Where limma validation comes from: scripts/06_discordance.py stores NO cameraPR cases. The validated cases are the
precomputed CAMERA_RESULTS shipped in MotrpacHumanPreSuspensionAnalysis 2.0.8 (data/raw/human_muscle_camera_v1.csv,
limma::cameraPR on the package's z.std with gene-symbol-collapsed feature ids), which scripts/05_pah_grid.py
'crosscheck' (lines ~689-715) reproduces with its port (deck key XC_camera_max_absdiff_port_vs_pkg = 2.4e-14).
test_port_reproduces_package_camera_results replays that crosscheck for every contrast in both layers. The optional
R test (marked `slow`, MPROBE_RUN_R=1) calls limma::cameraPR directly on a small stored case.
"""
import os
import subprocess

import numpy as np
import pandas as pd
import pytest

from motrpac_probe import core, legacy, paths, store

RSCRIPT = os.path.expanduser("~/miniconda3/envs/motrpac/bin/Rscript")


# ------------------------------------------------------------------------------------------------ (a) vec vs legacy
@pytest.mark.parametrize("G, m, seed", [(500, 2, 1), (2000, 9, 2), (5000, 60, 3), (300, 150, 4)])
def test_camera_vec_equals_legacy(G, m, seed):
    rng = np.random.default_rng(seed)
    stat = rng.standard_t(5, G) + rng.normal(0, 0.3)
    B = 25
    sets = np.array([rng.choice(G, m, replace=False) for _ in range(B)])
    signs = rng.choice([-1.0, 1.0], m)
    for sg in (signs, None):
        got = core.camera_vec(stat, sets, sg)
        for b in range(B):
            s = stat.copy()
            if sg is not None:
                s[sets[b]] = s[sets[b]] * sg
            mask = np.zeros(G, bool)
            mask[sets[b]] = True
            ref = legacy.camera_pr(s, mask, 0.01)[0]
            assert abs(got[b] - ref) < 1e-9 * max(1.0, abs(ref)), (b, got[b], ref)


def test_camera_vec_other_inter_gene_cor():
    rng = np.random.default_rng(11)
    stat = rng.normal(size=1000)
    sets = np.array([rng.choice(1000, 20, replace=False) for _ in range(5)])
    got = core.camera_vec(stat, sets, inter_gene_cor=0.05)
    for b in range(5):
        mask = np.zeros(1000, bool)
        mask[sets[b]] = True
        assert abs(got[b] - legacy.camera_pr(stat, mask, 0.05)[0]) < 1e-9


def test_camera_signed_matches_legacy_and_edge_cases():
    rng = np.random.default_rng(5)
    stat = rng.normal(size=400)
    member = np.zeros(400, bool)
    member[rng.choice(400, 12, replace=False)] = True
    signs = rng.choice([-1.0, 1.0], 12)
    t, p = core.camera_signed(stat, member, signs)
    s = stat.copy()
    s[member] *= signs
    rt, _, rp = legacy.camera_pr(s, member, 0.01)
    assert t == pytest.approx(rt, abs=1e-12) and p == pytest.approx(rp, abs=1e-12)
    one = np.zeros(400, bool)
    one[3] = True
    assert np.isnan(core.camera_signed(stat, one, np.ones(1))[0])
    assert np.isnan(core.camera_vec(stat, np.array([[1]]))[0])


# ------------------------------------------------------------------------------------------------ (b) vs package CAMERA
@pytest.fixture(scope="module")
def pkg_camera():
    cam = pd.read_csv(paths.RAW / "human_muscle_camera_v1.csv", usecols=["assay", "contrast_short", "set", "t"],
                      engine="pyarrow")
    return cam[cam.set == "GOBP_OXIDATIVE_PHOSPHORYLATION"]


@pytest.fixture(scope="module")
def gobp(S):
    g = pd.read_parquet(store.GENESETS)
    g = g[(g.gs_name == "GOBP_OXIDATIVE_PHOSPHORYLATION") & g.source.str.startswith("MotrpacHumanPreSuspensionAnalysis")]
    genes = set(g.gene_symbol)
    assert genes == S.gs["GOBP_OXIDATIVE_PHOSPHORYLATION"]  # no other source contributes members
    return genes


@pytest.mark.parametrize("layer, fn, assay", [("RNA", "rna", "transcript-rna-seq"), ("PROT", "prot", "prot-pr")])
def test_port_reproduces_package_camera_results(pkg_camera, gobp, deck, layer, fn, assay):
    f2g = pd.read_csv(paths.RAW / "human_feature_to_gene_muscle_v1.csv", usecols=["assay", "feature_id", "gene_symbol"],
                      low_memory=False)
    raw = pd.read_csv(paths.RAW / f"human_muscle_{fn}_DA_v1.csv", usecols=["contrast_short", "feature_id", "z.std"],
                      engine="pyarrow")
    raw = raw.merge(f2g[f2g.assay == assay].drop(columns="assay").drop_duplicates(), on="feature_id", how="left")
    raw["new_id"] = raw.gene_symbol.fillna(raw.feature_id)  # unmapped feature ids are kept as their own id
    pre = pkg_camera[pkg_camera.assay == assay].set_index("contrast_short").t
    diffs = {}
    for short, x in raw.groupby("contrast_short", sort=False):
        if short not in pre.index:
            continue
        x = x.reindex(x["z.std"].abs().sort_values(ascending=False).index).drop_duplicates("new_id")
        t = legacy.camera_pr(x["z.std"].to_numpy(), x.new_id.isin(gobp).to_numpy())[0]
        diffs[short] = abs(t - pre[short])
    assert len(diffs) >= 6
    assert max(diffs.values()) < 1e-8, diffs
    assert max(diffs.values()) <= max(1e-8, float(deck["XC_camera_max_absdiff_port_vs_pkg"]) * 1e3)


# ------------------------------------------------------------------------------------------------ (c) vs limma in R
R_CODE = r"""
suppressMessages(library(limma))
a <- commandArgs(TRUE)
x <- read.csv(a[1])
sets <- list()
for (k in setdiff(names(x), "stat")) sets[[k]] <- which(x[[k]] == 1)
r <- cameraPR(x$stat, sets, inter.gene.cor = 0.01, sort = FALSE)
write.csv(data.frame(set = rownames(r), PValue = r$PValue, Direction = r$Direction), a[2], row.names = FALSE)
"""


@pytest.mark.slow
@pytest.mark.skipif(not os.path.exists(RSCRIPT), reason="motrpac R env not available")
def test_legacy_port_equals_limma(tmp_path, S):
    """Small stored case (seeded synthetic stats + a real column subset) through limma::cameraPR."""
    rng = np.random.default_rng(20260926)
    cases = []
    G = 300
    stat = rng.standard_t(4, G)
    cases.append(("synthetic", stat, {f"s{m}": rng.choice(G, m, replace=False) for m in (2, 5, 17, 60)}))
    f = S.frame(store.column_id("human_acute", "VL", "RNA", "EE_vs_CON_24h"))
    st = f.stat.to_numpy()
    g = f.gene_symbol_human.astype(str).to_numpy()
    cases.append(("EE24_RNA", st, {"oxphos": np.flatnonzero(np.isin(g, list(S.gs["KEGG_OXIDATIVE_PHOSPHORYLATION"]))),
                                   "mito9": np.flatnonzero(np.isin(g, legacy.MITO9))}))
    for name, stat, sets in cases:
        df = pd.DataFrame({"stat": stat})
        for k, idx in sets.items():
            df[k] = np.isin(np.arange(len(stat)), idx).astype(int)
        inp, out = tmp_path / f"{name}.csv", tmp_path / f"{name}_out.csv"
        df.to_csv(inp, index=False)
        (tmp_path / "cam.R").write_text(R_CODE, encoding="utf-8")
        subprocess.run([RSCRIPT, str(tmp_path / "cam.R"), str(inp), str(out)], check=True, capture_output=True)
        res = pd.read_csv(out).set_index("set")
        for k in sets:
            mask = df[k].to_numpy() == 1
            t, _, p = legacy.camera_pr(stat, mask, 0.01)
            assert p == pytest.approx(res.loc[k, "PValue"], rel=1e-8, abs=1e-300), (name, k)
            assert ("Up" if t > 0 else "Down") == res.loc[k, "Direction"], (name, k)

