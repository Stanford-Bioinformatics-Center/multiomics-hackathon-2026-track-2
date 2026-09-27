"""Tests for motrpac_probe.discord_extra (protein-without-RNA and lag panels). Uses the real store."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402
from scipy.stats import false_discovery_control  # noqa: E402

from motrpac_probe import core  # noqa: E402
from motrpac_probe import discord_extra as DX  # noqa: E402

TISSUES = ["SKM-GN", "HEART"]
KEYS = {"tables", "figures", "numbers", "captions", "notes"}


@pytest.fixture(scope="module")
def S():
    return core.get_store()


@pytest.fixture(scope="module")
def prot(S):
    out = {t: DX.protein_without_rna(S, t, nboot=100, nnull=50) for t in TISSUES}
    yield out
    plt.close("all")


@pytest.fixture(scope="module")
def lag(S):
    out = {t: DX.lag_panel(S, t, nperm=50) for t in TISSUES}
    yield out
    plt.close("all")


def _check_structure(r, tables, figures):
    assert set(r) == KEYS
    assert set(r["tables"]) == set(tables)
    for k, df in r["tables"].items():
        assert isinstance(df, pd.DataFrame) and len(df) > 0, k
    assert set(r["figures"]) == set(figures)
    for f in r["figures"].values():
        assert isinstance(f, Figure)
    for k in list(tables) + list(figures):
        title, how = r["captions"][k]
        assert title and how
    assert r["notes"] and all(isinstance(n, str) for n in r["notes"])
    assert all(np.isscalar(v) for v in r["numbers"].values())


@pytest.mark.parametrize("tissue", TISSUES)
def test_protein_without_rna(prot, tissue):
    r = prot[tissue]
    _check_structure(r, ["classes", "nulls", "mirror", "stability", "illustrations"],
                     ["calls_vs_null", "trajectories"])
    nb = r["numbers"]
    assert nb["n_repl"] > 100
    assert nb["n_strict"] <= nb["n_lenient"] <= nb["n_repl_measurable"] <= nb["n_repl"]
    assert (nb["n_supported"] + nb["n_opposite"] + nb["n_lenient"] + nb["n_undetermined"]
            == nb["n_repl_measurable"])
    assert nb["n_every"] <= nb["n_lenient"]
    # label-shuffle analogue: chance produces essentially no protein-only calls
    assert nb["null_shuffle_mean_lenient"] < 5 and nb["null_shuffle_mean_strict"] < 5
    assert nb["null_shuffle_mean_repl"] < 5
    assert 0 < nb["protein_reliability"] <= 1
    ill = r["tables"]["illustrations"]
    assert 0 < len(ill) <= 12
    assert {"rna_F_1w", "prot_M_8w"} <= set(ill.columns)
    assert any("protonly/RESULTS.md" in n for n in r["notes"])


@pytest.mark.parametrize("tissue", TISSUES)
def test_lag_panel(lag, tissue):
    r = lag[tissue]
    _check_structure(r, ["models", "n_genes"], ["delta_r2"])
    m = r["tables"]["models"]
    assert len(m) == 3 * (3 + 4 + 3)
    # permutation null band (5-95%) of held-out dR2 contains 0 for every model x sex
    assert (m.null_heldout_q05 <= 0).all() and (m.null_heldout_q95 >= 0).all()
    assert m.perm_p_heldout.between(0, 1).all()
    assert any("lag/RESULTS.md" in n for n in r["notes"])


def test_lag_skmgn_reproduces_original(lag):
    """In-sample grid dR2 equals lag/RESULTS.md (pooled, raw) to the reported precision; the 8-week RNA reference
    adds more than 2-week RNA."""
    nb = lag["SKM-GN"]["numbers"]
    m = lag["SKM-GN"]["tables"]["models"]
    g = m[(m.family == "grid") & (m.sex == "pooled")].set_index("rna_week").dR2_in_sample
    assert g["2w"] == pytest.approx(0.0037, abs=5e-5)
    assert nb["fixed_8w_heldout_dR2"] > nb["fixed_2w_heldout_dR2"] > 0


def test_not_available(S):
    for fn in (DX.protein_without_rna, DX.lag_panel):
        r = fn(S, "BLOOD")          # rat blood has RNA only
        assert set(r) == KEYS and r["tables"] == {} and r["notes"][0].startswith("not available")


def test_bh_matches_scipy():
    p = np.random.default_rng(0).uniform(size=(3, 200)) ** 3
    q = DX._bh(p)
    for i in range(3):
        np.testing.assert_allclose(q[i], false_discovery_control(p[i], method="bh"))


def test_foldcv_matches_lstsq():
    r = np.random.default_rng(1)
    n = 500
    X = np.column_stack([np.ones(n), r.normal(size=(n, 2))])
    y, x = r.normal(size=n), r.normal(size=n)
    fold = r.permutation(n) % 5
    assert DX._FoldCV(X, y, fold, 5).sse(x) == pytest.approx(DX._oof_sse(np.column_stack([X, x]), y, fold, 5))
