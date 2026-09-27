"""Golden test: the PAH muscle example reproduces the numbers in the deck / NARRATIVE.

Reference values are read from ../deck_extracts/deck_numbers.csv (compared to the precision they were printed
with), ../fig/tables/pah_grid_summary.md, pah_grid_v2_summary.md and ../everywhere/pah_everywhere_ranking.md;
nothing is hard-coded where a reference file has the number.
"""
import re

import numpy as np
import pandas as pd
import pytest
from conftest import deck_close

from motrpac_probe import core, layers, legacy, nulls, paths, run, signature, store

PAH = paths.EXAMPLES / "pah_muscle_malenfant2015.csv"
SEED = 20260926
POOL_SETS = ("KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_CITRATE_CYCLE_TCA_CYCLE")
cid = store.column_id
H = lambda layer, c: cid("human_acute", "VL", layer, c)  # noqa: E731
EE24_RNA, EE24_PROT, RE24_PROT = H("RNA", "EE_vs_CON_24h"), H("PROT", "EE_vs_CON_24h"), H("PROT", "RE_vs_CON_24h")
F8_PROT, M8_PROT = cid("rat_train", "SKM-GN", "PROT", "F_8w"), cid("rat_train", "SKM-GN", "PROT", "M_8w")


# ------------------------------------------------------------------------------------------------ fixtures
@pytest.fixture(scope="module")
def R(S):
    """The full `mprobe run` computation for the PAH example (1000 draws, deck seed, fixed-pool null on)."""
    return run.compute(PAH, _path=str(PAH), quiet=True, nboot=1000, seed=SEED, pool_sets=POOL_SETS,
                       pool_groups=("OXPHOS-down", "TCA/transport-down"))


@pytest.fixture(scope="module")
def scores(R):
    return R["scores"].set_index("column_id")


@pytest.fixture(scope="module")
def sig(S):
    return signature.load(PAH)


def _md_table(path):
    lines = [ln for ln in path.read_text().splitlines() if ln.startswith("|")]
    head = [c.strip() for c in lines[0].strip("|").split("|")]
    rows = [[c.strip() for c in ln.strip("|").split("|")] for ln in lines[2:]]
    return pd.DataFrame(rows, columns=head)


# ------------------------------------------------------------------------------------------------ mapping
def test_pah_signature_mapping(sig):
    assert len(sig.genes) == 19
    assert (sig.dirs == -1).sum() == 9 and (sig.dirs == 1).sum() == 10
    assert list(np.array(sig.genes)[sig.dirs == -1]) == legacy.MITO9
    assert len(sig.shown) == 25
    assert (sig.table.status == "mapped").all()


# ------------------------------------------------------------------------------------------------ counts / sign test
def test_headline_counts(scores, deck):
    r = scores.loc[EE24_RNA]
    assert r.n_measured == int(float(deck["R1b_n_measured"])) == 18
    assert r.n_opposed == int(float(deck["R1b_obs_opposed"])) == 15
    r = scores.loc[RE24_PROT]
    assert (r.n_opposed, r.n_measured) == (2, 19)
    assert r.sign_p == pytest.approx(0.000729, rel=1e-3)
    r = scores.loc[F8_PROT]
    assert (r.n_opposed, r.n_measured) == (15, 19)
    r = scores.loc[EE24_PROT]
    assert (r.n_opposed, r.n_measured) == (4, 19)


@pytest.mark.parametrize("fname", ["pah_grid_summary.md", "pah_grid_v2_summary.md"])
def test_grid_summary_tables(scores, fname):
    """Every column of the deck's PAH grid: n measured / n opposed / sign-test p over the 19 counted genes."""
    t = _md_table(paths.FIGTAB / fname)
    assert len(t) >= 28
    for r in t.itertuples():
        s = scores.loc[cid(r.dataset, r.tissue, r.layer, r.contrast)]
        assert s.n_measured == int(r.n_measured_19), r
        assert s.n_opposed == int(r.n_opposed_19), r
        ref = float(r.sign_p_19)  # printed with 3 significant digits
        assert s.sign_p == pytest.approx(ref, rel=6e-3), r


# ------------------------------------------------------------------------------------------------ cameraPR
@pytest.mark.parametrize("col, tag", [(EE24_RNA, "hum_EE24_RNA"), (F8_PROT, "rat_GN_PROT_F8"),
                                      (M8_PROT, "rat_GN_PROT_M8")])
def test_down_half_camera_t(scores, deck, col, tag):
    obs = scores.loc[col, "camera_t_down"]
    assert abs(obs - float(deck[f"R1a_{tag}_obs_camera_t"])) < 1e-3


OX_KEYS = [("RNA", "EE_15to45min"), ("PROT", "EE_15to45min"), ("PROT", "RE_24h"), ("RNA", "EE_24h")]


@pytest.mark.parametrize("layer, tag", OX_KEYS)
def test_oxphos_full_vs_intrinsic_universe(S, deck, layer, tag):
    kegg = S.gs["KEGG_OXIDATIVE_PHOSPHORYLATION"]
    uni = run.intrinsic_universe(S)
    f = S.frame(H(layer, tag.replace("_", "_vs_CON_", 1)))
    g = f.gene_symbol_human.astype(str)
    t_full = legacy.camera_pr(f.stat.to_numpy(), g.isin(kegg).to_numpy())[0]
    keep = g.isin(uni).to_numpy()
    t_int = legacy.camera_pr(f.stat.to_numpy()[keep], g[keep].isin(kegg).to_numpy())[0]
    assert abs(t_full - float(deck[f"HY_{layer}_{tag}_t_OXPHOS_full"])) < 0.01
    assert abs(t_int - float(deck[f"HY_{layer}_{tag}_t_OXPHOS_intrinsic"])) < 0.01
    assert keep.sum() == int(float(deck[f"HY_n_universe_intrinsic_{layer}"]))


def test_oxphos_intrinsic_via_score_column(S, deck):
    """Same number through the tool's own `universe=` path (up-half-only signature of OXPHOS genes, sign +1)."""
    kegg = sorted(S.gs["KEGG_OXIDATIVE_PHOSPHORYLATION"])
    osig = core.Sig(name="ox", table=pd.DataFrame(), genes=kegg, dirs=-np.ones(len(kegg), int))
    r = core.score_column(S, H("RNA", "EE_vs_CON_24h"), osig, universe=run.intrinsic_universe(S))
    # dirs -1 -> member stats multiplied by +1 -> plain cameraPR t
    assert abs(r["camera_t"] - float(deck["HY_RNA_EE_24h_t_OXPHOS_intrinsic"])) < 0.01
    assert r["n_universe"] == int(float(deck["HY_n_universe_intrinsic_RNA"]))


# ------------------------------------------------------------------------------------------------ specificity nulls
@pytest.fixture(scope="module")
def r1(S, sig):
    """NARRATIVE R1a then R1b on ONE rng stream, exactly as scripts/05_pah_grid.py 'robust'."""
    rng = np.random.default_rng(SEED)
    pools = nulls.legacy_pools(S, list(POOL_SETS))
    cids = [cid(*k) for k in legacy.GRID_COLS_V2]
    pool, _ = nulls.pool_null(S, cids, legacy.MITO9, pools, rng)
    # R1b: 19-gene sets matched on AveExpr decile of the human EE 24 h RNA column, same direction vector
    pos = S.col(EE24_RNA)
    base = S.base[pos]
    assert np.isfinite(base).all()
    dec = pd.qcut(base, 10, labels=False)
    dmap = dict(zip(S.codes(sig.genes), sig.dirs))
    d = np.array([dmap.get(c, 0) for c in S.gcode[pos]])
    mi = np.flatnonzero(d != 0)
    draws = nulls.matched_draws(dec, mi, rng, 1000)
    lfc = S.lfc[pos]
    null_k = (np.sign(lfc[draws]) * d[mi] < 0).sum(1)
    obs_k = int((np.sign(lfc[mi]) * d[mi] < 0).sum())
    return dict(pools=pools, pool=pool.set_index("column_id"), n=len(mi), obs=obs_k, null=null_k,
                decile_median=float(np.median(dec[mi])))


def test_r1a_pools(r1, deck, S):
    assert len(r1["pools"]["human_acute"]) == int(float(deck["R1a_pool_human"]))
    assert len(r1["pools"]["rat_train"]) == int(float(deck["R1a_pool_rat"]))
    union = S.gs[POOL_SETS[0]] | S.gs[POOL_SETS[1]]
    assert len(set(legacy.MITO9) & union) == int(float(deck["R1a_mito9_in_pool"]))


@pytest.mark.parametrize("col, tag", [(EE24_RNA, "hum_EE24_RNA"), (F8_PROT, "rat_GN_PROT_F8"),
                                      (M8_PROT, "rat_GN_PROT_M8")])
def test_r1a_pool_null(r1, R, deck, col, tag):
    p = r1["pool"].loc[col]
    assert deck_close(p.pct_t, f"R1a_{tag}_pct_t", deck)
    assert deck_close(p.obs_t, f"R1a_{tag}_obs_camera_t", deck)
    assert deck_close(p.null_t_median, f"R1a_{tag}_null_t_median", deck)
    # the run pipeline (--pool-sets / --pool-groups) draws the same sets and gets the same percentile
    q = R["pool"].set_index("column_id").loc[col]
    assert q.pct_t == pytest.approx(p.pct_t)
    assert sorted(R["pool_genes"]) == sorted(legacy.MITO9)


def test_r1b_abundance_matched(r1, deck):
    assert r1["n"] == int(float(deck["R1b_n_measured"]))
    assert r1["obs"] == int(float(deck["R1b_obs_opposed"]))
    assert np.median(r1["null"]) == float(deck["R1b_null_median"])
    assert np.percentile(r1["null"], 95) == float(deck["R1b_null_q95"])
    p_emp = (1 + (r1["null"] >= r1["obs"]).sum()) / (len(r1["null"]) + 1)
    assert deck_close(p_emp, "R1b_p_empirical", deck)
    assert r1["decile_median"] == float(deck["R1b_decile_of_pah_genes_median"])


def test_effective_n(S, deck):
    e = core.effective_n(S, legacy.MITO9, [cid(*k) for k in legacy.GRID_COLS])
    assert e["n_genes"] == 9
    for k, key in [("meff_nyholt", "R3_meff_nyholt"), ("meff_liji", "R3_meff_liji"), ("mean_r", "R3_mean_offdiag_r"),
                   ("min_r", "R3_min_offdiag_r"), ("max_r", "R3_max_offdiag_r"), ("top_eigen_frac", "R3_top_eigen_frac")]:
        assert abs(e[k] - float(deck[key])) < 0.01, k
        assert deck_close(e[k], key, deck), k


# ------------------------------------------------------------------------------------------------ layers
def _rho_keys(deck):
    out = []
    for k in deck:
        m = re.fullmatch(r"rho_(hum|rat|ratVL)_(.+)_all", k)
        if m:
            out.append((k, m.group(1), m.group(2)))
    return out


def test_genome_wide_rna_prot_rho(R, deck):
    d = R["disc"]
    keys = _rho_keys(deck)
    assert len(keys) >= 12
    for key, kind, c in keys:
        if kind == "hum":
            m = d[(d.dataset == "human_acute") & (d.tissue == "VL") & (d.contrast == c) & ~d.cross]
        else:
            m = d[(d.dataset == "rat_train") & (d.tissue == "SKM-GN") & (d.contrast == c) & (d.cross == (kind == "ratVL"))]
        assert len(m) == 1, key
        assert abs(m.rho_all.iloc[0] - float(deck[key])) <= 0.005, (key, m.rho_all.iloc[0])


def test_layer_pairs_direct(S, sig, deck):
    pairs = layers.layer_pairs(S)
    assert pairs.cross.sum() == 8  # rat SKM-VL RNA vs SKM-GN PROT, 2 sexes x 4 weeks
    sel = pairs[pairs.prot.isin([EE24_PROT, H("PROT", "EE_vs_CON_15to45min"), F8_PROT,
                                 cid("rat_train", "SKM-GN", "PROT", "F_4w")])]
    summ, _ = layers.discordance(S, sig, sel)
    got = {(r.prot, r.cross): r.rho_all for r in summ.itertuples()}
    assert abs(got[(EE24_PROT, False)] - float(deck["rho_hum_EE_vs_CON_24h_all"])) <= 0.005
    assert abs(got[(H("PROT", "EE_vs_CON_15to45min"), False)] - float(deck["rho_hum_EE_vs_CON_15to45min_all"])) <= 0.005
    assert abs(got[(F8_PROT, False)] - float(deck["rho_rat_F_8w_all"])) <= 0.005
    assert abs(got[(cid("rat_train", "SKM-GN", "PROT", "F_4w"), True)] - float(deck["rho_ratVL_F_4w_all"])) <= 0.005


# ------------------------------------------------------------------------------------------------ everywhere
def test_everywhere_summary(R):
    txt = (paths.HACK / "everywhere" / "pah_everywhere_ranking.md").read_text()
    m = re.search(r"(\d+) ranked units\. Nominal sign_p < 0\.05: (\d+) opposed-majority and (\d+) same-majority"
                  r".*?After BH: (\d+) opposed-majority and (\d+) same-majority", txt, re.S)
    assert m, "reference sentence not found in pah_everywhere_ranking.md"
    ref = [int(x) for x in m.groups()]
    s = R["ev_summary"]
    assert [s["n_ranked"], s["nominal_opposed"], s["nominal_same"], s["q_opposed"], s["q_same"]] == ref
    assert ref == [258, 45, 19, 15, 3]
    n_units = re.search(r"BH over all (\d+) units", txt)
    assert s["n_units"] == int(n_units.group(1))


# ------------------------------------------------------------------------------------------------ up-10 table
UP10_COLS = [(c, layer) for c in ["EE_vs_CON_3.5to4h", "EE_vs_CON_24h", "RE_vs_CON_3.5to4h", "RE_vs_CON_24h"]
             for layer in ("RNA", "PROT")]  # scripts/06_discordance.py: hcol[1:3] + hcol[4:6], then layer


@pytest.fixture(scope="module")
def up10(S, sig):
    up = [g for g, d in zip(sig.genes, sig.dirs) if d == 1]
    s = core.Sig(name="up10", table=sig.table, genes=up, dirs=np.ones(len(up), int))
    return core.score_columns(S, [H(layer, c) for c, layer in UP10_COLS], s).set_index("column_id")


@pytest.mark.parametrize("k", range(8))
def test_up10_sign_table(up10, deck, k):
    c, layer = UP10_COLS[k]
    r = up10.loc[H(layer, c)]
    assert r.n_up == r.n_measured
    assert r.n_measured == int(float(deck[f"up10[{k}|n_measured]"]))
    assert r.n_opposed == int(float(deck[f"up10[{k}|n_opposed]"]))
    assert deck_close(r.sign_p, f"up10[{k}|sign_test_p]", deck)
    assert r.n_sig == int(float(deck[f"up10[{k}|n_fdrbh05]"]))
