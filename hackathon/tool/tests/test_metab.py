"""METAB layer: store rows, pathway map, metabolite signatures and the 3-layer pathway panel."""
import numpy as np
import pandas as pd
import pytest

from motrpac_probe import core, legacy, metab, store
from motrpac_probe.paths import EXAMPLES

TCA_DEMO = EXAMPLES / "tca_intermediates_demo.csv"


@pytest.fixture(scope="module")
def M():
    return store.load_metab()


@pytest.fixture(scope="module")
def MS():
    return metab.metab_store()


# ------------------------------------------------------------------------------------------------ store rows
def test_schema_equals_contrasts(M):
    rows, _ = M
    ref = pd.read_parquet(store.CONTRASTS).head(10)
    assert list(rows.columns) == list(ref.columns)
    kind = lambda d: "category" if isinstance(d, pd.CategoricalDtype) else str(d)  # noqa: E731
    assert [kind(d) for d in rows.dtypes] == [kind(d) for d in ref.dtypes]
    assert set(rows.layer) == {"METAB"}


def test_fdr_and_unique(M):
    rows, cols = M
    assert rows.fdr_bh.between(0, 1).all()
    assert not rows.duplicated(["column_id", "gene_symbol_human"]).any()
    assert set(rows.column_id) == set(cols.column_id)
    assert (rows.fdr_bh >= rows.p - 1e-12).all()


def test_rat_fdr_bh_before_collapse(M):
    """fdr_bh = BH over every tested platform feature of the tissue x contrast (incl. iSTDs), before collapsing."""
    rows, _ = M
    d = pd.read_csv(store.RAW_METAB / "rat_metab_DA_v1.csv.gz", low_memory=False,
                    usecols=["tissue", "dataset", "feature_ID", "sex", "comparison_group", "p_value"])
    d = d[(d.tissue == "SKM-GN") & (d.sex == "female") & (d.comparison_group == "8w")]
    d["fdr_bh"] = legacy.bh(d.p_value)
    ref = dict(zip(d.dataset + ":" + d.feature_ID, d.fdr_bh))
    got = rows[rows.column_id == "rat_train|SKM-GN|METAB|F_8w"]
    assert len(got) > 500
    np.testing.assert_allclose(got.fdr_bh.to_numpy(), got.feature_id.astype(str).map(ref).to_numpy(), rtol=1e-12)


def test_metab_store(MS):
    assert MS.rows.gene_symbol_human.dtype.name == "category"
    mc = metab.metab_columns(MS)
    assert mc and all(MS.cols.loc[c, "layer"] == "METAB" for c in mc)
    assert MS.cols.index.is_unique
    base_max = store.load_columns().col_order.max()
    assert (MS.cols.loc[mc, "col_order"] > base_max).all()
    assert len(MS.frame(mc[0])) == MS.cols.loc[mc[0], "n_genes"]
    assert metab.metab_columns(MS, ["SKM-GN"]) == [f"rat_train|SKM-GN|METAB|{s}_{w}w" for s in "FM"
                                                  for w in (1, 2, 4, 8)]


# ------------------------------------------------------------------------------------------------ pathway map
def test_pathway_map(MS):
    pm = store.load_pathway_map()
    assert {"pathway", "side", "member", "rule", "source"} <= set(pm.columns)
    assert set(pm.side) == {"gene", "metabolite"}
    assert not pm.duplicated(["pathway", "side", "member"]).any()
    measured = set(MS.genes)
    ok = 0
    for pw, d in pm.groupby("pathway"):
        ng = d[(d.side == "gene") & d.member.isin(measured)].member.nunique()
        nm = d[(d.side == "metabolite") & d.member.isin(measured)].member.nunique()
        ok += ng >= 3 and nm >= 3
    assert ok >= 5
    assert set(pm.pathway) == set(metab.PATHWAYS)


def test_pathway_panel_and_heatmap(MS):
    df = metab.pathway_panel(MS, tissue="SKM-GN", species="rat")
    assert set(df.layer) == {"RNA", "PROT", "METAB"}
    assert df.groupby("layer").contrast.nunique().eq(8).all()
    assert len(df) == 6 * 3 * 8
    assert df.t.notna().mean() > 0.9
    # sign: t > 0 = set up. Check one cell against legacy.camera_pr directly.
    r = df[(df.pathway == "TCA cycle") & (df.layer == "METAB")].iloc[0]
    f = MS.frame(r.column_id)
    mem = store.load_pathway_map().query("pathway == 'TCA cycle' and side == 'metabolite'").member
    t, _, _ = legacy.camera_pr(f.stat.to_numpy(), f.gene_symbol_human.isin(mem).to_numpy(), 0.01)
    assert np.isclose(t, r.t)
    h = metab.pathway_panel(MS, tissue="VL", species="human")
    assert sorted(h.contrast.unique()) == sorted(f"{g}_vs_CON_{t}" for g in ("EE", "RE")
                                                 for t in ("15to45min", "3.5to4h", "24h"))
    fig = metab.pathway_heatmap(df, "SKM-GN")
    assert len(fig.axes[0].texts) >= 100
    for w in ["cameraPR", "not flux", "many-to-many", "store/pathway_map.csv", "one-to-one"]:
        assert w in metab.PATHWAY_CAPTION


# ------------------------------------------------------------------------------------------------ signatures
def test_tca_demo(MS):
    assert metab.is_metab_signature(TCA_DEMO)
    sig = metab.load_metab_signature(TCA_DEMO)
    assert (sig.table.status == "mapped").sum() >= 8
    cids = metab.metab_columns(MS, ["SKM-GN", "VL"])
    assert any("|VL|" in c for c in cids) and any("|SKM-GN|" in c for c in cids)
    r = core.score_columns(MS, cids, sig)
    assert (r.n_measured >= 8).all()
    assert r.camera_t.notna().all()


def test_hmdb_and_unknown(tmp_path, MS):
    p = tmp_path / "hmdb.csv"
    p.write_text("hmdb,direction\nHMDB0000094,up\nHMDB00254,-1\nHMDB0000190,+1\nHMDB9999999,1\n", encoding="utf-8")
    assert metab.is_metab_signature(p)
    sig = metab.load_metab_signature(p)
    st = dict(zip(sig.table.hmdb, sig.table.status))
    assert st["HMDB9999999"].startswith("unmapped: HMDB id not in")
    assert {"Citric acid", "Succinic acid", "Lactic acid"} <= set(sig.genes)
    assert dict(zip(sig.genes, sig.dirs))["Succinic acid"] == -1
    q = tmp_path / "unknown.csv"
    q.write_text("refmet_name,direction\nNotAMetabolite xyz,1\nCitric acid,1\ncitric ACID,1\n", encoding="utf-8")
    s2 = metab.load_metab_signature(q)
    assert list(s2.table.status) == ["unmapped: RefMet name not measured in MoTrPAC metabolomics", "mapped",
                                     "dropped: duplicate of an earlier row"]
    assert s2.genes == ["Citric acid"]
    g = tmp_path / "genes.csv"
    g.write_text("gene_symbol,refmet_name,direction\nPPARGC1A,,1\n", encoding="utf-8")
    assert not metab.is_metab_signature(g)
