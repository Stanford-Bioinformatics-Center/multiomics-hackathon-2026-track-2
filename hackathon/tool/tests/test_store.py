"""Store schema and provenance: tool/store/*.parquet as documented in store.py (ROW_COLS / CAT_COLS; there is no
store/SCHEMA.md) and as recorded in provenance.parquet."""
import json

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import pytest

from motrpac_probe import paths, store

FLOAT64 = ["logFC", "stat", "p", "fdr", "fdr_bh"]
FLOAT32 = ["baseline_expr", "prot_n_missing", "training_fdr"]


@pytest.fixture(scope="module")
def rows():
    return store.load()


@pytest.fixture(scope="module")
def cols():
    return store.load_columns()


@pytest.fixture(scope="module")
def prov():
    return store.load_provenance()


def test_contrasts_schema(rows):
    schema_md = paths.STORE / "SCHEMA.md"
    if schema_md.exists():  # documented schema wins when present
        documented = [ln.split("|")[1].strip(" `") for ln in schema_md.read_text().splitlines()
                      if ln.startswith("|") and "`" in ln.split("|")[1]]
        assert set(store.ROW_COLS) <= set(documented) | set(store.ROW_COLS)
    assert list(rows.columns) == store.ROW_COLS
    for c in store.CAT_COLS:
        assert isinstance(rows[c].dtype, pd.CategoricalDtype), c
    for c in FLOAT64:
        assert rows[c].dtype == np.float64, c
    for c in FLOAT32:
        assert rows[c].dtype == np.float32, c
    assert rows.n_collapsed.dtype == np.int16
    assert (rows.n_collapsed >= 1).all()
    assert rows.stat.notna().all()
    assert rows.gene_symbol_human.notna().all()
    assert set(rows.species.cat.categories) == {"human", "rat"}
    assert set(rows.source.cat.categories) == {"join_table_v2", "alltissue_v1"}


def test_parquet_file_metadata():
    md = pq.ParquetFile(paths.CONTRASTS).schema_arrow
    assert md.names == store.ROW_COLS


def test_column_id_consistent(rows):
    sample = rows.iloc[:: 50_000]
    ids = [store.column_id(*t) for t in zip(sample.dataset.astype(str), sample.tissue.astype(str),
                                            sample.layer.astype(str), sample.contrast.astype(str))]
    assert ids == sample.column_id.astype(str).tolist()


def test_no_duplicate_gene_rows(rows):
    assert not rows[["column_id", "gene_symbol_human"]].duplicated().any()


def test_join_rows_match_csv(rows, prov):
    j = pd.read_csv(paths.JOIN, usecols=["dataset", "tissue", "layer", "contrast"], dtype=str, engine="pyarrow")
    assert len(j) == 542526 == int(prov["join_table_rows_csv"])
    jr = rows[rows.source == "join_table_v2"]
    assert len(jr) == len(j) == int(prov["join_table_rows_store"])
    csv_counts = (j.dataset + "|" + j.tissue + "|" + j.layer + "|" + j.contrast).value_counts().sort_index()
    st_counts = jr.column_id.value_counts()
    st_counts.index = st_counts.index.astype(str)
    st_counts = st_counts[st_counts > 0].sort_index()
    pd.testing.assert_series_equal(csv_counts, st_counts, check_names=False)
    # every column is entirely from one source
    src = rows.groupby("column_id", observed=True).source.nunique()
    assert (src == 1).all()


def test_columns_index(rows, cols, prov):
    assert len(cols) == 422 == int(prov["store_columns"])
    assert cols.column_id.is_unique
    assert list(cols.col_order) == list(range(len(cols)))
    n = rows.column_id.value_counts()
    n.index = n.index.astype(str)
    n = n[n > 0]
    assert (cols.set_index("column_id").n_genes == n.reindex(cols.column_id).to_numpy()).all()
    assert len(rows) == int(prov["store_rows"])
    assert set(cols.column_id) == set(n.index)
    for c in ["dataset", "tissue", "layer", "contrast", "kind", "label", "gtex_tissue"]:
        assert c in cols.columns
    assert set(cols.kind) <= {"exercise vs control", "exercise post−pre", "reference (non-exercise)",
                              "training vs sedentary"}
    early = cols[cols.early_human_muscle]
    assert len(early) and early.contrast.str.contains("15to45min").all() and (early.tissue == "VL").all()


def test_provenance(prov):
    assert prov["MotrpacRatTraining6moData"] == "2.0.0"
    assert prov["MotrpacHumanPreSuspensionAnalysis"] == "2.0.8"
    assert prov["join_table_version"].startswith("join_table_v2")
    xc = json.loads(prov["crosscheck_join_vs_rebuilt"])
    assert xc["n_join"] == xc["n_rebuilt"] == xc["n_both"] == xc["n_stat_identical"] == xc["n_fdr_bh_identical"]
    assert xc["n_join"] == 542526
    assert xc["n_join_only"] == xc["n_rebuilt_only"] == 0
    for k in ["built_at", "repo_git_sha", "gtex", "mitocarta", "legacy_script_sha256_16"]:
        assert prov.get(k), k
    inputs = [k for k in prov if k.startswith("input:")]
    assert "input:data/join_table_v2.csv" in inputs
    for k in inputs:
        v = json.loads(prov[k])
        assert v["bytes"] > 0


def test_annotations(rows):
    ann = store.load_annotations()
    assert ann.index.is_unique
    assert set(rows.gene_symbol_human.cat.categories) <= set(ann.index)
    for c in ["mitocarta3", "go_any_complex", "secreted_proxy", "go_contractile_fiber"]:
        assert ann[c].dtype == bool, c
    assert ann.loc["NDUFA9", "mitocarta3"] and not ann.loc["MYH7", "mitocarta3"]
    for t in {v for v in store.GTEX_TISSUE.values() if v}:
        assert f"gtex_tpm|{t}" in ann.columns


def test_genesets():
    gs = store.load_genesets()
    for k in ["KEGG_OXIDATIVE_PHOSPHORYLATION", "KEGG_CITRATE_CYCLE_TCA_CYCLE", "GOBP_OXIDATIVE_PHOSPHORYLATION",
              "MITOCARTA_ALL", "GOCC_CONTRACTILE_FIBER", "GOCC_SARCOPLASMIC_RETICULUM"]:
        assert len(gs[k]) > 10, k
    assert gs["MYH_SET"] == {"MYH1", "MYH2", "MYH4", "MYH7"}
    assert "NDUFA9" in gs["KEGG_OXIDATIVE_PHOSPHORYLATION"]


def test_id_map():
    m = store.load_idmap()
    assert list(m.columns) == ["id_type", "id", "gene_symbol_human"]
    assert set(m.id_type) == {"uniprot", "ensembl", "rat_symbol"}
    assert not m.duplicated().any()
    look = {(r.id_type, r.id): r.gene_symbol_human for r in m.itertuples()}
    assert look[("uniprot", "P04217")] == "A1BG"
    assert look[("ensembl", "ENSG00000139180")] == "NDUFA9"
    assert look[("rat_symbol", "Ndufa9")] == "NDUFA9"
