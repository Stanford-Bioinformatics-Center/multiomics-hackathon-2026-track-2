"""Robustness matrix assertions. Results come from tests/matrix/run_matrix.py (run it first; on the cluster:
sbatch the desktop job, ~1 min with 12 workers). Summary: tests/REPORT.md (tests/matrix/write_report.py)."""
from pathlib import Path

import pandas as pd
import pytest

M = Path(__file__).parent / "matrix"
pytestmark = pytest.mark.skipif(not (M / "results.csv").exists(), reason="run tests/matrix/run_matrix.py first")


@pytest.fixture(scope="module")
def res():
    return pd.read_csv(M / "results.csv").set_index("case")


def test_no_crash(res):
    bad = res[res.exit != 0]
    assert bad.empty, bad.stderr_tail.to_dict()


def test_reports_validate(res):
    bad = res[res.validation_problems.fillna("") != ""]
    assert bad.empty, bad.validation_problems.to_dict()


def test_runtime_under_90s(res):
    assert (res.runtime_s < 90).all(), res.runtime_s[res.runtime_s >= 90].to_dict()


def test_negative_controls_at_null(res):
    """Random gene sets (any size / direction pattern) and the matched negative control: median class-null
    percentile within [5, 95] and no majority of columns with a set-level verdict."""
    neg = res[res.index.str.contains("random") & ~res.index.str.contains("mito")]
    assert neg.median_pct_class.between(5, 95).all(), neg.median_pct_class.to_dict()
    assert (neg.n_opposed_fdr + neg.n_same_fdr <= 5).all()


def test_positive_control_sign(res):
    """Hostrup HIIT proteome (direction = up after training) agrees in the SAME direction with rat 8 wk protein."""
    assert res.loc["example_hostrup2022_hiit_proteome", "rat_gn_prot_f8_t"] < -3


def test_id_types_equivalent(res):
    """Symbols, UniProt, Ensembl, rat symbols, mixed case and direction words give identical results."""
    ref = res.loc["id_human_symbols", "rat_gn_prot_f8_t"]
    for c in ["id_uniprot_in_symbol", "id_ensembl", "id_rat_symbols", "id_mixed_case", "direction_words",
              "weights_present", "weights_absent"]:
        assert res.loc[c, "rat_gn_prot_f8_t"] == pytest.approx(ref), c


def test_unknown_and_duplicates_reported(res):
    assert res.loc["id_20pct_unknown", "n_unmapped"] == 6
    assert res.loc["id_duplicates", "n_unmapped"] == 11  # 5 duplicates + 3 conflicting pairs (both rows dropped)


def test_empty_overlap_is_a_cell_not_an_error(res):
    r = res.loc["empty_overlap_tissue"]
    assert r.exit == 0 and r.n_columns_no_genes > 0 and r.headline_has_no_genes_cell == 1


def test_seed_reproducibility():
    s = pd.read_csv(M / "seed_compare.csv")
    assert (s[s.kind == "non-null"].max_abs_diff == 0).all()
    nul = s[s.kind != "non-null"]
    assert (nul.median_abs_diff <= 3).all() and (nul.max_abs_diff <= 10).all()


def test_real_signatures_run(res):
    ex = res[res.group.str.startswith("example: disease")]
    assert len(ex) >= 6 and (ex.exit == 0).all()
