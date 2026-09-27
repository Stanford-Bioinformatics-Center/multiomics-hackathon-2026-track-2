"""motrpac_probe reproduces the team pipeline on main (`MoTrPAC Hackathon/`, collection c2.0).
Reference values are the team's own committed outputs, copied to tests/fixtures/team_main_c2/ (see SOURCE.md)."""
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from motrpac_probe import core, ranked, signature, store

FIX = Path(__file__).parent / "fixtures" / "team_main_c2"
EX = Path(__file__).parents[1] / "examples"
TIMES = {"during_20_min": "during20min", "during_40_min": "during40min", "post_10_min": "10min",
         "post_15_30_45_min": "15to45min", "post_3.5_4_hr": "3.5to4h", "post_24_hr": "24h"}
BLOOD = {tp: f"human_acute|BLOOD|RNA|EE_vs_CON_{c}" for tp, c in TIMES.items()}


@pytest.fixture(scope="module")
def blood():
    S = core.get_store()
    rk, sig, info = ranked.load_ranked(EX / "pah_blood_gse33463_ranked.csv.gz", exact_symbols=True)
    return S, rk, info


def test_rank_association_matches_team(blood):
    """Team script 08: Spearman rho, PAH PBMC t vs MoTrPAC blood RNA, six times (13,055 shared genes there;
    13,053 here because two genes' representative feature differs in the store's max-|t| collapse)."""
    S, rk, _ = blood
    ref = pd.read_csv(FIX / "pah_motrpac_blood_rank_association_by_time.csv").set_index("Timepoint")
    ra = ranked.rank_association(S, rk, list(BLOOD.values()), nperm=0).set_index("column_id")
    for tp, cid in BLOOD.items():
        assert ra.loc[cid, "rho"] == pytest.approx(ref.loc[tp, "spearman_rho"], abs=1e-3), tp
        assert abs(ra.loc[cid, "n_shared"] - ref.loc[tp, "shared_genes"]) <= 5


def test_disease_pathways_match_team(blood):
    """Team script 09: GO:BP cameraPR on the PAH ranking: 5,183 sets tested (10-500 genes), 63 at BH < 0.05."""
    _, _, info = blood
    qc = pd.read_csv(FIX / "pah_gobp_camera_qc.csv").set_index("metric").value
    dis = ranked.disease_pathways(info["full"])
    assert len(dis) == int(qc["gobp_sets_tested_10_to_500_genes"])
    assert int((dis.fdr < 0.05).sum()) == int(qc["gobp_sets_bh_below_0.05"])


def test_pathway_concordance_matches_team(blood):
    """Team scripts 10-11 (exercise side = MoTrPAC precomputed CAMERA): per-time shared sets, jointly significant,
    same / opposite direction; totals 40 = 27 + 13."""
    S, _, info = blood
    ref = pd.read_csv(FIX / "pah_motrpac_blood_gobp_time_summary.csv").set_index("Timepoint")
    dis = ranked.disease_pathways(info["full"])
    sm, _ = ranked.pathway_concordance(S, dis, list(BLOOD.values()), exercise_side="precomputed")
    sm = sm.set_index("column_id")
    for tp, cid in BLOOD.items():
        r = ref.loc[tp]
        assert sm.loc[cid, "shared_sets"] == r.shared_sets
        assert sm.loc[cid, "both_bh"] == r["both_bh_below_0.05"]
        assert sm.loc[cid, "same_direction"] == r.same_direction_both_bh
        assert sm.loc[cid, "opposite_direction"] == r.opposite_direction_both_bh
    assert (sm.both_bh.sum(), sm.same_direction.sum(), sm.opposite_direction.sum()) == (40, 27, 13)


def test_muscle_rna_timepoints_match_team(tmp_path):
    """Team script 00: the nine lower PAH proteins' RNA, EE vs control: positive logFC 0/9, 7/9, 9/9 and positive
    with published BH q < 0.05 0, 0, 7 at 15-45 min, 3.5-4 h, 24 h."""
    ref = pd.read_csv(FIX / "pah_muscle_rna_timepoint_summary.csv").set_index("motrpac_rna_timepoint")
    down9 = ["NDUFA9", "UQCRC2", "UQCRC1", "ATP5MG", "ATP5F1B", "IDH2", "OGDH", "SLC25A4", "ECH1"]
    rows = store.load(["column_id", "gene_symbol_human", "logFC", "fdr"])
    for tp, code in [("post_15_30_45_min", "15to45min"), ("post_3.5_4_hr", "3.5to4h"), ("post_24_hr", "24h")]:
        d = rows[(rows.column_id == f"human_acute|VL|RNA|EE_vs_CON_{code}") & rows.gene_symbol_human.isin(down9)]
        assert len(d) == ref.loc[tp, "matching_genes"]
        assert int((d.logFC > 0).sum()) == ref.loc[tp, "positive_logFC_count"]
        assert int(((d.logFC > 0) & (d.fdr < 0.05)).sum()) == ref.loc[tp, "positive_and_bh_q_lt_0_05_count"]
