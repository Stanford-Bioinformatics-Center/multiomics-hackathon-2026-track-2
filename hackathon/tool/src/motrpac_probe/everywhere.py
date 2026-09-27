"""Section F: the signature in all 422 MoTrPAC comparisons, ranked, with non-exercise contrasts as the
negative control.

Method details: docs/METHODS.md#everywherepy
"""
import numpy as np
import pandas as pd
from scipy.stats import binomtest, false_discovery_control

from .core import col_words, score_columns

MIN_N = 10


def unit_key(c):
    return c.label if c.dataset == "human_acute" else c.time


def run(S, sig, cutoff=0.05, min_n=MIN_N, cids=None):
    cols = S.cols.sort_values("col_order")
    if cids is not None:
        cols = cols[cols.column_id.isin(cids)]
    per = score_columns(S, cols.column_id, sig, cutoff=cutoff)
    per = cols[["column_id", "dataset", "species", "tissue", "layer", "contrast", "label", "kind", "sex", "time",
                "col_order", "early_human_muscle"]].reset_index(drop=True).merge(per, on="column_id")
    per["unit"] = [unit_key(c) for c in per.itertuples()]
    per["words"] = [col_words(c) for c in per.itertuples()]
    rows = []
    for (ds, t, l, u, kind), d in per.groupby(["dataset", "tissue", "layer", "unit", "kind"], sort=False):
        n, k, s = int(d.n_measured.sum()), int(d.n_opposed.sum()), int(d.n_same.sum())
        rows.append(dict(dataset=ds, tissue=t, layer=l, time=u, kind=kind, n_columns=len(d), n_measured=n,
                         n_opposed=k, n_same=s, frac_opposed=k / n if n else np.nan, frac_same=s / n if n else np.nan,
                         sign_p=binomtest(k, n, 0.5).pvalue if n else np.nan,
                         n_fdrbh05=int(d.n_sig.sum()), n_fdrbh05_opposed=int(d.n_sig_opposed.sum()),
                         n_fdrbh05_same=int(d.n_sig_same.sum()),
                         camera_t_mean=float(d.camera_t.mean()) if d.camera_t.notna().any() else np.nan,
                         early_human_muscle=bool(d.early_human_muscle.any())))
    rank = pd.DataFrame(rows)
    ok = rank.sign_p.notna()
    rank["sign_q_bh"] = np.nan
    if ok.any():
        rank.loc[ok, "sign_q_bh"] = false_discovery_control(rank.loc[ok, "sign_p"], method="bh")
    rank["ranked"] = rank.n_measured >= min_n
    rank["majority"] = np.where(rank.frac_opposed > 0.5, "opposed", "same")
    rank["unit_words"] = (rank.dataset.map({"rat_train": "Rat", "human_acute": "Human"}) + " " + rank.tissue + " "
                          + rank.layer + " " + rank.time.str.replace("−", "–"))
    view = rank[rank.ranked]
    summary = dict(
        n_columns=len(per), n_units=len(rank), n_ranked=len(view), min_n=min_n,
        nominal_opposed=int(((view.sign_p < 0.05) & (view.frac_opposed > 0.5)).sum()),
        nominal_same=int(((view.sign_p < 0.05) & (view.frac_same > 0.5)).sum()),
        q_opposed=int(((view.sign_q_bh < 0.05) & (view.frac_opposed > 0.5)).sum()),
        q_same=int(((view.sign_q_bh < 0.05) & (view.frac_same > 0.5)).sum()),
        expected_by_chance=round(0.05 * len(view), 1))
    bykind = view.groupby("kind").agg(units=("frac_opposed", "size"), median_frac_opposed=("frac_opposed", "median"),
                                      q_opposed=("sign_q_bh", lambda q: int(((q < 0.05) & (view.loc[q.index, "frac_opposed"] > 0.5)).sum())),
                                      q_same=("sign_q_bh", lambda q: int(((q < 0.05) & (view.loc[q.index, "frac_same"] > 0.5)).sum())))
    return per, rank, summary, bykind.reset_index()
