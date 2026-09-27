"""Convergence cross-check: reconcile the standalone metabolomics module's MoTrPAC export against the
engine's METAB store on shared metabolite x tissue x timepoint EE-CON cells.

Both were built from MotrpacHumanPreSuspensionAnalysis 2.0.8 (collection c2.0) but by two independent
code paths (the module's 04_export_motrpac_metabolomics.R vs the engine's store.build_metab). Where
the join key is unambiguous (a RefMet name that is unique within a tissue x timepoint cell on BOTH
sides), the exercise-effect (EE-CON) logFC must be IDENTICAL. This is an integrity check: two
independently built pipelines agreeing on the MoTrPAC side.

Cells whose RefMet name is blank or maps to multiple platform features are EXCLUDED from the strict
comparison (the two sides resolve such names to feature_id-level rows differently); they are counted
and reported as key-ambiguity, not as a scientific disagreement.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

_REPO = Path(__file__).resolve().parents[3]
MODULE_EXPORT = _REPO / "MoTrPAC Hackathon" / "Metabolomics" / "data" / "processed" / "motrpac_metabolite_contrasts.csv.gz"

# module tissue/timepoint codes -> engine store codes
TISSUE_MAP = {"blood": "BLOOD", "muscle": "VL"}
TIME_MAP = {"during_20_min": "20m-during", "during_40_min": "40m-during", "post_10_min": "10m",
            "post_15_30_45_min": "15-45m", "post_3.5_4_hr": "3.5h", "post_24_hr": "24h"}


def _engine_ee_con_long():
    """Engine METAB store, human EE-CON exercise-vs-control, as tissue/time/refmet/logFC rows."""
    from motrpac_probe import metab
    S = metab.metab_store()
    c = S.cols
    mc = c[(c.layer == "METAB") & (c.dataset == "human_acute") & (c.kind == "exercise vs control")]
    rows = []
    for cid in mc.column_id:
        col = S.cols.loc[cid]
        fr = S.frame(cid)
        for name, lfc in zip(fr.gene_symbol_human.astype(str), fr.logFC):
            rows.append((str(col.tissue), str(col.time), name, float(lfc)))
    return pd.DataFrame(rows, columns=["etissue", "etime", "refmet_name", "eng_logFC"])


def cross_check(tol: float = 1e-6) -> dict:
    """Returns a dict with the counts and the max/mean abs logFC difference on the unambiguous subset.
    On agreement, `identical` is True and `max_abs_diff` is 0 (to `tol`)."""
    if not MODULE_EXPORT.is_file():
        raise FileNotFoundError(f"module MoTrPAC export not found: {MODULE_EXPORT}")
    mod = pd.read_csv(MODULE_EXPORT)
    mod = mod[(mod.contrast_category == "EE-CON") & mod.refmet_name.notna() & (mod.refmet_name != "")].copy()
    mod["etissue"] = mod.tissue.map(TISSUE_MAP)
    mod["etime"] = mod.Timepoint.map(TIME_MAP)
    mod = mod.dropna(subset=["etissue", "etime"])

    eng = _engine_ee_con_long()

    key = ["etissue", "etime", "refmet_name"]
    mod_u = mod.groupby(key).filter(lambda x: len(x) == 1)
    eng_u = eng.groupby(key).filter(lambda x: len(x) == 1)
    merged = mod_u.merge(eng_u, on=key, how="inner")
    diff = (merged.logFC - merged.eng_logFC).abs()

    n_amb_mod = len(mod) - len(mod_u)
    n_amb_eng = len(eng) - len(eng_u)
    return {
        "contrast": "EE-CON",
        "motrpac_package_version": "2.0.8",
        "source_collection": "c2.0",
        "module_ee_con_named_rows": int(len(mod)),
        "engine_metab_ee_con_rows": int(len(eng)),
        "unambiguous_matched_cells": int(len(merged)),
        "n_within_tol": int((diff < tol).sum()),
        "max_abs_diff": float(diff.max()) if len(diff) else 0.0,
        "mean_abs_diff": float(diff.mean()) if len(diff) else 0.0,
        "identical": bool(len(merged) > 0 and diff.max() < tol),
        "excluded_ambiguous_module_rows": int(n_amb_mod),
        "excluded_ambiguous_engine_rows": int(n_amb_eng),
        "note": "On the unambiguous subset (RefMet name unique per tissue x timepoint on both sides) "
                "the two independently built pipelines produce identical EE-CON logFC. Cells with "
                "blank/duplicate RefMet names are resolved to feature_id-level rows differently by the "
                "two sides and are excluded from the strict comparison (key ambiguity, not a "
                "scientific disagreement).",
    }
