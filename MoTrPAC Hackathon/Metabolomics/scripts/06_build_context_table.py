"""Add healthy MoTrPAC exercise context to PAH metabolite lists.

Each list row is one metabolite from a PAH comparison. is_reported_hit = True
marks the PAH-altered metabolites; False marks the rest of what that study
measured, used only as background. MoTrPAC is not used to test PAH. For each
metabolite it answers:

    Does healthy endurance exercise change it (EE-CON), when, and for how long?
    Does it change in resting controls with no exercise (CON-CON drift)?
    Does muscle metabolomics change too (EE-CON, muscle)?

Lists use the standard columns in LIST_COLUMNS and come from:
    data/pah_metabolite_lists/*.csv   hand-curated, e.g. transcribed paper tables
    this script's adapter             ST000763 resting PAH vs healthy (script 03)
"""

import json

import numpy as np
import pandas as pd
from scipy import stats

from common import BLOOD_TIMES, DATA, OUTPUT, PROCESSED, Q, load_motrpac_one_feature_per_name, name_key

LIST_COLUMNS = [
    "list_id",          # short ID shared by every row of one list
    "source",           # paper or dataset
    "source_type",      # paper_reported or our_reanalysis
    "comparison",       # e.g. PAH rest vs healthy rest
    "metabolite_name",  # name as reported
    "refmet_name",      # RefMet name used for matching (blank: metabolite_name is used)
    "pah_log2_effect",  # PAH minus comparator, log2 scale; blank if not reported
    "pah_direction",    # up or down
    "pah_q_value",      # blank if not reported
    "is_reported_hit",  # True = PAH-altered; False = measured background
    "caveat",           # design limits that travel with the list
]
LATE_TIMES = ["post_3.5_4_hr", "post_24_hr"]
# A label needs BH q < 0.05 AND |log2 FC| >= 0.5 (about 1.4-fold). With
# MoTrPAC's sample size, q alone flags ~78% of metabolites as drifting in
# resting controls, which makes the label uninformative. Unthresholded values
# stay in the table.
MIN_LFC = 0.5
CURATED = DATA / "pah_metabolite_lists"
CURATED.mkdir(parents=True, exist_ok=True)


# 1. PAH lists.
def st000763_rest_list():
    """Adapter: script 03's resting PAH-versus-healthy comparison."""
    res = pd.read_csv(OUTPUT / "03_resting_feature_tests.csv")
    res = res[res.eligible_PAH_vs_Healthy & res.refmet_name.notna()].copy()
    res["n_rest"] = res[[c for c in res.columns if c.startswith("n_")]].sum(axis=1)
    # One row per RefMet name: smallest q, then the most-measured feature.
    res = (res.sort_values(["refmet_name", "q_welch_PAH_vs_Healthy", "n_rest"],
                           ascending=[True, True, False]).drop_duplicates("refmet_name"))
    out = pd.DataFrame({
        "list_id": "st000763_rest_pah_vs_healthy",
        "source": "Metabolomics Workbench ST000763 / PR000551",
        "source_type": "our_reanalysis",
        "comparison": "SSc-PAH rest vs healthy rest (script 03, Welch t, BH)",
        "metabolite_name": res.metabolite_name,
        "refmet_name": res.refmet_name,
        "pah_log2_effect": res.log2_difference_PAH_vs_Healthy,
        "pah_direction": np.where(res.log2_difference_PAH_vs_Healthy > 0, "up", "down"),
        "pah_q_value": res.q_welch_PAH_vs_Healthy,
        "is_reported_hit": res.q_welch_PAH_vs_Healthy < Q,
        "caveat": ("PAH samples all Cath, healthy all Non-invasive; PAH vs "
                   "normal-pressure SSc (same setting) had 0 hits"),
    })
    out.to_csv(PROCESSED / "list_st000763_rest_pah_vs_healthy.csv", index=False)
    return out


lists = [st000763_rest_list()]
for path in sorted(CURATED.glob("*.csv")):
    df = pd.read_csv(path)
    missing = set(LIST_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"{path.name} lacks columns {sorted(missing)}")
    lists.append(df[LIST_COLUMNS])
    print(f"Added curated list {path.name}: {len(df)} rows")
pah = pd.concat(lists, ignore_index=True)[LIST_COLUMNS]
pah["is_reported_hit"] = pah.is_reported_hit.astype(str).str.lower() == "true"
pah["match_key"] = pah.refmet_name.fillna(pah.metabolite_name).map(name_key)

# 2. MoTrPAC summaries per metabolite name.
mo, choice = load_motrpac_one_feature_per_name()


def wide(tissue, category):
    sub = mo[(mo.tissue == tissue) & (mo.contrast_category == category)]
    return sub.pivot_table(index="match_key", columns="Timepoint", values=["logFC", "adj_p_value"])


def summarise(tissue):
    ee, con, eee = wide(tissue, "EE-CON"), wide(tissue, "CON-CON"), wide(tissue, "EE-EE")
    times = [t for t in BLOOD_TIMES if t in ee["logFC"].columns]
    lfc, q = ee["logFC"][times], ee["adj_p_value"][times]
    sig = (q < Q) & (lfc.abs() >= MIN_LFC)
    peak_abs = lfc.where(sig).abs()
    peak_time = pd.Series([r.idxmax() if r.notna().any() else np.nan for _, r in peak_abs.iterrows()],
                          index=peak_abs.index)
    s = pd.DataFrame(index=lfc.index)
    s["exercise_sensitive"] = sig.any(axis=1)
    s["peak_time"] = peak_time
    s["peak_logFC"] = [lfc.loc[i, t] if isinstance(t, str) else np.nan for i, t in peak_time.items()]
    late = [t for t in LATE_TIMES if t in times]
    s["still_changed_late"] = sig[late].any(axis=1)
    s["control_drift"] = ((con["adj_p_value"] < Q) & (con["logFC"].abs() >= MIN_LFC)).any(axis=1)
    s["max_abs_control_logFC"] = con["logFC"].abs().max(axis=1)
    for t in times:
        s[f"ee_con_logFC_{t}"] = lfc[t]
        s[f"ee_con_q_{t}"] = q[t]
    for t in ["during_20_min", "during_40_min"]:
        if t in eee["logFC"].columns:
            s[f"ee_ee_logFC_{t}"] = eee["logFC"][t]
    return s


blood = summarise("blood").add_prefix("blood_")
muscle = summarise("muscle")[["exercise_sensitive", "peak_time", "peak_logFC"]].add_prefix("muscle_")
basis = (choice[choice.tissue == "blood"].set_index("match_key")[["feature_id", "platform", "match_basis"]]
         .add_prefix("motrpac_blood_"))

ctx = pah.join(basis, on="match_key").join(blood, on="match_key").join(muscle, on="match_key")
ctx["motrpac_blood_matched"] = ctx.motrpac_blood_feature_id.notna()
ctx["motrpac_muscle_matched"] = ctx.muscle_exercise_sensitive.notna()
same = np.sign(ctx.blood_peak_logFC) == np.where(ctx.pah_direction == "up", 1, -1)
ctx["exercise_vs_pah_direction"] = np.select(
    [~(ctx.blood_exercise_sensitive == True), same], ["no exercise change", "same"], "opposite")
# Different platforms and cohorts: rough scale comparisons only.
ctx["pah_effect_over_exercise_peak"] = ctx.pah_log2_effect.abs() / ctx.blood_peak_logFC.abs()
ctx["pah_effect_over_control_drift"] = ctx.pah_log2_effect.abs() / ctx.blood_max_abs_control_logFC


def label(r):
    if not r.motrpac_blood_matched:
        return "no MoTrPAC blood match"
    parts = []
    if r.blood_exercise_sensitive:
        parts.append("exercise-sensitive, persists >=3.5 h" if r.blood_still_changed_late
                     else "exercise-sensitive, transient")
    if r.blood_control_drift:
        parts.append("drifts without exercise")
    if r.motrpac_muscle_matched and r.muscle_exercise_sensitive:
        parts.append("also changes in muscle")
    return "; ".join(parts) if parts else "stable in MoTrPAC blood"


ctx["context_label"] = ctx.apply(label, axis=1)
ctx.drop(columns="match_key").to_csv(OUTPUT / "06_context_table_all.csv", index=False)
show = ["list_id", "refmet_name", "pah_direction", "pah_log2_effect", "pah_q_value", "context_label",
        "blood_ee_con_logFC_during_20_min", "blood_peak_time", "blood_peak_logFC",
        "blood_max_abs_control_logFC", "pah_effect_over_control_drift",
        "exercise_vs_pah_direction", "motrpac_blood_match_basis", "caveat"]
ctx[ctx.is_reported_hit].sort_values(["list_id", "pah_q_value"])[show].to_csv(
    OUTPUT / "06_context_table_hits.csv", index=False)

# 3. Are hits more exercise- or drift-sensitive than the study's other
#    matched metabolites? Only possible when a list includes background rows.
summary = {"motrpac_package_version": "2.0.8",
           "label_rule": f"BH q < {Q} and |log2 FC| >= {MIN_LFC}", "lists": {}}
for list_id, g in ctx.groupby("list_id"):
    m = g[g.motrpac_blood_matched]
    entry = {"n_hits": int(g.is_reported_hit.sum()),
             "n_hits_matched_blood": int(m.is_reported_hit.sum()),
             "n_background_matched_blood": int((~m.is_reported_hit).sum()),
             "hit_labels": g[g.is_reported_hit].context_label.value_counts().to_dict(),
             "median_pah_effect_over_control_drift_hits":
                 float(m.loc[m.is_reported_hit, "pah_effect_over_control_drift"].median())}
    if (~m.is_reported_hit).any() and m.is_reported_hit.any():
        for flag in ["blood_exercise_sensitive", "blood_control_drift"]:
            h, b = m.loc[m.is_reported_hit, flag].astype(bool), m.loc[~m.is_reported_hit, flag].astype(bool)
            table = [[int(h.sum()), int((~h).sum())], [int(b.sum()), int((~b).sum())]]
            entry[flag] = {"hits_fraction": round(float(h.mean()), 3),
                           "background_fraction": round(float(b.mean()), 3),
                           "fisher_p": round(float(stats.fisher_exact(table)[1]), 3)}
    summary["lists"][list_id] = entry
(OUTPUT / "06_context_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
