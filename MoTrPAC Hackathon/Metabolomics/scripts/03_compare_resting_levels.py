"""Compare resting metabolite levels between ST000763 groups.

Only rest samples are used, so no inferred pairing is needed. Every PAH,
borderline and normal-pressure sample was drawn in the catheterization
setting ("Cath"); every healthy and low-risk SSc sample was "Non-invasive".
PAH versus healthy therefore mixes PAH, scleroderma and sampling setting, so
same-setting comparisons and a setting check are run alongside it:

    PAH_vs_Healthy            requested contrast; disease + SSc + setting
    PAH_vs_Normal_Pressures   same setting, all SSc
    PAH_vs_SSc_no_PH_all      mixes Cath and Non-invasive comparators
    Normal_vs_LowRisk         setting check; both SSc without PH
    Cath_severity_trend       Normal = 0, Borderline = 1, PAH = 2, Cath only
"""

import gzip

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

from common import DATA, GROUPS, OUTPUT, PROCESSED, Q, bh, primary_st000763_features

CONTRASTS = {
    "PAH_vs_Healthy": (["PAH"], ["Healthy"]),
    "PAH_vs_Normal_Pressures": (["PAH"], ["Normal Pressures"]),
    "PAH_vs_SSc_no_PH_all": (["PAH"], ["Normal Pressures", "LowRisk"]),
    "Normal_vs_LowRisk": (["Normal Pressures"], ["LowRisk"]),
}
TREND_SCORE = {"Normal Pressures": 0, "Borderline Pressures": 1, "PAH": 2}
GMT = DATA / "refmet" / "metabolomics.workbench.refmet.2024.08.07.metabolites.gmt.gz"


def welch(a, b):
    return stats.ttest_ind(a, b, equal_var=False).pvalue


def enough(values, in_arm):
    """At least 70% of an arm measured, and at least seven samples."""
    return values[in_arm].notna().sum() >= max(7, int(np.ceil(0.7 * in_arm.sum())))


samples = pd.read_csv(OUTPUT / "01_st000763_samples.csv")
rest = samples[samples.exercise == "No"].set_index("sample_id")
features = pd.read_csv(PROCESSED / "st000763_features.csv").set_index("feature_key")
log2x = pd.read_csv(PROCESSED / "st000763_log2.csv.gz", index_col="feature_key")[rest.index]
group = rest.group
female = (rest.sex == "Female").astype(float)

# Sensitivity scale: subtract each sample's median within its assay run.
centered = log2x.copy()
for _, idx in features.groupby("analysis_id").groups.items():
    centered.loc[idx] = log2x.loc[idx] - log2x.loc[idx].median(axis=0)

# 1. Feature tests.
rows = []
for fid in log2x.index:
    y, yc = log2x.loc[fid], centered.loc[fid]
    row = {"feature_key": fid, "metabolite_name": features.loc[fid, "metabolite_name"],
           "refmet_name": features.loc[fid, "refmet_name"]}
    for g in GROUPS:
        v = y[group == g].dropna()
        row[f"n_{g}"] = len(v)
        row[f"mean_log2_{g}"] = v.mean() if len(v) else np.nan
    for name, (case, ctrl) in CONTRASTS.items():
        ok = enough(y, group.isin(case)) and enough(y, group.isin(ctrl))
        a, b = y[group.isin(case)].dropna(), y[group.isin(ctrl)].dropna()
        row[f"eligible_{name}"] = ok
        row[f"log2_difference_{name}"] = a.mean() - b.mean() if ok else np.nan
        row[f"p_welch_{name}"] = welch(a, b) if ok else np.nan
        row[f"p_mannwhitney_{name}"] = stats.mannwhitneyu(a, b).pvalue if ok else np.nan
        row[f"p_centered_{name}"] = (welch(yc[group.isin(case)].dropna(), yc[group.isin(ctrl)].dropna())
                                     if ok else np.nan)
        row[f"p_sex_adjusted_{name}"] = np.nan
        if ok:
            keep = group.isin(case + ctrl) & y.notna()
            X = sm.add_constant(pd.DataFrame({"case": group[keep].isin(case).astype(float),
                                              "female": female[keep]}), has_constant="add")
            row[f"p_sex_adjusted_{name}"] = sm.OLS(y[keep], X).fit().pvalues["case"]
    ok = all(enough(y, group == g) for g in TREND_SCORE)
    row["eligible_Cath_severity_trend"] = ok
    row["log2_slope_Cath_severity_trend"] = np.nan
    row["p_welch_Cath_severity_trend"] = np.nan
    if ok:
        yy = y[group.isin(TREND_SCORE)].dropna()
        fit = stats.linregress(group[yy.index].map(TREND_SCORE).astype(float), yy)
        row["log2_slope_Cath_severity_trend"] = fit.slope
        row["p_welch_Cath_severity_trend"] = fit.pvalue
    rows.append(row)

res = pd.DataFrame(rows)
for name in CONTRASTS:
    for method in ["welch", "mannwhitney", "centered", "sex_adjusted"]:
        res[f"q_{method}_{name}"] = bh(res[f"p_{method}_{name}"])
res["q_welch_Cath_severity_trend"] = bh(res.p_welch_Cath_severity_trend)
# Share of each PAH-versus-healthy difference already present in normal-pressure
# SSc (same Cath setting, no PAH). Near 1 means PAH is not needed to explain it.
res["normal_share_of_PAH_vs_Healthy"] = (
    (res["mean_log2_Normal Pressures"] - res.mean_log2_Healthy) / res.log2_difference_PAH_vs_Healthy)
res.to_csv(OUTPUT / "03_resting_feature_tests.csv", index=False)

# 2. RefMet class scores: z-scored primary features, median per class.
primary = primary_st000763_features(features.reset_index())
z = log2x.loc[primary.feature_key]
z.index = primary.refmet_name.to_numpy()
z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0)
classes = {}
with gzip.open(GMT, "rt") as handle:
    for line in handle:
        parts = line.rstrip("\n").split("\t")
        members = sorted(set(parts[2:]) & set(z.index))
        if len(members) >= 3:
            classes[parts[0]] = members
class_rows = []
for cname, members in classes.items():
    n_ok = z.loc[members].notna().sum()
    score = z.loc[members].median().where(n_ok >= max(3, int(np.ceil(0.5 * len(members)))))
    row = {"class": cname, "n_metabolites": len(members), "metabolites": "; ".join(members)}
    for g in GROUPS:
        row[f"mean_score_{g}"] = score[group == g].mean()
    for name, (case, ctrl) in CONTRASTS.items():
        a, b = score[group.isin(case)].dropna(), score[group.isin(ctrl)].dropna()
        row[f"difference_{name}"] = a.mean() - b.mean()
        row[f"p_{name}"] = welch(a, b) if len(a) >= 7 and len(b) >= 7 else np.nan
    class_rows.append(row)
cls = pd.DataFrame(class_rows)
for name in CONTRASTS:
    cls[f"q_{name}"] = bh(cls[f"p_{name}"])
cls.sort_values("class").to_csv(OUTPUT / "03_resting_class_tests.csv", index=False)

# 3. Summary.
summary = []
for name in list(CONTRASTS) + ["Cath_severity_trend"]:
    q = res[f"q_welch_{name}"]
    hits = res[q < Q]
    entry = {"test": name, "n_eligible": int(res[f"eligible_{name}"].sum()),
             "n_q_lt_0_05": len(hits), "n_unique_refmet": hits.refmet_name.nunique(),
             "smallest_q": q.min()}
    for method in ["mannwhitney", "centered", "sex_adjusted"]:
        if f"q_{method}_{name}" in res:
            entry[f"n_q_lt_0_05_{method}"] = int((res[f"q_{method}_{name}"] < Q).sum())
    if name != "Cath_severity_trend":
        entry["classes_q_lt_0_05"] = "; ".join(cls.loc[cls[f"q_{name}"] < Q, "class"])
    summary.append(entry)
summary = pd.DataFrame(summary)
summary.to_csv(OUTPUT / "03_resting_summary.csv", index=False)
hits = res[res.q_welch_PAH_vs_Healthy < Q]
print(summary.drop(columns="classes_q_lt_0_05").to_string(index=False))
print(f"\nPAH vs healthy hits: {(hits.log2_difference_PAH_vs_Healthy > 0).sum()} higher, "
      f"{(hits.log2_difference_PAH_vs_Healthy < 0).sum()} lower in PAH; median share already "
      f"in normal-pressure SSc = {hits.normal_share_of_PAH_vs_Healthy.median():.2f}")
