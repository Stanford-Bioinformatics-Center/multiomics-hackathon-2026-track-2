"""Compare the rest-to-peak change between ST000763 groups.

For each inferred pair and metabolite: change = log2(peak) - log2(rest).
Within a group, a one-sample t test asks whether the mean change differs from 0.
Between groups, a Welch t test compares mean changes, e.g.
(PAH peak - rest) - (healthy peak - rest). No MoTrPAC data are used here.
"""

import numpy as np
import pandas as pd
from scipy import stats

from common import GROUPS, OUTPUT, PROCESSED, Q, bh

CONTRASTS = {
    "PAH_vs_Healthy": (["PAH"], ["Healthy"]),
    "PAH_vs_Normal_Pressures": (["PAH"], ["Normal Pressures"]),
    "PAH_vs_Borderline_Pressures": (["PAH"], ["Borderline Pressures"]),
    "PAH_vs_SSc_no_PH_all": (["PAH"], ["Normal Pressures", "LowRisk"]),
}

features = pd.read_csv(PROCESSED / "st000763_features.csv").set_index("feature_key")
log2x = pd.read_csv(PROCESSED / "st000763_log2.csv.gz", index_col="feature_key")
pairs = pd.read_csv(OUTPUT / "01_inferred_pairs.csv")
change = pd.DataFrame({r.pair_id: log2x[r.peak_sample] - log2x[r.rest_sample]
                       for r in pairs.itertuples()}).T
group = pairs.set_index("pair_id").group.reindex(change.index)

rows = []
for fid in change.columns:
    d = change[fid]
    row = {"feature_key": fid, "metabolite_name": features.loc[fid, "metabolite_name"],
           "refmet_name": features.loc[fid, "refmet_name"]}
    for g in GROUPS:
        v = d[group == g].dropna()
        row[f"n_pairs_{g}"] = len(v)
        row[f"mean_change_{g}"] = v.mean() if len(v) else np.nan
        row[f"p_within_{g}"] = (stats.ttest_1samp(v, 0).pvalue
                                if len(v) >= 5 and v.std() > 0 else np.nan)
    for name, (case, ctrl) in CONTRASTS.items():
        a, b = d[group.isin(case)].dropna(), d[group.isin(ctrl)].dropna()
        # At least 70% of each arm's pairs measured, and at least seven pairs.
        ok = all(len(v) >= max(7, int(np.ceil(0.7 * group.isin(g).sum())))
                 for v, g in ((a, case), (b, ctrl)))
        row[f"eligible_{name}"] = ok
        row[f"change_difference_{name}"] = a.mean() - b.mean() if ok else np.nan
        row[f"p_{name}"] = stats.ttest_ind(a, b, equal_var=False).pvalue if ok else np.nan
    rows.append(row)

res = pd.DataFrame(rows)
for g in GROUPS:
    res[f"q_within_{g}"] = bh(res[f"p_within_{g}"])
for name in CONTRASTS:
    res[f"q_{name}"] = bh(res[f"p_{name}"])
res.to_csv(OUTPUT / "02_exercise_response_tests.csv", index=False)

summary = []
for name in CONTRASTS:
    q = res[f"q_{name}"]
    summary.append({"test": f"change difference, {name}", "n_eligible": int(res[f"eligible_{name}"].sum()),
                    "n_q_lt_0_05": int((q < Q).sum()), "n_raw_p_lt_0_05": int((res[f"p_{name}"] < 0.05).sum()),
                    "smallest_q": q.min()})
for g in GROUPS:
    q = res[f"q_within_{g}"]
    summary.append({"test": f"within-group change, {g}", "n_eligible": int(q.notna().sum()),
                    "n_q_lt_0_05": int((q < Q).sum()), "n_raw_p_lt_0_05": int((res[f"p_within_{g}"] < 0.05).sum()),
                    "smallest_q": q.min()})
summary = pd.DataFrame(summary)
summary.to_csv(OUTPUT / "02_exercise_response_summary.csv", index=False)
print(summary.to_string(index=False))
print("\nWithin-PAH changes at q < 0.05:")
print(res.loc[res["q_within_PAH"] < Q, ["refmet_name", "mean_change_PAH", "q_within_PAH"]].to_string(index=False))
