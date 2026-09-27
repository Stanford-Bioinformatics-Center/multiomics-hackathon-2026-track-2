"""Compare each ST000763 group's rest-to-peak change with MoTrPAC's changes.

Patients have peak - rest only (no non-exercising control group), so three
MoTrPAC contrasts are shown side by side at the two during-exercise times
closest to "peak":

    EE-CON   control-adjusted exercise effect
    EE-EE    exercisers' change only: the like-for-like contrast
    CON-CON  controls' change only: drift without exercise

Spearman rho across matched metabolites describes how similarly the two
studies order the changes. It is descriptive; no P value is reported because
metabolites are correlated with each other and the platforms differ.
"""

import pandas as pd
from scipy import stats

from common import GROUPS, OUTPUT, PROCESSED, name_key, load_motrpac_one_feature_per_name, primary_st000763_features

TIMES = ["during_20_min", "during_40_min"]
CONTRASTS = ["EE-CON", "EE-EE", "CON-CON"]

features = pd.read_csv(PROCESSED / "st000763_features.csv")
response = pd.read_csv(OUTPUT / "02_exercise_response_tests.csv")
primary = primary_st000763_features(features)[["feature_key"]]
patient = response.merge(primary, on="feature_key")
patient["match_key"] = patient.refmet_name.map(name_key)
patient["PAH_minus_Healthy"] = patient.mean_change_PAH - patient.mean_change_Healthy

mo, _ = load_motrpac_one_feature_per_name()
mo = mo[(mo.tissue == "blood") & mo.Timepoint.isin(TIMES)]
wide = mo.pivot_table(index=["match_key", "Timepoint"], columns="contrast_category",
                      values="logFC").reset_index()
basis = mo.drop_duplicates("match_key").set_index("match_key")[["feature_id", "match_basis"]]

matched = patient.merge(wide, on="match_key")
matched = matched.join(basis, on="match_key")
keep = (["refmet_name", "Timepoint", "feature_id", "match_basis"]
        + [f"mean_change_{g}" for g in GROUPS] + ["PAH_minus_Healthy"] + CONTRASTS)
matched[keep].sort_values(["Timepoint", "refmet_name"]).to_csv(
    OUTPUT / "05_matched_metabolites.csv", index=False)

# "all" uses every exact name match. "refmet_name_only" is a sensitivity check
# that drops names matched through MoTrPAC's feature_id (mostly lipids).
rows = []
match_sets = {"all": matched, "refmet_name_only": matched[matched.match_basis == "refmet_name"]}
for set_name, subset in match_sets.items():
    for time, g in subset.groupby("Timepoint"):
        for contrast in CONTRASTS:
            row = {"match_set": set_name, "motrpac_time": time, "motrpac_contrast": contrast}
            for col in [f"mean_change_{x}" for x in GROUPS] + ["PAH_minus_Healthy"]:
                ok = g[[col, contrast]].dropna()
                row[col.replace("mean_change_", "rho_")] = stats.spearmanr(ok[col], ok[contrast])[0]
            row["n_metabolites"] = len(g[["mean_change_PAH", contrast]].dropna())
            rows.append(row)
align = pd.DataFrame(rows)
align.to_csv(OUTPUT / "05_alignment_by_group.csv", index=False)
print(align.round(3).to_string(index=False))
print(f"\n{matched.match_key.nunique()} metabolite names matched; "
      f"{(matched.drop_duplicates('match_key').match_basis == 'feature_id').sum()} via MoTrPAC feature_id")
