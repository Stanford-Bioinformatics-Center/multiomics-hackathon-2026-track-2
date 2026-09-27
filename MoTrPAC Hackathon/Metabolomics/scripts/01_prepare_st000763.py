"""Read ST000763 sample labels and metabolite values; infer rest/peak pairs.

The Workbench files have no participant IDs. A rest sample is paired with the
next-numbered sample when that sample is a peak sample with the same group,
sex and sampling setting. This gives 108 *inferred* pairs, not verified ones.
"""

import json

import numpy as np
import pandas as pd

from common import OUTPUT, PROCESSED, WORKBENCH


def fields(text):
    z = dict(v.strip().split(":", 1) for v in text.split(" | "))
    return {k.lower(): v for k, v in z.items()}


# 1. Samples.
factor_json = json.loads((WORKBENCH / "factors.json").read_text(encoding="utf-8"))
samples = pd.DataFrame([{"sample_id": x["local_sample_id"], **fields(x["factors"])}
                        for x in factor_json.values()])
samples = samples.rename(columns={"type": "setting"})
samples["sample_num"] = samples.sample_id.str.replace("S", "", regex=False).astype(int)
samples = samples.sort_values("sample_num").reset_index(drop=True)
assert len(samples) == 218 and samples.sample_id.is_unique

# 2. Inferred rest/peak pairs.
by_num = samples.set_index("sample_num")
pairs = []
for row in samples.itertuples():
    if row.exercise != "No" or row.sample_num + 1 not in by_num.index:
        continue
    peak = by_num.loc[row.sample_num + 1]
    if peak.exercise == "Yes" and (row.group, row.sex, row.setting) == (peak.group, peak.sex, peak.setting):
        pairs.append({"pair_id": row.sample_id, "rest_sample": row.sample_id,
                      "peak_sample": peak.sample_id, "group": row.group})
pairs = pd.DataFrame(pairs)
assert len(pairs) == 108 and pairs.peak_sample.is_unique
samples["pair_status"] = np.where(samples.sample_id.isin(pairs.rest_sample), "rest in pair",
                         np.where(samples.sample_id.isin(pairs.peak_sample), "peak in pair", "unpaired"))

# 3. Metabolite values, log2 transformed ("Peak area normalized", all positive).
met_json = json.loads((WORKBENCH / "data.json").read_text(encoding="utf-8"))
meta, rows = [], []
for x in met_json.values():
    fid = x["analysis_id"] + ":" + x["metabolite_id"]
    meta.append({"feature_key": fid, "analysis_id": x["analysis_id"],
                 "analysis_summary": x.get("analysis_summary", ""),
                 "metabolite_name": x["metabolite_name"], "refmet_name": x["refmet_name"]})
    rows.append(pd.Series(x["DATA"], name=fid, dtype=float))
features = pd.DataFrame(meta)
intensity = pd.DataFrame(rows)[samples.sample_id]
assert features.feature_key.is_unique and intensity.min().min() > 0
features["n_measured"] = intensity.notna().sum(axis=1).to_numpy()
features["median_intensity"] = intensity.median(axis=1).to_numpy()
log2x = np.log2(intensity)
log2x.index.name = "feature_key"

# 4. Save.
samples.to_csv(OUTPUT / "01_st000763_samples.csv", index=False)
pairs.to_csv(OUTPUT / "01_inferred_pairs.csv", index=False)
counts = (samples.groupby(["group", "setting", "exercise"]).size()
          .unstack("exercise", fill_value=0).rename(columns={"No": "rest", "Yes": "peak"})
          .reset_index())
counts["inferred_pairs"] = counts.group.map(pairs.group.value_counts()).fillna(0).astype(int)
counts.to_csv(OUTPUT / "01_group_setting_counts.csv", index=False)
features.to_csv(PROCESSED / "st000763_features.csv", index=False)
log2x.to_csv(PROCESSED / "st000763_log2.csv.gz")

print(counts.to_string(index=False))
print(f"{len(features)} features ({features.refmet_name.nunique()} RefMet names); "
      f"{len(pairs)} inferred pairs")
