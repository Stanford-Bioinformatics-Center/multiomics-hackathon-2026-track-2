#!/usr/bin/env python3
"""Count phosphosite changes with an equivalently small mapped parent protein.

This is a descriptive summary-level audit, not a phosphosite-occupancy model.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from build_discordance import filter_and_audit, read_effects, representative_features


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--contrast-category", default="EE-CON")
    parser.add_argument("--contrast-type", default="exercise_with_controls")
    parser.add_argument("--times", nargs="+", required=True)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-effect", type=float, default=0.2)
    parser.add_argument("--equivalence-margin", type=float, default=0.2)
    args = parser.parse_args()
    if len(set(args.times)) != len(args.times):
        parser.error("Times must be distinct")
    # Reuse the catalog's mapping exclusions and contrast checks.
    args.target_time = args.times[0]
    args.earlier_times = args.times[1:]
    original = read_effects(args.input)
    rows, _ = filter_and_audit(original, args)
    protein = representative_features(rows, "protein")
    accession_count = (
        protein.groupby("gene_symbol")["uniprot"].nunique().rename("accession_count").reset_index()
    )
    protein = protein.merge(accession_count, on="gene_symbol", how="left")
    protein = protein.loc[
        protein["uniprot"].ne("")
        & protein["protein_feature_count"].eq(1)
        & protein["accession_count"].eq(1)
    ]
    sites = rows.loc[rows["layer"].eq("phosphosite") & rows["uniprot"].ne("")].copy()
    parent_columns = ["gene_symbol", "uniprot", "timepoint", "feature_id", "log2_fc", "ci_lower", "ci_upper", "q_value"]
    matched = sites.merge(
        protein[parent_columns], on=["gene_symbol", "uniprot", "timepoint"],
        how="inner", suffixes=("_site", "_parent"), validate="many_to_one",
    )
    matched["site_supported"] = (
        matched["log2_fc_site"].abs().ge(args.min_effect)
        & matched["q_value_site"].le(args.alpha)
        & (
            (matched["log2_fc_site"].gt(0) & matched["ci_lower_site"].gt(0))
            | (matched["log2_fc_site"].lt(0) & matched["ci_upper_site"].lt(0))
        )
    )
    matched["parent_equivalent"] = (
        matched["ci_lower_parent"].gt(-args.equivalence_margin)
        & matched["ci_upper_parent"].lt(args.equivalence_margin)
    )
    matched["candidate"] = matched["site_supported"] & matched["parent_equivalent"]

    summaries = []
    for time in args.times:
        source_sites = original.loc[
            original["species"].eq("human")
            & original["tissue"].eq(args.tissue)
            & original["contrast_category"].eq(args.contrast_category)
            & original["timepoint"].eq(time)
            & original["layer"].eq("phosphosite")
        ]
        kept_sites = sites.loc[sites["timepoint"].eq(time)]
        group = matched.loc[matched["timepoint"].eq(time)]
        candidate = group.loc[group["candidate"]]
        summaries.append({
            "tissue": args.tissue,
            "contrast_category": args.contrast_category,
            "timepoint": time,
            "raw_phosphosite_rows": len(source_sites),
            "after_mapping_checks_with_uniprot": len(kept_sites),
            "excluded_missing_or_ambiguous_site_mapping": len(source_sites) - len(kept_sites),
            "matched_to_single_gene_uniprot_single_feature_parent": len(group),
            "matched_with_both_site_and_parent_ci": int(
                group[["ci_lower_site", "ci_upper_site", "ci_lower_parent", "ci_upper_parent"]].notna().all(axis=1).sum()
            ),
            "supported_phosphosite_rows_among_matched": int(group["site_supported"].sum()),
            "parent_equivalent_rows_among_matched": int(group["parent_equivalent"].sum()),
            "supported_phosphosite_and_parent_equivalent_rows": len(candidate),
            "candidate_genes": int(candidate["gene_symbol"].nunique()),
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(summaries).to_csv(args.out_dir / "ptm_parent_summary.csv", index=False)
    matched.loc[matched["candidate"], [
        "gene_symbol", "uniprot", "timepoint", "feature_id_site", "feature_id_parent",
        "log2_fc_site", "ci_lower_site", "ci_upper_site", "q_value_site",
        "log2_fc_parent", "ci_lower_parent", "ci_upper_parent", "q_value_parent",
    ]].to_csv(args.out_dir / "ptm_parent_candidates.csv", index=False)
    summary_frame = pd.DataFrame(summaries)
    figure, axes = plt.subplots(1, 2, figsize=(11, 4), constrained_layout=True)
    locations = np.arange(len(summary_frame))
    axes[0].bar(locations - 0.18, summary_frame["supported_phosphosite_rows_among_matched"],
                width=0.36, label="Site supported", color="#4f7093")
    axes[0].bar(locations + 0.18, summary_frame["supported_phosphosite_and_parent_equivalent_rows"],
                width=0.36, label="Site supported + parent equivalent", color="#c36d4b")
    readable_times = {
        "post_15_30_45_min": "15–45 min",
        "post_3.5_4_hr": "3.5–4 h",
        "post_24_hr": "24 h",
    }
    display_times = [readable_times.get(str(t), str(t).replace("_", " "))
                     for t in summary_frame["timepoint"]]
    axes[0].set_xticks(locations, display_times)
    axes[0].set_ylabel("Phosphosite feature rows")
    axes[0].set_title("Evidence across time")
    axes[0].legend(frameon=False, fontsize=8)
    first = matched.loc[matched["timepoint"].eq(args.times[0])]
    axes[1].scatter(first["log2_fc_site"], first["log2_fc_parent"],
                    s=6, alpha=0.12, color="#4f7093", rasterized=True)
    candidates = first.loc[first["candidate"]]
    axes[1].scatter(candidates["log2_fc_site"], candidates["log2_fc_parent"],
                    s=10, alpha=0.55, color="#c36d4b", rasterized=True)
    axes[1].axhline(args.equivalence_margin, color="gray", lw=0.8, ls="--")
    axes[1].axhline(-args.equivalence_margin, color="gray", lw=0.8, ls="--")
    axes[1].axvline(0, color="gray", lw=0.7)
    axes[1].set(xlabel="Phosphosite log2 FC", ylabel="Parent protein log2 FC",
                title=readable_times.get(args.times[0], args.times[0].replace("_", " "))
                + " matched effects")
    figure.savefig(args.out_dir / "ptm_parent_overview.png", dpi=170)
    plt.close(figure)
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
