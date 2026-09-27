#!/usr/bin/env python3
"""Build a within-tissue MoTrPAC RNA/protein catalog and simple held-out model.

The input is a long table of *summary-level effects*, not individual samples.
See README.md for the accepted schema and interpretation boundaries.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


ALIASES = {
    "layer": ("assay",),
    "log2_fc": ("logFC",),
    "timepoint": ("Timepoint",),
    "ci_lower": ("CI.L_calculated", "CI.L", "ci_l"),
    "ci_upper": ("CI.R_calculated", "CI.R", "ci_r"),
    "q_value": ("adj_p_value", "padj", "fdr"),
    "ave_expr": ("AveExpr",),
    "feature_id": ("entity_id", "identifier"),
}
LAYERS = {
    "transcript-rna-seq": "rna",
    "prot-pr": "protein",
    "prot-ph": "phosphosite",
    "rna": "rna",
    "protein": "protein",
    "phosphosite": "phosphosite",
}
REQUIRED = ("tissue", "contrast_category", "timepoint", "layer", "log2_fc", "gene_symbol")


def read_effects(path: Path) -> pd.DataFrame:
    data = pd.read_csv(path, low_memory=False)
    for canonical, alternatives in ALIASES.items():
        if canonical not in data:
            for alternative in alternatives:
                if alternative in data:
                    data[canonical] = data[alternative]
                    break
    missing = [column for column in REQUIRED if column not in data]
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")
    if "source_feature_id" in data:
        # The canonical feature_id is typed for cross-study matching; use the
        # original assay feature to resolve within-gene duplicate measurements.
        source_ids = data["source_feature_id"].replace("", np.nan)
        data["feature_id"] = source_ids.fillna(data.get("feature_id", ""))
    data["layer"] = data["layer"].map(LAYERS).fillna(data["layer"])
    for column in ("tissue", "contrast_category", "timepoint", "layer", "gene_symbol"):
        data[column] = data[column].fillna("").astype(str).str.strip()
    if "species" not in data:
        data["species"] = "human"  # Native MoTrPAC human package omits this field.
    data["species"] = data["species"].fillna("").astype(str).str.lower().str.strip()
    if "uniprot" not in data:
        data["uniprot"] = ""
    data["uniprot"] = data["uniprot"].fillna("").astype(str).str.strip()
    if "feature_id" not in data:
        # Stable gene/accession IDs work for already-collapsed reference tables.
        # A missing phosphosite ID cannot establish site identity over time.
        data["feature_id"] = [
            f"gene:{gene}" if layer == "rna" else
            f"uniprot:{accession or gene}" if layer == "protein" else
            f"unidentified_site_row_{i}"
            for i, (layer, gene, accession) in enumerate(
                zip(data["layer"], data["gene_symbol"], data["uniprot"])
            )
        ]
    data["feature_id"] = data["feature_id"].fillna("").astype(str).str.strip()
    for column in ("log2_fc", "ci_lower", "ci_upper", "q_value", "ave_expr"):
        if column not in data:
            data[column] = np.nan
        data[column] = pd.to_numeric(data[column], errors="coerce")
    return data


def filter_and_audit(data: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, dict]:
    chosen_times = [args.target_time, *args.earlier_times]
    filtered = data.loc[
        (data["species"] == "human")
        & (data["tissue"] == args.tissue)
        & (data["contrast_category"] == args.contrast_category)
        & (data["timepoint"].isin(chosen_times))
        & (data["layer"].isin(("rna", "protein", "phosphosite")))
    ].copy()
    if "contrast_type" in filtered and args.contrast_type:
        filtered = filtered.loc[filtered["contrast_type"] == args.contrast_type].copy()
    if filtered.empty:
        raise ValueError("No human RNA/protein/phosphosite rows match the chosen tissue, contrast, and times")
    if "contrast_short" in filtered:
        per_time = filtered.groupby("timepoint")["contrast_short"].nunique(dropna=True)
        if (per_time > 1).any():
            raise ValueError("Multiple contrast_short values found for one timepoint; split contrasts first")
    if not (filtered["layer"].eq("rna") & filtered["timepoint"].eq(args.target_time)).any():
        raise ValueError("No target-time RNA rows")
    if not (filtered["layer"].eq("protein") & filtered["timepoint"].eq(args.target_time)).any():
        raise ValueError("No target-time total-protein rows")
    if "biospecimen" in filtered:
        specimens = {
            layer: set(
                filtered.loc[
                    filtered["layer"].eq(layer) & filtered["timepoint"].eq(args.target_time),
                    "biospecimen",
                ].dropna().astype(str).str.strip().str.lower()
            ) - {""}
            for layer in ("rna", "protein")
        }
        if specimens["rna"] and specimens["protein"] and specimens["rna"].isdisjoint(specimens["protein"]):
            raise ValueError(
                f"Target RNA and protein biospecimens differ: {sorted(specimens['rna'])} vs {sorted(specimens['protein'])}"
            )

    n_before = len(filtered)
    bad_gene = filtered["gene_symbol"].eq("") | filtered["gene_symbol"].str.contains(r"[,;|]", regex=True)
    bad_protein_id = filtered["layer"].isin(("protein", "phosphosite")) & filtered["uniprot"].str.contains(r"[,;|]", regex=True)
    filtered = filtered.loc[~(bad_gene | bad_protein_id)].copy()
    # The package annotation join can give a feature more than one gene or parent protein.
    key = ["layer", "timepoint", "feature_id"]
    gene_count = filtered.groupby(key)["gene_symbol"].transform("nunique")
    protein_count = filtered.groupby(key)["uniprot"].transform("nunique")
    ambiguous = (gene_count > 1) | (filtered["layer"].isin(("protein", "phosphosite")) & (protein_count > 1))
    ambiguous_rows = int(ambiguous.sum())
    filtered = filtered.loc[~ambiguous].copy()

    identity = ["layer", "timepoint", "feature_id", "gene_symbol", "uniprot"]
    duplicates = filtered.duplicated(identity, keep=False)
    if duplicates.any():
        conflict = filtered.loc[duplicates].groupby(identity)["log2_fc"].nunique(dropna=False)
        if (conflict > 1).any():
            raise ValueError("Conflicting effect sizes for the same layer/time/feature/mapping")
    filtered = filtered.drop_duplicates(identity, keep="first")
    audit = {
        "input_rows": int(len(data)),
        "selected_rows_before_mapping_checks": n_before,
        "rows_with_missing_or_multi_token_mapping": int((bad_gene | bad_protein_id).sum()),
        "ambiguous_feature_mapping_rows": ambiguous_rows,
        "exact_duplicate_rows_removed": int(n_before - int((bad_gene | bad_protein_id).sum()) - ambiguous_rows - len(filtered)),
        "usable_rows": int(len(filtered)),
        "chosen_times": chosen_times,
        "source_package_versions": sorted(filtered.get("source_package_version", pd.Series(dtype=str)).dropna().astype(str).unique().tolist()),
    }
    return filtered, audit


def representative_features(rows: pd.DataFrame, layer: str) -> pd.DataFrame:
    """Choose one assay feature per gene(/UniProt), without using response strength."""
    part = rows.loc[rows["layer"] == layer].copy()
    keys = ["gene_symbol"] + (["uniprot"] if layer == "protein" else [])
    if part.empty:
        return part
    features = (
        part.groupby(keys + ["feature_id"], dropna=False)
        .agg(time_coverage=("timepoint", "nunique"), mean_abundance=("ave_expr", "median"))
        .reset_index()
    )
    features["mean_abundance"] = features["mean_abundance"].fillna(-np.inf)
    counts = features.groupby(keys, dropna=False)["feature_id"].nunique().rename(f"{layer}_feature_count").reset_index()
    chosen = (
        features.sort_values(
            keys + ["time_coverage", "mean_abundance", "feature_id"],
            ascending=[True] * len(keys) + [False, False, True],
        )
        .drop_duplicates(keys, keep="first")
        [keys + ["feature_id"]]
    )
    chosen = chosen.merge(counts, on=keys, how="left")
    return part.merge(chosen, on=keys + ["feature_id"], how="inner", validate="many_to_one")


def supported_sign(row: pd.Series, prefix: str, alpha: float, min_effect: float) -> int:
    effect = row.get(f"{prefix}_log2_fc", np.nan)
    if not np.isfinite(effect) or abs(effect) < min_effect:
        return 0
    q = row.get(f"{prefix}_q_value", np.nan)
    lo = row.get(f"{prefix}_ci_lower", np.nan)
    hi = row.get(f"{prefix}_ci_upper", np.nan)
    has_q = np.isfinite(q)
    has_ci = np.isfinite(lo) and np.isfinite(hi)
    if not (has_q or has_ci):
        return 0
    if has_q and q > alpha:
        return 0
    if has_ci and not ((effect > 0 and lo > 0) or (effect < 0 and hi < 0)):
        return 0
    return int(np.sign(effect))


def equivalent(row: pd.Series, prefix: str, margin: float) -> bool:
    lo = row.get(f"{prefix}_ci_lower", np.nan)
    hi = row.get(f"{prefix}_ci_upper", np.nan)
    # A 95% interval inside the margin is conservative evidence of equivalence.
    return bool(np.isfinite(lo) and np.isfinite(hi) and lo > -margin and hi < margin)


def make_catalog(rna: pd.DataFrame, protein: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    r = rna.loc[rna["timepoint"] == args.target_time].copy()
    p = protein.loc[protein["timepoint"] == args.target_time].copy()
    base = ["gene_symbol", "uniprot", "feature_id", "log2_fc", "ci_lower", "ci_upper", "q_value", "ave_expr"]
    r = r[[column for column in base + ["rna_feature_count"] if column in r]].drop(columns="uniprot")
    p = p[[column for column in base + ["protein_feature_count"] if column in p]]
    r = r.rename(columns={column: f"rna_{column}" for column in r if column not in ("gene_symbol", "rna_feature_count")})
    p = p.rename(columns={column: f"protein_{column}" for column in p if column not in ("gene_symbol", "uniprot", "protein_feature_count")})
    catalog = p.merge(r, on="gene_symbol", how="outer", validate="many_to_one")
    catalog["uniprot"] = catalog["uniprot"].fillna("")
    accession_count = p.groupby("gene_symbol")["uniprot"].nunique().rename("protein_accession_count_for_gene")
    catalog = catalog.merge(accession_count, on="gene_symbol", how="left")
    catalog["mapping_status"] = np.select(
        [
            catalog["protein_log2_fc"].isna(),
            catalog["rna_log2_fc"].isna(),
            catalog["protein_accession_count_for_gene"].gt(1),
            catalog["uniprot"].eq(""),
            catalog["protein_feature_count"].fillna(0).gt(1),
            catalog["rna_feature_count"].fillna(0).gt(1),
        ],
        [
            "no_protein_measurement", "no_rna_measurement", "multiple_uniprot_accessions",
            "protein_uniprot_missing", "multiple_protein_features_representative_selected",
            "multiple_rna_features_representative_selected",
        ],
        default="one_to_one",
    )
    categories = []
    for _, row in catalog.iterrows():
        r_sign = supported_sign(row, "rna", args.alpha, args.min_effect)
        p_sign = supported_sign(row, "protein", args.alpha, args.min_effect)
        if not np.isfinite(row.get("rna_log2_fc", np.nan)):
            label = "no_rna_measurement"
        elif not np.isfinite(row.get("protein_log2_fc", np.nan)):
            label = "no_protein_measurement"
        elif r_sign and p_sign:
            label = "supported_concordant" if r_sign == p_sign else "supported_opposite"
        elif r_sign and equivalent(row, "protein", args.equivalence_margin):
            label = "rna_response_protein_equivalent"
        elif p_sign and equivalent(row, "rna", args.equivalence_margin):
            label = "protein_response_rna_equivalent"
        else:
            label = "indeterminate"
        categories.append(label)
    catalog["event_class"] = categories
    catalog["tissue"] = args.tissue
    catalog["contrast_category"] = args.contrast_category
    catalog["timepoint"] = args.target_time
    return catalog.sort_values(["gene_symbol", "uniprot"]).reset_index(drop=True)


def make_early_features(
    rows: pd.DataFrame, rna: pd.DataFrame, protein: pd.DataFrame, args: argparse.Namespace
) -> tuple[pd.DataFrame, dict]:
    out = pd.DataFrame({"gene_symbol": sorted(rows["gene_symbol"].unique())})
    coverage = {}
    for index, time in enumerate(args.earlier_times, start=1):
        early_rna = rna.loc[rna["timepoint"] == time, ["gene_symbol", "log2_fc"]]
        name = f"rna_early_{index}_log2_fc"
        out = out.merge(early_rna.rename(columns={"log2_fc": name}), on="gene_symbol", how="left", validate="one_to_one")
        coverage[name] = int(out[name].notna().sum())

        # A site-minus-parent response is descriptive; its uncertainty/covariance
        # cannot be reconstructed from summary statistics.
        sites = rows.loc[(rows["layer"] == "phosphosite") & (rows["timepoint"] == time)].copy()
        parent = protein.loc[protein["timepoint"] == time, ["gene_symbol", "uniprot", "log2_fc"]]
        matched = sites.merge(parent, on=["gene_symbol", "uniprot"], how="inner", suffixes=("_site", "_parent"))
        ptm_name = f"ptm_early_{index}_site_minus_parent_maxabs"
        count_name = f"ptm_early_{index}_matched_site_count"
        if not matched.empty:
            matched["site_minus_parent"] = matched["log2_fc_site"] - matched["log2_fc_parent"]
            matched["absolute"] = matched["site_minus_parent"].abs()
            counts = matched.groupby("gene_symbol").size().rename(count_name).reset_index()
            picked = matched.sort_values(["gene_symbol", "absolute", "feature_id"], ascending=[True, False, True])
            picked = picked.drop_duplicates("gene_symbol")[["gene_symbol", "site_minus_parent"]]
            picked = picked.rename(columns={"site_minus_parent": ptm_name})
            out = out.merge(picked, on="gene_symbol", how="left", validate="one_to_one")
            out = out.merge(counts, on="gene_symbol", how="left", validate="one_to_one")
        else:
            out[ptm_name] = np.nan
            out[count_name] = 0
        out[count_name] = out[count_name].fillna(0).astype(int)
        coverage[ptm_name] = int(out[ptm_name].notna().sum())
    return out, coverage


def fit_oof(data: pd.DataFrame, feature_columns: list[str], folds: list[np.ndarray]) -> np.ndarray:
    predictions = np.full(len(data), np.nan)
    for test in folds:
        train = np.setdiff1d(np.arange(len(data)), test)
        model = make_pipeline(
            SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
            StandardScaler(),
            Ridge(alpha=10.0),
        )
        model.fit(data.iloc[train][feature_columns], data.iloc[train]["protein_log2_fc"])
        predictions[test] = model.predict(data.iloc[test][feature_columns])
    return predictions


def train_models(catalog: pd.DataFrame, early: pd.DataFrame, args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    eligible = catalog.loc[
        catalog["rna_log2_fc"].notna()
        & catalog["protein_log2_fc"].notna()
        & catalog["protein_accession_count_for_gene"].eq(1)
    ].copy()
    eligible = eligible.merge(early, on="gene_symbol", how="left", validate="many_to_one")
    eligible = eligible.reset_index(drop=True)
    eligible["rna_responsive"] = eligible.apply(
        lambda row: supported_sign(row, "rna", args.alpha, args.min_effect) != 0, axis=1
    )
    summary = {"eligible_rows": len(eligible), "eligible_genes": int(eligible["gene_symbol"].nunique())}
    if eligible["gene_symbol"].nunique() < max(20, args.folds * 3):
        summary["model_status"] = "skipped_insufficient_unique_genes"
        return eligible.iloc[0:0], pd.DataFrame(), summary
    genes = eligible["gene_symbol"].drop_duplicates().to_numpy()
    rng = np.random.default_rng(args.seed)
    rng.shuffle(genes)
    gene_folds = np.array_split(genes, args.folds)
    folds = [np.flatnonzero(eligible["gene_symbol"].isin(group).to_numpy()) for group in gene_folds]
    eligible["cv_fold"] = -1
    for fold_index, indices in enumerate(folds, start=1):
        eligible.loc[indices, "cv_fold"] = fold_index
    eligible["prediction_zero"] = 0.0
    eligible["prediction_rna_only"] = fit_oof(eligible, ["rna_log2_fc"], folds)
    additional = [
        column for column in eligible
        if column.startswith(("rna_early_", "ptm_early_"))
        and not column.endswith("_matched_site_count")
        and eligible[column].notna().any()
    ]
    if additional:
        eligible["prediction_temporal"] = fit_oof(eligible, ["rna_log2_fc", *additional], folds)
        summary["model_status"] = "fitted"
    else:
        eligible["prediction_temporal"] = np.nan
        summary["model_status"] = "rna_only_no_earlier_predictors"
    summary["temporal_features"] = additional
    eligible["observed_discordance_residual"] = eligible["rna_log2_fc"] - eligible["protein_log2_fc"]
    for name in ("zero", "rna_only", "temporal"):
        eligible[f"predicted_discordance_residual_{name}"] = (
            eligible["rna_log2_fc"] - eligible[f"prediction_{name}"]
        )
    metrics = []
    responsive = eligible["rna_responsive"].to_numpy()
    subsets = [("all_mapped", np.ones(len(eligible), dtype=bool))]
    summary["rna_responsive_eligible_genes"] = int(eligible.loc[responsive, "gene_symbol"].nunique())
    if responsive.sum() >= 20 and eligible.loc[responsive, "cv_fold"].nunique() >= 2:
        subsets.append(("rna_responsive", responsive))
        summary["rna_responsive_evaluation"] = "reported_from_same_all-gene_out_of_fold_predictions"
    else:
        summary["rna_responsive_evaluation"] = "skipped_fewer_than_20_or_one_fold"
    for subset_name, mask in subsets:
        observed = eligible.loc[mask, "protein_log2_fc"].to_numpy()
        for name in ("zero", "rna_only", "temporal"):
            predicted = eligible.loc[mask, f"prediction_{name}"].to_numpy()
            if not np.isfinite(predicted).all():
                continue
            metrics.append({
                "subset": subset_name,
                "model": name,
                "n_rows": int(mask.sum()),
                "n_genes": int(eligible.loc[mask, "gene_symbol"].nunique()),
                "mae": float(mean_absolute_error(observed, predicted)),
                "rmse": float(np.sqrt(mean_squared_error(observed, predicted))),
                "r2": float(r2_score(observed, predicted)),
            })
    return eligible, pd.DataFrame(metrics), summary


def plot_outputs(catalog: pd.DataFrame, predictions: pd.DataFrame, output: Path) -> None:
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    counts = catalog["event_class"].value_counts().sort_values()
    bars = axes[0].barh(counts.index, counts.values, color="#4f7093")
    for bar, count in zip(bars, counts.values):
        axes[0].text(bar.get_width() + max(counts.values) * 0.01, bar.get_y() + bar.get_height() / 2,
                     f"{count:,}", va="center", fontsize=8)
    axes[0].set_xlim(0, max(counts.values) * 1.18)
    axes[0].set_title("Evidence classes")
    axes[0].set_xlabel("Mapped entries")
    matched = catalog.dropna(subset=["rna_log2_fc", "protein_log2_fc"])
    axes[1].scatter(matched["rna_log2_fc"], matched["protein_log2_fc"], s=9, alpha=0.4, color="#4f7093")
    axes[1].axhline(0, color="gray", lw=0.7)
    axes[1].axvline(0, color="gray", lw=0.7)
    axes[1].set(xlabel="RNA log2 FC", ylabel="Protein log2 FC", title="Same-time effects")
    if not predictions.empty:
        column = "prediction_temporal" if predictions["prediction_temporal"].notna().any() else "prediction_rna_only"
        axes[2].scatter(predictions["protein_log2_fc"], predictions[column], s=9, alpha=0.4, color="#ad6446")
        finite = np.concatenate((predictions["protein_log2_fc"].to_numpy(), predictions[column].dropna().to_numpy()))
        low, high = float(np.nanmin(finite)), float(np.nanmax(finite))
        padding = max((high - low) * 0.04, 0.02)
        low, high = low - padding, high + padding
        axes[2].plot([low, high], [low, high], color="black", lw=0.8)
        axes[2].set_xlim(low, high)
        axes[2].set_ylim(low, high)
        axes[2].set(xlabel="Observed protein log2 FC", ylabel="Out-of-fold prediction", title=column.replace("prediction_", "") + " model")
    else:
        axes[2].text(0.5, 0.5, "Model unavailable\n(see run_summary.json)", ha="center", va="center", transform=axes[2].transAxes)
        axes[2].set_axis_off()
    figure.savefig(output, dpi=170)
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Canonical or package-native long-format MoTrPAC effects CSV(.gz)")
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--tissue", required=True)
    parser.add_argument("--contrast-category", default="EE-CON")
    parser.add_argument("--contrast-type", default="exercise_with_controls")
    parser.add_argument("--target-time", required=True)
    parser.add_argument("--earlier-times", nargs="*", default=[])
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--min-effect", type=float, default=0.2, help="Minimum absolute log2 FC for supported directional event")
    parser.add_argument("--equivalence-margin", type=float, default=0.2, help="Absolute log2 FC margin; requires CI wholly inside")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=2026)
    args = parser.parse_args()
    if args.folds < 2 or not (0 < args.alpha < 1) or args.min_effect < 0 or args.equivalence_margin <= 0:
        parser.error("Require folds >= 2, 0 < alpha < 1, min-effect >= 0, equivalence-margin > 0")
    if args.target_time in args.earlier_times or len(set(args.earlier_times)) != len(args.earlier_times):
        parser.error("Target and earlier times must be distinct")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows, mapping_audit = filter_and_audit(read_effects(args.input), args)
    rna = representative_features(rows, "rna")
    protein = representative_features(rows, "protein")
    catalog = make_catalog(rna, protein, args)
    early, coverage = make_early_features(rows, rna, protein, args)
    predictions, metrics, model_summary = train_models(catalog, early, args)
    catalog.to_csv(args.out_dir / "catalog.csv", index=False)
    predictions.to_csv(args.out_dir / "model_predictions.csv", index=False)
    metrics.to_csv(args.out_dir / "model_metrics.csv", index=False)
    plot_outputs(catalog, predictions, args.out_dir / "discordance_overview.png")
    summary = {
        "input": str(args.input.resolve()),
        "tissue": args.tissue,
        "contrast_category": args.contrast_category,
        "contrast_type": args.contrast_type,
        "target_time": args.target_time,
        "earlier_times": args.earlier_times,
        "alpha": args.alpha,
        "min_effect": args.min_effect,
        "equivalence_margin": args.equivalence_margin,
        "catalog_classes": {k: int(v) for k, v in catalog["event_class"].value_counts().items()},
        "mapping_audit": mapping_audit,
        "earlier_feature_coverage": coverage,
        "model": model_summary,
        "limitations": [
            "Summary-level effects cannot establish participant-level RNA-protein covariance or causality.",
            "Same-time q-values are adjusted within source assays; catalog classes are descriptive, not joint-FDR claims.",
            "Site-minus-parent protein log2 FC is descriptive, not measured phosphosite occupancy.",
            "One representative RNA/protein feature per gene or gene/UniProt is chosen by time coverage and abundance, never by response strength.",
            "The temporal model is associative and its coefficients are not mechanistic rate estimates.",
        ],
    }
    (args.out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"catalog_rows": len(catalog), "model_status": model_summary["model_status"], "out_dir": str(args.out_dir.resolve())}))


if __name__ == "__main__":
    main()
