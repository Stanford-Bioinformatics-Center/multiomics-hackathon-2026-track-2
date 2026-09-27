"""Render auditable plots from query_core result tables.

Usage: python plot_results.py --results-dir DIR --out-dir DIR [--top N]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


GROUP = (
    "disease_study_id", "disease_contrast", "species", "reference_species",
    "tissue", "reference_study_id", "reference_contrast_category",
)
LAYER_ORDER = {"rna": 0, "protein": 1, "phosphosite": 2, "metabolite": 3, "pathway": 4}


def read_table(path: Path, required: tuple[str, ...]) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = set(required) - set(frame)
    if missing:
        raise ValueError(f"{path.name} lacks columns {sorted(missing)}")
    return frame


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def as_float(value: object) -> float | None:
    try:
        x = float(value)
    except (ValueError, TypeError):
        return None
    return x if math.isfinite(x) else None


def time_key(label: str, order: object, first_seen: int) -> tuple:
    declared = as_float(order)
    if declared is not None:
        return (0, declared, first_seen)
    name = str(label).lower()
    phase = 0 if "pre" in name or "baseline" in name else (
        1 if "during" in name else 2 if "post" in name else 3)
    match = re.search(r"(\d+(?:\.\d+)?)", name)
    if match:
        value = float(match.group(1))
        if re.search(r"(?:_hr|\bhr\b|\bhour|\d+h\b)", name):
            value *= 60
        return (1, phase, value, first_seen)
    return (2, phase, first_seen)


def group_token(values: tuple[str, ...], index: int) -> str:
    short = re.sub(r"[^A-Za-z0-9]+", "_", values[0]).strip("_")[:24] or "study"
    digest = hashlib.sha1("\x1f".join(values).encode()).hexdigest()[:8]
    return f"{index:03d}_{short}_{digest}"


def group_title(values: tuple[str, ...]) -> str:
    study, contrast, species, reference_species, tissue, ref_study, category = values
    return (f"{study}: {contrast}\n{tissue}; {species} disease vs "
            f"{reference_species} {ref_study} {category}")


def context_columns(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    for name in ("reference_assay", "reference_platform", "reference_timepoint_order",
                 "disease_gene_symbol"):
        if name not in result:
            result[name] = ""
    result["column_key"] = result[[
        "reference_layer", "reference_timepoint", "reference_assay", "reference_platform"
    ]].agg("\x1f".join, axis=1)
    result["source_order"] = np.arange(len(result))
    return result


def feature_label(layer: str, feature_id: str, gene_symbol: str, show_layer: bool) -> str:
    symbol = str(gene_symbol).strip()
    original = str(feature_id).strip()
    label = (f"{symbol} ({original})" if symbol and symbol.casefold() != original.casefold()
             else symbol or original)
    if len(label) > 54:
        label = label[:51] + "…"
    return f"{layer}: {label}" if show_layer else label


def select_bubbles(frame: pd.DataFrame, top: int, reference_field: str = "reference_log2_fc"
                   ) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Return plotted rows, every candidate with reason, and selection counts."""
    data = context_columns(frame)
    data["selection_reason"] = "eligible"
    data.loc[data.mapping_status != "unique", "selection_reason"] = "mapping_not_unique"
    data["reference_plot_value"] = numeric(data[reference_field])
    data.loc[(data.selection_reason == "eligible") &
             data.reference_plot_value.isna(), "selection_reason"] = f"{reference_field}_missing"

    # A unique query_core mapping is per reference group. Enforce a single
    # source and target identity within the exact plotted layer/time/assay.
    eligible = data[data.selection_reason == "eligible"]
    source_key = ["column_key", "disease_row_id"]
    target_key = ["column_key", "reference_row_id"]
    dup_source = eligible.duplicated(source_key, keep=False)
    dup_target = eligible.duplicated(target_key, keep=False)
    bad_index = eligible.index[dup_source | dup_target]
    data.loc[bad_index, "selection_reason"] = "not_one_to_one_in_plot_column"

    eligible = data[data.selection_reason == "eligible"].copy()
    eligible["disease_statistic_numeric"] = numeric(eligible.disease_statistic)
    eligible["disease_log2_fc_numeric"] = numeric(eligible.disease_log2_fc)
    eligible["selection_score"] = eligible.disease_statistic_numeric.abs().fillna(
        eligible.disease_log2_fc_numeric.abs())
    features = (eligible.sort_values("source_order")
                .drop_duplicates("disease_row_id")[[
                    "disease_row_id", "disease_layer", "disease_feature_id",
                    "disease_statistic_numeric", "disease_log2_fc_numeric", "selection_score",
                    "source_order",
                ]])
    features = features.sort_values(
        ["selection_score", "source_order"], ascending=[False, True], na_position="last")
    chosen = set(features.head(top).disease_row_id)
    data.loc[(data.selection_reason == "eligible") &
             ~data.disease_row_id.isin(chosen), "selection_reason"] = "outside_top_features"
    plotted = data[data.selection_reason == "eligible"].copy()
    plotted["reference_q_value_numeric"] = numeric(plotted.reference_q_value)
    summary = {
        "top_requested": top,
        "features_with_unique_numeric_matches": len(features),
        "features_displayed": len(chosen),
        "plotted_bubbles": len(plotted),
        "candidate_reason_counts": data.selection_reason.value_counts().to_dict(),
        "selection_rule": ("mapping_status=unique; one disease and reference identity per "
                           f"layer/time/assay/platform cell; finite {reference_field}; top "
                           "disease rows by |signed statistic|, falling back to |log2FC|"),
    }
    return plotted, data, summary


def bubble_size(q: float | None) -> float:
    if q is None or not math.isfinite(q):
        return 30.0
    return 55.0 + 24.0 * min(-math.log10(max(q, 1e-12)), 8.0)


def render_bubbles(frame: pd.DataFrame, out_dir: Path, title: str, top: int,
                   *, pathway: bool = False) -> dict:
    stem = "pathway_bubble" if pathway else "bubble"
    value_field = "reference_statistic" if pathway else "reference_log2_fc"
    plotted, selection, summary = select_bubbles(frame, top, value_field)
    selection.to_csv(out_dir / f"{stem}_selection.csv.gz", index=False)
    plotted.to_csv(out_dir / f"{stem}_plot_data.csv", index=False)
    if plotted.empty:
        summary["status"] = "no_eligible_bubbles"
        return summary

    # Keep source order for opaque time labels, but use a declared ordinal or
    # the physical numeric time when it is available.
    columns = plotted.drop_duplicates("column_key").copy()
    columns["sort_key"] = columns.apply(lambda r: (
        LAYER_ORDER.get(r.reference_layer, 99), r.reference_assay,
        r.reference_platform,
        time_key(r.reference_timepoint, r.reference_timepoint_order, int(r.source_order))), axis=1)
    columns = columns.sort_values("sort_key")
    xmap = {key: i for i, key in enumerate(columns.column_key)}
    duplicate_layer_time = plotted.groupby(["reference_layer", "reference_timepoint"]).column_key.nunique()
    labels = []
    for row in columns.itertuples():
        extra = (f" [{row.reference_platform or row.reference_assay}]"
                 if duplicate_layer_time[(row.reference_layer, row.reference_timepoint)] > 1 else "")
        labels.append(f"{row.reference_layer}{extra}\n{row.reference_timepoint}")
    disease_rows = plotted.drop_duplicates("disease_row_id")[[
        "disease_row_id", "disease_layer", "disease_feature_id", "disease_gene_symbol",
        "disease_statistic",
        "disease_log2_fc", "source_order",
    ]].copy()
    disease_rows["score"] = numeric(disease_rows.disease_statistic).abs().fillna(
        numeric(disease_rows.disease_log2_fc).abs())
    disease_rows = disease_rows.sort_values(["score", "source_order"],
                                            ascending=[False, True], na_position="last")
    ymap = {row.disease_row_id: i for i, row in enumerate(disease_rows.itertuples())}
    show_layer = disease_rows.disease_layer.nunique() > 1
    ylabels = [feature_label(r.disease_layer, r.disease_feature_id,
                             r.disease_gene_symbol, show_layer)
               for r in disease_rows.itertuples()]

    color = plotted.reference_plot_value.to_numpy(dtype=float)
    max_abs = max(float(np.nanmax(np.abs(color))), 0.1)
    norm = TwoSlopeNorm(vcenter=0, vmin=-max_abs, vmax=max_abs)
    q = plotted.reference_q_value_numeric.to_numpy(dtype=float)
    sizes = [bubble_size(v) for v in q]
    edge = ["#202733" if math.isfinite(v) and v < 0.05 else "#8c939c"
            if not math.isfinite(v) else "none" for v in q]
    fig, ax = plt.subplots(figsize=(max(8, len(columns) * 0.9 + 4),
                                    max(6.8 if pathway else 5.5,
                                        len(disease_rows) * 0.36 + 3.5)))
    scatter = ax.scatter(
        [xmap[k] for k in plotted.column_key],
        [ymap[k] for k in plotted.disease_row_id],
        c=color, cmap="coolwarm", norm=norm, s=sizes, edgecolors=edge,
        linewidths=0.8, alpha=0.9,
    )
    ax.set_xticks(range(len(columns)), labels, rotation=55, ha="right", fontsize=8)
    ax.set_yticks(range(len(disease_rows)), ylabels, fontsize=8)
    ax.set_xlim(-0.65, len(columns) - 0.35)
    ax.set_ylim(len(disease_rows) - 0.5, -0.6)
    ax.grid(axis="both", color="#e6e9ed", linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_xlabel("MoTrPAC reference layer and time")
    ax.set_ylabel("Top disease features")
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=11)
    cbar = fig.colorbar(scatter, ax=ax, shrink=0.63, pad=0.02)
    if pathway:
        ref_methods = sorted(set(plotted.get("reference_statistic_type", pd.Series(dtype=str))
                                 .dropna().astype(str)) - {""})
        disease_methods = sorted(set(plotted.get("disease_statistic_type", pd.Series(dtype=str))
                                     .dropna().astype(str)) - {""})
        color_label = ("MoTrPAC CAMERA Z" if ref_methods and
                       all("CAMERA" in method.upper() for method in ref_methods) else
                       f"MoTrPAC {ref_methods[0]}" if len(ref_methods) == 1 else
                       "MoTrPAC signed pathway statistic")
        methods_note = (f"Disease method: {', '.join(disease_methods) or 'unspecified'}; "
                        f"reference method: {', '.join(ref_methods) or 'unspecified'}. "
                        "Scores use distinct methods and are not fold changes.")
        reference_scopes = sorted(set(plotted.get("reference_scope", pd.Series(dtype=str))
                                      .dropna().astype(str)) - {""})
        if reference_scopes:
            methods_note += "\nReference scope: " + "; ".join(reference_scopes)
    else:
        color_label = "MoTrPAC log2 fold change"
        methods_note = "Separate cohorts; descriptive context only."
    cbar.set_label(color_label)
    handles = [ax.scatter([], [], s=bubble_size(qv), color="#8894a5",
                          edgecolors="#202733" if qv < 0.05 else "none",
                          linewidths=0.8, label=f"q={qv:g}")
               for qv in (0.1, 0.05, 0.01)]
    handles.append(ax.scatter([], [], s=bubble_size(None), color="#8894a5",
                              edgecolors="#8c939c", linewidths=0.8,
                              label="q unavailable"))
    ax.legend(handles=handles, title="Reference q / point area", frameon=False,
              loc="upper left", bbox_to_anchor=(1.15, 1))
    fig.text(0.01, 0.035,
             "Unique one-to-one mappings; top disease features by |signed statistic|, else |log2FC|.\n"
             "Dot area uses source MoTrPAC q; small grey-edged dots mean q unavailable.\n"
             + methods_note,
             fontsize=7.5, color="#4b5563")
    fig.subplots_adjust(bottom=0.38 if pathway else 0.28, top=0.82, right=0.78,
                        left=min(0.47, max(0.16, 0.06 + max(map(len, ylabels)) * 0.007)))
    fig.savefig(out_dir / f"{stem}_heatmap.png", dpi=180)
    fig.savefig(out_dir / f"{stem}_heatmap.svg")
    plt.close(fig)
    summary["status"] = "rendered"
    summary["column_order"] = labels
    summary["feature_order"] = ylabels
    summary["color_field"] = value_field
    summary["color_label"] = color_label
    if pathway:
        summary["disease_statistic_methods"] = disease_methods
        summary["reference_statistic_methods"] = ref_methods
        summary["reference_scopes"] = reference_scopes
    return summary


def _time_order(frame: pd.DataFrame, label_col: str) -> list[str]:
    first = frame.drop_duplicates(label_col).reset_index(drop=True)
    if "reference_timepoint_order" not in first:
        first["reference_timepoint_order"] = ""
    first["key"] = first.apply(lambda r: time_key(
        r[label_col], r.reference_timepoint_order, int(r.name)), axis=1)
    return first.sort_values("key")[label_col].tolist()


def render_rank(frame: pd.DataFrame, out_dir: Path, title: str) -> dict:
    data = frame.copy()
    data["spearman_rho_numeric"] = numeric(data.spearman_rho)
    data["plot_status"] = np.where(
        (data.status == "descriptive") & data.spearman_rho_numeric.notna(),
        "eligible", "not_descriptive_or_missing_rho")
    eligible = data[data.plot_status == "eligible"]
    repeated = eligible.duplicated(["reference_layer", "reference_timepoint"], keep=False)
    data.loc[eligible.index[repeated], "plot_status"] = "multiple_rows_same_layer_time"
    data.to_csv(out_dir / "rank_plot_data.csv", index=False)
    plotted = data[data.plot_status == "eligible"]
    if plotted.empty:
        return {"status": "no_descriptive_rank_points", "plotted_points": 0}
    order = _time_order(plotted, "reference_timepoint")
    xmap = {t: i for i, t in enumerate(order)}
    fig, ax = plt.subplots(figsize=(max(8, len(order) * 1.05 + 2), 5.5))
    colors = plt.get_cmap("tab10")
    for i, (layer, rows) in enumerate(sorted(plotted.groupby("reference_layer"),
                                              key=lambda item: LAYER_ORDER.get(item[0], 99))):
        rows = rows.sort_values("reference_timepoint", key=lambda s: s.map(xmap))
        ax.plot([xmap[t] for t in rows.reference_timepoint], rows.spearman_rho_numeric,
                marker="o", linewidth=1.7, markersize=5, label=layer, color=colors(i % 10))
    ax.axhline(0, color="#545d69", linewidth=0.9)
    ax.set_xticks(range(len(order)), order, rotation=40, ha="right")
    max_abs = min(1, max(0.2, float(plotted.spearman_rho_numeric.abs().max()) * 1.2))
    ax.set_ylim(-max_abs, max_abs)
    ax.set_ylabel("Spearman rho: disease vs MoTrPAC feature ranks")
    ax.set_xlabel("MoTrPAC reference time")
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=11)
    ax.grid(axis="y", color="#e6e9ed")
    ax.legend(title="Reference layer", frameon=False)
    rank_field = ", ".join(sorted(set(plotted.rank_field)))
    fig.text(0.01, 0.035,
             f"Rank field: {rank_field}. Unique one-to-one identity pairs only; "
             "rho is descriptive and has no inferential P value. Distinct cohorts.",
             fontsize=8, color="#4b5563")
    fig.subplots_adjust(bottom=0.27, top=0.79)
    fig.savefig(out_dir / "rank_timecourse.png", dpi=180)
    fig.savefig(out_dir / "rank_timecourse.svg")
    plt.close(fig)
    return {"status": "rendered", "plotted_points": len(plotted), "time_order": order,
            "layers": sorted(set(plotted.reference_layer))}


def render_coverage(frame: pd.DataFrame, out_dir: Path, title: str) -> dict:
    data = frame.copy()
    for col in ("n_unique", "n_ambiguous", "n_unmatched"):
        data[col] = numeric(data[col])
    data["plot_status"] = "eligible"
    repeated = data.duplicated(["reference_layer", "reference_timepoint"], keep=False)
    data.loc[repeated, "plot_status"] = "multiple_rows_same_layer_time"
    data.to_csv(out_dir / "coverage_plot_data.csv", index=False)
    plotted = data[data.plot_status == "eligible"].copy()
    if plotted.empty:
        return {"status": "no_unique_coverage_groups", "plotted_bars": 0}
    plotted["source_order"] = np.arange(len(plotted))
    if "reference_timepoint_order" not in plotted:
        plotted["reference_timepoint_order"] = ""
    plotted["sort_key"] = plotted.apply(lambda r: (
        LAYER_ORDER.get(r.reference_layer, 99),
        time_key(r.reference_timepoint, r.reference_timepoint_order, int(r.source_order))), axis=1)
    plotted = plotted.sort_values("sort_key")
    labels = [f"{r.reference_layer}\n{r.reference_timepoint}" for r in plotted.itertuples()]
    x = np.arange(len(plotted))
    fig, ax = plt.subplots(figsize=(max(8, len(plotted) * 0.7 + 2), 5.5))
    bottom = np.zeros(len(plotted))
    for col, label, color in (
        ("n_unique", "Unique", "#497e9b"),
        ("n_ambiguous", "Ambiguous", "#d7a55a"),
        ("n_unmatched", "Unmatched", "#c4cbd1"),
    ):
        vals = plotted[col].fillna(0).to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=color, label=label, width=0.75)
        bottom += vals
    ax.set_xticks(x, labels, rotation=50, ha="right", fontsize=8)
    ax.set_ylabel("Disease feature rows")
    fig.suptitle(title, x=0.08, y=0.985, ha="left", fontsize=11)
    ax.legend(frameon=False)
    ax.grid(axis="y", color="#e6e9ed")
    ax.set_axisbelow(True)
    fig.text(0.01, 0.035,
             "Coverage is counted separately for each MoTrPAC layer/time. "
             "Ambiguous candidates are not used in the bubble plot or rank correlation.",
             fontsize=8, color="#4b5563")
    fig.subplots_adjust(bottom=0.27, top=0.79)
    fig.savefig(out_dir / "coverage.png", dpi=180)
    fig.savefig(out_dir / "coverage.svg")
    plt.close(fig)
    return {"status": "rendered", "plotted_bars": len(plotted)}


def render(results_dir: Path, out_dir: Path, top: int = 25) -> dict:
    if top < 1:
        raise ValueError("--top must be positive")
    results_dir, out_dir = Path(results_dir), Path(out_dir)
    matches = read_table(results_dir / "matches.csv.gz", GROUP + (
        "mapping_status", "disease_row_id", "reference_row_id", "disease_layer",
        "disease_feature_id", "reference_layer", "reference_timepoint",
        "disease_statistic", "disease_log2_fc", "reference_log2_fc", "reference_q_value",
    ))
    ranks = read_table(results_dir / "rank_correlation.csv", GROUP + (
        "reference_layer", "reference_timepoint", "rank_field", "spearman_rho", "status",
    ))
    coverage = read_table(results_dir / "summary_by_group.csv", GROUP + (
        "reference_layer", "reference_timepoint", "n_unique", "n_ambiguous", "n_unmatched",
    ))
    group_frames = [frame[list(GROUP)].drop_duplicates() for frame in (matches, ranks, coverage)
                    if not frame.empty]
    if not group_frames:
        raise ValueError("No disease/reference groups in query_core outputs")
    group_values = sorted(set(tuple(row) for row in pd.concat(group_frames).itertuples(index=False, name=None)))
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"results_dir": str(results_dir.resolve()), "top": top, "groups": []}
    for i, values in enumerate(group_values, 1):
        token = group_token(values, i)
        dest = out_dir / token
        dest.mkdir(parents=True, exist_ok=True)
        mask = lambda frame: np.logical_and.reduce([
            frame[col].to_numpy() == val for col, val in zip(GROUP, values)])
        title = group_title(values)
        m = matches.loc[mask(matches)]
        r = ranks.loc[mask(ranks)]
        c = coverage.loc[mask(coverage)]
        feature_matches = m[m.reference_layer != "pathway"]
        pathway_matches = m[m.reference_layer == "pathway"]
        result = {"group": dict(zip(GROUP, values)), "directory": token,
                  "bubble": render_bubbles(feature_matches, dest, title, top)
                            if not feature_matches.empty else {"status": "no_matches"},
                  "pathway_bubble": render_bubbles(pathway_matches, dest, title, top,
                                                    pathway=True)
                            if not pathway_matches.empty else {"status": "no_matches"},
                  "rank": render_rank(r, dest, title) if not r.empty else
                          {"status": "no_rank_rows"},
                  "coverage": render_coverage(c, dest, title) if not c.empty else
                              {"status": "no_coverage_rows"}}
        (dest / "plot_metadata.json").write_text(json.dumps(result, indent=2) + "\n",
                                                 encoding="utf-8")
        manifest["groups"].append(result)
    (out_dir / "plot_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n",
                                                encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-dir", type=Path, required=True,
                        help="query_core output directory")
    parser.add_argument("--out-dir", type=Path, required=True,
                        help="Output directory for plots and plotted-data CSVs")
    parser.add_argument("--top", type=int, default=25,
                        help="Maximum disease features in each bubble heatmap (default 25)")
    args = parser.parse_args()
    manifest = render(args.results_dir, args.out_dir, args.top)
    print(f"Rendered {len(manifest['groups'])} disease/reference plot groups in {args.out_dir}")


if __name__ == "__main__":
    main()
