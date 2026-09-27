#!/usr/bin/env python3
"""Test whether a source-condition model transfers to a second condition.

This is an exploratory stress test. It trains only on source responses and
evaluates on the other condition's measured protein effects.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


def train_predict(train: pd.DataFrame, test: pd.DataFrame, features: list[str]) -> np.ndarray:
    model = make_pipeline(
        SimpleImputer(strategy="median", add_indicator=True, keep_empty_features=True),
        StandardScaler(),
        Ridge(alpha=10.0),
    )
    model.fit(train[features], train["protein_log2_fc"])
    return model.predict(test[features])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--train-dir", type=Path, required=True, help="Source build_discordance.py output directory")
    parser.add_argument("--test-dir", type=Path, required=True, help="Target build_discordance.py output directory")
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    source_info = json.loads((args.train_dir / "run_summary.json").read_text(encoding="utf-8"))
    target_info = json.loads((args.test_dir / "run_summary.json").read_text(encoding="utf-8"))
    for field in ("tissue", "target_time", "earlier_times"):
        if source_info[field] != target_info[field]:
            parser.error(f"Source and target have different {field}")
    if source_info["contrast_category"] == target_info["contrast_category"]:
        parser.error("Source and target contrast categories must differ")
    source = pd.read_csv(args.train_dir / "model_predictions.csv")
    target = pd.read_csv(args.test_dir / "model_predictions.csv")
    if source.empty or target.empty:
        parser.error("A source or target model prediction table is empty")
    key = ["gene_symbol", "uniprot"]
    if source.duplicated(key).any() or target.duplicated(key).any():
        parser.error("Model tables must have unique gene/UniProt pairs")
    common = source[key].merge(target[key], on=key, how="inner")
    if len(common) < 20:
        parser.error("Fewer than 20 shared gene/UniProt pairs")
    source = source.merge(common, on=key, how="inner", validate="one_to_one")
    target = target.merge(common, on=key, how="inner", validate="one_to_one")
    source = source.sort_values(key).reset_index(drop=True)
    target = target.sort_values(key).reset_index(drop=True)
    temporal_features = source_info["model"].get("temporal_features", [])
    if not temporal_features:
        parser.error("Source run has no earlier temporal predictors")
    missing = [feature for feature in temporal_features if feature not in target]
    if missing:
        parser.error(f"Target lacks temporal features: {missing}")
    predictions = target[key + ["rna_log2_fc", "protein_log2_fc"]].copy()
    predictions["prediction_zero"] = 0.0
    predictions["prediction_rna_only"] = train_predict(source, target, ["rna_log2_fc"])
    predictions["prediction_temporal"] = train_predict(source, target, ["rna_log2_fc", *temporal_features])
    predictions["observed_discordance_residual"] = (
        predictions["rna_log2_fc"] - predictions["protein_log2_fc"]
    )
    for model in ("zero", "rna_only", "temporal"):
        predictions[f"predicted_discordance_residual_{model}"] = (
            predictions["rna_log2_fc"] - predictions[f"prediction_{model}"]
        )
    metrics = []
    observed = predictions["protein_log2_fc"].to_numpy()
    for model in ("zero", "rna_only", "temporal"):
        forecast = predictions[f"prediction_{model}"].to_numpy()
        metrics.append({
            "model": model,
            "n_shared_genes": int(predictions["gene_symbol"].nunique()),
            "mae": float(mean_absolute_error(observed, forecast)),
            "rmse": float(np.sqrt(mean_squared_error(observed, forecast))),
            "r2": float(r2_score(observed, forecast)),
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    predictions.to_csv(args.out_dir / "transfer_predictions.csv", index=False)
    pd.DataFrame(metrics).to_csv(args.out_dir / "transfer_metrics.csv", index=False)
    summary = {
        "train_contrast_category": source_info["contrast_category"],
        "test_contrast_category": target_info["contrast_category"],
        "tissue": source_info["tissue"],
        "target_time": source_info["target_time"],
        "earlier_times": source_info["earlier_times"],
        "shared_gene_uniprot_pairs": len(common),
        "temporal_features": temporal_features,
        "limitation": "Different exercise arms may share the same control group and all data come from the same release; this tests modality transfer, not external replication or causality.",
    }
    (args.out_dir / "transfer_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps({"shared_pairs": len(common), "out_dir": str(args.out_dir.resolve())}))


if __name__ == "__main__":
    main()
