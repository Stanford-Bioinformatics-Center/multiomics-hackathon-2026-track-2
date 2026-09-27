"""Turn a named Metabolomics Workbench study into a signed feature table.

This is a local-file adapter, not a study-design inference engine. The caller
must supply the contrast, scale, and (for repeated measures) real subject IDs.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


REQUIRED_OUTPUT = [
    "study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
    "log2_fc", "statistic", "p_value", "q_value", "contrast",
]
DESIGN_TYPES = {"unpaired", "paired_change", "difference_in_changes"}


def read_records(path: Path, kind: str) -> list[dict]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    records = list(payload.values()) if isinstance(payload, dict) else payload
    if not isinstance(records, list) or not records or not all(isinstance(r, dict) for r in records):
        raise ValueError(f"{kind} must be a nonempty JSON object or array of objects")
    return records


def parse_factors(value: object) -> dict[str, str]:
    if isinstance(value, dict):
        parts = value.items()
    elif isinstance(value, str):
        parts = []
        for token in re.split(r"\s+\|\s+", value.strip()):
            if ":" not in token:
                raise ValueError(f"Cannot parse sample factor {token!r}")
            parts.append(token.split(":", 1))
    else:
        raise ValueError("Each sample must have a factors string or object")
    result = {}
    for key, val in parts:
        key = str(key).strip().casefold()
        if not key or key in result:
            raise ValueError(f"Missing or repeated factor name {key!r}")
        result[key] = str(val).strip()
    return result


def required_text(obj: dict, key: str) -> str:
    value = obj.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Required nonempty string: {key}")
    return value.strip()


def clean_annotation(value: object) -> str | None:
    if value is None:
        return None
    result = str(value).strip()
    return None if result.casefold() in {"", "-", "na", "n/a", "null", "unknown"} else result


def selected_values(obj: dict, key: str) -> set[str]:
    value = obj.get(key)
    if not isinstance(value, list) or not value or not all(isinstance(v, str) and v for v in value):
        raise ValueError(f"{key} must be a nonempty list of group values")
    return set(value)


def factor_value(factors: dict, name: str, sample_id: str) -> str:
    key = name.strip().casefold()
    if key not in factors or not factors[key]:
        raise ValueError(f"Sample {sample_id} lacks required factor {name!r}")
    return factors[key]


def load_samples(records: list[dict], study_id: str) -> pd.DataFrame:
    rows = []
    for rec in records:
        if rec.get("study_id") != study_id:
            raise ValueError("Factors contain a different study_id")
        sample_id = required_text(rec, "local_sample_id")
        rows.append({"sample_id": sample_id, "factors": parse_factors(rec.get("factors"))})
    samples = pd.DataFrame(rows).set_index("sample_id", verify_integrity=True)
    return samples


def load_subject_ids(samples: pd.DataFrame, design: dict, design_path: Path) -> pd.Series:
    factor = design.get("subject_id_factor")
    csv_name = design.get("subject_map_csv")
    if bool(factor) == bool(csv_name):
        raise ValueError("Paired designs require exactly one subject_id_factor or subject_map_csv")
    if factor:
        values = {sid: factor_value(row.factors, factor, sid) for sid, row in samples.iterrows()}
        return pd.Series(values, name="subject_id")
    path = Path(csv_name)
    if not path.is_absolute():
        path = design_path.parent / path
    mapping = pd.read_csv(path, dtype=str, keep_default_na=False)
    if not {"sample_id", "subject_id"}.issubset(mapping.columns):
        raise ValueError("subject_map_csv needs sample_id and subject_id columns")
    if mapping.sample_id.duplicated().any() or (mapping.subject_id.str.strip() == "").any():
        raise ValueError("Subject map has repeated sample IDs or blank subject IDs")
    subject = mapping.set_index("sample_id").subject_id.reindex(samples.index)
    if subject.isna().any():
        raise ValueError(f"Subject IDs missing for {int(subject.isna().sum())} samples")
    return subject


def choose_samples(samples: pd.DataFrame, cfg: dict, design_path: Path) -> tuple[pd.DataFrame, str]:
    design = cfg.get("design")
    if not isinstance(design, dict):
        raise ValueError("design must be an object")
    kind = required_text(design, "type")
    if kind not in DESIGN_TYPES:
        raise ValueError(f"Unsupported design {kind!r}; supported: {sorted(DESIGN_TYPES)}")
    group_factor = required_text(design, "group_factor")
    case = selected_values(design, "case_values")
    control = selected_values(design, "control_values") if kind != "paired_change" else set()
    if case & control:
        raise ValueError("Case and control group values overlap")
    time_factor = design.get("time_factor")
    if kind != "unpaired" and not time_factor:
        raise ValueError("Paired designs require time_factor")
    if kind == "unpaired" and ("subject_id_factor" in design or "subject_map_csv" in design):
        raise ValueError("Use a paired design when subject IDs affect the contrast")
    rows = []
    for sid, row in samples.iterrows():
        group = factor_value(row.factors, group_factor, sid)
        arm = "case" if group in case else ("control" if group in control else "exclude")
        time = factor_value(row.factors, time_factor, sid) if time_factor and arm != "exclude" else None
        rows.append({"sample_id": sid, "arm": arm, "time": time})
    selected = pd.DataFrame(rows).set_index("sample_id")
    if kind == "unpaired":
        if time_factor:
            time_value = required_text(design, "time_value")
            selected = selected[selected.time == time_value]
        elif "time_value" in design:
            raise ValueError("time_value requires time_factor")
    else:
        baseline = required_text(design, "baseline_value")
        endpoint = required_text(design, "endpoint_value")
        if baseline == endpoint:
            raise ValueError("baseline_value and endpoint_value must differ")
        selected = selected[(selected.arm != "exclude") & selected.time.isin([baseline, endpoint])]
        selected["subject_id"] = load_subject_ids(samples.loc[selected.index], design, design_path)
        selected["time"] = selected.time.map({baseline: "baseline", endpoint: "endpoint"})
        pairs = selected[selected.arm != "exclude"].reset_index()
        if pairs.duplicated(["subject_id", "time"]).any():
            raise ValueError("A subject has multiple samples at one selected time")
        if (pairs.groupby("subject_id").arm.nunique() > 1).any():
            raise ValueError("A subject changes case/control group between times")
    selected = selected[selected.arm != "exclude"]
    if selected.empty or not (selected.arm == "case").any():
        raise ValueError("No case samples match the design")
    if kind != "paired_change" and not (selected.arm == "control").any():
        raise ValueError("No control samples match the design")
    return selected, kind


def parse_number(value: object, feature_id: str, sample_id: str) -> float:
    if value is None or value == "" or (isinstance(value, str) and value.strip().casefold() in {"na", "nan", "null"}):
        return math.nan
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Non-numeric value for {feature_id} / {sample_id}: {value!r}") from exc
    if not math.isfinite(number):
        return math.nan
    return number


def load_features(records: list[dict], study_id: str, samples: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata, values = [], []
    known = set(samples.index)
    for rec in records:
        if rec.get("study_id") != study_id:
            raise ValueError("Data contain a different study_id")
        analysis_id = required_text(rec, "analysis_id")
        metabolite_id = required_text(rec, "metabolite_id")
        feature_id = f"{analysis_id}:{metabolite_id}"
        data = rec.get("DATA")
        if not isinstance(data, dict):
            raise ValueError(f"Feature {feature_id} lacks sample-value DATA object")
        unknown = set(data) - known
        if unknown:
            raise ValueError(f"Feature {feature_id} has sample IDs absent from factors: {sorted(unknown)[:3]}")
        metadata.append({
            "feature_id": feature_id, "analysis_id": analysis_id,
            "refmet_name": clean_annotation(rec.get("refmet_name")),
            "refmet_id": clean_annotation(rec.get("refmet_id")),
            "metabolite_name": clean_annotation(rec.get("metabolite_name")),
            "units": clean_annotation(rec.get("units")),
        })
        values.append(pd.Series({sid: parse_number(v, feature_id, sid) for sid, v in data.items()}, name=feature_id))
    feature = pd.DataFrame(metadata).set_index("feature_id", verify_integrity=True)
    matrix = pd.DataFrame(values).reindex(columns=samples.index)
    matrix.index = feature.index
    return feature, matrix


def normalize(matrix: pd.DataFrame, feature: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    norm = cfg.get("normalization")
    if not isinstance(norm, dict):
        raise ValueError("normalization must be an object")
    scale = required_text(norm, "input_scale")
    if scale not in {"linear", "log2"}:
        raise ValueError("input_scale must be linear or log2")
    center = norm.get("sample_center", "none")
    if center not in {"none", "median_within_analysis"}:
        raise ValueError("Unsupported sample_center")
    pseudo = float(norm.get("pseudocount", 0))
    if not math.isfinite(pseudo) or pseudo < 0:
        raise ValueError("pseudocount must be finite and nonnegative")
    if scale == "log2" and pseudo != 0:
        raise ValueError("Do not add a pseudocount to already log2 data")
    result = matrix.copy()
    if scale == "linear":
        if ((result + pseudo) <= 0).any().any():
            raise ValueError("Linear measurements must be positive after the chosen pseudocount")
        result = np.log2(result + pseudo)
    if center == "median_within_analysis":
        for analysis_id, idx in feature.groupby("analysis_id").groups.items():
            block = result.loc[idx]
            result.loc[idx] = block.sub(block.median(axis=0), axis=1)
    return result


def bh(p_values: pd.Series) -> pd.Series:
    q = pd.Series(np.nan, index=p_values.index, dtype=float)
    valid = p_values.dropna().sort_values()
    if valid.empty:
        return q
    m = len(valid)
    adjusted = (valid.to_numpy() * m / np.arange(1, m + 1))[::-1]
    adjusted = np.minimum.accumulate(adjusted)[::-1]
    q.loc[valid.index] = np.minimum(adjusted, 1)
    return q


def paired_deltas(row: pd.Series, selected: pd.DataFrame, arm: str) -> tuple[pd.Series, int]:
    part = selected[selected.arm == arm].copy()
    if part.empty:
        return pd.Series(dtype=float), 0
    # Build design pairs from verified subject IDs before considering whether a
    # feature was measured. A lone baseline or endpoint is not a paired unit.
    wide = part.reset_index().pivot(index="subject_id", columns="time", values="sample_id")
    for time in ("baseline", "endpoint"):
        if time not in wide:
            return pd.Series(dtype=float), 0
    complete = wide.dropna(subset=["baseline", "endpoint"])
    delta = pd.Series(
        row.reindex(complete.endpoint.to_numpy()).to_numpy()
        - row.reindex(complete.baseline.to_numpy()).to_numpy(),
        index=complete.index,
    )
    return delta.dropna(), len(complete)


def estimate(row: pd.Series, selected: pd.DataFrame, kind: str, min_n: int, min_fraction: float) -> tuple:
    if kind == "unpaired":
        case = row.reindex(selected.index[selected.arm == "case"])
        control = row.reindex(selected.index[selected.arm == "control"])
        case_total, ctrl_total = len(case), len(control)
    else:
        case, case_total = paired_deltas(row, selected, "case")
        control, ctrl_total = (paired_deltas(row, selected, "control")
                               if kind == "difference_in_changes" else (pd.Series(dtype=float), 0))
    case, control = case.dropna(), control.dropna()
    enough_case = len(case) >= max(min_n, math.ceil(min_fraction * case_total))
    enough_ctrl = kind == "paired_change" or len(control) >= max(min_n, math.ceil(min_fraction * ctrl_total))
    if not enough_case or not enough_ctrl:
        return (math.nan, math.nan, math.nan, len(case), len(control),
                "insufficient_coverage", case_total, ctrl_total)
    if kind == "paired_change":
        effect = float(case.mean())
        test = stats.ttest_1samp(case, 0)
    else:
        effect = float(case.mean() - control.mean())
        test = stats.ttest_ind(case, control, equal_var=False)
    if not math.isfinite(float(test.statistic)) or not math.isfinite(float(test.pvalue)):
        return (effect, math.nan, math.nan, len(case), len(control),
                "not_estimable", case_total, ctrl_total)
    return (effect, float(test.statistic), float(test.pvalue), len(case), len(control),
            "tested", case_total, ctrl_total)


def run(factors_path: Path, data_path: Path, design_path: Path) -> pd.DataFrame:
    cfg = json.loads(design_path.read_text(encoding="utf-8"))
    study_id = required_text(cfg, "study_id")
    species = required_text(cfg, "species")
    if species not in {"human", "rat"}:
        raise ValueError("species must be 'human' or 'rat' for the current MoTrPAC reference")
    tissue = required_text(cfg, "tissue")
    biospecimen = required_text(cfg, "biospecimen") if "biospecimen" in cfg else None
    contrast = required_text(cfg, "contrast")
    min_n = cfg.get("min_samples_per_arm", 3)
    min_fraction = cfg.get("min_fraction_per_arm", 0.7)
    if not isinstance(min_n, int) or min_n < 2 or not isinstance(min_fraction, (float, int)) or not 0 < min_fraction <= 1:
        raise ValueError("min_samples_per_arm >= 2 and 0 < min_fraction_per_arm <= 1 are required")
    samples = load_samples(read_records(factors_path, "factors"), study_id)
    selected, kind = choose_samples(samples, cfg, design_path)
    data_records = read_records(data_path, "data")
    analysis_ids = cfg.get("analysis_ids")
    if analysis_ids is not None:
        if not isinstance(analysis_ids, list) or not analysis_ids or not all(isinstance(x, str) for x in analysis_ids):
            raise ValueError("analysis_ids must be a nonempty list of strings")
        data_records = [r for r in data_records if r.get("analysis_id") in analysis_ids]
        if not data_records:
            raise ValueError("No named features match analysis_ids")
    units = {clean_annotation(r.get("units")) for r in data_records}
    units.discard(None)
    if len(units) > 1 and not cfg.get("allow_mixed_units", False):
        raise ValueError(f"Different measurement units {sorted(units)}; select compatible analysis_ids")
    feature, raw = load_features(data_records, study_id, samples)
    matrix = normalize(raw, feature, cfg)
    rows = []
    for fid, values in matrix.iterrows():
        effect, statistic, p, n_case, n_control, status, n_design_case, n_design_control = (
            estimate(values, selected, kind, min_n, min_fraction))
        canonical_id = feature.at[fid, "refmet_id"] or feature.at[fid, "refmet_name"]
        namespace = "RefMet" if feature.at[fid, "refmet_id"] else (
            "RefMetName" if feature.at[fid, "refmet_name"] else None)
        rows.append({
            "study_id": study_id, "species": species, "tissue": tissue,
            "biospecimen": biospecimen,
            "layer": "metabolite", "feature_id": canonical_id,
            "id_namespace": namespace, "source_feature_id": fid,
            "log2_fc": effect, "statistic": statistic, "p_value": p,
            "contrast": contrast, "refmet_name": feature.at[fid, "refmet_name"],
            "refmet_id": feature.at[fid, "refmet_id"],
            "metabolite_name": feature.at[fid, "metabolite_name"],
            "analysis_id": feature.at[fid, "analysis_id"], "units": feature.at[fid, "units"],
            "n_case": n_case, "n_control": n_control,
            "n_design_case": n_design_case, "n_design_control": n_design_control,
            "status": status,
        })
    output = pd.DataFrame(rows)
    output["q_value"] = bh(output.p_value)
    key = output[["id_namespace", "feature_id"]].astype(str).agg("|".join, axis=1)
    counts = key.map(key.value_counts())
    output["duplicate_canonical_id"] = output.id_namespace.notna() & (counts > 1)
    return output[REQUIRED_OUTPUT + [c for c in output if c not in REQUIRED_OUTPUT]]


def split_output(result: pd.DataFrame, duplicate_policy: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create one-row-per-ID signature when explicitly requested.

    The policy uses coverage only; statistical magnitude and P values never
    determine which assay feature represents a repeated RefMet annotation.
    """
    matched = result[result.id_namespace.notna()].copy()
    unmatched = result[result.id_namespace.isna()].copy()
    duplicates = matched[matched.duplicate_canonical_id].copy()
    if duplicate_policy not in {"error", "keep_all", "highest_coverage"}:
        raise ValueError("duplicate_policy must be error, keep_all, or highest_coverage")
    if not duplicates.empty and duplicate_policy == "error":
        raise ValueError(f"{len(duplicates)} rows have duplicate canonical IDs; choose a reviewed duplicate_policy")
    if duplicate_policy == "highest_coverage":
        matched["coverage_n"] = matched.n_case + matched.n_control
        matched["eligible_for_test"] = matched.status == "tested"
        matched = (matched.sort_values(
            ["id_namespace", "feature_id", "eligible_for_test", "coverage_n", "source_feature_id"],
            ascending=[True, True, False, False, True])
            .drop_duplicates(["id_namespace", "feature_id"])
            .drop(columns=["coverage_n", "eligible_for_test"]))
        chosen = set(matched.source_feature_id)
        duplicates["chosen_for_signature"] = duplicates.source_feature_id.isin(chosen)
    else:
        duplicates["chosen_for_signature"] = duplicate_policy == "keep_all"
    return matched, unmatched, duplicates


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--factors", type=Path, required=True, help="Local Workbench /factors JSON")
    parser.add_argument("--data", type=Path, required=True, help="Local named-metabolite /data JSON")
    parser.add_argument("--design", type=Path, required=True, help="User-reviewed study design JSON")
    parser.add_argument("--output", type=Path, required=True, help="Output disease-signature CSV")
    args = parser.parse_args()
    result = run(args.factors, args.data, args.design)
    config = json.loads(args.design.read_text(encoding="utf-8"))
    matched, unmatched, duplicates = split_output(result, config.get("duplicate_policy", "error"))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    matched.to_csv(args.output, index=False)
    unmatched_path = args.output.with_name(args.output.stem + "_unmatched.csv")
    unmatched.to_csv(unmatched_path, index=False)
    duplicate_path = args.output.with_name(args.output.stem + "_duplicates.csv")
    duplicates.to_csv(duplicate_path, index=False)
    print(f"Wrote {len(matched)} mapped signature rows to {args.output}; "
          f"{len(unmatched)} unmapped rows to {unmatched_path}; "
          f"{len(duplicates)} duplicate-ID source rows to {duplicate_path}. "
          f"Tested {int((result.status == 'tested').sum())} source features.")


if __name__ == "__main__":
    main()
