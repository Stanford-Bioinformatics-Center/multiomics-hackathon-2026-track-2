"""Preview a user metabolite list against the committed c2.0 MoTrPAC export.

The input format is the Metabolomics module's documented 11-column list format.
This adapter reads precomputed contrast estimates; it fits no statistical model and
never changes the source output or its BH correction family.
"""
from __future__ import annotations

import csv
import gzip
import io
import math
import re
from collections import defaultdict
from functools import lru_cache
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3] / "MoTrPAC Hackathon" / "Metabolomics"
EXPORT = ROOT / "data" / "processed" / "motrpac_metabolite_contrasts.csv.gz"
REQUIRED = (
    "list_id", "source", "source_type", "comparison", "metabolite_name",
    "refmet_name", "pah_log2_effect", "pah_direction", "pah_q_value",
    "is_reported_hit", "caveat",
)
CONTRASTS = ("EE-CON", "EE-EE", "CON-CON")


def _key(name: str) -> str:
    """Same case/whitespace normalization as Metabolomics/scripts/common.py."""
    return re.sub(r"\s+", " ", name.strip().casefold())


def _present(value: str | None) -> bool:
    return value is not None and value.strip().upper() not in {"", "NA", "NAN", "NULL"}


def _number(value: str | None) -> float | None:
    try:
        number = float(value) if value not in (None, "") else float("nan")
        return number if math.isfinite(number) else None
    except ValueError:
        return None


@lru_cache(maxsize=1)
def _reference() -> dict[tuple[str, str], dict]:
    """Choose one assay feature per name/tissue without using exercise results.

    This repeats the published module's deterministic rule: targeted platform
    first, then highest average abundance, then feature ID as a tie-breaker.
    """
    choices: dict[tuple[str, str], tuple[tuple[int, float, str], str, str]] = {}
    feature_abundance: dict[tuple[str, str, str, str, str], list[float]] = defaultdict(list)
    source_rows: list[tuple[tuple[str, str], dict]] = []
    with gzip.open(EXPORT, "rt", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            if row["source_package_version"] != "2.0.8" or row["source_collection"] != "c2.0":
                raise ValueError("MoTrPAC metabolite export has an unexpected data version")
            name = row["refmet_name"] if _present(row["refmet_name"]) else row["feature_id"]
            key = (row["tissue"], _key(name))
            feature = row["feature_id"]
            basis = "refmet_name" if _present(row["refmet_name"]) else "feature_id"
            abundance = _number(row["AveExpr"])
            feature_key = (key[0], key[1], feature, row["platform"], basis)
            feature_abundance.setdefault(feature_key, [])
            if abundance is not None:
                feature_abundance[feature_key].append(abundance)
            source_rows.append((key, row))
    for (tissue, name_key, feature, platform, basis), abundances in feature_abundance.items():
        key = (tissue, name_key)
        rank = (1 if platform.startswith("metab-t-") else 0,
                sum(abundances) / len(abundances) if abundances else float("-inf"), feature)
        prior = choices.get(key)
        if prior is None or rank[:2] > prior[0][:2] or (rank[:2] == prior[0][:2] and feature < prior[0][2]):
            choices[key] = (rank, feature, basis)
    result = {key: {"feature_id": selected, "match_basis": basis, "contrasts": []}
              for key, (_, selected, basis) in choices.items()}
    for key, row in source_rows:
        selected = result[key]
        if row["feature_id"] != selected["feature_id"] or row["contrast_category"] not in CONTRASTS:
            continue
        selected["contrasts"].append({
            "contrast": row["contrast_category"], "timepoint": row["Timepoint"],
            "logFC": _number(row["logFC"]), "p_value": _number(row["p_value"]),
            "bh_q": _number(row["adj_p_value"]),
        })
    return result


def preview_metabolite_list(csv_text: str) -> dict:
    if len(csv_text.encode("utf-8")) > 2_000_000:
        raise ValueError("Metabolite CSV exceeds the 2 MB preview limit")
    reader = csv.DictReader(io.StringIO(csv_text, newline=""))
    if reader.fieldnames is None:
        raise ValueError("Metabolite CSV is empty")
    missing = [field for field in REQUIRED if field not in reader.fieldnames]
    if missing:
        raise ValueError(f"Metabolite CSV lacks required columns: {', '.join(missing)}")
    inputs = list(reader)
    if not inputs:
        raise ValueError("Metabolite CSV has no data rows")
    if len(inputs) > 1000:
        raise ValueError("Preview accepts at most 1,000 metabolite rows")

    reference = _reference()
    output = []
    unmapped = []
    for number, row in enumerate(inputs, start=2):
        direction = row["pah_direction"].strip().lower()
        if direction not in {"up", "down"}:
            raise ValueError(f"Row {number}: pah_direction must be up or down")
        name = row["refmet_name"].strip() or row["metabolite_name"].strip()
        if not name:
            raise ValueError(f"Row {number}: refmet_name or metabolite_name is required")
        matches = {tissue: reference.get((tissue, _key(name))) for tissue in ("blood", "muscle")}
        blood = matches["blood"]
        if blood is None:
            unmapped.append(name)
        observations = []
        for tissue, entry in matches.items():
            if entry is None:
                continue
            for contrast in entry["contrasts"]:
                observations.append({"tissue": tissue, "feature_id": entry["feature_id"],
                                     "match_basis": entry["match_basis"], **contrast})
        def changed(o: dict, tissue: str, contrast: str) -> bool:
            return (o["tissue"] == tissue and o["contrast"] == contrast
                    and o["bh_q"] is not None and o["bh_q"] < 0.05
                    and o["logFC"] is not None and abs(o["logFC"]) >= 0.5)

        exercise_sensitive = any(changed(o, "blood", "EE-CON") for o in observations)
        late = any(changed(o, "blood", "EE-CON") and
                   o["timepoint"] in {"post_3.5_4_hr", "post_24_hr"} for o in observations)
        control_drift = any(changed(o, "blood", "CON-CON") for o in observations)
        muscle_sensitive = any(changed(o, "muscle", "EE-CON") for o in observations)
        parts = []
        if exercise_sensitive:
            parts.append("exercise-sensitive, persists >=3.5 h" if late else
                         "exercise-sensitive, transient")
        if control_drift:
            parts.append("drifts without exercise")
        if matches["muscle"] is not None and muscle_sensitive:
            parts.append("also changes in muscle")
        label = ("no MoTrPAC blood match" if blood is None else
                 "; ".join(parts) or "stable in MoTrPAC blood")
        output.append({
            "input_row": number, "list_id": row["list_id"], "name": name,
            "pah_direction": direction, "pah_log2_effect": _number(row["pah_log2_effect"]),
            "pah_q_value": _number(row["pah_q_value"]),
            "is_reported_hit": row["is_reported_hit"].strip().lower() == "true",
            "blood_matched": blood is not None, "muscle_matched": matches["muscle"] is not None,
            "context_label": label, "observations": observations,
            "source": "Metabolomics/data/processed/motrpac_metabolite_contrasts.csv.gz",
            "caveat": row["caveat"],
        })
    return {
        "n_input_rows": len(inputs), "n_blood_matched": len(inputs) - len(unmapped),
        "unmapped_ids": unmapped, "rows": output,
        "rule": "BH q < 0.05 and |log2 FC| >= 0.5 in the committed c2.0 contrast export",
        "note": "Reference annotation only; uploaded disease effects are not tested against MoTrPAC or re-adjusted.",
    }
