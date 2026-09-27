"""Merge normalized MoTrPAC reference exports with release and duplicate audits."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from .adapter_legacy import TIMEPOINT_ORDER
from .schema import REFERENCE_COLUMNS, open_text

AUDIT_COLUMNS = (
    "key", "previous_input", "replacement_input", "same_numeric_values",
    "previous_log2_fc", "replacement_log2_fc", "previous_statistic",
    "replacement_statistic", "previous_q_value", "replacement_q_value",
)


def _key(row: dict[str, str]) -> tuple[str, ...]:
    """Source assay feature and contrast identify repeated exported results."""
    return tuple((row.get(field) or "").strip() for field in (
        "study_id", "species", "tissue", "layer", "assay", "platform",
        "source_feature_id", "timepoint", "contrast_category",
    ))


def merge(inputs: list[Path], output: Path, audit: Path, *,
          expected_package_version: str, expected_collection: str) -> dict:
    if len(inputs) < 2:
        raise ValueError("Pass at least two --input files")
    inputs = [Path(p) for p in inputs]
    output, audit = Path(output), Path(audit)
    rows: dict[tuple[str, ...], tuple[dict[str, str], str]] = {}
    columns: list[str] = []
    duplicates: list[dict[str, str]] = []
    counts = Counter()
    for path in inputs:
        with open_text(path, "rt") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or not set(REFERENCE_COLUMNS).issubset(reader.fieldnames):
                raise ValueError(f"{path}: not a normalized reference CSV")
            for field in reader.fieldnames:
                if field not in columns:
                    columns.append(field)
            for line, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(f"{path}:{line}: malformed CSV row")
                if row["study_id"] != "MoTrPAC" or row["species"] != "human":
                    raise ValueError(f"{path}:{line}: this merger accepts human MoTrPAC exports only")
                version = row.get("source_package_version", "")
                collection = row.get("source_collection", "")
                if (version != expected_package_version or
                        collection != expected_collection):
                    raise ValueError(
                        f"{path}:{line}: source release {version!r}/{collection!r} "
                        f"does not match {expected_package_version!r}/{expected_collection!r}")
                if not row.get("source_feature_id"):
                    raise ValueError(f"{path}:{line}: source_feature_id required for safe de-duplication")
                if not row.get("timepoint_order"):
                    row["timepoint_order"] = str(TIMEPOINT_ORDER.get(row["timepoint"], ""))
                key = _key(row)
                if key in rows:
                    previous, previous_path = rows[key]
                    numeric = ("log2_fc", "statistic", "p_value", "q_value")
                    same = all((previous.get(field) or "") == (row.get(field) or "")
                               for field in numeric)
                    duplicates.append({
                        "key": "|".join(key), "previous_input": previous_path,
                        "replacement_input": str(path),
                        "same_numeric_values": str(same).lower(),
                        "previous_log2_fc": previous.get("log2_fc", ""),
                        "replacement_log2_fc": row.get("log2_fc", ""),
                        "previous_statistic": previous.get("statistic", ""),
                        "replacement_statistic": row.get("statistic", ""),
                        "previous_q_value": previous.get("q_value", ""),
                        "replacement_q_value": row.get("q_value", ""),
                    })
                rows[key] = (row, str(path))  # Last input is the explicit winner.
                counts[path.name] += 1
    output.parent.mkdir(parents=True, exist_ok=True)
    audit.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output, "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row, _ in rows.values())
    with open_text(audit, "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=AUDIT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(duplicates)
    return {"input_rows": sum(counts.values()), "output_rows": len(rows),
            "duplicates_replaced": len(duplicates),
            "duplicates_with_different_numeric_values": sum(
                x["same_numeric_values"] == "false" for x in duplicates),
            "source_package_version": expected_package_version,
            "source_collection": expected_collection,
            "inputs": {str(path): counts[path.name] for path in inputs},
            "output": str(output.resolve()), "duplicate_audit": str(audit.resolve())}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, action="append", required=True,
                        help="Repeated in priority order; the last duplicate wins")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--expected-package-version", required=True)
    parser.add_argument("--expected-collection", required=True)
    args = parser.parse_args(argv)
    result = merge(args.input, args.output, args.audit,
                   expected_package_version=args.expected_package_version,
                   expected_collection=args.expected_collection)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
