"""Exploratory competitive rank-sum enrichment for a declared measured universe."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import numpy as np
from scipy.stats import norm, rankdata


REQUIRED_SIGNATURE = (
    "study_id", "species", "tissue", "layer", "feature_id", "id_namespace",
    "log2_fc", "statistic", "p_value", "q_value", "contrast",
)
REQUIRED_PROVENANCE = (
    "source", "inclusion_rule", "feature_id_policy", "gmt_source", "gmt_version",
)
OUTPUT_FIELDS = REQUIRED_SIGNATURE + (
    "biospecimen", "pathway_membership_version",
    "pathway_collection", "effect_scale", "statistic_type", "score_field",
    "rank_biserial", "n_set_members", "n_overlap", "n_universe",
    "gmt_description", "q_value_scope", "test_status",
)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if path.suffix == ".gz" else path.open(
        "rt", encoding="utf-8-sig", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def member_key(value: str) -> str:
    """Case/whitespace normalization only; this does not resolve synonyms."""
    return re.sub(r"\s+", " ", value.strip()).casefold()


def read_signature(path: Path, score_field: str) -> tuple[list[dict[str, str]], dict[str, str]]:
    with open_text(path) as stream:
        reader = csv.DictReader(stream)
        absent = set(REQUIRED_SIGNATURE) - set(reader.fieldnames or ())
        if absent:
            raise ValueError(f"Signature lacks canonical columns: {sorted(absent)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Signature is empty")
    context_fields = ("study_id", "species", "tissue", "layer", "id_namespace", "contrast")
    context = {}
    for field in context_fields:
        values = {row[field].strip() for row in rows}
        if len(values) != 1 or not next(iter(values)):
            raise ValueError(f"One nonempty {field} is required per run; observed {sorted(values)}")
        context[field] = next(iter(values))
    if "biospecimen" in rows[0]:
        biospecimens = {row["biospecimen"].strip() for row in rows}
        if len(biospecimens) != 1:
            raise ValueError("One biospecimen is required per ranked signature")
        context["biospecimen"] = next(iter(biospecimens))
    if context["layer"] not in {"rna", "metabolite"}:
        raise ValueError("Only full-rank RNA or metabolite signatures are supported")
    allowed = ({"hgnc", "ratgenesymbol", "rgd"} if context["layer"] == "rna" else
               {"refmet", "refmetname"})
    if context["id_namespace"].casefold() not in allowed:
        raise ValueError(f"Unsupported {context['layer']} feature namespace: {context['id_namespace']}")
    seen = set()
    for number, row in enumerate(rows, 2):
        key = member_key(row["feature_id"])
        if not key or key in seen:
            raise ValueError(f"Missing or duplicated feature_id at CSV line {number}")
        seen.add(key)
        try:
            score = float(row[score_field])
        except (TypeError, ValueError):
            raise ValueError(f"Missing numeric {score_field} at CSV line {number}") from None
        if not math.isfinite(score):
            raise ValueError(f"Nonfinite {score_field} at CSV line {number}")
    return rows, context


def read_provenance(path: Path, n_rows: int) -> dict:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError("Universe provenance must be a JSON object")
    absent = [field for field in REQUIRED_PROVENANCE
              if data.get(field) is None or not str(data.get(field, "")).strip()]
    if absent:
        raise ValueError(f"Universe provenance lacks: {absent}")
    if data.get("complete_measured_universe") is not True:
        raise ValueError("Universe provenance must declare complete_measured_universe=true")
    if "expected_feature_count" in data and data["expected_feature_count"] != n_rows:
        raise ValueError("expected_feature_count does not match signature rows")
    return data


def read_gmt(path: Path) -> list[tuple[str, str, set[str]]]:
    sets = []
    names = set()
    with open_text(path) as stream:
        for line_number, line in enumerate(stream, 1):
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) < 3 or not fields[0].strip():
                raise ValueError(f"Malformed GMT line {line_number}")
            name = fields[0].strip()
            if member_key(name) in names:
                raise ValueError(f"Repeated pathway name at GMT line {line_number}: {name}")
            names.add(member_key(name))
            members = {member_key(value) for value in fields[2:] if member_key(value)}
            if not members:
                raise ValueError(f"Pathway has no members at GMT line {line_number}")
            sets.append((name, fields[1].strip(), members))
    if not sets:
        raise ValueError("GMT contains no pathways")
    return sets


def bh_adjust(p_values: list[float]) -> list[float]:
    if not p_values:
        return []
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values))
    tail = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        position = int(order[rank_index])
        tail = min(tail, float(p_values[position]) * len(order) / (rank_index + 1))
        adjusted[position] = tail
    return adjusted.tolist()


def write_csv(path: Path, rows: list[dict], fields: tuple[str, ...]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(signature: Path, gmt: Path, provenance: Path, out_dir: Path, *,
        gmt_namespace: str, pathway_namespace: str, collection: str,
        score_field: str = "statistic", min_size: int = 10,
        max_size: int = 500) -> dict:
    signature, gmt, provenance, out_dir = map(Path, (signature, gmt, provenance, out_dir))
    if score_field not in {"statistic", "log2_fc"}:
        raise ValueError("score_field must be statistic or log2_fc")
    if min_size < 2 or max_size < min_size:
        raise ValueError("Require 2 <= min_size <= max_size")
    if pathway_namespace not in {"GO", "BIOCARTA", "Pathway"} or not collection.strip():
        raise ValueError("Declare pathway namespace (GO/BIOCARTA/Pathway) and collection")
    rows, context = read_signature(signature, score_field)
    if gmt_namespace.casefold() != context["id_namespace"].casefold():
        raise ValueError("GMT member namespace must equal signature id_namespace; crosswalk upstream")
    source = read_provenance(provenance, len(rows))
    sets = read_gmt(gmt)
    values = np.asarray([float(row[score_field]) for row in rows])
    ranked = rankdata(values, method="average")
    n = len(ranked)
    if n < min_size + 2:
        raise ValueError("Measured universe must leave at least two out-of-set features")
    _, tie_counts = np.unique(values, return_counts=True)
    tie_correction = 1 - float(np.sum(tie_counts.astype(float) ** 3 - tie_counts)) / (n ** 3 - n)
    if tie_correction <= 0:
        raise ValueError("All full-universe scores are tied")
    positions = {member_key(row["feature_id"]): i for i, row in enumerate(rows)}
    results = []
    p_values = []
    eligible_indices = []
    for name, description, members in sets:
        indices = [positions[member] for member in members if member in positions]
        m = len(indices)
        outside = n - m
        status = ("below_min_size" if m < min_size else
                  "above_max_size" if m > max_size else
                  "outside_universe_too_small" if outside < 2 else "eligible")
        item = dict.fromkeys(OUTPUT_FIELDS, "")
        item.update({
            "study_id": context["study_id"], "species": context["species"],
            "tissue": context["tissue"], "layer": "pathway", "feature_id": name,
            "id_namespace": pathway_namespace, "contrast": context["contrast"],
            "biospecimen": context.get("biospecimen", ""),
            "pathway_membership_version": source["gmt_version"],
            "pathway_collection": collection,
            "effect_scale": "competitive within-universe rank shift",
            "statistic_type": "asymptotic Wilcoxon rank-sum Z",
            "score_field": score_field, "n_set_members": len(members),
            "n_overlap": m, "n_universe": n, "gmt_description": description,
            "q_value_scope": "BH over eligible sets in this GMT and contrast",
            "test_status": status,
        })
        if status == "eligible":
            rank_sum = float(np.sum(ranked[indices]))
            u = rank_sum - m * (m + 1) / 2
            expected_u = m * outside / 2
            variance_u = m * outside * (n + 1) * tie_correction / 12
            z = (u - expected_u) / math.sqrt(variance_u)
            p = min(1.0, 2 * float(norm.sf(abs(z))))
            item.update(statistic=z, p_value=p,
                        rank_biserial=2 * u / (m * outside) - 1)
            eligible_indices.append(len(results))
            p_values.append(p)
        results.append(item)
    for index, q in zip(eligible_indices, bh_adjust(p_values)):
        results[index]["q_value"] = q
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "pathway_tests.csv", results, OUTPUT_FIELDS)
    write_csv(out_dir / "pathway_signature.csv",
              [results[index] for index in eligible_indices], OUTPUT_FIELDS)
    summary = {
        "signature": str(signature.resolve()), "signature_sha256": sha256(signature),
        "gmt": str(gmt.resolve()), "gmt_sha256": sha256(gmt),
        "universe_provenance": str(provenance.resolve()),
        "universe_provenance_sha256": sha256(provenance),
        "user_universe_declaration": source, "universe_verification":
            "declared by user; completeness cannot be proved from a rank table",
        "source_context": context, "gmt_member_namespace": gmt_namespace,
        "pathway_namespace": pathway_namespace, "pathway_collection": collection,
        "score_field": score_field, "n_measured_features": n,
        "n_gmt_sets": len(sets), "n_eligible_sets": len(eligible_indices),
        "min_size": min_size, "max_size": max_size,
        "test": "two-sided competitive Wilcoxon rank-sum normal approximation with tie correction",
        "null": "pathway feature identities exchangeable within the supplied measured universe",
        "q_value_scope": "Benjamini-Hochberg over all eligible sets in this GMT and contrast",
        "limits": [
            "Gene or metabolite correlation is not modeled; p and q values are exploratory.",
            "Feature ID matching is case/whitespace-normalized only; no synonym or RefMet crosswalk.",
            "The set statistic is a rank shift, not a log2 fold change or CAMERA/PAGE statistic.",
        ],
    }
    (out_dir / "manifest.json").write_text(json.dumps(summary, indent=2) + "\n",
                                            encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--gmt", type=Path, required=True)
    parser.add_argument("--universe-provenance", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--gmt-namespace", required=True,
                        help="Exact identifier namespace of GMT members, e.g. HGNC or RefMetName")
    parser.add_argument("--pathway-namespace", choices=("GO", "BIOCARTA", "Pathway"),
                        default="Pathway")
    parser.add_argument("--collection", required=True,
                        help="Set collection identity, e.g. C5/GOBP or RefMet2024")
    parser.add_argument("--score-field", choices=("statistic", "log2_fc"),
                        default="statistic")
    parser.add_argument("--min-size", type=int, default=10)
    parser.add_argument("--max-size", type=int, default=500)
    args = parser.parse_args()
    result = run(args.signature, args.gmt, args.universe_provenance, args.out_dir,
                 gmt_namespace=args.gmt_namespace,
                 pathway_namespace=args.pathway_namespace, collection=args.collection,
                 score_field=args.score_field, min_size=args.min_size,
                 max_size=args.max_size)
    print(f"Tested {result['n_eligible_sets']} of {result['n_gmt_sets']} sets "
          f"against {result['n_measured_features']} measured features")


if __name__ == "__main__":
    main()
