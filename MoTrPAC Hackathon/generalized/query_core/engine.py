"""Pure-Python query engine for normalized disease and MoTrPAC contrasts."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

from .ortholog import OrthologMap
from .schema import identifiers, load_csv, number, open_text, sign_basis, sign_relation

GROUP_FIELDS = ("study_id", "species", "tissue", "layer", "timepoint",
                "contrast_category", "contrast")
MATCH_FIELDS = (
    "disease_row_id", "reference_row_id", "disease_study_id", "reference_study_id",
    "species", "reference_species", "tissue", "disease_contrast", "reference_contrast_category",
    "reference_contrast", "reference_timepoint", "reference_timepoint_order",
    "disease_layer", "reference_layer",
    "disease_feature_id", "disease_id_namespace", "disease_gene_symbol",
    "reference_feature_id", "reference_id_namespace", "reference_gene_symbol",
    "reference_source_feature_id", "reference_assay",
    "reference_platform", "reference_match_basis", "matched_on",
    "match_confidence", "match_relation", "mapping_status",
    "ortholog_mapping",
    "disease_biospecimen", "reference_biospecimen", "biospecimen_compatibility",
    "disease_log2_fc", "reference_log2_fc", "disease_statistic",
    "reference_statistic", "disease_p_value", "reference_p_value",
    "disease_q_value", "reference_q_value", "disease_q_value_scope",
    "disease_source_q_value", "disease_source_q_scope",
    "sign_relation", "evidence_status", "supported_direction",
    "disease_effect_scale", "reference_effect_scale", "disease_statistic_type",
    "reference_statistic_type", "disease_pathway_collection",
    "reference_pathway_collection", "disease_pathway_membership_version",
    "reference_pathway_membership_version", "sign_basis", "reference_scope",
    "reference_source_collection", "reference_source_package_version",
)
COVERAGE_FIELDS = (
    "disease_row_id", "study_id", "species", "tissue", "layer", "feature_id",
    "id_namespace", "contrast", "biospecimen", "supported_identifiers",
    "unique_groups", "ambiguous_groups", "candidate_rows", "coverage_status",
)
GROUP_SUMMARY_FIELDS = (
    "disease_study_id", "disease_contrast", "reference_study_id", "species",
    "reference_species", "tissue", "reference_layer", "reference_timepoint",
    "reference_contrast_category", "reference_contrast", "n_disease_rows",
    "n_unique", "n_ambiguous", "n_unmatched", "n_reference_rows",
)
RANK_FIELDS = GROUP_SUMMARY_FIELDS[:10] + (
    "rank_field", "n_unique_pairs", "n_used", "n_excluded_duplicate_identity",
    "n_excluded_missing_value", "spearman_rho", "status",
)


def group_key(row: dict[str, str]) -> tuple[str, ...]:
    return tuple(row[field] for field in GROUP_FIELDS)


def hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, columns: tuple[str, ...], rows) -> None:
    with open_text(path, "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _rank(values: list[float]) -> list[float]:
    """Average ranks for ties, including ties at zero."""
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    pos = 0
    while pos < len(order):
        end = pos + 1
        while end < len(order) and values[order[end]] == values[order[pos]]:
            end += 1
        average = (pos + 1 + end) / 2
        for i in order[pos:end]:
            ranks[i] = average
        pos = end
    return ranks


def spearman(x: list[float], y: list[float]) -> float | None:
    if len(x) < 3:
        return None
    a, b = _rank(x), _rank(y)
    ma, mb = sum(a) / len(a), sum(b) / len(b)
    va = sum((v - ma) ** 2 for v in a)
    vb = sum((v - mb) ** 2 for v in b)
    if va == 0 or vb == 0:
        return None
    return sum((u - ma) * (v - mb) for u, v in zip(a, b)) / math.sqrt(va * vb)


def _compatibility(disease: dict[str, str], reference: dict[str, str]) -> str:
    d = disease.get("biospecimen", "").strip().casefold()
    r = reference.get("biospecimen", "").strip().casefold()
    if not d or not r:
        return "unknown"
    return "same" if d == r else "different"


def _evidence(disease: dict[str, str], reference: dict[str, str],
              q_threshold: float, min_abs_log2_fc: float) -> tuple[str, str]:
    dq, rq = number(disease, "q_value"), number(reference, "q_value")
    dlfc, rlfc = number(disease, "log2_fc"), number(reference, "log2_fc")
    if dq is None or rq is None:
        return "missing_q", ""
    if dlfc is None or rlfc is None:
        if (disease["layer"] == reference["layer"] == "pathway" and
                min_abs_log2_fc == 0 and sign_basis(disease, reference) == "signed_statistic"):
            dlfc, rlfc = number(disease, "statistic"), number(reference, "statistic")
        else:
            return "missing_effect", ""
    if (dq < q_threshold and rq < q_threshold and
            abs(dlfc) > min_abs_log2_fc and abs(rlfc) > min_abs_log2_fc):
        relation = sign_relation(disease, reference)
        return "jointly_supported", ("concordant" if relation == "same" else
                                      "discordant" if relation == "opposite" else "")
    return "not_jointly_supported", ""


def _match_record(disease: dict[str, str], reference: dict[str, str],
                  identity_types: set[str], confidence: str, status: str,
                  q_threshold: float, min_abs_log2_fc: float) -> dict[str, str]:
    evidence_status, supported_direction = _evidence(
        disease, reference, q_threshold, min_abs_log2_fc)
    relation = ("ortholog_gene_link" if confidence == "ortholog" else
                ("pathway_identity" if disease.get("pathway_membership_version") and
                 disease.get("pathway_membership_version") ==
                 reference.get("pathway_membership_version") else
                 "pathway_name_overlap") if disease["layer"] == "pathway" else
                "metabolite_name_only" if disease["layer"] == "metabolite" and
                confidence == "name_only" else
                "protein_site_link" if {disease["layer"], reference["layer"]} ==
                {"protein", "phosphosite"} and "uniprot" in identity_types else
                "gene_linked_cross_layer" if disease["layer"] != reference["layer"]
                else "protein_linked_site" if disease["layer"] == "phosphosite" and
                "uniprot" in identity_types else
                "gene_linked_same_layer" if disease["layer"] in {"protein", "phosphosite"}
                else "same_analyte")
    return {
        "disease_row_id": disease["row_id"], "reference_row_id": reference["row_id"],
        "disease_study_id": disease["study_id"],
        "reference_study_id": reference["study_id"],
        "species": disease["species"], "reference_species": reference["species"],
        "tissue": disease["tissue"],
        "disease_contrast": disease["contrast"],
        "reference_contrast_category": reference["contrast_category"],
        "reference_contrast": reference["contrast"],
        "reference_timepoint": reference["timepoint"],
        "reference_timepoint_order": reference.get("timepoint_order", ""),
        "disease_layer": disease["layer"], "reference_layer": reference["layer"],
        "disease_feature_id": disease["feature_id"],
        "disease_id_namespace": disease["id_namespace"],
        "disease_gene_symbol": disease.get("gene_symbol", ""),
        "reference_feature_id": reference["feature_id"],
        "reference_id_namespace": reference["id_namespace"],
        "reference_gene_symbol": reference.get("gene_symbol", ""),
        "reference_source_feature_id": reference.get("source_feature_id", ""),
        "reference_assay": reference.get("assay", ""),
        "reference_platform": reference.get("platform", ""),
        "reference_match_basis": reference.get("match_basis", ""),
        "matched_on": ";".join(sorted(identity_types)),
        "match_confidence": confidence, "match_relation": relation,
        "mapping_status": status,
        "ortholog_mapping": (disease["species"] + ":" + (disease.get("gene_symbol") or disease["feature_id"]) +
                              " -> " + reference["species"] + ":" + reference.get("gene_symbol", "")
                              if confidence == "ortholog" else ""),
        "disease_biospecimen": disease.get("biospecimen", ""),
        "reference_biospecimen": reference.get("biospecimen", ""),
        "biospecimen_compatibility": _compatibility(disease, reference),
        "disease_log2_fc": disease["log2_fc"],
        "reference_log2_fc": reference["log2_fc"],
        "disease_statistic": disease["statistic"],
        "reference_statistic": reference["statistic"],
        "disease_p_value": disease["p_value"],
        "reference_p_value": reference["p_value"],
        "disease_q_value": disease["q_value"],
        "reference_q_value": reference["q_value"],
        "disease_q_value_scope": disease.get("q_value_scope", ""),
        "disease_source_q_value": disease.get("source_q_value", ""),
        "disease_source_q_scope": disease.get("source_q_scope", ""),
        "sign_relation": sign_relation(disease, reference),
        "evidence_status": evidence_status,
        "supported_direction": supported_direction,
        "disease_effect_scale": disease.get("effect_scale", ""),
        "reference_effect_scale": reference.get("effect_scale", ""),
        "disease_statistic_type": disease.get("statistic_type", ""),
        "reference_statistic_type": reference.get("statistic_type", ""),
        "disease_pathway_collection": disease.get("pathway_collection", ""),
        "reference_pathway_collection": reference.get("pathway_collection", ""),
        "disease_pathway_membership_version": disease.get("pathway_membership_version", ""),
        "reference_pathway_membership_version": reference.get("pathway_membership_version", ""),
        "sign_basis": sign_basis(disease, reference),
        "reference_scope": reference.get("reference_scope", ""),
        "reference_source_collection": reference.get("source_collection", ""),
        "reference_source_package_version": reference.get("source_package_version", ""),
    }


def run_query(disease_path: Path, reference_path: Path, out_dir: Path, *,
              species: str | None = None, reference_species: str | None = None,
              tissue: str | None = None,
              disease_contrast: str | None = None,
              reference_contrast_category: str | None = None,
              reference_contrast: str | None = None,
              timepoints: tuple[str, ...] = (), rank_field: str = "auto",
              q_threshold: float = 0.05,
              min_abs_log2_fc: float = 0.0,
              ortholog_map: Path | None = None) -> dict:
    if rank_field not in {"auto", "log2_fc", "statistic"}:
        raise ValueError("rank_field must be auto, log2_fc or statistic")
    if not 0 < q_threshold <= 1 or min_abs_log2_fc < 0:
        raise ValueError("q_threshold must be in (0,1], min_abs_log2_fc >= 0")
    disease_path, reference_path, out_dir = map(Path,
                                                 (disease_path, reference_path, out_dir))
    orthologs = OrthologMap(ortholog_map) if ortholog_map else None
    disease = load_csv(disease_path, reference=False)
    reference = load_csv(reference_path, reference=True)
    if species:
        disease = [r for r in disease if r["species"] == species.lower()]
    if reference_species:
        reference = [r for r in reference if r["species"] == reference_species.lower()]
    if tissue:
        disease = [r for r in disease if r["tissue"] == tissue.lower()]
        reference = [r for r in reference if r["tissue"] == tissue.lower()]
    if disease_contrast:
        disease = [r for r in disease if r["contrast"] == disease_contrast]
    if reference_contrast_category:
        reference = [r for r in reference if
                     r["contrast_category"] == reference_contrast_category]
    if reference_contrast:
        reference = [r for r in reference if r["contrast"] == reference_contrast]
    if timepoints:
        reference = [r for r in reference if r["timepoint"] in timepoints]
    if not disease:
        raise ValueError("No disease rows remain after filters")
    if not reference:
        raise ValueError("No reference rows remain after filters")

    # A query may combine human and rat releases, but must not pool different
    # known releases of the same reference study/species as one comparison.
    release_groups: dict[tuple[str, str], dict[str, object]] = defaultdict(
        lambda: {"package_versions": set(), "collections": set(),
                 "rows_without_package_version": 0, "rows_without_collection": 0}
    )
    for r in reference:
        group = release_groups[(r["study_id"], r["species"])]
        version = r.get("source_package_version", "")
        collection = r.get("source_collection", "")
        if version:
            group["package_versions"].add(version)
        else:
            group["rows_without_package_version"] += 1
        if collection:
            group["collections"].add(collection)
        else:
            group["rows_without_collection"] += 1
    release_audit = []
    for (study_id, ref_species), group in sorted(release_groups.items()):
        versions = sorted(group["package_versions"])
        collections = sorted(group["collections"])
        if len(versions) > 1 or len(collections) > 1:
            raise ValueError(
                f"Mixed known reference releases for {study_id}/{ref_species}: "
                f"versions={versions}, collections={collections}"
            )
        release_audit.append({
            "study_id": study_id, "species": ref_species,
            "package_versions": versions, "collections": collections,
            "rows_without_package_version": group["rows_without_package_version"],
            "rows_without_collection": group["rows_without_collection"],
        })

    ref_by_id = {r["row_id"]: r for r in reference}
    index: dict[tuple[str, str, str, str], list[str]] = defaultdict(list)
    groups: dict[tuple[str, ...], list[str]] = defaultdict(list)
    for r in reference:
        groups[group_key(r)].append(r["row_id"])
        high, low = identifiers(r)
        for kind, key in high | low:
            index[(r["species"], r["tissue"], kind, key)].append(r["row_id"])

    matches: list[dict[str, str]] = []
    coverage: list[dict[str, str]] = []
    unmatched: list[dict[str, str]] = []
    ambiguous: list[dict[str, str]] = []
    unique_by_group: dict[tuple[str, str, tuple[str, ...]], list[tuple[dict, dict]]] = defaultdict(list)
    unique_ids: dict[tuple[str, str, tuple[str, ...]], set[str]] = defaultdict(set)
    ambiguous_ids: dict[tuple[str, str, tuple[str, ...]], set[str]] = defaultdict(set)

    for d in disease:
        high, low = identifiers(d)
        found: dict[str, dict[tuple[str, ...], dict[str, set[str]]]] = {
            "high": defaultdict(lambda: defaultdict(set)),
            "ortholog": defaultdict(lambda: defaultdict(set)),
            "name_only": defaultdict(lambda: defaultdict(set)),
        }
        for confidence, keys in (("high", high), ("name_only", low)):
            for kind, key in keys:
                for ref_id in index.get((d["species"], d["tissue"], kind, key), ()):
                    rg = group_key(ref_by_id[ref_id])
                    found[confidence][rg][ref_id].add(kind)
        if orthologs:
            for kind, key in high:
                if kind != "gene_symbol":
                    continue
                for target_species, target_gene in orthologs.targets(d["species"], key):
                    for ref_id in index.get((target_species, d["tissue"],
                                             "gene_symbol", target_gene), ()):
                        rg = group_key(ref_by_id[ref_id])
                        found["ortholog"][rg][ref_id].add("explicit_ortholog")
        candidate_groups = (set(found["high"]) | set(found["ortholog"]) |
                            set(found["name_only"]))
        n_unique = n_ambiguous = n_candidates = 0
        for rg in sorted(candidate_groups):
            confidence = ("high" if rg in found["high"] else
                          "ortholog" if rg in found["ortholog"] else "name_only")
            candidates = found[confidence][rg]
            status = "unique" if len(candidates) == 1 else "ambiguous"
            n_candidates += len(candidates)
            if status == "unique":
                n_unique += 1
                unique_ids[(d["study_id"], d["contrast"], rg)].add(d["row_id"])
            else:
                n_ambiguous += 1
                ambiguous_ids[(d["study_id"], d["contrast"], rg)].add(d["row_id"])
            for ref_id, kinds in sorted(candidates.items(), key=lambda item: int(item[0])):
                r = ref_by_id[ref_id]
                record = _match_record(d, r, kinds, confidence, status,
                                       q_threshold, min_abs_log2_fc)
                matches.append(record)
                if status == "unique":
                    unique_by_group[(d["study_id"], d["contrast"], rg)].append((d, r))
                else:
                    ambiguous.append(record)
        row = {
            "disease_row_id": d["row_id"], "study_id": d["study_id"],
            "species": d["species"], "tissue": d["tissue"],
            "layer": d["layer"], "feature_id": d["feature_id"],
            "id_namespace": d["id_namespace"], "contrast": d["contrast"],
            "biospecimen": d.get("biospecimen", ""),
            "supported_identifiers": ";".join(sorted(
                f"{kind}:{key}" for kind, key in high | low)),
            "unique_groups": n_unique, "ambiguous_groups": n_ambiguous,
            "candidate_rows": n_candidates,
            "coverage_status": ("unmatched" if not candidate_groups else
                                "has_ambiguous" if n_ambiguous else "unique_only"),
        }
        coverage.append(row)
        if not candidate_groups:
            row = dict(row)
            row["unmatched_reason"] = ("no_supported_identifier" if not high and not low
                                       else "no_compatible_reference_identity")
            unmatched.append(row)

    # Summaries include zero-coverage rows, which are invisible in matches.csv.gz.
    disease_by_axis: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for d in disease:
        disease_by_axis[(d["study_id"], d["species"], d["tissue"], d["contrast"])].append(d)
    group_summary: list[dict[str, str]] = []
    ranks: list[dict[str, str]] = []
    for rg in sorted(groups):
        ref_study, ref_species, ref_tissue, ref_layer, ref_time, ref_category, ref_contrast = rg
        for (d_study, d_species, d_tissue, d_contrast), source_rows in sorted(disease_by_axis.items()):
            if (d_tissue != ref_tissue or
                    not (d_species == ref_species or
                         orthologs and orthologs.allows(d_species, ref_species))):
                continue
            key = (d_study, d_contrast, rg)
            n_unique = len(unique_ids[key])
            n_ambiguous = len(ambiguous_ids[key])
            base = {
                "disease_study_id": d_study,
                "disease_contrast": d_contrast,
                "reference_study_id": ref_study,
                "species": d_species, "reference_species": ref_species,
                "tissue": ref_tissue,
                "reference_layer": ref_layer,
                "reference_timepoint": ref_time,
                "reference_contrast_category": ref_category,
                "reference_contrast": ref_contrast,
            }
            group_summary.append({**base, "n_disease_rows": len(source_rows),
                                  "n_unique": n_unique, "n_ambiguous": n_ambiguous,
                                  "n_unmatched": len(source_rows) - n_unique - n_ambiguous,
                                  "n_reference_rows": len(groups[rg])})
            pairs = unique_by_group[key]
            source_counts = Counter((d["id_namespace"], d["feature_id"])
                                    for d, _ in pairs)
            reference_counts = Counter(r["row_id"] for _, r in pairs)
            valid_pairs = [(d, r) for d, r in pairs
                           if source_counts[(d["id_namespace"], d["feature_id"])] == 1
                           and reference_counts[r["row_id"]] == 1]
            excluded_duplicate = len(pairs) - len(valid_pairs)
            chosen_field = rank_field
            if rank_field == "auto":
                all_signed_stats = len(valid_pairs) >= 3 and all(
                    number(d, "statistic") is not None and
                    number(r, "statistic") is not None for d, r in valid_pairs)
                chosen_field = "statistic" if all_signed_stats else "log2_fc"
            numeric = [(number(d, chosen_field), number(r, chosen_field))
                       for d, r in valid_pairs]
            numeric = [(x, y) for x, y in numeric if x is not None and y is not None]
            rho = spearman([x for x, _ in numeric], [y for _, y in numeric])
            ranks.append({**base, "rank_field": chosen_field,
                          "n_unique_pairs": len(pairs), "n_used": len(numeric),
                          "n_excluded_duplicate_identity": excluded_duplicate,
                          "n_excluded_missing_value": len(valid_pairs) - len(numeric),
                          "spearman_rho": "" if rho is None else f"{rho:.12g}",
                          "status": ("fewer_than_3" if len(numeric) < 3 else
                                     "constant_rank" if rho is None else "descriptive")})

    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "matches.csv.gz", MATCH_FIELDS, matches)
    write_csv(out_dir / "coverage.csv", COVERAGE_FIELDS, coverage)
    write_csv(out_dir / "unmatched.csv", COVERAGE_FIELDS + ("unmatched_reason",), unmatched)
    write_csv(out_dir / "ambiguous.csv.gz", MATCH_FIELDS, ambiguous)
    write_csv(out_dir / "summary_by_group.csv", GROUP_SUMMARY_FIELDS, group_summary)
    write_csv(out_dir / "rank_correlation.csv", RANK_FIELDS, ranks)
    summary = {
        "disease_input": str(disease_path.resolve()),
        "disease_sha256": hash_file(disease_path),
        "reference_input": str(reference_path.resolve()),
        "reference_sha256": hash_file(reference_path),
        "reference_releases": release_audit,
        "filters": {"species": species, "reference_species": reference_species,
                    "tissue": tissue,
                    "disease_contrast": disease_contrast,
                    "reference_contrast_category": reference_contrast_category,
                    "reference_contrast": reference_contrast,
                    "timepoints": list(timepoints),
                    "ortholog_map": str(Path(ortholog_map).resolve()) if ortholog_map else None},
        "thresholds": {"q_threshold": q_threshold,
                       "min_abs_log2_fc": min_abs_log2_fc},
        "rank_field": rank_field,
        "counts": {"disease_rows": len(disease), "reference_rows": len(reference),
                   "matched_candidate_rows": len(matches),
                   "disease_rows_unmatched_everywhere": len(unmatched),
                   "ambiguous_candidate_rows": len(ambiguous),
                   "different_biospecimen_candidate_rows": sum(
                       m["biospecimen_compatibility"] == "different" for m in matches)},
        "interpretation": (
            "Signs compare distinct disease and reference contrasts. "
            "Joint support requires both source q values and effect thresholds; "
            "q above threshold is not evidence of no change. "
            "Rank correlations are descriptive and receive no inferential p value. "
            "A broad tissue match with different biospecimens is flagged per row. " +
            ("Cross-species matches use only the supplied explicit one-to-one ortholog map. "
             if orthologs else "No cross-species orthology is performed. ") +
            "No gene-metabolite pathway bridge or PTM occupancy inference is performed."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n",
                                           encoding="utf-8")
    return summary
