"""Convert the bundled legacy MoTrPAC exports into a typed reference CSV.

This adapter reads published summary tables already shipped with the legacy
repository. It does not download data, fit models, or require R.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

from .schema import REFERENCE_COLUMNS, clean_optional, open_text

OUTPUT_COLUMNS = REFERENCE_COLUMNS + (
    "gene_symbol", "uniprot", "refmet_id", "refmet_name", "pathway_collection",
    "effect_scale", "statistic_type", "biospecimen", "assay", "platform",
    "source_feature_id", "source_package_version", "source_collection",
    "reference_scope", "timepoint_order", "match_basis",
)

TIMEPOINT_ORDER = {
    "during_20_min": 1, "during_40_min": 2, "post_10_min": 3,
    "post_15_30_45_min": 4, "post_3.5_4_hr": 5, "post_24_hr": 6,
}

SOURCES = (
    ("muscle_rna", "data/processed/motrpac_muscle_rna_ee_con_24h.csv.gz"),
    ("muscle_protein", "data/processed/motrpac_muscle_protein_ee_con.csv.gz"),
    ("blood_genes", "data/processed/blood_gene_ranked.csv.gz"),
    ("metabolites", "Metabolomics/data/processed/motrpac_metabolite_contrasts.csv.gz"),
    ("gobp", "outputs/PAH blood/motrpac_blood_gobp_camera.csv.gz"),
    ("biocarta", "Transcriptomics/outputs/motrpac_biocarta_camera.csv"),
)


def _preferred_id(source: dict[str, str], layer: str) -> tuple[str, str]:
    if layer in {"rna", "protein"}:
        if layer == "protein" and source.get("uniprot"):
            return source["uniprot"], "UniProt"
        if source.get("gene_symbol"):
            return source["gene_symbol"], "HGNC"
    elif layer == "metabolite":
        if source.get("refmet_id"):
            return source["refmet_id"], "RefMet"
        if source.get("refmet_name"):
            return source["refmet_name"], "RefMetName"
    return source["feature_id"], "MoTrPAC"


def _read(path: Path):
    with open_text(path, "rt") as stream:
        yield from csv.DictReader(stream)


def _molecular(source: dict[str, str], kind: str,
               metabolite_feature_name_fallback: bool) -> dict[str, str]:
    assay = source.get("assay", "")
    if kind == "muscle_rna" or assay == "transcript-rna-seq":
        layer = "rna"
    elif kind in {"muscle_protein", "blood_genes"}:
        layer = "protein"
    else:
        layer = "metabolite"
    source = dict(source)
    for field in ("gene_symbol", "uniprot", "refmet_id", "refmet_name"):
        source[field] = clean_optional(source.get(field))
    if (layer == "metabolite" and metabolite_feature_name_fallback and
            not source.get("refmet_id") and not source.get("refmet_name")):
        source["refmet_name"] = source["feature_id"]
        match_basis = "feature_id_name_fallback"
    else:
        match_basis = ("uniprot" if layer == "protein" and source.get("uniprot") else
                       "gene_symbol" if layer in {"rna", "protein"} and
                       source.get("gene_symbol") else
                       "refmet_id" if layer == "metabolite" and source.get("refmet_id") else
                       "refmet_name" if layer == "metabolite" and
                       source.get("refmet_name") else "unmapped_source_feature")
    feature_id, namespace = _preferred_id(source, layer)
    category = source.get("contrast_category") or "EE-CON"
    return {
        "study_id": "MoTrPAC", "species": "human",
        "tissue": source["tissue"].lower(), "layer": layer,
        "feature_id": feature_id, "id_namespace": namespace,
        "log2_fc": source.get("logFC", ""),
        "statistic": source.get("z.std", ""),
        "p_value": source.get("p_value", ""),
        "q_value": source.get("adj_p_value", ""),
        "contrast": source.get("contrast_short") or category,
        "timepoint": source["Timepoint"],
        "contrast_category": category,
        "gene_symbol": source.get("gene_symbol", ""),
        "uniprot": source.get("uniprot", ""),
        "refmet_id": source.get("refmet_id", ""),
        "refmet_name": source.get("refmet_name", ""),
        "pathway_collection": "", "effect_scale": "log2 fold change",
        "statistic_type": "MoTrPAC z.std" if source.get("z.std") else "",
        "biospecimen": "whole blood" if layer == "rna" and
                       source["tissue"].lower() == "blood" else
                       "skeletal muscle" if source["tissue"].lower() == "muscle" else "",
        "assay": assay or "metab", "platform": source.get("platform", ""),
        "source_feature_id": source["feature_id"],
        "source_package_version": source.get("source_package_version", ""),
        "source_collection": source.get("source_collection", ""),
        "reference_scope": "published MoTrPAC summary export",
        "timepoint_order": TIMEPOINT_ORDER.get(source["Timepoint"], ""),
        "match_basis": match_basis,
    }


def _pathway(source: dict[str, str], kind: str) -> dict[str, str]:
    collection = source["collection"] + "/" + source["database"]
    is_go = kind == "gobp"
    return {
        "study_id": "MoTrPAC", "species": "human",
        "tissue": source["tissue"].lower(), "layer": "pathway",
        "feature_id": source["set"], "id_namespace": "GO" if is_go else "BIOCARTA",
        "log2_fc": "", "statistic": source["z.std"],
        "p_value": source["p_value"], "q_value": source["adj_p_value"],
        "contrast": source["contrast_short"],
        "timepoint": source.get("Timepoint") or source.get("timepoint", ""),
        "contrast_category": source.get("contrast_category") or
                             ("EE-CON" if source.get("modality") == "EE" else "RE-CON"),
        "gene_symbol": "", "uniprot": "", "refmet_id": "", "refmet_name": "",
        "pathway_collection": collection,
        "effect_scale": "CAMERA pathway Z", "statistic_type": "signed CAMERA z.std",
        "biospecimen": "whole blood", "assay": source.get("assay", ""),
        "platform": "", "source_feature_id": source["set"],
        "source_package_version": source.get("source_package_version", ""),
        "source_collection": source.get("source_collection") or (
            "c2.0" if source.get("source_package_version") == "2.0.8" else ""),
        "reference_scope": ("full bundled MoTrPAC GOBP table" if is_go else
                            "partial BioCarta table: six PAH-selected sets only"),
        "timepoint_order": TIMEPOINT_ORDER.get(
            source.get("Timepoint") or source.get("timepoint", ""), ""),
        "match_basis": "collection_and_set_label",
    }


def convert(legacy_root: Path, output: Path, *, include_biocarta: bool = False,
            skip_missing: bool = False,
            metabolite_feature_name_fallback: bool = False) -> Counter:
    legacy_root, output = Path(legacy_root), Path(output)
    if not legacy_root.is_dir():
        raise FileNotFoundError(legacy_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter = Counter()
    with open_text(output, "wt") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, extrasaction="ignore",
                                lineterminator="\n")
        writer.writeheader()
        for kind, relative in SOURCES:
            if kind == "biocarta" and not include_biocarta:
                continue
            path = legacy_root / relative
            if not path.exists():
                if skip_missing:
                    counts[f"skipped_missing_{kind}"] += 1
                    continue
                raise FileNotFoundError(path)
            for source in _read(path):
                normalized = (_pathway(source, kind) if kind in {"gobp", "biocarta"}
                              else _molecular(source, kind,
                                              metabolite_feature_name_fallback))
                writer.writerow(normalized)
                counts[kind] += 1
    return counts


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-root", type=Path, required=True,
                        help="Path to the bundled 'MoTrPAC Hackathon' folder")
    parser.add_argument("--output", type=Path, required=True,
                        help="Normalized reference CSV or CSV.gz")
    parser.add_argument("--include-biocarta", action="store_true",
                        help="Include the partial six-set BioCarta table")
    parser.add_argument("--skip-missing", action="store_true",
                        help="Build from available exports and report missing ones")
    parser.add_argument("--metabolite-feature-name-fallback", action="store_true",
                        help="Use a blank-RefMet metabolite feature_id as a lower-confidence RefMetName")
    args = parser.parse_args(argv)
    counts = convert(args.legacy_root, args.output,
                     include_biocarta=args.include_biocarta,
                     skip_missing=args.skip_missing,
                     metabolite_feature_name_fallback=args.metabolite_feature_name_fallback)
    for kind, n in sorted(counts.items()):
        print(f"{kind}: {n}")
    print(f"Total normalized rows: {sum(n for kind, n in counts.items() if not kind.startswith('skipped_'))}")
    print(f"Reference: {args.output.resolve()}")


if __name__ == "__main__":
    main()
