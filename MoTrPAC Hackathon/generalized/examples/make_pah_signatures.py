"""Convert the existing PAH examples to the shared processed-signature schema.

This is an example adapter, not a general disease differential-analysis fitter.
It reads published/reanalysed results already present in the legacy repository.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
from pathlib import Path


FIELDS = [
    "study_id", "species", "tissue", "biospecimen", "layer", "feature_id",
    "id_namespace", "log2_fc", "statistic", "p_value", "q_value",
    "contrast", "gene_symbol", "uniprot", "refmet_id", "refmet_name",
    "source_feature_id", "source_gene_symbol", "alias_mapping_basis",
    "source", "statistic_type", "source_q_value",
    "source_q_scope", "pathway_collection", "evidence_note",
]


def read_csv(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle)


def write_csv(path: Path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        count = 0
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in FIELDS})
            count += 1
    return count


def blood_rows(root: Path):
    path = root / "outputs/PAH blood/pah_blood_gene_ranked.csv.gz"
    for row in read_csv(path):
        symbol = row["gene_symbol"].strip()
        if not symbol:
            continue
        yield {
            "study_id": "GSE33463", "species": "human", "tissue": "blood",
            "biospecimen": "PBMC",
            "layer": "rna", "feature_id": symbol, "id_namespace": "HGNC",
            "log2_fc": row["logFC"], "statistic": row["t"],
            "p_value": row["p_value"], "q_value": "",
            "contrast": "IPAH minus healthy (resting PBMC)",
            "gene_symbol": symbol, "source_feature_id": row["probe_id"],
            "source": "GSE33463 reanalysis, legacy script 07",
            "source_q_value": row["adj_p_value"],
            "source_q_scope": "probe-level BH carried by outcome-selected gene probe",
            "evidence_note": (
                "Gene represented by outcome-selected probe; formal gene-level "
                "q_value unavailable. PBMC versus MoTrPAC whole blood."
            ),
        }


def load_aliases(path: Path):
    aliases = {}
    for row in read_csv(path):
        old = row["source_symbol"].strip()
        if not old or old in aliases:
            raise ValueError(f"Duplicate or blank alias in {path}: {old!r}")
        aliases[old] = (row["current_symbol"].strip(), row["mapping_basis"].strip())
    return aliases


def muscle_rows(root: Path, aliases):
    path = root / "data/from paper/pah_lower_proteins_malenfant2015.csv"
    for row in read_csv(path):
        original_symbol = row["paper_symbol"].strip()
        current_symbol, alias_basis = aliases.get(original_symbol, (original_symbol, ""))
        ratio = float(row["pah_to_control_ratio"])
        if not math.isfinite(ratio) or ratio <= 0:
            raise ValueError(f"Invalid PAH/control ratio for {row['paper_symbol']}")
        reported_p = row["reported_p_value_text"].strip()
        p_value = reported_p if reported_p and not reported_p.startswith("<") else ""
        yield {
            "study_id": "Malenfant2015", "species": "human", "tissue": "muscle",
            "biospecimen": "skeletal muscle",
            "layer": "protein", "feature_id": row["uniprot_accession"],
            "id_namespace": "UniProt", "log2_fc": math.log2(ratio),
            "p_value": p_value, "contrast": "resting PAH minus control",
            "gene_symbol": current_symbol,
            "uniprot": row["uniprot_accession"],
            "source_feature_id": row["uniprot_accession"],
            "source_gene_symbol": original_symbol,
            "alias_mapping_basis": alias_basis,
            "source": row["source_doi"],
            "evidence_note": (
                "Exploratory iTRAQ subset (four PAH, four controls); "
                f"reported P text={reported_p}; no multiple-testing q reported."
            ),
        }


def metabolite_rows(root: Path, comparison: str):
    processed = root / "Metabolomics"
    features = list(read_csv(processed / "data/processed/st000763_features.csv"))
    # Select duplicate assay features without inspecting their disease effect.
    features.sort(key=lambda row: (
        row["refmet_name"],
        -int(row["n_measured"] or 0),
        -float(row["median_intensity"] or 0),
        row["feature_key"],
    ))
    selected = {}
    for row in features:
        name = row["refmet_name"].strip()
        if name and name not in selected:
            selected[name] = row["feature_key"]
    tests = {
        row["feature_key"]: row
        for row in read_csv(processed / "output/03_resting_feature_tests.csv")
    }
    config = {
        "healthy": (
            "PAH_vs_Healthy", "SSc-PAH minus healthy at rest",
            "Confounded by systemic sclerosis and catheterization setting."
        ),
        "normal_pressure": (
            "PAH_vs_Normal_Pressures",
            "SSc-PAH minus SSc with normal pressures at rest",
            "Same catheterization setting; smaller comparator group."
        ),
    }
    prefix, contrast, caveat = config[comparison]
    for name, feature_key in selected.items():
        row = tests.get(feature_key)
        if not row or row[f"eligible_{prefix}"].lower() != "true":
            continue
        yield {
            "study_id": "ST000763", "species": "human", "tissue": "blood",
            "biospecimen": "plasma",
            "layer": "metabolite", "feature_id": name,
            "id_namespace": "RefMetName", "log2_fc": row[f"log2_difference_{prefix}"],
            "p_value": row[f"p_welch_{prefix}"],
            "q_value": row[f"q_welch_{prefix}"],
            "contrast": contrast, "refmet_name": name,
            "source_feature_id": feature_key,
            "source": "Metabolomics Workbench ST000763 reanalysis, legacy script 03",
            "source_q_scope": "feature-level BH before result-independent name collapse",
            "evidence_note": (
                caveat + " RefMet name match is lower confidence than an ID; "
                "q_value is feature-level BH before result-independent name collapse."
            ),
        }


def biocarta_rows(root: Path):
    path = root / "Transcriptomics/data/paper/table2_selected_pathways.csv"
    for row in read_csv(path):
        if row["annotation"] != "BioCarta" or not row["set"].strip():
            continue
        yield {
            "study_id": "Cheadle2012", "species": "human", "tissue": "blood",
            "biospecimen": "PBMC", "layer": "pathway",
            "feature_id": row["set"], "id_namespace": "BIOCARTA",
            "statistic": row["ipah_page_z"],
            "contrast": "IPAH minus healthy (resting PBMC)",
            "source_feature_id": row["paper_label"],
            "source": "10.1371/journal.pone.0034951 Table 2",
            "statistic_type": "PAGE_Z",
            "pathway_collection": "C2/BIOCARTA",
            "evidence_note": (
                "Published selected pathway result, not a full tested pathway "
                "universe. No pathway log2 fold change, P value, or q value "
                "is supplied. PAGE and MoTrPAC CAMERA Z are different tests."
            ),
        }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--legacy-root", required=True, type=Path,
                        help="Path to the existing MoTrPAC Hackathon folder")
    parser.add_argument("--alias-map", type=Path,
                        default=Path(__file__).with_name("pah_gene_aliases.csv"),
                        help="Reviewed old-to-current symbol map for this example")
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    root = args.legacy_root.resolve()
    out = args.out.resolve()
    aliases = load_aliases(args.alias_map)
    products = {
        "pah_blood_rna.csv.gz": blood_rows(root),
        "pah_muscle_protein.csv.gz": muscle_rows(root, aliases),
        "st000763_rest_pah_vs_healthy.csv.gz": metabolite_rows(root, "healthy"),
        "st000763_rest_pah_vs_normal_pressure.csv.gz":
            metabolite_rows(root, "normal_pressure"),
        "cheadle2012_biocarta_ipah.csv.gz": biocarta_rows(root),
    }
    for name, rows in products.items():
        path = out / name
        print(f"{path}: {write_csv(path, rows)} rows")


if __name__ == "__main__":
    main()
