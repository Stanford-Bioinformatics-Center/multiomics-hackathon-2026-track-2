"""Match the PAH paper's nine lower proteins to MoTrPAC muscle protein results.

Run from any directory:
    python scripts/03_compare_pah_motrpac_protein.py

Inputs are a manually curated paper table and a precomputed MoTrPAC export.
This script does not fit either study's statistical model.
"""

import argparse
import csv
import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TIMES = ("post_15_30_45_min", "post_3.5_4_hr", "post_24_hr")


def under_root(path: Path) -> Path:
    return path if path.is_absolute() else ROOT / path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pah",
        type=Path,
        default=Path("data/from paper/pah_lower_proteins_malenfant2015.csv"),
        help="PAH paper Table 2 transcription, relative to the project root",
    )
    args = parser.parse_args()

    pah_path = under_root(args.pah)
    protein_path = ROOT / "data/processed/motrpac_muscle_protein_ee_con.csv.gz"
    output_path = ROOT / "outputs/PAH protein/pah_motrpac_protein_matches.csv"

    with pah_path.open("r", encoding="utf-8", newline="") as handle:
        paper = list(csv.DictReader(handle))
    if len(paper) != 9:
        raise ValueError(f"Expected nine PAH paper proteins; found {len(paper)}")
    accessions = [row["uniprot_accession"] for row in paper]
    if len(set(accessions)) != len(accessions):
        raise ValueError("The PAH paper table has duplicate UniProt accessions")

    # Only retain the accession/time pairs needed for this comparison.
    protein = {}
    with gzip.open(protein_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            if row["uniprot"] not in accessions or row["Timepoint"] not in TIMES:
                continue
            if (row["tissue"], row["assay"], row["contrast_category"]) != (
                "muscle", "prot-pr", "EE-CON"
            ):
                raise ValueError("A matched MoTrPAC row has an unexpected assay or contrast")
            key = (row["uniprot"], row["Timepoint"])
            if key in protein:
                raise ValueError(f"Duplicate MoTrPAC accession/time match: {key}")
            protein[key] = row

    output = []
    for paper_row in paper:
        accession = paper_row["uniprot_accession"]
        for time in TIMES:
            key = (accession, time)
            if key not in protein:
                raise ValueError(f"Missing MoTrPAC accession/time match: {key}")
            mo = protein[key]
            q = float(mo["adj_p_value"])
            output.append(
                {
                    "paper_symbol": paper_row["paper_symbol"],
                    "uniprot_accession": accession,
                    "pah_to_control_ratio": paper_row["pah_to_control_ratio"],
                    "paper_reported_p_value_text": paper_row["reported_p_value_text"],
                    "motrpac_symbol": mo["gene_symbol"],
                    "motrpac_timepoint": time,
                    "motrpac_logFC": mo["logFC"],
                    "motrpac_bh_adj_p_value": mo["adj_p_value"],
                    "motrpac_bh_q_lt_0_05": "yes" if q < 0.05 else "no",
                    "source_doi": paper_row["source_doi"],
                    "motrpac_package_version": mo["source_package_version"],
                }
            )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    significant = sum(row["motrpac_bh_q_lt_0_05"] == "yes" for row in output)
    print(f"Matched {len(paper)} PAH proteins at {len(TIMES)} times: {len(output)} rows")
    print(f"MoTrPAC protein results with BH q < 0.05: {significant}/{len(output)}")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
