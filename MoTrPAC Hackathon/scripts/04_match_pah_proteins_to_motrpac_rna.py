"""Look up 24-hour healthy MoTrPAC muscle RNA for nine PAH paper proteins.

Run from any directory:
    python scripts/04_match_pah_proteins_to_motrpac_rna.py

The preceding protein-match output supplies the paper UniProt accession to
current MoTrPAC gene-symbol mapping. No statistical model is fitted here.
"""

import argparse
import csv
import gzip
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


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
    protein_matches_path = ROOT / "outputs/PAH protein/pah_motrpac_protein_matches.csv"
    rna_path = ROOT / "data/processed/motrpac_muscle_rna_ee_con_24h.csv.gz"
    output_path = ROOT / "outputs/PAH protein/pah_protein_motrpac_rna_matches.csv"

    with pah_path.open("r", encoding="utf-8", newline="") as handle:
        paper = list(csv.DictReader(handle))
    if len(paper) != 9 or len({r["uniprot_accession"] for r in paper}) != 9:
        raise ValueError("Expected nine distinct PAH protein accessions")

    # Each accession appears at three protein times. All three annotations
    # must agree before one current gene symbol is used for the RNA lookup.
    symbols_by_accession = {}
    rows_by_accession = {}
    with protein_matches_path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            symbols_by_accession.setdefault(row["uniprot_accession"], set()).add(
                row["motrpac_symbol"]
            )
            rows_by_accession[row["uniprot_accession"]] = (
                rows_by_accession.get(row["uniprot_accession"], 0) + 1
            )
    paper_accessions = {r["uniprot_accession"] for r in paper}
    if set(symbols_by_accession) != paper_accessions:
        raise ValueError("Protein matches do not cover exactly the PAH paper accessions")
    if any(count != 3 for count in rows_by_accession.values()):
        raise ValueError("Each protein accession must have three MoTrPAC time-point matches")
    if any(len(symbols) != 1 or "" in symbols for symbols in symbols_by_accession.values()):
        raise ValueError("A protein accession has missing or conflicting gene symbols")
    current_symbol = {key: next(iter(value)) for key, value in symbols_by_accession.items()}

    wanted_symbols = set(current_symbol.values())
    rna_by_symbol = {}
    with gzip.open(rna_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            symbol = row["gene_symbol"]
            if symbol not in wanted_symbols:
                continue
            if (row["tissue"], row["assay"], row["contrast_category"], row["Timepoint"]) != (
                "muscle", "transcript-rna-seq", "EE-CON", "post_24_hr"
            ):
                raise ValueError("A matched RNA row has an unexpected assay or contrast")
            if symbol in rna_by_symbol:
                raise ValueError(f"More than one RNA feature maps to {symbol}")
            rna_by_symbol[symbol] = row
    if set(rna_by_symbol) != wanted_symbols:
        raise ValueError(f"Missing RNA symbols: {sorted(wanted_symbols - set(rna_by_symbol))}")

    output = []
    for paper_row in paper:
        accession = paper_row["uniprot_accession"]
        symbol = current_symbol[accession]
        rna = rna_by_symbol[symbol]
        logfc = float(rna["logFC"])
        q = float(rna["adj_p_value"])
        output.append(
            {
                "paper_symbol": paper_row["paper_symbol"],
                "uniprot_accession": accession,
                "pah_to_control_protein_ratio": paper_row["pah_to_control_ratio"],
                "motrpac_current_symbol": symbol,
                "motrpac_rna_feature_id": rna["feature_id"],
                "motrpac_rna_timepoint": rna["Timepoint"],
                "motrpac_rna_logFC": rna["logFC"],
                "motrpac_rna_bh_adj_p_value": rna["adj_p_value"],
                "motrpac_rna_positive": "yes" if logfc > 0 else "no",
                "motrpac_rna_positive_and_bh_q_lt_0_05": "yes"
                if logfc > 0 and q < 0.05
                else "no",
                "source_doi": paper_row["source_doi"],
                "motrpac_package_version": rna["source_package_version"],
            }
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)

    positive = sum(row["motrpac_rna_positive"] == "yes" for row in output)
    positive_q = sum(
        row["motrpac_rna_positive_and_bh_q_lt_0_05"] == "yes" for row in output
    )
    print(f"Matched {len(output)} PAH paper proteins to 24-hour MoTrPAC muscle RNA")
    print(f"Positive RNA estimates: {positive}/{len(output)}")
    print(f"Positive RNA estimates with BH q < 0.05: {positive_q}/{len(output)}")
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
