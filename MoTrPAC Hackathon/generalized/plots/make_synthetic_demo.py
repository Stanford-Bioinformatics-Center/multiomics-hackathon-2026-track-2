"""Create a small fictional query result and render it for visual QA."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from query_core.engine import run_query

from plot_results import render


def write_csv(path: Path, rows: list[dict]) -> None:
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    root = Path(__file__).resolve().parent / "demo_synthetic"
    root.mkdir(parents=True, exist_ok=True)
    genes = [("GENEA", 1.2, 4.8), ("GENEB", -0.8, -3.5),
             ("GENEC", 0.6, 2.5), ("GENED", -0.3, -1.3),
             ("GENEE", 0.2, 0.9)]
    disease = [{
        "study_id": "SyntheticDisease", "species": "human", "tissue": "muscle",
        "layer": "rna", "feature_id": gene, "id_namespace": "HGNC",
        "log2_fc": fc, "statistic": stat, "p_value": 0.01, "q_value": 0.04,
        "contrast": "case minus control",
    } for gene, fc, stat in genes]
    reference = []
    for layer in ("rna", "protein"):
        for time, multiplier in (("post_24_hr", 1), ("during_20_min", -0.6)):
            for gene, fc, _ in genes:
                ref_fc = round(fc * multiplier * (0.8 if layer == "protein" else 1), 3)
                reference.append({
                    "study_id": "MoTrPAC", "species": "human", "tissue": "muscle",
                    "layer": layer, "feature_id": gene, "id_namespace": "HGNC",
                    "log2_fc": ref_fc, "statistic": ref_fc, "p_value": 0.01,
                    "q_value": "" if gene == "GENEC" else 0.03,
                    "contrast": f"exercise minus control {time}",
                    "timepoint": time, "contrast_category": "EE-CON",
                })
    # The extra protein mapping shows how an ambiguous cell is excluded.
    reference.append({**reference[-2], "feature_id": "GENED_alt",
                      "gene_symbol": "GENED"})
    write_csv(root / "disease.csv", disease)
    write_csv(root / "reference.csv", reference)
    run_query(root / "disease.csv", root / "reference.csv", root / "query")
    manifest = render(root / "query", root / "plots", top=5)
    print(root / "plots" / manifest["groups"][0]["directory"])


if __name__ == "__main__":
    main()
