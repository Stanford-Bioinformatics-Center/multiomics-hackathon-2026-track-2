"""Command-line entry point for typed disease-to-MoTrPAC queries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .engine import run_query


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--disease", type=Path, required=True,
                        help="Normalized disease signature CSV or CSV.gz")
    parser.add_argument("--reference", type=Path, required=True,
                        help="Normalized MoTrPAC long-format CSV or CSV.gz")
    parser.add_argument("--out", type=Path, required=True,
                        help="Output directory for this query")
    parser.add_argument("--species", choices=("human", "rat"))
    parser.add_argument("--reference-species", choices=("human", "rat"),
                        help="Optional independent reference species filter")
    parser.add_argument("--tissue")
    parser.add_argument("--disease-contrast")
    parser.add_argument("--reference-contrast-category")
    parser.add_argument("--reference-contrast")
    parser.add_argument("--timepoint", action="append", default=[],
                        help="May be repeated; omit to use all available times")
    parser.add_argument("--rank-field", choices=("auto", "log2_fc", "statistic"),
                        default="auto", help="auto uses signed statistics only when all unique pairs have them; otherwise signed log2_fc")
    parser.add_argument("--q-threshold", type=float, default=0.05)
    parser.add_argument("--min-abs-log2-fc", type=float, default=0.0)
    parser.add_argument("--ortholog-map", type=Path,
                        help="Explicit one-to-one source/target gene-symbol CSV")
    args = parser.parse_args(argv)
    summary = run_query(
        args.disease, args.reference, args.out, species=args.species,
        reference_species=args.reference_species,
        tissue=args.tissue, disease_contrast=args.disease_contrast,
        reference_contrast_category=args.reference_contrast_category,
        reference_contrast=args.reference_contrast, timepoints=tuple(args.timepoint),
        rank_field=args.rank_field, q_threshold=args.q_threshold,
        min_abs_log2_fc=args.min_abs_log2_fc, ortholog_map=args.ortholog_map,
    )
    print(json.dumps(summary["counts"], indent=2))
    print(f"Results: {args.out.resolve()}")


if __name__ == "__main__":
    main()
