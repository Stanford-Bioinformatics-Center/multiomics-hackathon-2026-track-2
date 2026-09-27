"""Verify bundled inputs and run the rat protein comparison from any directory.

Requires Python 3.10+ and Rscript. All analysis scripts use base R.
"""

import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPTS = (
    "01_map_paper_to_rat.R",
    "02_prepare_rat_samples.R",
    "03_estimate_training_course.R",
    "04_compare_pooling_choices.R",
    "05_plot_training_course.R",
)


def verify_inputs() -> None:
    paper = ROOT / "data" / "paper" / "pah_lower_proteins_malenfant2015.csv"
    with paper.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if len(rows) != 9:
        raise ValueError(f"Expected nine downregulated paper proteins; found {len(rows)}")
    if len({row["paper_symbol"] for row in rows}) != 9:
        raise ValueError("Paper symbols are not unique")
    if len({row["uniprot_accession"] for row in rows}) != 9:
        raise ValueError("Paper UniProt accessions are not unique")
    for row in rows:
        ratio = float(row["pah_to_control_ratio"])
        if not 0 < ratio < 1:
            raise ValueError(f"Unexpected PAH/control ratio for {row['paper_symbol']}")
        if row["source_doi"] != "10.1007/s00109-014-1244-0":
            raise ValueError(f"Unexpected paper DOI for {row['paper_symbol']}")

    rat_dir = ROOT / "data" / "rat"
    manifest = json.loads((rat_dir / "source_manifest.json").read_text(encoding="utf-8"))
    for name, expected in manifest["files"].items():
        path = rat_dir / name
        digest = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b""):
                digest.update(block)
        if digest.hexdigest() != expected:
            raise ValueError(f"SHA-256 mismatch for {name}")
        print(f"Verified {name}", flush=True)
    print("Verified the nine paper rows and all rat data files.", flush=True)


def find_rscript(supplied: str | None) -> str:
    if supplied:
        path = Path(supplied).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Rscript not found: {path}")
        return str(path)
    on_path = shutil.which("Rscript")
    if on_path:
        return on_path
    if sys.platform == "win32":
        base = Path("C:/Program Files/R")
        candidates = sorted(base.glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    raise FileNotFoundError("Rscript not found. Install R or pass --rscript PATH.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rscript", help="Path to Rscript when it is not on PATH")
    parser.add_argument("--verify-only", action="store_true", help="Check inputs without running R")
    args = parser.parse_args()
    verify_inputs()
    if args.verify_only:
        return
    rscript = find_rscript(args.rscript)
    for name in SCRIPTS:
        print(f"Running {name}", flush=True)
        subprocess.run([rscript, str(ROOT / "scripts" / name)], cwd=ROOT, check=True)
    print("Pipeline complete. Results are in output/.", flush=True)


if __name__ == "__main__":
    main()
