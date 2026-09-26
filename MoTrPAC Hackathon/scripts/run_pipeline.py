"""Run the modular PAH muscle and blood analyses in dependency order.

Examples:
    python scripts/run_pipeline.py
    python scripts/run_pipeline.py --module blood --rscript /path/to/Rscript
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
MUSCLE = (
    "00_audit_pah_muscle_rna_timepoints.R",
    "01_export_motrpac_muscle_rna.R",
    "02_export_motrpac_muscle_protein.R",
    "03_compare_pah_motrpac_protein.py",
    "04_match_pah_proteins_to_motrpac_rna.py",
    "04.5_audit_pah_oxphos_overlap.R",
    "05_export_motrpac_oxphos_pathway.R",
)
BLOOD = (
    "06_prepare_gse33463_blood.R",
    "07_fit_gse33463_blood_limma.R",
    "07.5_export_motrpac_blood_gene_ranks.R",
    "08_compare_pah_motrpac_blood_ranks.R",
    "09_test_pah_gobp_pathways.R",
    "10_export_motrpac_blood_gobp.R",
    "11_match_pah_motrpac_blood_gobp.R",
    "12_check_pah_motrpac_blood_rank_sensitivity.R",
)


def find_rscript(explicit: str | None) -> str:
    if explicit:
        path = Path(explicit).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Rscript was not found at {path}")
        return str(path)
    on_path = shutil.which("Rscript")
    if on_path:
        return on_path
    if sys.platform == "win32":
        candidates = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"))
        if candidates:
            return str(candidates[-1])
    raise FileNotFoundError("Rscript is not on PATH; pass --rscript /path/to/Rscript")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", choices=("all", "muscle", "blood"),
                        default="all")
    parser.add_argument("--rscript", help="Path to Rscript, if it is not on PATH")
    args = parser.parse_args()
    rscript = find_rscript(args.rscript)
    steps = (MUSCLE if args.module == "muscle" else
             BLOOD if args.module == "blood" else MUSCLE + BLOOD)
    for index, name in enumerate(steps, 1):
        script = SCRIPTS / name
        executable = sys.executable if script.suffix == ".py" else rscript
        print(f"[{index}/{len(steps)}] {name}", flush=True)
        subprocess.run((executable, str(script)), cwd=ROOT, check=True)
    print("Pipeline completed successfully.", flush=True)


if __name__ == "__main__":
    main()
