"""Run the Metabolomics module from any folder.

    python run_pipeline.py                       all steps
    python run_pipeline.py --skip-r              reuse the bundled MoTrPAC export
    python run_pipeline.py --r-library PATH      MoTrPAC package in its own library
    python run_pipeline.py --rscript PATH        Rscript not on PATH

Stops at the first error. Never downloads data or installs packages.
"""

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STEPS = [
    "00_verify_inputs.py",
    "01_prepare_st000763.py",
    "02_compare_exercise_response.py",
    "03_compare_resting_levels.py",
    "04_export_motrpac_metabolomics.R",
    "05_align_exercise_with_motrpac.py",
    "06_build_context_table.py",
    "07_plot_resting_setting_effect.py",
]


def find_rscript(supplied):
    if supplied:
        return supplied
    found = shutil.which("Rscript")
    if found:
        return found
    if sys.platform == "win32":
        candidates = sorted(Path("C:/Program Files/R").glob("R-*/bin/Rscript.exe"), reverse=True)
        if candidates:
            return str(candidates[0])
    raise FileNotFoundError("Rscript not found. Install R, pass --rscript PATH, or use --skip-r.")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--rscript")
    parser.add_argument("--r-library", help="Folder holding MotrpacHumanPreSuspensionAnalysis 2.0.8")
    parser.add_argument("--skip-r", action="store_true",
                        help="Skip step 04 and use data/processed/motrpac_metabolite_contrasts.csv.gz")
    args = parser.parse_args()

    env = os.environ.copy()
    if args.r_library:
        env["R_LIBS"] = str(Path(args.r_library).resolve())
    for step in STEPS:
        path = ROOT / "scripts" / step
        if step.endswith(".R"):
            if args.skip_r:
                print(f"Skipping {step}; using the bundled MoTrPAC export.", flush=True)
                continue
            command = [find_rscript(args.rscript), str(path)]
        else:
            command = [sys.executable, str(path)]
        print(f"\n=== {step} ===", flush=True)
        subprocess.run(command, cwd=ROOT / "scripts", env=env, check=True)
    print("\nPipeline complete. Results are in output/.")


if __name__ == "__main__":
    main()
