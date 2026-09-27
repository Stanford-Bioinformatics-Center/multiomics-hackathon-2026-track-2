#!/bin/bash
# Slurm job: build tool/store/*.parquet from the hackathon CSVs (runs on danilogin where the mprobe env lives).
#SBATCH -J mprobe-store
set -euo pipefail
cd /mnt/nfs/projects/multiomics-hackathon-2026-track-2/hackathon/tool
~/miniconda3/envs/mprobe/bin/pip install -q -e . >/dev/null
~/miniconda3/envs/mprobe/bin/mprobe store build
