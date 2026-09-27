#!/bin/bash
# rebuild example runs, compare, library, discord and site on desktop.
set -euo pipefail
T=/mnt/nfs/projects/multiomics-hackathon-2026-track-2/hackathon/tool
mkdir -p /scratch/mprobe/store
rsync -a --exclude raw_metab --exclude library $T/store/ /scratch/mprobe/store/
export MPROBE_STORE=/scratch/mprobe/store MPLBACKEND=Agg
cd $T
~/miniconda3/envs/mprobe/bin/python scripts/build_all.py --skip-store --only examples
~/miniconda3/envs/mprobe/bin/python -m motrpac_probe.cli site build
