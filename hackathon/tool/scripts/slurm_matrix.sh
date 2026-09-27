#!/bin/bash
# robustness matrix on desktop. Stages the store to /scratch.
set -euo pipefail
T=/mnt/nfs/projects/multiomics-hackathon-2026-track-2/hackathon/tool
mkdir -p /scratch/mprobe/store /scratch/mprobe/matrix
rsync -a --exclude raw_metab --exclude library $T/store/ /scratch/mprobe/store/
export MPROBE_STORE=/scratch/mprobe/store MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd $T
~/miniconda3/envs/mprobe/bin/python tests/matrix/run_matrix.py --workers 12 --out /scratch/mprobe/matrix
