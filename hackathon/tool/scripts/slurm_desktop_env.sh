#!/bin/bash
# create the mprobe env on desktop, same pins as danilogin.
set -euo pipefail
cd /tmp
~/miniconda3/bin/conda create -y -q -n mprobe python=3.11 pip
~/miniconda3/envs/mprobe/bin/pip install -q -r /mnt/nfs/projects/multiomics-hackathon-2026-track-2/hackathon/tool/requirements.txt
~/miniconda3/envs/mprobe/bin/pip install -q -e /mnt/nfs/projects/multiomics-hackathon-2026-track-2/hackathon/tool
~/miniconda3/envs/mprobe/bin/python -c "import motrpac_probe, pandas, pyarrow; print('ok', pandas.__version__)"
