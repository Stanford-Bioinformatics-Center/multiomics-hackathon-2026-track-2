# motrpac_probe. `make all` rebuilds store -> examples -> compare -> library -> discord -> site.
# On danilogin heavy steps go through Slurm: `make slurm-all`.
PY ?= python
.PHONY: all store examples site test slurm-all clean-out

all:
	$(PY) scripts/build_all.py

store:
	$(PY) -m motrpac_probe.cli store build

examples:
	$(PY) scripts/build_all.py --only examples

site:
	$(PY) -m motrpac_probe.cli site build

test:
	$(PY) -m pytest -q

slurm-all:   # desktop node: examples + compare + library + discord + site (~2 min); matrix: scripts/slurm_matrix.sh
	cd /tmp && sbatch -p cpu -w desktop -c 8 --mem=16G --time=00:40:00 -J mprobe-rebuild \
	  -o /mnt/nfs/shared/%x-%j.out $(CURDIR)/scripts/slurm_rebuild.sh

clean-out:
	rm -rf out/*
