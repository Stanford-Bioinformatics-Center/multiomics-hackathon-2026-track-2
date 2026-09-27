# Run the app locally

The app has **two processes**: a Python **API** (FastAPI) and a **web** dev server (Vite/React).
You open the **web** URL in your browser — not the API. They run on different ports.

> Common mistake: `127.0.0.1:8000` is NOT this app. If you see an app called "aegis" there, that is a
> different project occupying port 8000. This app's web UI runs on **8443** (or 5173) and the API on a
> port you pick below — choose a free one to avoid clashing with other servers.

## Prerequisites (one time)

For a fresh clone, the one-command bootstrap documented in [BOOTSTRAP.md](BOOTSTRAP.md)
runs all of the steps below (editable installs + `mprobe store fetch` + a smoke check):

```
bash scripts/bootstrap.sh
```

Or do the same steps manually:

```
cd "multiomics-hackathon-2026-track-2/hackathon/tool"
python -m venv .venv
source .venv/bin/activate
pip install -e ".[test]"
pip install -e "../../apps/api[api]"
mprobe store fetch
```
Node 22 for the web app (pinned via `apps/web/.mise.toml`).

## One-off: pathway table for the explorer

The explorer's *Pathways* input reads MoTrPAC's precomputed CAMERA results. Build the table once (about a minute, needs R) from the repository root:

```
Rscript apps/api/scripts_build_pathways.R
```

It writes `apps/api/data/motrpac_camera_pathways.csv.gz` (not committed; about 32 MB). Without it, gene/protein and metabolite lists still work and the Pathways tab explains how to build it.

## Terminal A — the API

```
cd "multiomics-hackathon-2026-track-2/apps/api"
source ../../hackathon/tool/.venv/bin/activate
uvicorn motrpac_probe_service.app:app --host 127.0.0.1 --port 8765
```
Leave it running. Sanity check: `curl http://127.0.0.1:8765/api/health` -> `{"status":"ok",...}`.
Port 8765 is an arbitrary free choice; use anything not already taken.

## Terminal B — the web app (this is the URL you open)

```
cd "multiomics-hackathon-2026-track-2/apps/web"
pnpm install --frozen-lockfile
VITE_API_BASE_URL=http://127.0.0.1:8765/api PORT=8443 pnpm run dev
```
Vite prints a URL like `http://localhost:8443/`. Open that in your browser. `VITE_API_BASE_URL` tells
the frontend where the API is; without it the live views cannot fetch data.

If 8443 (or 5173) is busy, pick another: `PORT=5199 ... pnpm run dev`. If an old Vite is stuck, find and
stop it: `lsof -nP -iTCP:5173 -sTCP:LISTEN` then `kill <PID>`.

## What you'll see (8 views in the top nav)

| # | View | Data source |
|---|---|---|
| 01-04, 08 | Architecture / Technical flow / Workflow / Results dashboard / Judge slide | static design (mock) |
| 05 | **Explorer** | paste or upload any gene/protein, metabolite or pathway list (or pick an example); returns every matching MoTrPAC comparison by omic layer (human acute and rat training, all tissues, contrasts, times, sexes) with raw P / BH / Bonferroni |
| 06 | **Metabolomics (ST000763)** | standalone metabolomics module (read-only case study) |
| 07 | **Generalized query** | `query_core` engine — pick a signature, run vs MoTrPAC; PAH examples reproduce the legacy result |

Views 05-07 fetch real data from the API. If the API is not running, they show a
"backend not reachable" state with a retry button.

## CLI (no server) — the generalized pipeline directly

Run these from inside `MoTrPAC Hackathon/generalized/` (that is where `query_core` is importable;
from the repo root you get `ModuleNotFoundError: No module named 'query_core'`).

```
cd "multiomics-hackathon-2026-track-2/MoTrPAC Hackathon/generalized"
source ../../hackathon/tool/.venv/bin/activate

python -m query_core.cli --disease examples/signatures/pah_blood_rna.csv.gz \
  --reference query_core/motrpac_reference.csv.gz --out out_blood \
  --tissue blood --reference-contrast-category EE-CON

python -m query_core.cli --disease examples/signatures/pah_muscle_protein.csv.gz \
  --reference query_core/motrpac_reference_full.csv.gz --out out_muscle \
  --tissue muscle --reference-contrast-category EE-CON

python run_config.py --config examples/pah_jobs.example.json --out demo_all
```
For the blood run, the `rna` rows of `out_blood/rank_correlation.csv` are the six Spearman values that
reproduce the legacy hardcoded pipeline.

## Ports summary
- API (uvicorn): you choose — examples use 8765. JSON only, not a browser page.
- Web (Vite): 8443 by default here (or 5173) — this is the browser URL.
- 8000: unrelated ("aegis") — do not use it for this app.


## CORS (cross-origin)
The web app (Vite) and API run on different origins, so the browser sends a preflight OPTIONS before
each request. The API enables CORS for localhost/127.0.0.1 on any port by default, so local dev works
out of the box. For a non-localhost origin, set MPROBE_CORS_ORIGINS before starting the API server,
e.g. MPROBE_CORS_ORIGINS="https://my.host".
