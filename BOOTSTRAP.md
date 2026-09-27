# Fresh-clone bootstrap

This document describes the one-command bootstrap for a fresh clone (R8) and the
`parents[3]` path-coupling constraint that ties the API adapters to the repository
layout.

## One-command invocation

From anywhere inside the clone (the script resolves the repository root itself):

```
bash scripts/bootstrap.sh
```

On a fresh clone this single command:

1. Creates (or reuses) a Python virtual environment at `hackathon/tool/.venv`.
2. Editable-installs the engine — `pip install -e "hackathon/tool[test]"`.
3. Editable-installs the service — `pip install -e "apps/api[api,test]"`.
4. Fetches the prebuilt public store — `mprobe store fetch`.
5. Runs a smoke check that imports the engine and the service, confirms the store
   is present, and confirms the FastAPI app builds its catalog.

The bootstrap is designed to finish within **600 seconds**; it records its own
elapsed time and fails if it exceeds that deadline.

To run the local app after bootstrapping, follow [RUN_LOCAL.md](RUN_LOCAL.md) — it
starts the two processes (FastAPI API + Vite web) and documents the port wiring.

### Exit-code contract

The script uses `set -euo pipefail` fail-fast semantics and emits a targeted
error message for every failure mode:

| Situation | Behavior | Exit code |
| --- | --- | --- |
| Everything succeeds | prints `bootstrap succeeded` | `0` |
| Editable install of `hackathon/tool` fails | halts, names `hackathon/tool`, runs no later step | non-zero |
| Editable install of `apps/api` fails | halts, names `apps/api`, runs no later step | non-zero |
| `mprobe store fetch` fails or the store is unavailable | halts, emits a store-fetch error | non-zero |
| A smoke check fails | halts, names the specific failed check | non-zero |
| Bootstrap exceeds 600 s | halts, reports the elapsed time vs the deadline | non-zero |

The engine install runs first because the service (`apps/api`) depends on the
engine as an editable package rather than a pinned dependency; a failure in the
engine install therefore stops before the service install is attempted.

## Package manager

The project uses a single declared package manager. The Python packages install
editable via `pip` as shown above; the web app uses **pnpm** (honoring the
committed `apps/web/pnpm-lock.yaml`), and CI uses pnpm as well. No `npm`
invocation remains in the setup instructions or CI.

## The `parents[3]` path-coupling constraint

The API adapters locate the committed demo/output files by climbing from their own
file location to the repository root:

```python
_REPO = Path(__file__).resolve().parents[3]
```

`parents[3]` means "go up **four** directory levels from the adapter file"
(`parents[0]` is the file's own directory, `parents[1]` its parent, and so on).
This only resolves to the repository root when each adapter sits **exactly three
directories below the repo root**:

```
<repo-root>/            <- parents[3]  (the value _REPO must equal)
  apps/                 <- parents[2]
    api/                <- parents[1]
      motrpac_probe_service/   <- parents[0]  (the adapter's own directory)
        discordance_casestudy.py   <- Path(__file__)
        metab_casestudy.py
        convergence.py
        generalized_query.py
```

So the coupling is: **an adapter file must live at
`apps/api/motrpac_probe_service/<name>.py`** — three directories (`apps` → `api` →
`motrpac_probe_service`) beneath the repository root. From `_REPO`, the adapters
then reference repo-root paths such as `_REPO / "MoTrPAC Hackathon" / ...`.

### How a reader verifies a given file location

Count the directory segments between the repository root and the file:

- Let `depth` = the number of directories from the repo root down to (and
  including) the directory that contains the file, **not** counting the file
  itself.
- The constraint is satisfied when `depth == 3` (the file is
  `<root>/a/b/c/file.py`), because `Path(file).resolve().parents[3]` is then the
  repo root.
- If the file is moved shallower (e.g. to `apps/api/<name>.py`, depth 2) or
  deeper (depth 4), `parents[3]` no longer points at the repo root and the
  adapter fails to find the committed data files.

Equivalently, to check the current adapter without counting by hand:

```
cd "<repo-root>"
python3 -c "import pathlib, motrpac_probe_service.discordance_casestudy as m, sys; \
print('OK' if pathlib.Path(m.__file__).resolve().parents[3] == pathlib.Path('.').resolve() else 'MISMATCH')"
```

`OK` means the adapter's `parents[3]` resolves to the repository root and the
path-coupling constraint holds. The `apps/api` test suite pins this with
`test_parents3_path_coupling` in `apps/api/tests/test_discordance.py`.
