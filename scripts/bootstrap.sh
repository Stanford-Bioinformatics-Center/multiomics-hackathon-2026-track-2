#!/usr/bin/env bash
# =============================================================================
# One-command fresh-clone bootstrap (R8 AC1-AC4)
#
# Usage (from anywhere; the script cd's to the repo root itself):
#
#     bash scripts/bootstrap.sh
#
# What it does, in order, on a fresh clone:
#   1. Create/reuse a Python virtual environment under hackathon/tool/.venv.
#   2. Editable-install the engine   (hackathon/tool      -> pip install -e ".[test]").
#   3. Editable-install the service  (apps/api            -> pip install -e ".[api,test]").
#   4. Fetch the prebuilt mprobe store (mprobe store fetch).
#   5. Run a smoke check that imports the engine + service and confirms the
#      store is present and the FastAPI app builds its catalog.
#
# Contract (R8):
#   - AC1: single command performs the installs + store fetch + smoke check,
#          designed to complete within 600 seconds (see BOOTSTRAP_DEADLINE_S).
#   - AC2: any failed editable install halts immediately, exits non-zero, and
#          names the failed target; no subsequent step runs.
#   - AC3: a failed/unavailable store fetch halts, exits non-zero, and emits a
#          store-fetch error.
#   - AC4: the smoke check exits zero iff every check passes; otherwise it
#          exits non-zero and names the failed check.
#
# Fail-fast semantics: `set -euo pipefail` plus explicit per-step guards so the
# emitted message always names the specific failed target/check.
# =============================================================================

set -euo pipefail

# ----------------------------------------------------------------------------- #
# Configuration
# ----------------------------------------------------------------------------- #
# The bootstrap must finish within 600 seconds (R8 AC1). We record the start
# time and, on success, assert we stayed under the deadline.
BOOTSTRAP_DEADLINE_S=600
START_EPOCH="$(date +%s)"

# Resolve the repo root as the parent of this script's directory
# (repo-root/scripts/bootstrap.sh -> repo root).
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

TOOL_DIR="${REPO_ROOT}/hackathon/tool"
API_DIR="${REPO_ROOT}/apps/api"
VENV_DIR="${TOOL_DIR}/.venv"

# ----------------------------------------------------------------------------- #
# Logging helpers
# ----------------------------------------------------------------------------- #
log()  { printf '[bootstrap] %s\n' "$*"; }
fail() { printf '[bootstrap] ERROR: %s\n' "$*" >&2; exit 1; }

# ----------------------------------------------------------------------------- #
# Step 0 - Python virtual environment
# ----------------------------------------------------------------------------- #
PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
  fail "python interpreter '${PYTHON_BIN}' not found on PATH (set PYTHON_BIN to override)"
fi

if [ ! -d "${VENV_DIR}" ]; then
  log "creating virtual environment at ${VENV_DIR}"
  "${PYTHON_BIN}" -m venv "${VENV_DIR}" \
    || fail "failed to create virtual environment at ${VENV_DIR}"
fi

# Use the venv's interpreter directly so we never depend on shell activation.
VENV_PY="${VENV_DIR}/bin/python"
[ -x "${VENV_PY}" ] || fail "virtual environment python not found at ${VENV_PY}"

log "upgrading pip"
"${VENV_PY}" -m pip install --upgrade pip >/dev/null \
  || fail "failed to upgrade pip inside the virtual environment"

# ----------------------------------------------------------------------------- #
# Step 1 - editable install: engine (hackathon/tool)   [R8 AC2]
# ----------------------------------------------------------------------------- #
# The engine must install first; the service depends on it (editable, not pinned).
log "editable install: engine (hackathon/tool)"
if ! "${VENV_PY}" -m pip install -e "${TOOL_DIR}[test]"; then
  # AC2: name the failed target, run no subsequent step.
  fail "editable install failed for target: hackathon/tool"
fi

# ----------------------------------------------------------------------------- #
# Step 2 - editable install: service (apps/api)        [R8 AC2]
# ----------------------------------------------------------------------------- #
log "editable install: service (apps/api)"
if ! "${VENV_PY}" -m pip install -e "${API_DIR}[api,test]"; then
  # AC2: name the failed target, run no subsequent step.
  fail "editable install failed for target: apps/api"
fi

# ----------------------------------------------------------------------------- #
# Step 3 - fetch the mprobe store                      [R8 AC3]
# ----------------------------------------------------------------------------- #
# `mprobe store fetch` downloads the prebuilt public store bundle. If the fetch
# fails or the store is unavailable, halt with a store-fetch error.
log "fetching the mprobe store (mprobe store fetch)"
if ! "${VENV_DIR}/bin/mprobe" store fetch; then
  # AC3: store-fetch error, non-zero exit, no smoke check.
  fail "mprobe store could not be fetched (store fetch failed or store unavailable)"
fi

# ----------------------------------------------------------------------------- #
# Step 4 - smoke check                                 [R8 AC4]
# ----------------------------------------------------------------------------- #
# Each check is named; the first failure exits non-zero naming that check. The
# smoke check imports the engine + service, confirms the store is present, and
# confirms the FastAPI app can build its catalog (i.e. the wiring is intact).
log "running smoke check"
if ! "${VENV_PY}" - <<'PYEOF'
import sys

def check(name, fn):
    try:
        fn()
    except Exception as exc:  # noqa: BLE001 - surface any failure by name
        print(f"[bootstrap] smoke check FAILED: {name}: {exc}", file=sys.stderr)
        sys.exit(1)
    print(f"[bootstrap] smoke check ok: {name}")

def _import_engine():
    import motrpac_probe  # noqa: F401

def _store_present():
    from motrpac_probe import store
    if not store.CONTRASTS.exists():
        raise RuntimeError(f"engine store not found at {store.CONTRASTS}")

def _import_service():
    from motrpac_probe_service import app as _app  # noqa: F401

def _catalog_builds():
    # Building the catalog exercises the app wiring without starting a server.
    from motrpac_probe_service import catalog
    result = catalog.build_catalog()
    if not result:
        raise RuntimeError("catalog.build_catalog() returned empty result")

check("import motrpac_probe (engine)", _import_engine)
check("engine store present", _store_present)
check("import motrpac_probe_service (service)", _import_service)
check("service catalog builds", _catalog_builds)

print("[bootstrap] all smoke checks passed")
PYEOF
then
  # AC4: the specific failed check is named by the Python block above.
  fail "smoke check failed (see the named check above)"
fi

# ----------------------------------------------------------------------------- #
# Timing guard (R8 AC1)
# ----------------------------------------------------------------------------- #
END_EPOCH="$(date +%s)"
ELAPSED_S=$(( END_EPOCH - START_EPOCH ))
log "completed in ${ELAPSED_S}s (deadline ${BOOTSTRAP_DEADLINE_S}s)"
if [ "${ELAPSED_S}" -gt "${BOOTSTRAP_DEADLINE_S}" ]; then
  fail "bootstrap exceeded the ${BOOTSTRAP_DEADLINE_S}s deadline (took ${ELAPSED_S}s)"
fi

# AC4: exit zero only when every step and every smoke check passed.
log "bootstrap succeeded"
exit 0
