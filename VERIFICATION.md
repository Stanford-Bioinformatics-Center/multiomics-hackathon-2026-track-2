# VERIFICATION — clean-clone / CI reproduction log

Purpose: prove the integrated project reproduces from a clean state. Records exact
commands, versions, checksums, and passed/skipped tests. Append per gate.

## Environment baseline
- OS: macOS (darwin, arm64)
- git: 2.47.0
- git-lfs: 3.7.0 (present but NOT required — see below)

## Repo state at start of integration
- Integration branch base: `origin/t3code/build-motrpac-probe-tool` @ `9962a43`
- Merged: `origin/main` @ `57e01e7` (merge commit `09ae0d7`)
- Merge base of tool branch and main: `993f746`
- React app source repo: `github.com/StanchPillow55/ExerciseSignatureExplorerHandoff` @ `d9fd564`

## LFS finding (Gate 1)
- The React repo's `.gitattributes` declares LFS rules, but `git lfs ls-files` is EMPTY.
- The only tracked binary, `src/imports/Screenshot_2026-09-26_at_6.48.13_PM.png`, is a
  real 152,756-byte PNG git blob (`git cat-file -p HEAD:...` begins with PNG magic bytes,
  not a `version https://git-lfs...` pointer).
- Conclusion: there is NO LFS object to materialize; a subtree import carries the real
  blob. The parent repo has no LFS dependency.

---

## Gate checkpoints
(Appended as each gate completes.)

### Gate 1 checkpoint — Repository & integration (PASS)
Files changed / added:
- `DECISIONS.md`, `VERIFICATION.md`, `REQUIREMENTS_TRACEABILITY.md` (new; traceability trail)
- `apps/web/**` (React app imported via `git subtree add --prefix=apps/web`, no `--squash`)
- `apps/web/.gitattributes` (LFS rules replaced with plain `binary` attributes)

Why each change was necessary:
- Branch off tool branch + merge main → single integration line with engine + analysis
  modules (ADR-0002).
- Subtree import (not submodule, not rm-rf) → self-contained repo, preserved history
  (ADR-0003).
- `.gitattributes` neutralized → fresh clone works without Git LFS.

Commands run (key):
- `git switch -c integration/exercise-signature-explorer origin/t3code/build-motrpac-probe-tool`
- `git merge --no-ff origin/main`  → merge commit `09ae0d7`
- `git subtree add --prefix=apps/web ese-src main`  (source `d9fd564`)  → commit `7863879`
- `GIT_LFS_SKIP_SMUDGE=1 git clone --branch integration/... <repo>.git /tmp/ese-freshclone`

Results:
- Integration branch HEAD after Gate 1: `728f36f`
- `hackathon/tool/` (mprobe) and `MoTrPAC Hackathon/Rat Comparison PAH Proteins/` coexist.
- React app at `apps/web`, tracked as real files (not a gitlink).
- PNG `apps/web/src/imports/Screenshot_2026-09-26_at_6.48.13_PM.png`: real PNG blob
  (152,756 bytes), NOT an LFS pointer, verified in a fresh clone with LFS smudge skipped.
- In-repo build: `npm run build` OK (JS 269.10 kB / gzip 81.42 kB), `npm test` 16/16 pass.
- Fresh clone (no original repo, LFS skipped): install OK, build OK, tests 16/16 pass.

Skips / assumptions / limitations:
- The React app ships both `pnpm-lock.yaml` (original) and now installs cleanly with npm;
  the stray `package-lock.json` from an npm install was removed to avoid dual lockfiles.
  Toolchain (`.mise.toml`) prefers pnpm 10.34.3; CI should standardize on one. (Open item.)
- Source repo URL: `github.com/StanchPillow55/ExerciseSignatureExplorerHandoff` @ `d9fd564`.

### Gate 2 checkpoint — Fixtures & mapping (PASS)
Files added / changed:
- `hackathon/tool/examples/pah_muscle_malenfant2015.provenance.csv` (new; lossless superset
  schema: paper vs canonical ids, ratio, disease_log2fc, reported p-value text, correction notes)
- `hackathon/tool/examples/pah_muscle_malenfant2015.mapping_manifest.{csv,json}` (new; committed
  pre-approved mapping for the built-in fixture)
- `hackathon/tool/examples/pah_muscle_sensitivity_views.json` (new; full19/lower9/upper10/fibre)
- `hackathon/tool/src/motrpac_probe/signature.py` (added `candidates_for()`, `mapping_audit()`,
  `_all_maps()`; existing `map_genes`/`load` UNCHANGED so golden outputs are preserved)
- `hackathon/tool/tests/test_pah_fixture_table2.py` (new; 8 tests)
- `hackathon/tool/tests/test_mapping_candidates.py` (new; 4 tests)
- `hackathon/tool/.gitignore` (added `.venv/`)
- NOT changed: `pah_muscle_malenfant2015.csv`, `pah_muscle_lower9_malenfant2015.csv` (ADR-0005).

Why each change was necessary:
- Enrich via companion table, not by editing the fixture → preserves the live input + golden
  outputs while adding full paper↔canonical provenance (ADR-0005).
- `candidates_for`/`mapping_audit` → surface all candidates (176 Ensembl ids map to >1 human
  symbol) for the mandatory confirmation gate; never `distinct()` away duplication; conserve rows.
- Sensitivity views → planned subgroup analyses (incl. fibre/myosin competing explanation), not
  post-hoc gene removal.

Environment / setup:
- Created `hackathon/tool/.venv` (gitignored), `pip install -e ".[test]"`, `mprobe store fetch`
  (store bundle `mprobe-store-v1`, 10 files verified, store hash `d615c2f0b9aec182`).

Commands run (key):
- `python analysis.py` before edits → `/tmp/golden_full19_baseline.txt`
- `python -m pytest -q -rs` before edits → 110 passed, 34 skipped (allowlist baseline)
- `python -m pytest tests/test_pah_fixture_table2.py tests/test_mapping_candidates.py -q` → 12 passed
- `python analysis.py` after edits → byte-identical to baseline (diff empty) ⇒ NO scientific drift
- `python -m pytest -q -rs` after edits → 122 passed, 34 skipped

Golden numbers reproduced (full-19, unchanged):
- 19 genes counted, 25/25 rows mapped.
- Rat gastrocnemius male protein 8 wk: 13/19 opposed, cameraPR t +2.23 (matches QA).
- Rat gastrocnemius female protein 8 wk: 15/19 opposed, cameraPR t +4.62.

Expected-skip allowlist (34): test_golden.py (30), test_camera.py (3: 2 deck + 1 R env),
test_metab.py (1), test_store.py (1) — all "needs local hackathon/ data & deck numbers (not in
git)" or "R env not available". CI must use `pytest -q -rs` and allowlist exactly these.

Fixture provenance / integrity notes:
- Fixture `pah_muscle_malenfant2015.csv` sha256_16 = `0c80d82e5546ed9b`.
- All 25 rows map unambiguously to a single measured MoTrPAC symbol (n_ambiguous=0) under the
  released ID map; new uploads still require interactive confirmation.
- MYH1 correction recorded: paper prints Q9UKX2 (which is MYH2); canonical MYH1 = P12882. Both
  preserved. Isoform suffixes dropped for GLO1 (Q04760-2→Q04760) and ATP2A1 (O14983-2→O14983).
