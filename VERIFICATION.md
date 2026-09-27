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
