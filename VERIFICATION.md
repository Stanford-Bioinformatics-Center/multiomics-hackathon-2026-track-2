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
