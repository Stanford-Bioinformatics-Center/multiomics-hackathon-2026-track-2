# Merge-to-Main Plan and Pre-Merge Checklist (Documentation Only)

> [!IMPORTANT]
> **This document is DOCUMENTATION ONLY. It executes nothing.** Authoring or reading
> this plan does NOT trigger, invoke, or execute any `git`, `gh`, or shell command, and
> creating it as part of the *Finalize Discordance MVP* feature performs no version-control
> action of any kind (R10 AC1). Every command shown below is presented as literal text
> inside a fenced code block so it can be *read and copied*, never run automatically. A
> human operator runs these commands deliberately, later, only after the preconditions
> below are satisfied. The preconditions are **documented as text**; this plan does not
> itself check git state, run the verification battery, or query any remote.

This plan describes how the `integration/exercise-signature-explorer` branch is merged into
`main` once the work is complete. It records the preconditions that gate the merge, the
block rules that stop an unsafe merge, the rebase constraint that protects shared history,
the exact PR-creation command, and the required PR description template.

---

## 1. Pre-merge preconditions (ALL must hold before merge is permitted)

The merge is **permitted only when all three** of the following preconditions hold
simultaneously (R10 AC2). None of these is checked by this document — the operator verifies
each one manually before proceeding.

1. **All verification gates report a passing status.**
   Every gate in the verification battery (see `VERIFICATION.md` and the R9 CI gate in
   `.github/workflows/ci.yml`) reports **passing**:
   - Engine `pytest` battery: **exactly 123 passed / 34 skipped** (R9 AC7, updated to 123 by the motrpac-explorer merge — ADR-0023).
   - `apps/api` tests: **zero failures**.
   - Web `pnpm test` **and** `pnpm run build`: **zero failures**.
   - httpx API smoke check: passing (headline `q = 0.0584 ± 0.0001` and not-significant;
     discordance per-class counts `supported_concordant=1`, `supported_opposite=0`,
     `rna_response_protein_equivalent=285`, `indeterminate=5642`; metabolomics null result).
   - Playwright end-to-end battery: passing (every live view renders real data; the
     browser-level point→evidence→source chain holds).
   - Golden-hash guard: `analysis.py` is byte-identical to
     `hackathon/tool/analysis.py.sha256`.

2. **The working tree is clean.**
   `git status` reports **zero modified, zero staged, zero untracked, and zero deleted**
   files (R10 AC2). The command an operator uses to confirm this (documented, not executed):

   ```bash
   git status --porcelain
   ```

   A **clean** tree produces **no output** from that command. Any line of output means the
   precondition is NOT satisfied.

3. **The local branch is not ahead of its upstream remote.**
   The local `integration/exercise-signature-explorer` branch is **zero commits ahead** of
   its upstream remote (R10 AC2). The command an operator uses to confirm this (documented,
   not executed):

   ```bash
   git rev-list --count @{upstream}..HEAD
   ```

   A value of **`0`** means the branch is not ahead. Any value **≥ 1** means the precondition
   is NOT satisfied.

---

## 2. Block rules (when the merge is BLOCKED)

The merge is **BLOCKED** if **any** of the following is true. When blocked, this plan
requires that the operator **identify which specific precondition failed** before taking any
further action (R10 AC3):

- **Any verification gate reports a non-passing status** → BLOCKED.
  Identify the failing gate by name (e.g. "engine battery reported 121 passed / 34 skipped",
  "web build failed", "httpx smoke value mismatch on headline q", "golden-hash guard: mismatch").

- **The working tree contains one or more uncommitted or untracked files** → BLOCKED.
  I.e. `git status --porcelain` produced any output. Identify the offending files (modified,
  staged, untracked, or deleted).

- **The local branch is one or more commits ahead of its upstream remote** → BLOCKED.
  I.e. `git rev-list --count @{upstream}..HEAD` returned a value ≥ 1. Identify the count of
  unpushed commits.

| Failing condition | What to record | Resolution before retry |
|-------------------|----------------|-------------------------|
| A verification gate is non-passing | Which gate + the observed vs. expected value | Fix the underlying failure, re-run the full battery to green |
| Working tree not clean | The offending files from `git status --porcelain` | Commit or discard changes; ensure zero untracked files |
| Branch ahead of upstream | The number of unpushed commits | Push the branch so it is zero commits ahead of upstream |

Only when **none** of these block conditions holds may the operator proceed to Section 4.

---

## 3. Rebase constraint (protect shared history)

**No rebase may be performed that replays commits over
`origin/t3code/build-motrpac-probe-tool`** (R10 AC4).

`integration/exercise-signature-explorer` was branched off the tool branch, and
`origin/t3code/build-motrpac-probe-tool` is shared upstream history. Rebasing the
integration branch onto (or replaying its commits over) that branch would rewrite shared
commits and break provenance. The merge to `main` is therefore performed as a **merge /
pull request only** — never via a rebase that replays commits over
`origin/t3code/build-motrpac-probe-tool`. Do **not** run any `git rebase` whose base or
`--onto` target is `origin/t3code/build-motrpac-probe-tool`.

---

## 4. Create the pull request

Once every precondition in Section 1 holds and no block rule in Section 2 applies, the
operator creates the pull request with the **exact** command below. This is documented as
literal text — it is **not** executed by this plan (R10 AC1, R10 AC5):

```bash
gh pr create --base main --head integration/exercise-signature-explorer --title "Finalize Discordance MVP" --body-file PR_BODY.md
```

- **Base branch:** `main` (R10 AC5).
- **Head branch:** `integration/exercise-signature-explorer` (R10 AC5).
- The `--body-file` uses the PR description template in Section 5. (An operator may instead
  paste the template inline via `--body`; the three labeled sections are what matters.)

---

## 5. Pull request description template

The pull request description **must contain exactly three labeled sections** — a summary
section, a what-was-tested section, and a deferred-items section (R10 AC6). Copy the
template below into the PR body and fill each section:

```markdown
## Summary
<!-- What this PR does and why. The Track 2 "Omic Discordance Explained" work: the MoTrPAC
     Explorer is the default view (any molecule list in → every matching MoTrPAC comparison
     out, with a fused RNA↔protein co-view); alongside it, a read-only discordance API over
     committed demo outputs, a live query builder with a mapping-preview confirmation gate,
     generalized upload + shared results table, mock views retired, README rewritten to the
     eight-section scaffold anchored to Track 2, fresh-clone bootstrap, and automated
     verification wired into CI. Honest framing: cross-cohort association, NOT a PAH
     treatment claim; cohorts kept separate; headline male rat SKM-GN protein 8-week result
     is q = 0.0584 (not significant) under the frozen 52-column BH family. -->

## What was tested
<!-- The verification battery and its results: engine pytest (exactly 123 passed / 34
     skipped), apps/api tests (zero failures), web pnpm test + pnpm run build (zero
     failures), httpx API smoke (headline q = 0.0584 ± 0.0001, discordance per-class counts
     1/0/285/5642, metabolomics null result), Playwright e2e (every live view renders real
     data; browser point→evidence→source chain), and the golden-hash guard on analysis.py.
     Reference the VERIFICATION.md checkpoints. -->

## Deferred items
<!-- Out-of-scope follow-ons explicitly NOT included: live recompute of the discordance
     catalog (the adapter serves committed demo outputs only), Option-2 engine live
     metabolite scoring, and any other deferred work carried in HANDOFF_NEXT_FEATURES.md. -->
```

---

## Scope note (honest)

This plan is documentation only (task 13, R10 AC1). It records the merge procedure,
preconditions, block rules, rebase constraint, PR-create command, and PR template as **text**.
It does not run the verification battery, does not inspect the git working tree, does not
query any remote, and does not create a pull request. Every `git` and `gh` command above is
shown inside a fenced code block for a human operator to run deliberately, later, only after
the preconditions are satisfied. Engine, `apps/api`, and `apps/web` behavior is unchanged
(docs only).
