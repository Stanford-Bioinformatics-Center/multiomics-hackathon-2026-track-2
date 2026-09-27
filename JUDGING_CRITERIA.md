# Judging criteria — how this project maps to the Track 2 rubric

This document maps the Stanford Multi-omics Hackathon 2026 judging criteria (slide 36) to
concrete evidence in this repository and in the presentation. It also records the
minimum-defensible-submission checklist with an honest current status.

The judging criteria are evaluated through the **presentation**, the **GitHub repository**,
or **both**, as noted per row. There are no published weights, point values, or tie-breaks;
this map is for orientation, not scoring.

> **Framing note (from the slides).** Technical complexity is only one of seven criteria, and
> the deck explicitly warns that complexity alone does not make a project stronger. This
> project is deliberately a **coherent, reproducible analysis path with honest framing** — a
> useful prototype, not a polished product — which is the posture the slides reward.

---

## Criteria → evidence

| # | Criterion | Evaluated through | Where it is evidenced |
| --- | --- | --- | --- |
| 1 | **Scientific novelty** | Presentation | The **MoTrPAC Explorer**: any molecule list (genes / proteins / metabolites / pathways) in → every matching MoTrPAC comparison out, across layers/species/tissues/contrasts, with a fused RNA↔protein co-view for reading discordance directly. Plus the within-tissue discordance catalog + out-of-fold ridge model. See [README.md](README.md) §1–§3 and [WORKFLOWS.md](WORKFLOWS.md) §0. |
| 2 | **Impact of findings** | Presentation | A reusable paper→demo lookup for MoTrPAC discordance that any exercise-biology researcher can run on their own set. Honest headline (directed q = 0.0584, **not significant**) and a near-zero-R² weak-prediction model reported as a result, not a success claim. See [README.md](README.md) §1 (main result) and §7. |
| 3 | **Content** | Presentation | The presentation is expected to cover question → why → attempt → data/methods → result → meaning → what the data supports → what didn't work → next steps → contributions. The repo backs each with the section noted in the checklist below. **Gap:** the slide deck / demo recording is a presentation artifact, not in this repo. |
| 4 | **Collaboration** | Presentation + GitHub | Multiple contributors across engineering, science, and docs (see [README.md](README.md) §8 and the git history). Merge discipline preserving both authors' history is documented in [DECISIONS.md](DECISIONS.md) (ADR-0023) and [MERGE_PLAN.md](MERGE_PLAN.md). **Gap:** explicit per-person role assignments are marked to-be-confirmed in the README. |
| 5 | **Methods & approach** | Presentation + GitHub | Source-of-truth boundary (engines compute; API/UI serve/shape/render), the frozen 52-column BH family, and every design trade-off recorded as ADRs in [DECISIONS.md](DECISIONS.md) (ADR-0001…0024). Full method trace per analysis path in [WORKFLOWS.md](WORKFLOWS.md). |
| 6 | **Documentation & integrity** | GitHub | Eight-section README scaffold ([README.md](README.md)); reproducible clean-clone log with commands, counts, and hashes ([VERIFICATION.md](VERIFICATION.md)); honest framing throughout (cross-cohort association, cohorts never merged, PTM occupancy caveat); AI usage disclosed with the integrity boundary; golden-hash guard on `analysis.py`. |
| 7 | **Technical complexity** | GitHub | Engine (`motrpac_probe`) + FastAPI service + React app + `query_core` + standalone discordance/metabolomics modules, with property-based tests (Hypothesis + fast-check), httpx smoke, and Playwright e2e wired into CI. Framed per the slide caveat: complexity in service of a coherent, reproducible result — not complexity for its own sake. See [VERIFICATION.md](VERIFICATION.md) and `.github/workflows/ci.yml`. |

---

## Reproducibility anchors (the defensible result)

Verified in [VERIFICATION.md](VERIFICATION.md) on a clean clone:

- Directed headline **q = 0.0584**, labeled **not significant** (frozen 52-column BH family).
- Discordance per-class counts: `supported_concordant` = 1, `supported_opposite` = 0,
  `rna_response_protein_equivalent` = 285, `indeterminate` = 5642.
- Metabolomics null result: `hits_fraction = 0.0`, Fisher `p = 0.609`.
- Battery: engine **123 passed / 34 skipped**, `apps/api` **81 passed**, web **38 passed**
  + build OK, golden-hash on `analysis.py` **OK**
  (`1c1089723ec7c9d67c4c0f7ab51262717f80d95322c7da8df60187a06cab7372`).

---

## Minimum-defensible-submission checklist

Status legend: **[x]** met in-repo · **[~]** partially met / needs a team action ·
**[ ]** presentation artifact, tracked here but produced outside the repo.

- [x] GitHub repository is accessible to judges — public repo with root docs.
- [x] README identifies the question, team, users, and significance — [README.md](README.md) §1–§2, §8.
- [x] Methods and data provenance are documented — [README.md](README.md) §6, [WORKFLOWS.md](WORKFLOWS.md), [DECISIONS.md](DECISIONS.md).
- [x] Setup and execution commands documented and tested — [BOOTSTRAP.md](BOOTSTRAP.md), [RUN_LOCAL.md](RUN_LOCAL.md), [VERIFICATION.md](VERIFICATION.md).
- [x] Input files and required columns are explained — [README.md](README.md) §5, [`hackathon/tool/store/SCHEMA.md`](hackathon/tool/store/SCHEMA.md).
- [x] At least one result is reproducible — the headline battery in [VERIFICATION.md](VERIFICATION.md).
- [x] Expected outputs and their locations are documented — [README.md](README.md) §5, [WORKFLOWS.md](WORKFLOWS.md) contract table.
- [~] Figures or screenshots demonstrate that the project works — committed demo plots exist; **Explorer/Discordance UI screenshots to be inserted** (placeholder in [README.md](README.md) §7).
- [x] Known limitations and failure modes stated honestly — [README.md](README.md) §7.
- [~] Contributor roles are recorded — names recorded; **per-person roles to be confirmed** ([README.md](README.md) §8).
- [x] External code and AI assistance are disclosed — [README.md](README.md) §6 (AI usage + vendor licenses).
- [x] License and citations are included — MIT [LICENSE](LICENSE); citations in [README.md](README.md) §6.
- [ ] Presentation tells the same story as the repository — **presentation artifact** (ensure it leads with the Explorer paper→demo trace and the honest headline).
- [ ] Demo path is short and reliable — **presentation artifact** (recommended: bootstrap → open Explorer → run a cited example).
- [x] Next steps clearly separated from completed work — [README.md](README.md) §8 + [HANDOFF_NEXT_FEATURES.md](HANDOFF_NEXT_FEATURES.md).

**Open items for the team before submission:** (1) insert Explorer/Discordance screenshots
into README §7; (2) confirm per-person contributor roles in README §8; (3) prepare the
short demo/presentation so it tells the same Explorer-first story as this repository.
