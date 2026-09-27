# Metabolomics in the app

Two complementary metabolomics assets live in this repo. They are NOT competing answers; they are
different layers, mirroring the protein/RNA split already in the project.

| Asset | What it is | Analog on the protein/RNA side |
|---|---|---|
| Engine METAB layer (`hackathon/tool/src/motrpac_probe/metab.py`, 217 store columns) | the reusable metabolite **scorer**: score any directed metabolite signature against the MoTrPAC METAB store (cameraPR, pathway panel) | `mprobe` core scorer |
| Standalone module (`MoTrPAC Hackathon/Metabolomics/`) | a specific PAH **case study** (ST000763 SSc-PAH plasma) analysed end-to-end and mapped onto MoTrPAC blood+muscle metabolomics | the `MoTrPAC Hackathon/` R scripts (case study) |

"Source of truth" (ADR-0001) means the code that computes a given number. The two assets compute
**different comparisons on different data**, so there is no single number for them to disagree on —
this was never a "double truth" problem.

## What was wired into the app (this pass)

A **separate, read-only analysis type**: `metabolomics_case_study`.

- **API (Option-1 adapter).** `apps/api/motrpac_probe_service/metab_casestudy.py` reads the module's
  COMMITTED outputs (`06_context_table_hits.csv`, `06_context_summary.json`, `source_manifest.json`)
  and normalizes them to the service's evidence/provenance schema. No statistics, no R at request
  time. Endpoints: `GET /api/metabolomics/casestudy`, `.../casestudy/export`,
  `GET /api/metabolomics/convergence`. The catalog gains an `analysis_types` registry so the UI can
  tell the directed-signature run apart from this case study.
- **UI.** `apps/web/src/components/MetabolomicsCaseStudy.tsx` is its own top-level view
  ("Metabolomics (ST000763)"), never a tab inside the muscle-protein dashboard, so the two PAH
  stories are never blended. It renders the honest conclusion, the 41 resting hits with MoTrPAC
  context labels, the null-result stats, both contrasts + the EE-EE like-for-like caveat, the
  convergence badge, and an export link.

The engine's live directed-signature metabolite scoring (Option 2) remains a follow-on: the
capability matrix keeps that path `api:false / react:false`; only the case-study analysis type is
exposed here.

## Convergence cross-check (integrity win)

`apps/api/motrpac_probe_service/convergence.py` reconciles the module's `04` MoTrPAC export against
the engine METAB store on shared **metabolite × tissue × timepoint** EE-CON cells. Both were built
from `MotrpacHumanPreSuspensionAnalysis` 2.0.8 (c2.0) by two independent code paths.

Result: on the **unambiguous subset** (RefMet name unique within a tissue×timepoint cell on both
sides), the EE-CON logFC values are **identical — 1,156/1,156 cells, max abs diff 0.0**. Cells whose
RefMet name is blank or maps to multiple platform features are resolved to feature_id-level rows
differently by the two sides and are excluded from the strict comparison (key ambiguity, not a
scientific disagreement). Test: `apps/api/tests/test_metabolomics.py::test_convergence_identical_on_unambiguous_subset`.

## The scientific caveat (carried verbatim into the UI and exports)

In ST000763, the acute exercise response does not differ between SSc-PAH and healthy people (0 hits).
At rest PAH differs strongly from healthy, but normal-pressure SSc sampled in the same catheterization
setting shows almost the same differences (PAH vs that group: 0 hits). The leading explanation is
**sampling setting and/or scleroderma, not PAH**. This is a cross-study description, not proof of
mechanism, and not evidence that exercise treats PAH. All PAH samples are catheterization-setting and
all healthy samples are non-invasive; rest/peak pairs are inferred (108, unverified). PTM-style
occupancy caveats do not apply here, but the "not modification occupancy / enzyme activity /
functional consequence" caution is retained across the app for consistency.
