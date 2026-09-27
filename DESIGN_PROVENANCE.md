# Design Provenance — Migrated Mock Views

This document preserves the content of the original **Figma-Make** system-design mock
views (numbered 01–04 and 08) that were built during early ideation and have since been
retired from the live web application in favor of the working analysis views
(Discordance, Live, Generalized, Metabolomics). It exists so the design story, the
narrative content, and the Figma-Make design provenance are not lost when the mock views
are removed from `apps/web/src/App.tsx` (mock-view removal is handled in a separate task).

> **Purpose and status.** Everything below is **design mock content**, not a live result.
> The mock views used illustrative, hardcoded placeholder values (for example the
> gene identifiers `PPARGC1A`, `SOD2`, and `COL1A1`) purely to communicate layout and
> intent. None of these numbers are computed by the engines and none should be read as an
> actual analysis result. The authoritative, computed workflows live in the live views and
> are documented in the analysis workflow reference (see below).

## Figma-Make provenance (applies to all migrated views)

- **Design origin:** These views were authored in **Figma Make**. The Figma-Make tooling
  and generated site scaffolding remain visible in the repository under
  `apps/web/.figma/` and in `apps/web/vite.config.ts` (the `figma*` Vite plugins and
  `.figma/make/site.json`).
- Each migrated view below additionally carries its own **Figma-Make design origin**
  attribution line so provenance is preserved per view, per Requirement 6.5.

## Relationship to the analysis workflow diagrams

The authoritative mermaid workflow diagrams live in the root `WORKFLOWS.md` (authored under
Requirement 1). To avoid divergent copies, this document **does not duplicate** those
diagrams; it **references them by name/anchor** instead. The diagram sources in
`WORKFLOWS.md` are the single source of truth and are not modified by this migration.

Referenced `WORKFLOWS.md` diagrams (one per analysis type):

- `directed_signature`
- `metabolomics_case_study`
- `generalized_query`
- discordance catalog + predictive model
- PTM-parent audit

When reading a migrated view below, cross-reference the corresponding `WORKFLOWS.md`
diagram for the executable, computed workflow.

---

## View 01 — System architecture

> _Figma-Make design origin. Migrated from the retired `Architecture` mock view (nav
> eyebrow "01", "System architecture"). Mock content — illustrative only._

**Kicker:** Technical handoff
**Title:** A deterministic evidence pipeline
**Copy:** Five inspectable layers carry a directed disease signature into a cross-cohort,
multi-omic comparison. Every transform remains visible.

The architecture was presented as five stacked layers. Nodes marked _(planned)_ were
future/aspirational capabilities in the mock, not implemented behavior.

**Layer 01 — User + frontend**
Researcher · Curated PAH example · Structured CSV upload · Query builder · Context filters
(badge: _source species · target species · derived design · tissue · sex · time point · omics_) ·
Results dashboard · Export · _(planned)_ Natural-language query

**Layer 02 — Application API**
Signature validation · Dataset availability resolver · Query orchestration · Results API ·
Provenance API · Export service

**Layer 03 — Scientific processing**
Identifier normalization · Alias resolution · Human → rat ortholog mapping ·
Omic-specific adapters · Disease ↔ exercise comparison · Cross-omic classifier ·
Pathway / context annotation · Evidence + limitation labels

**Layer 04 — Data products**
Directed signature registry · MoTrPAC transcriptomics · MoTrPAC proteomics ·
Identifier + ortholog tables · Pathway resources · Versioned manifests ·
Study-context catalog · Layer capability registry ·
_(planned)_ Genomics adapter · _(planned)_ Epigenomics adapter · _(planned)_ Metabolomics adapter

**Layer 05 — Outputs**
Disease ↔ exercise plot · RNA ↔ protein plot · Gene-by-layer heatmap · Evidence summary cards ·
Filterable results · Mapping + provenance report · CSV / JSON bundle · _(planned)_ Predictive modeling

**Optional LLM assistant _(planned)_.** A convenience layer isolated from the scientific
calculation path: (01) natural language → structured query, (02) plain-language result
summaries.

**System key.** Available in PAH demo · Planned capability · Assumption / caution.

**Deterministic by design.** All identifiers, contrasts, mappings, and calculations are
deterministically validated.

**Assumption strip.** Human → rat mapping and RNA → protein comparison are explicit
analytical assumptions, preserved in every provenance report.

_Cross-reference: `WORKFLOWS.md` → `directed_signature` diagram._

---

## View 02 — Technical flow

> _Figma-Make design origin. Migrated from the retired `TechnicalFlow` mock view (nav
> eyebrow "02", "Technical flow"). Mock content — illustrative only._

**Kicker:** Deployment and data flow
**Title:** Every interaction has a named contract
**Copy:** The target species resolves a read-only study design from the catalog before
availability, mapping, or scientific processing begins.

This view embedded the accessible SVG technical flow diagram
(`apps/web/src/components/TechnicalFlowDiagram.tsx`), whose accessible summary reads: a
researcher submits a directed disease signature, its source species, a target MoTrPAC
species, and one or more omic layers; the catalog derives study design from species; query
validation and dataset availability precede identifier mapping; compatible contexts continue
through independent omic adapters into multi-omic evidence integration; incompatible contexts
produce a scientific empty state and alternatives; the optional language model remains
outside the calculation path.

**API contract (mock illustrative endpoint list):**

- `GET /api/v1/catalog`
- `GET /api/v1/catalog/species/:species`
- `POST /api/v1/signatures/validate`
- `POST /api/v1/availability`
- `POST /api/v1/comparisons`
- `GET /api/v1/comparisons/:id/features`
- `GET /api/v1/comparisons/:id/provenance`
- `GET /api/v1/comparisons/:id/export`

> These are the mock's illustrative endpoint names. The **authoritative** request→response
> contract for the real API lives in the `WORKFLOWS.md` endpoint contract table.

**Validation contract.** All identifiers, contexts, mappings, and calculations are
deterministic. The optional LLM may propose a query, explain an empty state, or summarize
validated output. It cannot alter species, omics, identifiers, effects, significance, or
classifications.

_Cross-reference: `WORKFLOWS.md` → endpoint contract table and the `directed_signature`
diagram._

---

## View 03 — Researcher workflow

> _Figma-Make design origin. Migrated from the retired `Workflow` mock view (nav eyebrow
> "03", "Researcher workflow"). Mock content — illustrative only._

**Kicker:** Researcher workflow
**Title:** No feature disappears in the pipeline
**Copy:** A deliberate confirmation checkpoint separates computational normalization from
researcher-approved analysis.

**Current comparison contract (mock).** PAH skeletal muscle vs healthy control →
Exercise response: Chronically trained rat vs sex-matched sedentary rat.

**Workflow steps (ordered):**

1. **Choose** — Select the curated PAH example or upload a directed signature.
2. **Review** — Inspect recognized identifiers, directions, effects, and source.
3. **Resolve** — Resolve aliases and human-to-rat species mappings.
4. **Confirm** _(required confirmation checkpoint)_ — Approve, edit, or exclude every
   normalization decision.
5. **Resolve context** — Choose the target MoTrPAC species; study design is derived from the
   catalog.
6. **Compare** — Choose compatible omics, tissue, sex, and time point, then run precomputed
   results.
7. **Classify** — Label disease opposition and RNA/protein agreement.
8. **Explore** — Review plots, pathways, provenance, and limitations.
9. **Export** — Download a reproducible CSV/JSON result bundle.

**Explicit mapping contract.** The application must never silently drop or remap a feature.
Unmapped and ambiguous records remain visible, downloadable, and attributable to the source
row.

**Catalog-derived context — the researcher selects species, not study design:**

- Target species **Human** → Acute exercise study
- Target species **Rat** → Chronic exercise-training study
- Context mismatch → Scientific empty state (alternatives remain visible)
- Planned layer → Unsupported omic (run supported layers instead)

**Independent classification A — Disease relationship:** Opposes disease direction · Matches
disease direction · Near-zero / uncertain · Missing / unmapped.

**Independent classification B — Cross-omic relationship (MVP):** Concordant up / concordant
down · Higher measured RNA · lower measured protein · Lower measured RNA · higher measured
protein · RNA only · Protein only · Near-zero · Insufficient data.

**Strongest exploratory state.** RNA and protein agree — and both oppose the disease
direction. A prioritized signal for hypothesis generation and follow-up, **not** evidence
that exercise treats PAH. RNA and protein remain independent measurements; neither is treated
as a proxy for the other.

_Cross-reference: `WORKFLOWS.md` → `directed_signature` and `generalized_query` diagrams._

---

## View 04 — Results dashboard

> _Figma-Make design origin. Migrated from the retired `Dashboard` mock view (nav eyebrow
> "04", "Results dashboard"). Mock content — illustrative only; all figures below are
> hardcoded placeholders, not computed results._

**Kicker:** Desktop wireframe
**Title:** Evidence first, claims constrained
**Copy:** The dashboard keeps contrast definitions, assumptions, missingness, and provenance
beside the scientific signal.

**Top bar (mock).** Brand "Exercise Signature Explorer — Cross-cohort molecular comparison";
current signature "PAH skeletal muscle v1.2"; dataset/version derived from study context;
Methodology link; Export bundle action.

**Contrast banner (mock).** Disease numerator: Human PAH skeletal muscle vs denominator:
Healthy human control. Exercise numerator: Chronically trained rat (rat) / Human post-acute
exercise (human) vs denominator: Sex-matched sedentary rat / Matched human baseline.

**Summary cards (mock placeholder counts):**

| Value | Label | Detail |
| --- | --- | --- |
| 128 | Input features | 100% retained |
| 112 | Successfully mapped | 87.5% |
| 47 | Opposite direction | 42% of mapped |
| 31 | Adjusted-significant | cross-cohort support |
| 26 | RNA / protein agree | of 39 dual-measured |
| 16 | Missing / ambiguous | requires review |

**Dashboard tabs:** Disease vs Exercise · RNA vs Protein · Heatmap · Evidence Table ·
Methods & Limitations.

- **Disease vs Exercise / RNA vs Protein** — quadrant scatter plots. With a single omic
  layer selected, the RNA vs Protein tab showed an empty state ("One layer selected —
  disease-versus-exercise analysis remains available"); with two layers a quadrant plot; the
  matrix mode showed a heatmap plus pairwise concordance matrix.
- **Heatmap** — a gene-by-layer signed-effect grid (Disease / RNA / Protein / Metabolite
  columns; Metabolite shown as a planned layer).
- **Pairwise matrix** — RNA/Protein/Metabolite concordance matrix (Metabolite entries
  labeled planned).
- **Evidence Table (mock placeholder rows):**

| Feature | Mapping | Disease | RNA | Protein | Classification |
| --- | --- | --- | --- | --- | --- |
| PPARGC1A | PPARGC1A → Ppargc1a | −1.24 | +0.88 | +0.61 | Opposes · concordant |
| SOD2 | SOD2 → Sod2 | −0.76 | +0.54 | +0.42 | Opposes · concordant |
| COL1A1 | COL1A1 → Col1a1 | +1.45 | −0.63 | −0.38 | Opposes · concordant |
| STAT3 | STAT3 → Stat3 | +0.71 | +0.39 | +0.26 | Matches · concordant |
| VEGFA | VEGFA → Vegfa | +0.32 | +0.05 | — | Near-zero · RNA only |

  Expandable provenance rows were described as including source accession, numerator,
  denominator, q-values, ortholog method, mapping confidence, and pipeline version.

- **Methods & Limitations** — four method blocks (cross-cohort source studies; reproducible
  version manifests; visible mapping assumptions; and a warning block: no treatment or causal
  claim — directional opposition is exploratory evidence for follow-up and does not show that
  exercise treats PAH). Plus a scientific-state gallery: no matching tissue/time point; no
  mapped features; no adjusted-significant support (measurements retained); technical error
  only.

_Cross-reference: `WORKFLOWS.md` → `directed_signature` diagram (the live counterpart is the
Live view)._

---

## View 08 — Judge slide

> _Figma-Make design origin. Migrated from the retired `JudgeSlide` mock view (nav eyebrow
> "08", "Judge slide"). Mock content — illustrative only._

**Kicker:** Presentation slide
**Title:** One contract, extensible evidence
**Copy:** A concise version of the system story for review panels and judges.

**Headline.** Exercise Signature Explorer — From a directed signature to inspectable
multi-omic evidence. Thesis: cross-cohort comparison for hypothesis generation.

**Flow (ordered steps):**

1. Directed disease signature — Signed effects + source
2. Species + omic selection — Study design derived from catalog
3. Availability + mapping validation — No silent reinterpretation
4. Compatible MoTrPAC data — Independent RNA + protein
5. Two-axis classification — Opposition + concordance
6. Interactive evidence dashboard — Plots + provenance + export

**Result framing.**

- **Prioritize** — Concordant RNA + protein that oppose disease direction.
- **Follow up** — Generate hypotheses; do not infer treatment or causality.
- **Reproduce** — Export mappings, versions, assumptions, and calculations.

**Current MVP.** Human → acute context from catalog · Rat → chronic context from catalog ·
Transcriptomics + proteomics implemented · Additional omics use the same planned adapter
contract.

**Caption.** "PAH skeletal muscle is our worked example; the same contract can accept
additional directed signatures and omic-specific adapters."

_Cross-reference: `WORKFLOWS.md` → all five analysis-type diagrams (this slide summarizes the
whole system story)._

---

## Notes

- **No diagram sources were modified** by this migration. The mermaid diagrams remain
  authoritative in `WORKFLOWS.md` and are referenced above by name/anchor only. If
  `WORKFLOWS.md` is not yet present at the time you read this, the anchor names above are the
  agreed diagram identifiers to link to once it lands (per the R1 task).
- The retired mock views used illustrative placeholder identifiers (`PPARGC1A`, `SOD2`,
  `COL1A1`) that must not appear in any navigable live view. They are retained here only as
  historical mock content.
- The About/Methods navigation entry links to this document (see the R6 nav task).
