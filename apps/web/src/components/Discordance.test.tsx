// Discordance component tests (Vitest + Testing Library) — task 4.6.
//
// These example-based tests complement the fast-check suite in Discordance.property.test.tsx.
// They render <Discordance /> against a mocked successful catalog/model/PTM payload and assert:
//   1. Three-stage order (R3 AC2): the Stage 1 (QUERY) → Stage 2 (VISUALIZE) → Stage 3 (INTERPRET)
//      markers appear in the DOM in that document order.
//   2. PTM-parent divergence summary present (R3 AC6) when PTM data is provided.
//   3. PTM occupancy caveat present (R3 AC7) as visible text.
//   4. Weak-prediction statement present (R3 AC5) alongside the model metrics.
//
// Requirements: 3.2, 3.5, 3.6, 3.7.

import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { cleanup, render, screen, within } from "@testing-library/react";
import type {
  DiscordanceCatalogResponse,
  DiscordanceModelResponse,
  PtmParentResponse,
} from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
// Discordance.tsx calls api.discordanceCatalog() then, best-effort, api.discordanceModel() and
// api.discordancePtmParent(). Back each with a mutable holder so each test can install its own
// resolved/rejected payload before rendering. ApiError stays a real class because the component's
// errText does `e instanceof ApiError` (mirrors the mocking approach in Discordance.property.test.tsx).
type CatalogResult = { ok: true; value: DiscordanceCatalogResponse } | { ok: false };
type ModelResult = { ok: true; value: DiscordanceModelResponse } | { ok: false };
type PtmResult = { ok: true; value: PtmParentResponse } | { ok: false };

const holder: { catalog: CatalogResult; model: ModelResult; ptm: PtmResult } = {
  catalog: { ok: false },
  model: { ok: false },
  ptm: { ok: false },
};

vi.mock("../api/client", async () => {
  const actual = await vi.importActual<typeof import("../api/client")>("../api/client");
  return {
    ApiError: actual.ApiError,
    api: {
      discordanceCatalog: () =>
        holder.catalog.ok
          ? Promise.resolve(holder.catalog.value)
          : Promise.reject(new actual.ApiError(500, "catalog unavailable")),
      discordanceModel: () =>
        holder.model.ok
          ? Promise.resolve(holder.model.value)
          : Promise.reject(new actual.ApiError(404, "model unavailable")),
      discordancePtmParent: () =>
        holder.ptm.ok
          ? Promise.resolve(holder.ptm.value)
          : Promise.reject(new actual.ApiError(404, "ptm unavailable")),
    },
  };
});

// Import AFTER vi.mock so the component binds to the mocked client.
import Discordance from "./Discordance";

// ─────────────────────────── Realistic mock payloads ───────────────────────────
// Consistent with the client interfaces in api/client.ts and the committed demo counts documented
// in the spec (supported_concordant=1, supported_opposite=0, rna_response_protein_equivalent=285,
// indeterminate=5642, no_protein_measurement=9228, no_rna_measurement=255).
const WEAK_PREDICTION_NOTE =
  "Out-of-fold R² is near zero: the temporal model does not meaningfully predict protein change. " +
  "Treat this as a weak prediction, a result rather than a success claim.";
const OCCUPANCY_CAVEAT =
  "PTM signal is not automatically a measure of modification occupancy, enzyme activity, or " +
  "functional consequence.";

function makeCatalog(): DiscordanceCatalogResponse {
  return {
    run_id: "run-discordance-demo",
    schema_version: "1",
    analysis_type: "discordance_catalog",
    case_study_id: "cs-muscle-ee",
    cohort_label: "Rat SKM-GN · endurance exercise vs control",
    title: "Omic discordance catalog",
    tissue: "SKM-GN",
    contrast_category: "EE-CON",
    total_events: 1 + 0 + 285 + 5642,
    classification_counts: {
      supported_concordant: 1,
      supported_opposite: 0,
      rna_response_protein_equivalent: 285,
      indeterminate: 5642,
    },
    coverage_counts: {
      no_protein_measurement: 9228,
      no_rna_measurement: 255,
    },
    occupancy_caveat: OCCUPANCY_CAVEAT,
    limitations: ["Cross-cohort association only; not a treatment claim."],
    provenance: { source: "committed demo outputs" },
  };
}

function makeModel(): DiscordanceModelResponse {
  return {
    run_id: "run-discordance-model",
    schema_version: "1",
    analysis_type: "discordance_model",
    case_study_id: "cs-muscle-ee",
    tissue: "SKM-GN",
    contrast_category: "EE-CON",
    metrics: [
      { subset: "all_mapped", model: "zero", n_rows: 5928, n_genes: 5928, mae: 0.412, rmse: 0.601, r2: 0.0 },
      { subset: "all_mapped", model: "rna_only", n_rows: 5928, n_genes: 5928, mae: 0.408, rmse: 0.598, r2: 0.01 },
      { subset: "all_mapped", model: "temporal", n_rows: 5928, n_genes: 5928, mae: 0.405, rmse: 0.595, r2: 0.02 },
      { subset: "rna_responsive", model: "zero", n_rows: 286, n_genes: 286, mae: 0.52, rmse: 0.71, r2: -0.03 },
      { subset: "rna_responsive", model: "rna_only", n_rows: 286, n_genes: 286, mae: 0.51, rmse: 0.7, r2: 0.0 },
      { subset: "rna_responsive", model: "temporal", n_rows: 286, n_genes: 286, mae: 0.5, rmse: 0.69, r2: 0.04 },
    ],
    weak_prediction_note: WEAK_PREDICTION_NOTE,
    limitations: ["Out-of-fold R² near zero across all variants."],
    provenance: { source: "committed demo outputs" },
  };
}

function makePtm(): PtmParentResponse {
  return {
    run_id: "run-ptm-parent",
    schema_version: "1",
    analysis_type: "ptm_parent_audit",
    summary_rows: [
      {
        tissue: "SKM-GN",
        contrast_category: "EE-CON",
        timepoint: "8w",
        raw_phosphosite_rows: 1200,
        after_mapping_checks_with_uniprot: 980,
        excluded_missing_or_ambiguous_site_mapping: 220,
        matched_to_single_gene_uniprot_single_feature_parent: 640,
        matched_with_both_site_and_parent_ci: 512,
        supported_phosphosite_rows_among_matched: 88,
        parent_equivalent_rows_among_matched: 401,
        supported_phosphosite_and_parent_equivalent_rows: 37,
        candidate_genes: 44,
      },
    ],
    candidate_count: 44,
    occupancy_caveat: OCCUPANCY_CAVEAT,
    provenance: { source: "committed demo outputs" },
  };
}

beforeEach(() => {
  holder.catalog = { ok: true, value: makeCatalog() };
  holder.model = { ok: true, value: makeModel() };
  holder.ptm = { ok: true, value: makePtm() };
});

afterEach(() => {
  cleanup();
  holder.catalog = { ok: false };
  holder.model = { ok: false };
  holder.ptm = { ok: false };
});

describe("Discordance view — three-stage flow and honest framing (R3)", () => {
  test("renders the three stage markers in QUERY → VISUALIZE → INTERPRET document order (R3 AC2)", async () => {
    const { container } = render(<Discordance />);

    // Wait for the ready state.
    await screen.findByText("STAGE 1 · QUERY");

    // The stage markers are the `.status-mark--assumption` strips. Query them in document order and
    // filter to the three stage strips (the "WEAK PREDICTION" marker shares the same class).
    const marks = Array.from(
      container.querySelectorAll<HTMLElement>(".status-mark--assumption"),
    ).map((el) => el.textContent?.trim() ?? "");

    const stageMarks = marks.filter((t) => t.startsWith("STAGE "));
    expect(stageMarks).toEqual(["STAGE 1 · QUERY", "STAGE 2 · VISUALIZE", "STAGE 3 · INTERPRET"]);

    // Cross-check with compareDocumentPosition: stage 1 precedes stage 2 precedes stage 3.
    const stage1 = screen.getByText("STAGE 1 · QUERY");
    const stage2 = screen.getByText("STAGE 2 · VISUALIZE");
    const stage3 = screen.getByText("STAGE 3 · INTERPRET");
    expect(stage1.compareDocumentPosition(stage2) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
    expect(stage2.compareDocumentPosition(stage3) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
  });

  test("renders the PTM-parent divergence summary when PTM data is provided (R3 AC6)", async () => {
    render(<Discordance />);

    // The Stage-3 PTM-parent divergence summary header includes the candidate-gene count.
    const summaryHeader = await screen.findByText(/PTM-parent divergence summary/);
    expect(summaryHeader).toBeTruthy();
    expect(summaryHeader.textContent).toContain("44 candidate genes");

    // The summary table shows the timepoint row from the payload.
    const table = summaryHeader.closest(".visual-card");
    expect(table).not.toBeNull();
    expect(within(table as HTMLElement).getByText("8w")).toBeTruthy();
  });

  test("renders the PTM occupancy caveat as visible text (R3 AC7)", async () => {
    render(<Discordance />);

    const caveat = await screen.findByText(OCCUPANCY_CAVEAT);
    expect(caveat).toBeTruthy();

    // It sits under the visible "Occupancy caveat" heading.
    expect(screen.getByText("Occupancy caveat")).toBeTruthy();
  });

  test("renders a visible weak-prediction statement alongside the model metrics (R3 AC5)", async () => {
    render(<Discordance />);

    // The weak-prediction note renders inside the model-metrics card.
    const note = await screen.findByText(WEAK_PREDICTION_NOTE);
    expect(note).toBeTruthy();

    // A visible "WEAK PREDICTION" marker labels it.
    expect(screen.getByText("WEAK PREDICTION")).toBeTruthy();

    // The note appears alongside the metrics: within the same model-metrics card that carries the
    // "Model metrics (R² / MAE / RMSE)" header.
    const metricsHeader = screen.getByText("Model metrics (R² / MAE / RMSE)");
    const metricsCard = metricsHeader.closest(".visual-card");
    expect(metricsCard).not.toBeNull();
    expect(within(metricsCard as HTMLElement).getByText(WEAK_PREDICTION_NOTE)).toBeTruthy();
  });
});
