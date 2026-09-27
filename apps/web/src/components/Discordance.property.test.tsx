// # Feature: finalize-discordance-mvp, Property 14
//
// Property 14 — result completeness + empty-state exclusivity (Requirements 3.3, 3.4, 3.8, 3.9).
//
// On SUCCESS (a non-empty discordance catalog is delivered):
//   - all four classification-class counts render, INCLUDING any class whose count is zero, and
//   - for each model variant present in the payload (zero / rna_only / temporal) its
//     out-of-fold R² / MAE / RMSE render, and
//   - the empty/error block is ABSENT.
//
// On FAILURE / EMPTY (the catalog fetch rejects, or the delivered catalog carries no events —
// total_events === 0, or all four classes sum to zero):
//   - the empty/error block is PRESENT, and
//   - NONE of the four class counts, the model metrics, or the PTM-parent summary render.
//
// These two branches are mutually exclusive: exactly one of them holds for any generated payload.
// fast-check drives ≥100 iterations; the API client module is mocked so no network is touched.

import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, render, screen, waitFor, within } from "@testing-library/react";
import fc from "fast-check";
import type {
  DiscordanceCatalogResponse,
  DiscordanceModelResponse,
  DiscordanceModelMetricRow,
  PtmParentResponse,
} from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
// Discordance.tsx calls api.discordanceCatalog() then, best-effort, api.discordanceModel() and
// api.discordancePtmParent(). We back each with a mutable holder so every fast-check iteration can
// install its own resolved/rejected payloads before rendering. ApiError must remain a real class
// (the component's errText does `e instanceof ApiError`).
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

// ─────────────────────────── Arbitraries ───────────────────────────
const MODELS: Array<DiscordanceModelMetricRow["model"]> = ["zero", "rna_only", "temporal"];
const SUBSETS: Array<DiscordanceModelMetricRow["subset"]> = ["all_mapped", "rna_responsive"];

const finiteNum = fc.float({ min: -50, max: 50, noNaN: true, noDefaultInfinity: true });
const count = fc.integer({ min: 0, max: 10000 });

/**
 * Four classification counts where at least ONE is zero and the total is > 0, so the catalog is a
 * genuine non-empty result whose rendering must still surface a zero class (Requirement 3.3).
 */
const nonEmptyClassCounts = fc
  .tuple(count, count, count, count)
  .chain(([a, b, c, d]) => {
    // Force at least one zero, and guarantee at least one positive so the catalog is non-empty.
    const arr = [a, b, c, d];
    return fc.integer({ min: 0, max: 3 }).chain((zeroAt) => {
      arr[zeroAt] = 0;
      return fc.integer({ min: 0, max: 3 }).map((posAt) => {
        const at = posAt === zeroAt ? (posAt + 1) % 4 : posAt;
        if (arr[at] === 0) arr[at] = 1; // ensure a positive so total > 0
        return {
          supported_concordant: arr[0],
          supported_opposite: arr[1],
          rna_response_protein_equivalent: arr[2],
          indeterminate: arr[3],
        };
      });
    });
  });

function catalogArb(
  classification: DiscordanceCatalogResponse["classification_counts"],
): fc.Arbitrary<DiscordanceCatalogResponse> {
  const total =
    classification.supported_concordant +
    classification.supported_opposite +
    classification.rna_response_protein_equivalent +
    classification.indeterminate;
  return fc.record({
    run_id: fc.string(),
    schema_version: fc.constant("1"),
    analysis_type: fc.constant("discordance_catalog" as const),
    case_study_id: fc.string({ minLength: 1 }).map((s) => `cs-${s}`),
    cohort_label: fc.string(),
    title: fc.string(),
    tissue: fc.constantFrom("SKM-GN", "SKM-VL", "muscle"),
    contrast_category: fc.constantFrom("EE-CON", "acute", "chronic"),
    total_events: fc.constant(total),
    classification_counts: fc.constant(classification),
    coverage_counts: fc.record({
      no_protein_measurement: count,
      no_rna_measurement: count,
    }),
    occupancy_caveat: fc.string(),
    limitations: fc.array(fc.string(), { maxLength: 3 }),
    provenance: fc.constant({}),
  });
}

function modelArb(models: Array<DiscordanceModelMetricRow["model"]>): fc.Arbitrary<DiscordanceModelResponse> {
  const metricRows: fc.Arbitrary<DiscordanceModelMetricRow[]> = fc
    .tuple(
      ...models.map((m) =>
        fc.record({
          subset: fc.constantFrom(...SUBSETS),
          model: fc.constant(m),
          n_rows: fc.integer({ min: 1, max: 9999 }),
          n_genes: fc.integer({ min: 1, max: 9999 }),
          mae: finiteNum,
          rmse: finiteNum.map((x) => Math.abs(x)),
          r2: finiteNum,
        }),
      ),
    )
    .map((rows) => rows as DiscordanceModelMetricRow[]);
  return metricRows.chain((metrics) =>
    fc.record({
      run_id: fc.string(),
      schema_version: fc.constant("1"),
      analysis_type: fc.constant("discordance_model" as const),
      case_study_id: fc.string(),
      tissue: fc.constantFrom("SKM-GN", "muscle"),
      contrast_category: fc.constantFrom("EE-CON", "chronic"),
      metrics: fc.constant(metrics),
      weak_prediction_note: fc.constant("Out-of-fold R2 is near zero: this is a weak prediction."),
      limitations: fc.array(fc.string(), { maxLength: 2 }),
      provenance: fc.constant({}),
    }),
  );
}

const ptmArb: fc.Arbitrary<PtmParentResponse> = fc.record({
  run_id: fc.string(),
  schema_version: fc.constant("1"),
  analysis_type: fc.constant("ptm_parent_audit" as const),
  summary_rows: fc.array(
    fc.record({
      tissue: fc.constant("SKM-GN"),
      contrast_category: fc.constant("EE-CON"),
      timepoint: fc.constantFrom("1w", "2w", "4w", "8w"),
      raw_phosphosite_rows: fc.integer({ min: 0, max: 999 }),
      after_mapping_checks_with_uniprot: fc.integer({ min: 0, max: 999 }),
      excluded_missing_or_ambiguous_site_mapping: fc.integer({ min: 0, max: 999 }),
      matched_to_single_gene_uniprot_single_feature_parent: fc.integer({ min: 0, max: 999 }),
      matched_with_both_site_and_parent_ci: fc.integer({ min: 0, max: 999 }),
      supported_phosphosite_rows_among_matched: fc.integer({ min: 0, max: 999 }),
      parent_equivalent_rows_among_matched: fc.integer({ min: 0, max: 999 }),
      supported_phosphosite_and_parent_equivalent_rows: fc.integer({ min: 0, max: 999 }),
      candidate_genes: fc.integer({ min: 0, max: 999 }),
    }),
    { minLength: 1, maxLength: 4 },
  ),
  candidate_count: fc.integer({ min: 0, max: 999 }),
  occupancy_caveat: fc.string(),
  provenance: fc.constant({}),
});

// Empty-but-delivered catalog: either total_events === 0 or all four classes sum to zero.
const emptyCatalogArb: fc.Arbitrary<DiscordanceCatalogResponse> = fc
  .oneof(
    // all classes zero (sum === 0) — total may be anything, the component treats sum 0 as empty
    catalogArb({
      supported_concordant: 0,
      supported_opposite: 0,
      rna_response_protein_equivalent: 0,
      indeterminate: 0,
    }),
  );

// Model variant labels the component renders as group heads.
const MODEL_LABELS: Record<DiscordanceModelMetricRow["model"], string> = {
  zero: "Zero-change baseline",
  rna_only: "RNA-only",
  temporal: "Temporal",
};

const CLASS_LABELS = [
  "Supported concordant",
  "Supported opposite",
  "RNA response · protein equivalent",
  "Indeterminate",
];

afterEach(() => {
  cleanup();
  holder.catalog = { ok: false };
  holder.model = { ok: false };
  holder.ptm = { ok: false };
});

describe("Property 14: result completeness + empty-state exclusivity", () => {
  test("on SUCCESS all four class counts (incl. zero) + each present variant's R²/MAE/RMSE render; empty/error absent", async () => {
    await fc.assert(
      fc.asyncProperty(
        nonEmptyClassCounts.chain((cc) => catalogArb(cc)),
        // 0..3 model variants present; empty subset => model still present but with no metric rows.
        fc.subarray(MODELS, { minLength: 0, maxLength: 3 }).chain((variants) =>
          fc.record({
            variants: fc.constant(variants),
            model: fc.oneof(fc.constant(null), modelArb(variants)),
          }),
        ),
        fc.oneof(fc.constant(null), ptmArb), // ptm resource (null => unavailable)
        async (catalog, modelChoice, ptm) => {
          cleanup();
          const { variants, model } = modelChoice;
          const modelAvailable = model !== null;
          holder.catalog = { ok: true, value: catalog };
          holder.model = modelAvailable ? { ok: true, value: model } : { ok: false };
          holder.ptm = ptm !== null ? { ok: true, value: ptm } : { ok: false };

          render(<Discordance />);

          // Wait for the ready state. The first class label appears in BOTH the Stage-2 chart and
          // the Stage-3 summary cards, so use findAllByText to await render without a multi-match
          // throw.
          await screen.findAllByText(CLASS_LABELS[0]);

          // Empty/error block must be ABSENT on the success path.
          expect(screen.queryByTestId("discordance-empty-or-error")).toBeNull();

          // All four classification classes render, including any that are zero. Each class label
          // sits in a summary-card next to a <strong> count; assert the count text is present in the
          // Stage-3 summary card for that class.
          const cc = catalog.classification_counts;
          const expectedByLabel: Record<string, number> = {
            [CLASS_LABELS[0]]: cc.supported_concordant,
            [CLASS_LABELS[1]]: cc.supported_opposite,
            [CLASS_LABELS[2]]: cc.rna_response_protein_equivalent,
            [CLASS_LABELS[3]]: cc.indeterminate,
          };
          for (const label of CLASS_LABELS) {
            const nodes = screen.getAllByText(label);
            const cardNode = nodes
              .map((n) => n.closest(".summary-card"))
              .find((c): c is Element => c !== null);
            expect(cardNode).toBeTruthy();
            const strong = cardNode!.querySelector("strong");
            expect(strong).not.toBeNull();
            expect(strong!.textContent).toBe(String(expectedByLabel[label]));
          }
          // At least one class is guaranteed zero by the generator.
          expect(Object.values(expectedByLabel).some((v) => v === 0)).toBe(true);

          // For each model variant PRESENT in the payload, its group renders R²/MAE/RMSE headers.
          if (modelAvailable) {
            const presentVariants = new Set(variants);
            for (const variant of MODELS) {
              const head = screen.queryByText(MODEL_LABELS[variant]);
              if (presentVariants.has(variant)) {
                expect(head).not.toBeNull();
                // The variant's table head row carries R²/MAE/RMSE column labels.
                const headRow = head!.closest(".table-row");
                expect(headRow).not.toBeNull();
                const scope = within(headRow as HTMLElement);
                expect(scope.getByText("R²")).toBeTruthy();
                expect(scope.getByText("MAE")).toBeTruthy();
                expect(scope.getByText("RMSE")).toBeTruthy();
              } else {
                // Absent variant: no group head for it.
                expect(head).toBeNull();
              }
            }
          }
        },
      ),
      { numRuns: 100 },
    );
  }, 60000);

  test("on FAILURE or EMPTY the empty/error block shows and NO counts / metrics / PTM render", async () => {
    await fc.assert(
      fc.asyncProperty(
        // failureMode: 0 => catalog fetch rejects; 1 => catalog delivered but empty.
        fc.integer({ min: 0, max: 1 }),
        emptyCatalogArb,
        // model / ptm resources that COULD resolve — but must never render on the empty/error path.
        fc.oneof(fc.constant(null), modelArb(["zero", "rna_only", "temporal"])),
        fc.oneof(fc.constant(null), ptmArb),
        async (failureMode, emptyCatalog, model, ptm) => {
          cleanup();
          if (failureMode === 0) {
            holder.catalog = { ok: false };
          } else {
            holder.catalog = { ok: true, value: emptyCatalog };
          }
          holder.model = model !== null ? { ok: true, value: model } : { ok: false };
          holder.ptm = ptm !== null ? { ok: true, value: ptm } : { ok: false };

          render(<Discordance />);

          // The empty/error block must appear.
          const block = await screen.findByTestId("discordance-empty-or-error");
          expect(block).toBeTruthy();

          // Nothing from the success path renders:
          await waitFor(() => {
            // No classification-class labels.
            for (const label of CLASS_LABELS) {
              expect(screen.queryByText(label)).toBeNull();
            }
            // No model-metrics section header.
            expect(screen.queryByText("Model metrics (R² / MAE / RMSE)")).toBeNull();
            // No PTM-parent divergence summary header (matched via substring).
            expect(screen.queryByText(/PTM-parent divergence summary/)).toBeNull();
          });
        },
      ),
      { numRuns: 100 },
    );
  }, 60000);
});
