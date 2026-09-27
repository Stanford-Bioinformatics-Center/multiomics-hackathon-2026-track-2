// # Feature: finalize-discordance-mvp, Property 9
//
// Property 9 — catalog-driven selector enablement (Requirements 4.6, 4.7).
//
// QueryBuilder sources its target selectors (species, omic layers, tissue, sex, timepoint) from the
// live `/api/catalog`. The enablement rules under test:
//
//   - PARTIAL DEGRADATION (Req 4.7): when the catalog loads, each per-dimension selector is enabled
//     iff that dimension has ≥1 value for the selected species context, and disabled iff the
//     dimension is empty. An empty dimension disables ONLY its own selector — dimensions that DO
//     carry values stay enabled. This holds independently across tissue / sex / timepoint / omics.
//
//   - TOTAL FAILURE (Req 4.6): when `api.catalog()` rejects, ALL target selectors are disabled and a
//     catalog-unavailable error is surfaced.
//
// The generator produces a Catalog with a single context for the selected species ("rat", the
// component's initial target) whose tissue entries yield arbitrary presence/absence of omics, sexes,
// timepoints (and tissues). fast-check drives ≥100 iterations; the API client module is mocked so no
// network is touched. `api.mappingsPreview` is mocked to resolve so the component mounts cleanly.

import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, render, screen, waitFor } from "@testing-library/react";
import fc from "fast-check";
import type { Catalog, CatalogContext, CatalogTissue, MappingPreview } from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
// QueryBuilder.tsx calls api.catalog() on mount and, later, api.mappingsPreview() behind a button.
// A mutable holder lets each fast-check iteration install its own catalog payload (resolved for the
// partial-degradation cases, rejected for the total-failure case) before rendering. ApiError stays a
// real class because the component does `error instanceof ApiError`.
type CatalogResult = { ok: true; value: Catalog } | { ok: false };

const holder: { catalog: CatalogResult } = { catalog: { ok: false } };

const emptyPreview: MappingPreview = { counts: {}, rows: [], requires_confirmation: false };

vi.mock("../api/client", async () => {
  const actual = await vi.importActual<typeof import("../api/client")>("../api/client");
  return {
    ApiError: actual.ApiError,
    api: {
      catalog: () =>
        holder.catalog.ok
          ? Promise.resolve(holder.catalog.value)
          : Promise.reject(new actual.ApiError(503, "catalog unavailable")),
      // Resolve preview so the component never throws on any incidental interaction.
      mappingsPreview: () => Promise.resolve(emptyPreview),
    },
  };
});

// Import AFTER vi.mock so the component binds to the mocked client.
import QueryBuilder from "./QueryBuilder";

// ─────────────────────────── Arbitraries ───────────────────────────
// The component's initial target species is "rat"; the matching context drives the dimensions.
const TARGET_SPECIES = "rat";

// Catalog omic strings the component maps onto domain OmicLayers (see omicToLayer in QueryBuilder).
// Using recognizable strings keeps the "omics dimension non-empty" branch meaningful.
const OMIC_VALUES = ["transcriptomics", "proteomics", "rna", "protein"] as const;
const TISSUE_VALUES = ["SKM-GN", "SKM-VL", "muscle", "heart", "liver"] as const;
const SEX_VALUES = ["male", "female"] as const;
const TIMEPOINT_VALUES = ["1w", "2w", "4w", "8w", "post_24_hr"] as const;

/**
 * Build a CatalogContext for `rat` whose aggregated dimensions (across its tissue entries) match the
 * requested presence/absence per dimension. Because `dimensionsForContext` unions every tissue's
 * omics/sexes/timepoints and collects each tissue's name, we can control each aggregate dimension by
 * putting the chosen values on a single tissue entry.
 *
 * Note the coupling: an empty `tissues` array forces ALL of tissues/omics/sexes/timepoints empty
 * (there is nothing to aggregate over). So the "tissues non-empty but omics empty" style cases are
 * expressed by keeping ≥1 tissue entry whose omics list is empty.
 */
function contextArb(): fc.Arbitrary<{
  context: CatalogContext;
  present: { tissues: boolean; omics: boolean; sexes: boolean; timepoints: boolean };
}> {
  return fc
    .record({
      hasTissue: fc.boolean(),
      omics: fc.subarray([...OMIC_VALUES]),
      sexes: fc.subarray([...SEX_VALUES]),
      timepoints: fc.subarray([...TIMEPOINT_VALUES]),
      tissueName: fc.constantFrom(...TISSUE_VALUES),
    })
    .map(({ hasTissue, omics, sexes, timepoints, tissueName }) => {
      const tissues: CatalogTissue[] = hasTissue
        ? [
            {
              tissue: tissueName,
              tissue_label: tissueName,
              layers: [...omics],
              omics: [...omics],
              sexes: [...sexes],
              timepoints: [...timepoints],
            },
          ]
        : [];
      const context: CatalogContext = {
        species: TARGET_SPECIES,
        dataset: "motrpac-rat",
        study_design: "chronic",
        study_design_label: "Chronic (rat)",
        tissues,
      };
      // With no tissue entry, every aggregate dimension is empty.
      const present = hasTissue
        ? {
            tissues: true,
            omics: omics.length > 0,
            sexes: sexes.length > 0,
            timepoints: timepoints.length > 0,
          }
        : { tissues: false, omics: false, sexes: false, timepoints: false };
      return { context, present };
    });
}

function catalogArb(context: CatalogContext): fc.Arbitrary<Catalog> {
  return fc.record({
    contexts: fc.constant([context]),
    capability_matrix: fc.constant({}),
    // Keep the source-species set broad so it never interferes; it is a separate selector.
    supported_source_species: fc.constant(["human", "rat", "mouse", "other"]),
    unsupported_source_species: fc.constant([]),
    store_hash: fc.hexaString({ minLength: 8, maxLength: 12 }),
  });
}

// ─────────────────────────── DOM query helpers ───────────────────────────
// Each Field renders a ".query-field" wrapper with a ".query-field__label" span, then the option
// buttons inside ".query-options". We locate a dimension's option buttons by its label text and
// inspect their `disabled` attribute.
function fieldByLabel(label: string): HTMLElement | null {
  const spans = Array.from(document.querySelectorAll(".query-field__label"));
  // Match on a stable prefix so label decorations (unit hints, the "· multi-select" suffix, which
  // QueryBuilder renders with literal backslash escapes) don't break lookup.
  const match = spans.find((s) => (s.textContent ?? "").startsWith(label));
  return match ? (match.closest(".query-field") as HTMLElement | null) : null;
}

function optionButtons(label: string): HTMLButtonElement[] {
  const field = fieldByLabel(label);
  if (!field) return [];
  const options = field.querySelector(".query-options");
  if (!options) return [];
  return Array.from(options.querySelectorAll("button"));
}

// Labels exactly as rendered by QueryBuilder's <Field label=...> props.
const TISSUE_LABEL = "Tissue";
const SEX_LABEL = "Sex";
const TIMEPOINT_LABEL = "Time point";
// Prefix only — QueryBuilder renders this field label as "Omic layers \u00b7 multi-select" with
// literal backslash escapes (JSX string attributes do not process \u escapes), so match the prefix.
const OMIC_LABEL = "Omic layers";

afterEach(() => {
  cleanup();
  holder.catalog = { ok: false };
});

describe("Property 9: catalog-driven selector enablement", () => {
  test("each per-dimension selector is enabled iff its dimension has ≥1 value; an empty dimension disables only its own selector", async () => {
    await fc.assert(
      fc.asyncProperty(
        contextArb().chain(({ context, present }) =>
          catalogArb(context).map((catalog) => ({ catalog, present })),
        ),
        async ({ catalog, present }) => {
          cleanup();
          holder.catalog = { ok: true, value: catalog };

          render(<QueryBuilder onQueryChange={() => {}} />);

          // Wait for the catalog to load: the total-failure banner must be absent and the derived
          // study-context block (always rendered once past loading) must be present.
          await waitFor(() => {
            expect(screen.queryByText(/analysis catalog is unavailable/i)).toBeNull();
          });

          // ── Tissue selector ──
          const tissueButtons = optionButtons(TISSUE_LABEL);
          expect(tissueButtons.length).toBeGreaterThan(0);
          if (present.tissues) {
            // Dimension has values → its buttons are enabled.
            expect(tissueButtons.every((b) => !b.disabled)).toBe(true);
          } else {
            // Empty dimension → disabled.
            expect(tissueButtons.every((b) => b.disabled)).toBe(true);
          }

          // ── Sex selector ──
          const sexButtons = optionButtons(SEX_LABEL);
          expect(sexButtons.length).toBeGreaterThan(0);
          if (present.sexes) {
            expect(sexButtons.every((b) => !b.disabled)).toBe(true);
          } else {
            expect(sexButtons.every((b) => b.disabled)).toBe(true);
          }

          // ── Time point selector ──
          const tpButtons = optionButtons(TIMEPOINT_LABEL);
          expect(tpButtons.length).toBeGreaterThan(0);
          if (present.timepoints) {
            expect(tpButtons.every((b) => !b.disabled)).toBe(true);
          } else {
            expect(tpButtons.every((b) => b.disabled)).toBe(true);
          }

          // ── Omic layers selector ──
          // The omic buttons enumerate the fixed capability list; each is disabled when the catalog
          // does not expose that layer for the context. When the omics dimension is EMPTY, none of
          // the layers are exposed, so every omic button is disabled. When the omics dimension has
          // values, at least one omic button is enabled (the recognizable strings map to a layer).
          const omicButtons = optionButtons(OMIC_LABEL);
          expect(omicButtons.length).toBeGreaterThan(0);
          if (present.omics) {
            expect(omicButtons.some((b) => !b.disabled)).toBe(true);
          } else {
            expect(omicButtons.every((b) => b.disabled)).toBe(true);
          }

          // ── "Empty dimension disables only its own selector" (Req 4.7) ──
          // For every dimension that DOES carry values, its selector stays enabled regardless of
          // whether OTHER dimensions are empty. Assert this cross-dimension independence directly.
          if (present.tissues) expect(tissueButtons.some((b) => !b.disabled)).toBe(true);
          if (present.sexes) expect(sexButtons.some((b) => !b.disabled)).toBe(true);
          if (present.timepoints) expect(tpButtons.some((b) => !b.disabled)).toBe(true);
          if (present.omics) expect(omicButtons.some((b) => !b.disabled)).toBe(true);

          cleanup();
        },
      ),
      { numRuns: 100 },
    );
  }, 60000);

  test("total catalog failure disables ALL target selectors and surfaces a catalog-unavailable error", async () => {
    await fc.assert(
      fc.asyncProperty(fc.constant(null), async () => {
        cleanup();
        // api.catalog() rejects → total failure (Req 4.6).
        holder.catalog = { ok: false };

        render(<QueryBuilder onQueryChange={() => {}} />);

        // Catalog-unavailable error is surfaced.
        await screen.findByText(/analysis catalog is unavailable/i);

        // Every target selector across all dimensions is disabled. Species buttons live in the
        // "Target MoTrPAC species" field; omics/tissue/sex/timepoint in their own fields.
        const speciesButtons = optionButtons("Target MoTrPAC species");
        const tissueButtons = optionButtons(TISSUE_LABEL);
        const sexButtons = optionButtons(SEX_LABEL);
        const tpButtons = optionButtons(TIMEPOINT_LABEL);
        const omicButtons = optionButtons(OMIC_LABEL);

        for (const group of [speciesButtons, tissueButtons, sexButtons, tpButtons, omicButtons]) {
          expect(group.length).toBeGreaterThan(0);
          expect(group.every((b) => b.disabled)).toBe(true);
        }

        cleanup();
      }),
      { numRuns: 100 },
    );
  }, 60000);
});
