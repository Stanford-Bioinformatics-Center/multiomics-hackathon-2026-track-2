// # Feature: finalize-discordance-mvp, Property 10
//
// Property 10 — Mapping-preview confirmation gate (Requirements 4.10, 4.11, 4.12).
//
// QueryBuilder never calls /api/comparisons itself: it calls api.mappingsPreview on
// "Preview mappings", and only after the user clicks the confirm button does it emit
// onRun({ request, preview }). Task 5.2's LiveDashboard is the sole caller of api.runComparison,
// so the gate lives in one place. We therefore model the /api/comparisons gate through the onRun
// spy: onRun firing IS the comparison being requested.
//
// This suite drives ≥100 fast-check iterations. The API client module is mocked so no network is
// touched; api.catalog resolves and api.mappingsPreview either resolves (with a generated set of
// rows, some ambiguous / multi-candidate) or rejects.
//
// The three asserted branches:
//   1. GATE (R4 AC10): before "Preview mappings", and after previewing but BEFORE confirming, the
//      onRun spy has NOT fired. Only after clicking confirm is onRun called exactly once, with the
//      resolved preview attached.
//   2. SURFACING (R4 AC12): every ambiguous / multi-candidate preview row is surfaced in the DOM
//      (its candidates listed) and no candidate is pre-selected.
//   3. FAILURE (R4 AC11): when api.mappingsPreview rejects, onRun never fires, the preview error is
//      shown, and the query selections are retained (the run request the component would build is
//      unchanged from the default, since nothing was toggled).

import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import fc from "fast-check";
import type { Catalog, MappingAuditRow, MappingPreview } from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
// The builder calls api.catalog() on mount and api.mappingsPreview() on "Preview mappings". It must
// NOT call api.runComparison — we spy on it so any accidental call fails the gate. Each fast-check
// iteration installs its own resolved/rejected preview through a mutable holder before rendering.
// ApiError stays a real class (the component does `error instanceof ApiError`).
type PreviewResult = { ok: true; value: MappingPreview } | { ok: false };

// Hoisted so the vi.mock factory (which is itself hoisted to the top of the module) can safely
// reference them without a temporal-dead-zone error.
const { holder, runComparisonSpy } = vi.hoisted(() => ({
  holder: { preview: { ok: false } as PreviewResult },
  runComparisonSpy: vi.fn(),
}));

// A minimal but shape-valid live catalog so the target selectors enable (human + rat contexts with
// the omics/tissue/sex/timepoint values the default query uses). Kept constant across iterations —
// the gate under test is independent of catalog contents.
const CATALOG: Catalog = {
  contexts: [
    {
      species: "rat",
      dataset: "motrpac-pass1b",
      study_design: "chronic",
      study_design_label: "Chronic training",
      tissues: [
        {
          tissue: "SKM-VL",
          tissue_label: "Vastus lateralis",
          layers: ["transcriptomics", "proteomics"],
          omics: ["transcriptomics", "proteomics"],
          sexes: ["female", "male"],
          timepoints: ["8w"],
        },
      ],
    },
    {
      species: "human",
      dataset: "motrpac-human",
      study_design: "acute",
      study_design_label: "Acute exercise",
      tissues: [
        {
          tissue: "SKM-VL",
          tissue_label: "Vastus lateralis",
          layers: ["transcriptomics", "proteomics"],
          omics: ["transcriptomics", "proteomics"],
          sexes: ["female", "male"],
          timepoints: ["8w"],
        },
      ],
    },
  ],
  capability_matrix: {},
  supported_source_species: ["human", "rat", "mouse", "other"],
  unsupported_source_species: [],
  store_hash: "deadbeefcafef00d",
};

vi.mock("../api/client", async () => {
  const actual = await vi.importActual<typeof import("../api/client")>("../api/client");
  return {
    ...actual,
    ApiError: actual.ApiError,
    api: {
      catalog: () => Promise.resolve(CATALOG),
      mappingsPreview: () =>
        holder.preview.ok
          ? Promise.resolve(holder.preview.value)
          : Promise.reject(new actual.ApiError(502, "mapping preview unavailable")),
      // The builder must never invoke this. The spy records any call so the gate assertions catch it.
      runComparison: runComparisonSpy,
    },
  };
});

// Import AFTER vi.mock so the component binds to the mocked client.
import QueryBuilder, { initialAnalysisQuery } from "./QueryBuilder";

// ─────────────────────────── Arbitraries ───────────────────────────

// A single mapping-audit row. `flavour` decides whether it is ambiguous / multi-candidate.
function rowArb(flavour: "single" | "ambiguous" | "multi"): fc.Arbitrary<MappingAuditRow> {
  const idArb = fc.string({ minLength: 1, maxLength: 8 }).map((s) => `GENE_${s.replace(/\s/g, "_")}`);
  const candListArb = fc
    .array(fc.string({ minLength: 1, maxLength: 6 }).map((s) => `C_${s.replace(/[\s,]/g, "_")}`), {
      minLength: 2,
      maxLength: 4,
    })
    .map((cs) => Array.from(new Set(cs)))
    .filter((cs) => cs.length >= 2);
  const singleCandArb = fc.string({ minLength: 1, maxLength: 6 }).map((s) => `C_${s.replace(/[\s,]/g, "_")}`);

  return fc.record({ inputRow: fc.integer({ min: 0, max: 500 }), inputId: idArb }).chain(({ inputRow, inputId }) => {
    if (flavour === "single") {
      return singleCandArb.map<MappingAuditRow>((cand) => ({
        input_row: inputRow,
        input_id: inputId,
        selected_symbol: cand,
        selected_method: "exact",
        status: "mapped",
        n_candidates: 1,
        ambiguous: false,
        candidates: cand,
      }));
    }
    // ambiguous or multi: >1 candidates, joined into the single `candidates` string the UI renders.
    return candListArb.map<MappingAuditRow>((cands) => ({
      input_row: inputRow,
      input_id: inputId,
      // AC12: none auto-selected — the component leaves selection to the user. We model "not
      // auto-selected" as an empty selected_symbol so the DOM shows no chosen candidate.
      selected_symbol: "",
      selected_method: "ambiguous",
      status: "ambiguous",
      n_candidates: cands.length,
      ambiguous: flavour === "ambiguous",
      candidates: cands.join(", "),
    }));
  });
}

// A preview whose rows are a mix — guaranteed to contain at least one ambiguous / multi-candidate
// row, plus optional single-candidate rows. Distinct input ids keep the DOM keys unique.
const previewArb: fc.Arbitrary<MappingPreview> = fc
  .record({
    ambig: fc.array(fc.oneof(rowArb("ambiguous"), rowArb("multi")), { minLength: 1, maxLength: 4 }),
    singles: fc.array(rowArb("single"), { minLength: 0, maxLength: 3 }),
  })
  .map(({ ambig, singles }) => {
    // Re-key every row's input_row/input_id to be unique so React keys `${input_row}-${input_id}`
    // never collide, which would otherwise drop DOM nodes and break the surfacing assertion.
    const rows: MappingAuditRow[] = [...ambig, ...singles].map((row, index) => ({
      ...row,
      input_row: index,
      input_id: `${row.input_id}_${index}`,
    }));
    const total = rows.length;
    const ambiguousCount = rows.filter((r) => r.ambiguous || r.n_candidates > 1).length;
    return {
      counts: { total, ambiguous: ambiguousCount, mapped: total - ambiguousCount },
      rows,
      requires_confirmation: ambiguousCount > 0,
    } satisfies MappingPreview;
  });

// ─────────────────────────── Helpers ───────────────────────────

/** The "Preview mappings" trigger (initial state, before any preview exists). */
function getPreviewButton(): HTMLButtonElement {
  return screen.getByRole("button", { name: /Preview mappings/ });
}

/**
 * Wait until the catalog has loaded and the "Preview mappings" trigger is enabled. jest-dom
 * matchers are not configured in this suite, so we assert on the native `disabled` property.
 */
async function waitForEnabledPreview(): Promise<void> {
  await waitFor(() => {
    const button = getPreviewButton();
    if (button.disabled) throw new Error("preview button still disabled");
  });
}

afterEach(() => {
  cleanup();
  holder.preview = { ok: false };
  runComparisonSpy.mockClear();
});

describe("Property 10: mapping-preview confirmation gate", () => {
  test("comparison is gated: onRun fires exactly once, only after confirm, with the preview attached", async () => {
    await fc.assert(
      fc.asyncProperty(previewArb, async (preview) => {
        cleanup();
        holder.preview = { ok: true, value: preview };
        runComparisonSpy.mockClear();
        const onRun = vi.fn();

        render(<QueryBuilder onQueryChange={() => {}} onRun={onRun} />);

        // Wait until the catalog has loaded and the preview trigger is enabled.
        await waitForEnabledPreview();

        // BEFORE previewing: comparison is not requested.
        expect(onRun).not.toHaveBeenCalled();
        expect(runComparisonSpy).not.toHaveBeenCalled();

        // Click "Preview mappings" — this calls api.mappingsPreview, NOT the comparison.
        fireEvent.click(getPreviewButton());

        // The preview resolves and its confirm button appears.
        const confirmButton = await screen.findByRole("button", {
          name: /(Confirm mappings & run comparison|Run comparison)/,
        });

        // AFTER previewing but BEFORE confirming: still gated.
        expect(onRun).not.toHaveBeenCalled();
        expect(runComparisonSpy).not.toHaveBeenCalled();

        // Click confirm — the gate opens and the builder emits the contract exactly once.
        fireEvent.click(confirmButton);

        expect(onRun).toHaveBeenCalledTimes(1);
        // The builder itself never calls /api/comparisons.
        expect(runComparisonSpy).not.toHaveBeenCalled();

        // The emitted contract carries the confirmed preview and a run request.
        const contract = onRun.mock.calls[0][0] as { request: unknown; preview: MappingPreview };
        expect(contract.preview).toBe(preview);
        expect(contract.request).toBeDefined();
      }),
      { numRuns: 100 },
    );
  }, 60000);

  test("every ambiguous / multi-candidate row is surfaced with all candidates and none pre-selected", async () => {
    await fc.assert(
      fc.asyncProperty(previewArb, async (preview) => {
        cleanup();
        holder.preview = { ok: true, value: preview };
        const onRun = vi.fn();

        render(<QueryBuilder onQueryChange={() => {}} onRun={onRun} />);
        await waitForEnabledPreview();

        fireEvent.click(getPreviewButton());
        // Wait for the confirm button as a signal the preview has rendered.
        await screen.findByRole("button", {
          name: /(Confirm mappings & run comparison|Run comparison)/,
        });

        const ambiguousRows = preview.rows.filter((r) => r.ambiguous || r.n_candidates > 1);
        // The generator guarantees at least one ambiguous / multi-candidate row.
        expect(ambiguousRows.length).toBeGreaterThan(0);

        // The ambiguous section must be present and list every such row's candidates. Each list item
        // renders `${input_id}: ${candidates}` and there is no selected-candidate control (no radio /
        // checkbox is checked), so nothing is auto-selected.
        for (const row of ambiguousRows) {
          // input_id is surfaced.
          const idNodes = screen.getAllByText(row.input_id);
          expect(idNodes.length).toBeGreaterThan(0);
          // Every candidate token in the row's candidate string is present in the DOM.
          const listItem = idNodes
            .map((n) => n.closest("li"))
            .find((li): li is HTMLLIElement => li !== null);
          expect(listItem).toBeTruthy();
          const itemText = listItem!.textContent ?? "";
          for (const candidate of row.candidates.split(",").map((c) => c.trim())) {
            expect(itemText).toContain(candidate);
          }
        }

        // No candidate is auto-selected: there are no checked radios/checkboxes anywhere in the
        // preview, and no element carries an aria-selected / "selected" state for a candidate.
        const checkedInputs = document.querySelectorAll(
          'input[type="radio"]:checked, input[type="checkbox"]:checked',
        );
        expect(checkedInputs.length).toBe(0);
        expect(document.querySelectorAll('[aria-selected="true"]').length).toBe(0);

        // Nothing has run merely from previewing.
        expect(onRun).not.toHaveBeenCalled();
        expect(runComparisonSpy).not.toHaveBeenCalled();
      }),
      { numRuns: 100 },
    );
  }, 60000);

  test("preview failure blocks the comparison, shows an error, and retains selections", async () => {
    await fc.assert(
      fc.asyncProperty(fc.constant(null), async () => {
        cleanup();
        holder.preview = { ok: false }; // api.mappingsPreview rejects
        runComparisonSpy.mockClear();
        const onRun = vi.fn();
        const onQueryChange = vi.fn();

        render(<QueryBuilder onQueryChange={onQueryChange} onRun={onRun} />);
        await waitForEnabledPreview();

        fireEvent.click(getPreviewButton());

        // The mapping-preview error is shown and the comparison is blocked.
        const errorNode = await screen.findByRole("alert");
        expect(errorNode.textContent).toMatch(/Mapping preview unavailable/i);

        // No preview confirm button exists — there is nothing to confirm.
        expect(
          screen.queryByRole("button", {
            name: /(Confirm mappings & run comparison|Run comparison)/,
          }),
        ).toBeNull();

        // The comparison never runs.
        expect(onRun).not.toHaveBeenCalled();
        expect(runComparisonSpy).not.toHaveBeenCalled();

        // Selections are retained: the query was never mutated by the failed preview, so no
        // onQueryChange fired, and the preview trigger is still enabled to retry with the same
        // selections (default query unchanged).
        expect(onQueryChange).not.toHaveBeenCalled();
        await waitForEnabledPreview();
        // Sanity: the default selections that would build the run request are still in place.
        expect(initialAnalysisQuery.selectedOmics.length).toBeGreaterThan(0);
      }),
      { numRuns: 100 },
    );
  }, 60000);
});
