// # Feature: finalize-discordance-mvp, Property 8
//
// Property 8 — "never a bare value + ortholog link" (Requirements 5.6, 5.7, 5.8, 5.9, 11.5, 12.5).
//
// The shared ResultsTable is the single presentation surface for every view. Its labeling invariant
// is that a value is NEVER rendered bare: every rendered data row carries all three of species,
// dataset, and contrast (R5 AC7/AC8, R11.5 — "show species, dataset, and contrast on every result").
// Additionally, whenever a row maps a query feature to a DIFFERENT-species target (an ortholog),
// the cross-species ortholog-link relation must render for that row (R5 AC9).
//
// This property drives ≥100 fast-check iterations over arbitrary arrays of fully-qualified rows —
// some of which carry an ortholog link with targetSpecies !== querySpecies — and asserts:
//   1. Every rendered data row shows non-empty species, dataset, and contrast cells (no bare value).
//   2. The number of rendered data rows equals the number of input rows (no row silently dropped —
//      cohorts stay whole and distinct, R12.5).
//   3. Every cross-species row (orthologLink.targetSpecies !== orthologLink.querySpecies) renders an
//      ortholog-link relation within that row (R5 AC9).

import { afterEach, describe, expect, test } from "vitest";
import { cleanup, render, screen, within } from "@testing-library/react";
import fc from "fast-check";
import ResultsTable from "./ResultsTable";
import type { OrthologLink, ResultRow } from "./ResultsTable";

// ─────────────────────────── Arbitraries ───────────────────────────

// Non-empty text with no leading/trailing whitespace ambiguity: we constrain to visible content so
// the "non-empty cell" assertions are meaningful (a cell full of whitespace would fail textContent
// trimming, which is exactly the bare-value condition we forbid).
const nonEmptyText = fc
  .string({ minLength: 1, maxLength: 24 })
  .map((s) => s.replace(/\s+/g, " ").trim())
  .filter((s) => s.length > 0);

const species = fc.constantFrom(
  "Homo sapiens",
  "Rattus norvegicus",
  "Mus musculus",
  "Sus scrofa",
);

// A statistic cell may legitimately be null (the engine may not produce a value); that is orthogonal
// to the species/dataset/contrast invariant, which must always hold.
const statValue = fc.option(fc.float({ min: 0, max: 1, noNaN: true, noDefaultInfinity: true }), {
  nil: null,
});

// An ortholog link whose target species is guaranteed to DIFFER from the query species — this is the
// true cross-species relation that must render (R5 AC9).
function crossSpeciesLink(querySpecies: string): fc.Arbitrary<OrthologLink> {
  return fc.record({
    querySymbol: nonEmptyText,
    querySpecies: fc.constant(querySpecies),
    targetSymbol: nonEmptyText,
    targetSpecies: species.filter((s) => s !== querySpecies),
    method: fc.option(fc.constantFrom("HomoloGene", "UniProt", "manual"), { nil: undefined }),
  });
}

// A same-species link — targetSpecies === querySpecies — which does NOT constitute a cross-species
// relation. It still renders the ortholog markup in this component, but our cross-species assertion
// only requires the relation for rows where the species genuinely differ, so we track that.
function sameSpeciesLink(querySpecies: string): fc.Arbitrary<OrthologLink> {
  return fc.record({
    querySymbol: nonEmptyText,
    querySpecies: fc.constant(querySpecies),
    targetSymbol: nonEmptyText,
    targetSpecies: fc.constant(querySpecies),
    method: fc.option(fc.constantFrom("HomoloGene", "UniProt"), { nil: undefined }),
  });
}

const rowArb: fc.Arbitrary<ResultRow> = species.chain((sp) =>
  fc
    .record({
      id: fc.uuid(),
      species: fc.constant(sp),
      dataset: nonEmptyText,
      contrast: nonEmptyText,
      label: nonEmptyText,
      camera_p: statValue,
      camera_fdr: statValue,
      verdict: fc.option(fc.constantFrom("concordant", "opposite", "equivalent"), {
        nil: undefined,
      }),
      significant: fc.option(fc.boolean(), { nil: undefined }),
      // ~1/3 no link, ~1/3 cross-species link, ~1/3 same-species link.
      orthologLink: fc.oneof(
        fc.constant<OrthologLink | undefined>(undefined),
        crossSpeciesLink(sp),
        sameSpeciesLink(sp),
      ),
    })
    .map((r) => r as ResultRow),
);

// Rows must have unique ids (React keys); dedupe by id after generation.
const rowsArb: fc.Arbitrary<ResultRow[]> = fc
  .array(rowArb, { minLength: 1, maxLength: 12 })
  .map((rows) => {
    const seen = new Set<string>();
    return rows.filter((r) => (seen.has(r.id) ? false : (seen.add(r.id), true)));
  });

afterEach(() => {
  cleanup();
});

describe("Property 8: never a bare value + ortholog link", () => {
  test("every rendered row shows species/dataset/contrast; no row dropped; cross-species rows show the ortholog link", () => {
    fc.assert(
      fc.property(rowsArb, (rows) => {
        cleanup();
        render(<ResultsTable rows={rows} />);

        const table = screen.getByTestId("results-table");
        // Data rows are `.table-row` elements EXCLUDING the header row (`.table-row--head`).
        const dataRows = Array.from(
          table.querySelectorAll<HTMLElement>(".table-row:not(.table-row--head)"),
        );

        // (2) No row silently dropped — one rendered data row per input row.
        expect(dataRows.length).toBe(rows.length);

        // Pair each rendered row with its source row (render order preserves input order).
        rows.forEach((sourceRow, i) => {
          const rowEl = dataRows[i];
          const scope = within(rowEl);

          // (1) species + dataset + contrast cells present and non-empty — never a bare value.
          const speciesCell = rowEl.querySelector<HTMLElement>('[data-field="species"]');
          const datasetCell = rowEl.querySelector<HTMLElement>('[data-field="dataset"]');
          const contrastCell = rowEl.querySelector<HTMLElement>('[data-field="contrast"]');

          expect(speciesCell).not.toBeNull();
          expect(datasetCell).not.toBeNull();
          expect(contrastCell).not.toBeNull();

          expect((speciesCell!.textContent ?? "").trim().length).toBeGreaterThan(0);
          expect((datasetCell!.textContent ?? "").trim().length).toBeGreaterThan(0);
          // The contrast cell also hosts the ortholog markup; the contrast text itself is still there.
          expect((contrastCell!.textContent ?? "").trim().length).toBeGreaterThan(0);
          expect(contrastCell!.textContent).toContain(sourceRow.contrast);
          expect(speciesCell!.textContent).toContain(sourceRow.species);
          expect(datasetCell!.textContent).toContain(sourceRow.dataset);

          // (3) Cross-species rows render the ortholog-link relation within that row.
          const link = sourceRow.orthologLink;
          const isCrossSpecies = link !== undefined && link.targetSpecies !== link.querySpecies;
          if (isCrossSpecies) {
            const rel = scope.getByTestId("ortholog-link");
            expect(rel).toBeTruthy();
            expect(rel.textContent).toContain(link!.querySymbol);
            expect(rel.textContent).toContain(link!.targetSymbol);
            expect(rel.textContent).toContain(link!.targetSpecies);
          }
        });

        // The set of cross-species rows in the input equals the count of ortholog-link relations
        // rendered across the table — no cross-species relation is missing and none is fabricated for
        // a row without a differing-species link.
        const expectedCrossSpecies = rows.filter(
          (r) => r.orthologLink !== undefined && r.orthologLink.targetSpecies !== r.orthologLink.querySpecies,
        ).length;
        // Same-species links also render the markup, so the rendered count is >= the cross-species
        // count. What we require is that EVERY cross-species row has one (checked per-row above); here
        // we assert the rendered relations at least cover them.
        const renderedRelations = screen.queryAllByTestId("ortholog-link").length;
        expect(renderedRelations).toBeGreaterThanOrEqual(expectedCrossSpecies);
      }),
      { numRuns: 100 },
    );
  });
});
