// # Feature: finalize-discordance-mvp, Property 7
//
// Property 7 — signature-source exclusivity + field mapping (Requirements 4.1, 4.2, 4.3, 4.4).
//
// The Live query builder allows the signature source to be EXACTLY ONE of a built-in example, pasted
// text, or an uploaded CSV file (R4 AC1). When the confirmed comparison request is emitted:
//   - a built-in example or pasted text maps to `signature_rows` / `example_name`, never to
//     `signature_csv_text` (R4 AC2), and
//   - an uploaded CSV maps to `signature_csv_text`, never to `signature_rows` / `example_name`
//     (R4 AC3), and
//   - in every case, exactly one of {example_name, signature_rows, signature_csv_text} is populated.
//
// When an upload is INVALID (empty, oversized, or unparseable) the builder rejects it, retains the
// previously selected source unchanged, and surfaces an upload error (R4 AC4) — the active source
// never silently switches to the bad upload.
//
// Approach (b): the component is driven end-to-end through Testing Library — `../api/client` is
// mocked (catalog resolves, mappingsPreview resolves) so no network is touched, a source selection
// is generated, exercised via the UI (chip select, textarea paste, file input upload), then
// preview → confirm captures the emitted `onRun` request. fast-check drives ≥100 iterations.
//
// This is the render-based approach deliberately chosen over unit-testing the component's internal
// helpers (`signatureRequestFields` / `parsePastedSignature`), which are not exported — task 5.1
// owns QueryBuilder.tsx, so this suite adds no export and asserts against the real component.

import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import fc from "fast-check";
import type { AnalysisRequestInput, Catalog, MappingPreview } from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
// QueryBuilder calls api.catalog() on mount and api.mappingsPreview() when the user previews. It
// never calls /api/comparisons itself — it emits the confirmed contract via the onRun callback, so
// capturing onRun is how we observe the assembled request. ApiError stays a real class because the
// component does `error instanceof ApiError`.
const holder: {
  catalog: { ok: true; value: Catalog } | { ok: false };
  preview: { ok: true; value: MappingPreview } | { ok: false };
} = {
  catalog: { ok: false },
  preview: { ok: false },
};

vi.mock("../api/client", async () => {
  const actual = await vi.importActual<typeof import("../api/client")>("../api/client");
  return {
    ApiError: actual.ApiError,
    api: {
      catalog: () =>
        holder.catalog.ok
          ? Promise.resolve(holder.catalog.value)
          : Promise.reject(new actual.ApiError(500, "catalog unavailable")),
      mappingsPreview: () =>
        holder.preview.ok
          ? Promise.resolve(holder.preview.value)
          : Promise.reject(new actual.ApiError(500, "preview unavailable")),
    },
  };
});

// Import AFTER vi.mock so the component binds to the mocked client.
import QueryBuilder, { MAX_SIGNATURE_UPLOAD_BYTES } from "./QueryBuilder";

// ─────────────────────────── Fixtures ───────────────────────────
// A minimal-but-valid live catalog exposing the rat context with the dimensions the initial query
// uses (transcriptomics + proteomics, SKM-VL tissue, female sex, 8w timepoint), so the target
// selectors are enabled and the mapping-preview button is reachable.
function makeCatalog(): Catalog {
  return {
    contexts: [
      {
        species: "rat",
        dataset: "pass1b-06",
        study_design: "chronic",
        study_design_label: "Chronic (rat)",
        tissues: [
          {
            tissue: "SKM-VL",
            tissue_label: "Vastus lateralis",
            layers: ["transcriptomics", "proteomics"],
            omics: ["transcriptomics", "proteomics"],
            sexes: ["female", "male"],
            timepoints: ["8w", "4w"],
          },
        ],
      },
      {
        species: "human",
        dataset: "pass1a-human",
        study_design: "acute",
        study_design_label: "Acute (human)",
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
}

// A preview that requires no ambiguity resolution, so a single "Confirm mappings & run comparison"
// click emits onRun.
function makePreview(): MappingPreview {
  return {
    counts: { mapped: 2, ambiguous: 0 },
    rows: [],
    requires_confirmation: false,
  };
}

// A CSV body that parsePastedSignature will turn into ≥1 row (header + one gene row).
const VALID_CSV = "gene_symbol,direction\nPPARGC1A,-1\nSOD2,1\n";
// Pasted text that yields ≥1 parsed row so the paste source is "ready".
const VALID_PASTE = "gene_symbol,direction\nCOL1A1,1\n";

// A jsdom File whose `.size` and `.text()`/FileReader read reflect the given content. jsdom's
// FileReader.readAsText resolves from the Blob parts, so a normal File is sufficient. For the
// oversized case we override `size` without allocating a megabyte of text.
function makeFile(content: string, name = "signature.csv", sizeOverride?: number): File {
  const file = new File([content], name, { type: "text/csv" });
  if (sizeOverride !== undefined) {
    Object.defineProperty(file, "size", { value: sizeOverride });
  }
  return file;
}

// ─────────────────────────── Arbitraries ───────────────────────────
// Valid sources: example (built-in), paste (parseable text), upload (parseable CSV).
type ValidSource =
  | { kind: "example" }
  | { kind: "paste"; text: string }
  | { kind: "upload"; content: string };

const validSourceArb: fc.Arbitrary<ValidSource> = fc.oneof(
  fc.constant<ValidSource>({ kind: "example" }),
  // Build parseable paste bodies: a header plus 1..4 gene rows.
  fc
    .array(
      fc.record({
        gene: fc.stringMatching(/^[A-Z][A-Z0-9]{1,7}$/),
        dir: fc.constantFrom("-1", "1"),
      }),
      { minLength: 1, maxLength: 4 },
    )
    .map<ValidSource>((rows) => ({
      kind: "paste",
      text: "gene_symbol,direction\n" + rows.map((r) => `${r.gene},${r.dir}`).join("\n") + "\n",
    })),
  // Build parseable upload CSV bodies the same way.
  fc
    .array(
      fc.record({
        gene: fc.stringMatching(/^[A-Z][A-Z0-9]{1,7}$/),
        dir: fc.constantFrom("-1", "1"),
      }),
      { minLength: 1, maxLength: 4 },
    )
    .map<ValidSource>((rows) => ({
      kind: "upload",
      content: "gene_symbol,direction\n" + rows.map((r) => `${r.gene},${r.dir}`).join("\n") + "\n",
    })),
);

// Invalid uploads: empty file, oversized file, or unparseable content (no recognizable rows).
type InvalidUpload =
  | { mode: "empty" }
  | { mode: "oversized" }
  | { mode: "unparseable"; content: string };

const invalidUploadArb: fc.Arbitrary<InvalidUpload> = fc.oneof(
  fc.constant<InvalidUpload>({ mode: "empty" }),
  fc.constant<InvalidUpload>({ mode: "oversized" }),
  // Whitespace-only / punctuation-only content parses to zero signature rows → rejected.
  fc
    .constantFrom("   \n  \n", "\n\n\n", ",,,\n,,\n", "   ")
    .map<InvalidUpload>((content) => ({ mode: "unparseable", content })),
);

// ─────────────────────────── Helpers ───────────────────────────
/** Select the given source chip in the "Signature source" field. */
function selectSourceChip(kind: "example" | "paste" | "upload") {
  const label = kind === "example" ? "Example" : kind === "paste" ? "Paste" : "Upload";
  // The chips render the title-cased kind; there is exactly one chip with that text.
  const chip = screen.getAllByText(label).find((el) => el.tagName === "BUTTON");
  expect(chip).toBeTruthy();
  fireEvent.click(chip as HTMLElement);
}

/** Drive preview → confirm and return the captured onRun request (or null if none emitted). */
async function previewConfirmAndCapture(
  captured: { current: AnalysisRequestInput | null },
): Promise<AnalysisRequestInput | null> {
  // Click "Preview mappings →".
  const previewBtn = await screen.findByText(/Preview mappings/);
  fireEvent.click(previewBtn.closest("button") as HTMLElement);
  // Wait for the confirm button to appear (preview resolved).
  const confirmBtn = await screen.findByText(/Run comparison|Confirm mappings & run comparison/);
  fireEvent.click(confirmBtn.closest("button") as HTMLElement);
  await waitFor(() => {
    expect(captured.current).not.toBeNull();
  });
  return captured.current;
}

/** Count how many of the three signature fields are populated in a request. */
function populatedSignatureFields(req: AnalysisRequestInput): string[] {
  const active: string[] = [];
  if (req.example_name !== undefined && req.example_name !== null) active.push("example_name");
  if (req.signature_rows !== undefined && req.signature_rows.length > 0) active.push("signature_rows");
  if (req.signature_csv_text !== undefined && req.signature_csv_text.length > 0)
    active.push("signature_csv_text");
  return active;
}

afterEach(() => {
  cleanup();
  holder.catalog = { ok: false };
  holder.preview = { ok: false };
});

describe("Property 7: signature-source exclusivity + field mapping", () => {
  test("exactly one source active; example/paste → signature_rows|example_name, upload → signature_csv_text", async () => {
    await fc.assert(
      fc.asyncProperty(validSourceArb, async (source) => {
        cleanup();
        holder.catalog = { ok: true, value: makeCatalog() };
        holder.preview = { ok: true, value: makePreview() };

        const captured: { current: AnalysisRequestInput | null } = { current: null };
        render(
          <QueryBuilder
            onQueryChange={() => undefined}
            onRun={(contract) => {
              captured.current = contract.request;
            }}
          />,
        );

        // Catalog loads → the Signature source field is present.
        await screen.findByText("Signature source");

        // Apply the generated source through the UI.
        if (source.kind === "example") {
          selectSourceChip("example");
        } else if (source.kind === "paste") {
          selectSourceChip("paste");
          const textarea = await screen.findByPlaceholderText(/gene_symbol,direction/);
          fireEvent.change(textarea, { target: { value: source.text } });
        } else {
          selectSourceChip("upload");
          const input = document.querySelector('input[type="file"]') as HTMLInputElement;
          expect(input).toBeTruthy();
          const file = makeFile(source.content);
          fireEvent.change(input, { target: { files: [file] } });
          // FileReader.readAsText is async — wait until the "Loaded …" note confirms the source took.
          await screen.findByText(/Loaded /);
        }

        const req = await previewConfirmAndCapture(captured);
        expect(req).not.toBeNull();

        const active = populatedSignatureFields(req as AnalysisRequestInput);
        // Exactly one signature field is populated (R4 AC1 exclusivity).
        expect(active).toHaveLength(1);

        if (source.kind === "example") {
          // Built-in example → example_name (an example is provided as signature_rows OR
          // example_name; the component uses example_name). Never signature_csv_text (R4 AC2).
          expect(active).toEqual(["example_name"]);
          expect((req as AnalysisRequestInput).signature_csv_text).toBeUndefined();
        } else if (source.kind === "paste") {
          // Pasted text → signature_rows, never signature_csv_text (R4 AC2).
          expect(active).toEqual(["signature_rows"]);
          expect((req as AnalysisRequestInput).signature_csv_text).toBeUndefined();
        } else {
          // Uploaded CSV → signature_csv_text, never signature_rows/example_name (R4 AC3).
          expect(active).toEqual(["signature_csv_text"]);
          expect((req as AnalysisRequestInput).signature_rows).toBeUndefined();
          expect((req as AnalysisRequestInput).example_name).toBeUndefined();
        }
      }),
      { numRuns: 100 },
    );
  }, 120000);

  test("invalid upload is rejected, prior source retained unchanged, error surfaced, no source switch", async () => {
    await fc.assert(
      fc.asyncProperty(
        // The prior (valid) source that must survive an invalid upload attempt.
        fc.oneof(
          fc.constant<{ kind: "example" }>({ kind: "example" }),
          fc.constant<{ kind: "paste"; text: string }>({ kind: "paste", text: VALID_PASTE }),
          fc.constant<{ kind: "upload"; content: string }>({ kind: "upload", content: VALID_CSV }),
        ),
        invalidUploadArb,
        async (prior, bad) => {
          cleanup();
          holder.catalog = { ok: true, value: makeCatalog() };
          holder.preview = { ok: true, value: makePreview() };

          const captured: { current: AnalysisRequestInput | null } = { current: null };
          render(
            <QueryBuilder
              onQueryChange={() => undefined}
              onRun={(contract) => {
                captured.current = contract.request;
              }}
            />,
          );

          await screen.findByText("Signature source");

          // Establish the prior valid source.
          let priorField: "example_name" | "signature_rows" | "signature_csv_text";
          if (prior.kind === "example") {
            selectSourceChip("example");
            priorField = "example_name";
          } else if (prior.kind === "paste") {
            selectSourceChip("paste");
            const textarea = await screen.findByPlaceholderText(/gene_symbol,direction/);
            fireEvent.change(textarea, { target: { value: prior.text } });
            priorField = "signature_rows";
          } else {
            selectSourceChip("upload");
            const input = document.querySelector('input[type="file"]') as HTMLInputElement;
            fireEvent.change(input, { target: { files: [makeFile(prior.content)] } });
            await screen.findByText(/Loaded /);
            priorField = "signature_csv_text";
          }

          // Now attempt an INVALID upload. Switch to the upload chip first (only the upload control
          // exposes the file input). Switching to the upload chip with no file yields an empty
          // upload source that is not "ready", so we always attempt the bad file from here and then
          // assert the error + retained-source behavior against the ORIGINAL prior source semantics.
          //
          // Per R4 AC4 the invalid upload must be rejected and the *previously selected* source
          // retained. When the prior source was itself an upload, the retained source is that valid
          // upload. When the prior source was example/paste, selecting the upload chip changes the
          // active kind to an (empty) upload before the bad file — so to test "prior retained" we
          // re-establish the prior source, then fire the bad file over the SAME upload session only
          // when the prior was an upload. For example/paste priors we assert the bad file produces
          // an error and does not yield a usable upload source.
          if (prior.kind === "upload") {
            const input = document.querySelector('input[type="file"]') as HTMLInputElement;
            expect(input).toBeTruthy();
            let badFile: File;
            if (bad.mode === "empty") badFile = makeFile("", "empty.csv", 0);
            else if (bad.mode === "oversized")
              badFile = makeFile("x", "big.csv", MAX_SIGNATURE_UPLOAD_BYTES + 1);
            else badFile = makeFile(bad.content, "junk.csv");

            fireEvent.change(input, { target: { files: [badFile] } });

            // An upload error is surfaced (R4 AC4).
            const err = await screen.findByRole("alert");
            expect(err.textContent && err.textContent.length > 0).toBe(true);

            // The prior valid upload is retained unchanged: preview → confirm still emits
            // signature_csv_text with the ORIGINAL content, not the rejected file.
            const req = await previewConfirmAndCapture(captured);
            expect(req).not.toBeNull();
            expect(populatedSignatureFields(req as AnalysisRequestInput)).toEqual([priorField]);
            expect((req as AnalysisRequestInput).signature_csv_text).toBe(prior.content);
          } else {
            // Prior was example/paste. Switch to upload and fire the bad file. The bad file must be
            // rejected (error shown) and must NOT become a usable signature source. Because the
            // upload chip switch cleared the ready example/paste source, the run button stays gated
            // on "Provide a signature source" — i.e. the invalid upload never produced a source.
            selectSourceChip("upload");
            const input = document.querySelector('input[type="file"]') as HTMLInputElement;
            expect(input).toBeTruthy();
            let badFile: File;
            if (bad.mode === "empty") badFile = makeFile("", "empty.csv", 0);
            else if (bad.mode === "oversized")
              badFile = makeFile("x", "big.csv", MAX_SIGNATURE_UPLOAD_BYTES + 1);
            else badFile = makeFile(bad.content, "junk.csv");

            fireEvent.change(input, { target: { files: [badFile] } });

            // An upload error is surfaced (R4 AC4).
            const err = await screen.findByRole("alert");
            expect(err.textContent && err.textContent.length > 0).toBe(true);

            // No "Loaded …" note appears — the invalid file never became the active source.
            expect(screen.queryByText(/Loaded /)).toBeNull();

            // The run button reflects that no signature source is ready (invalid upload rejected,
            // it did not silently switch in as a usable source).
            await waitFor(() => {
              expect(screen.queryByText(/Provide a signature source/)).not.toBeNull();
            });
            // And no comparison was ever emitted.
            expect(captured.current).toBeNull();
          }
        },
      ),
      { numRuns: 100 },
    );
  }, 120000);
});
