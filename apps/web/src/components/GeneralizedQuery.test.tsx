// GeneralizedQuery component tests (Vitest + Testing Library) — task 7.5.
//
// These example-based tests render <GeneralizedQuery /> against a mocked API client and assert the
// R5 AC3/AC4 behaviors for the paste + upload controls added alongside the bundled picklist:
//   1. All three source modes (bundled / paste / upload) are offered.
//   2. A pasted signature that cannot be parsed (client-side pre-check) is rejected with a visible
//      parse-failure indication and the analyst's pasted text is NOT cleared (R5 AC4).
//   3. A backend 422 (naming the failing schema requirement) is surfaced to the user, and the
//      pasted text is still NOT cleared (R5 AC2/AC4).
//   4. A conforming paste submits via generalizedQueryUpload() and renders into the shared result
//      panel (summary cards + interpretation) — same panel as the bundled path.
//   5. The bundled picklist still runs via the GET route into the same result panel.
//
// Requirements: 5.3, 5.4.

import { afterEach, beforeEach, describe, expect, test, vi } from "vitest";
import { cleanup, fireEvent, render, screen, waitFor } from "@testing-library/react";
import type {
  GeneralizedQueryResult,
  GeneralizedQueryUploadInput,
  GeneralizedSignature,
} from "../api/client";

// ─────────────────────────── Mock the API client module ───────────────────────────
type UploadResult = { ok: true; value: GeneralizedQueryResult } | { ok: false; status: number; detail: string };
type BundledResult = { ok: true; value: GeneralizedQueryResult } | { ok: false; status: number; detail: string };

const holder: {
  signatures: GeneralizedSignature[];
  upload: UploadResult;
  bundled: BundledResult;
  lastUploadBody: GeneralizedQueryUploadInput | null;
} = {
  signatures: [],
  upload: { ok: false, status: 422, detail: "" },
  bundled: { ok: false, status: 500, detail: "" },
  lastUploadBody: null,
};

vi.mock("../api/client", async () => {
  const actual = await vi.importActual<typeof import("../api/client")>("../api/client");
  return {
    ApiError: actual.ApiError,
    api: {
      generalizedSignatures: () =>
        Promise.resolve({ signatures: holder.signatures, note: "test" }),
      generalizedQuery: (_id: string) =>
        holder.bundled.ok
          ? Promise.resolve(holder.bundled.value)
          : Promise.reject(new actual.ApiError(holder.bundled.status, holder.bundled.detail)),
      generalizedQueryUpload: (body: GeneralizedQueryUploadInput) => {
        holder.lastUploadBody = body;
        return holder.upload.ok
          ? Promise.resolve(holder.upload.value)
          : Promise.reject(new actual.ApiError(holder.upload.status, holder.upload.detail));
      },
    },
  };
});

import GeneralizedQuery from "./GeneralizedQuery";

// ─────────────────────────── Fixtures ───────────────────────────
const SIGS: GeneralizedSignature[] = [
  { id: "pah_blood_rna", label: "PAH blood RNA", tissue: "blood", reference: "motrpac_blood", reference_contrast_category: "EE-CON" },
];

function makeResult(id: string): GeneralizedQueryResult {
  return {
    analysis_type: "generalized_query",
    schema_version: "1",
    signature_id: id,
    signature_label: id,
    tissue: "muscle",
    reference: "motrpac_muscle",
    counts: {
      disease_rows: 42,
      matched_candidate_rows: 30,
      disease_rows_unmatched_everywhere: 5,
      ambiguous_candidate_rows: 7,
    },
    filters: {},
    thresholds: {},
    interpretation: "Descriptive rank correlation only; no causal claim.",
    rank_correlation: [],
    coverage_rows: 30,
    provenance: { disease_sha256: "abcdef0123456789", reference_sha256: "0123456789abcdef" },
  };
}

const VALID_CSV = "feature_id,log2fc,q_value\nGENE1,1.2,0.01\nGENE2,-0.5,0.04\n";

beforeEach(() => {
  holder.signatures = SIGS;
  holder.upload = { ok: false, status: 422, detail: "" };
  holder.bundled = { ok: false, status: 500, detail: "" };
  holder.lastUploadBody = null;
});

afterEach(() => {
  cleanup();
});

async function switchTo(tabLabel: string) {
  const tab = await screen.findByRole("tab", { name: tabLabel });
  fireEvent.click(tab);
  return tab;
}

describe("GeneralizedQuery — upload + paste controls (R5 AC3/AC4)", () => {
  test("offers all three signature source modes", async () => {
    render(<GeneralizedQuery />);
    expect(await screen.findByRole("tab", { name: "Bundled example" })).toBeTruthy();
    expect(screen.getByRole("tab", { name: "Paste CSV" })).toBeTruthy();
    expect(screen.getByRole("tab", { name: "Upload CSV" })).toBeTruthy();
  });

  test("rejects an unparseable pasted signature with a visible parse-failure indication and does NOT clear the pasted text (R5 AC4)", async () => {
    render(<GeneralizedQuery />);
    await switchTo("Paste CSV");

    const textarea = screen.getByLabelText(/Paste query_core-format signature CSV/) as HTMLTextAreaElement;
    // Header only, no data row → client-side pre-check fails.
    fireEvent.change(textarea, { target: { value: "just-a-header-no-rows" } });
    fireEvent.click(screen.getByRole("button", { name: /Run query/ }));

    // Parse-failure indication is visible.
    const alert = await screen.findByRole("alert");
    expect(alert.textContent).toMatch(/[Pp]arse failure/);

    // The upload endpoint must NOT have been called (rejected before running).
    expect(holder.lastUploadBody).toBeNull();

    // The analyst's pasted text is preserved.
    expect((screen.getByLabelText(/Paste query_core-format signature CSV/) as HTMLTextAreaElement).value)
      .toBe("just-a-header-no-rows");
  });

  test("surfaces the backend 422 (naming the failing requirement) and keeps the pasted text (R5 AC2/AC4)", async () => {
    holder.upload = { ok: false, status: 422, detail: "signature_csv_text: missing required column 'log2fc'" };
    render(<GeneralizedQuery />);
    await switchTo("Paste CSV");

    const textarea = screen.getByLabelText(/Paste query_core-format signature CSV/) as HTMLTextAreaElement;
    fireEvent.change(textarea, { target: { value: VALID_CSV } });
    fireEvent.click(screen.getByRole("button", { name: /Run query/ }));

    const alert = await screen.findByRole("alert");
    expect(alert.textContent).toContain("422");
    expect(alert.textContent).toMatch(/log2fc/);

    // The endpoint was actually called (passes the client pre-check), with the pasted text.
    await waitFor(() => expect(holder.lastUploadBody).not.toBeNull());
    expect(holder.lastUploadBody?.signature_csv_text).toBe(VALID_CSV);

    // The pasted text is NOT cleared after the 422.
    expect((screen.getByLabelText(/Paste query_core-format signature CSV/) as HTMLTextAreaElement).value)
      .toBe(VALID_CSV);
  });

  test("a conforming paste submits via generalizedQueryUpload() and renders in the shared result panel", async () => {
    holder.upload = { ok: true, value: makeResult("pasted") };
    render(<GeneralizedQuery />);
    await switchTo("Paste CSV");

    fireEvent.change(screen.getByLabelText(/Paste query_core-format signature CSV/), {
      target: { value: VALID_CSV },
    });
    fireEvent.click(screen.getByRole("button", { name: /Run query/ }));

    // Shared result panel renders (summary card + interpretation).
    expect(await screen.findByText("Descriptive rank correlation only; no causal claim.")).toBeTruthy();
    expect(screen.getByText("Disease rows")).toBeTruthy();

    // The submitted body carried the required query_core fields with sensible defaults.
    expect(holder.lastUploadBody?.tissue).toBeTruthy();
    expect(holder.lastUploadBody?.reference_contrast_category).toBe("EE-CON");
  });

  test("the bundled picklist still runs into the same result panel", async () => {
    holder.bundled = { ok: true, value: makeResult("pah_blood_rna") };
    render(<GeneralizedQuery />);

    // Default mode is bundled; pick the example and run.
    fireEvent.click(await screen.findByRole("button", { name: "PAH blood RNA" }));
    fireEvent.click(screen.getByRole("button", { name: /Run query/ }));

    expect(await screen.findByText("Descriptive rank correlation only; no causal claim.")).toBeTruthy();
    expect(screen.getByText("Disease rows")).toBeTruthy();
    // Upload endpoint untouched by the bundled path.
    expect(holder.lastUploadBody).toBeNull();
  });
});
