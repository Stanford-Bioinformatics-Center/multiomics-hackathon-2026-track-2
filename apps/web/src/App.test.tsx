// App navigation tests + banned-identifier grep — task 8.4.
//
// Two concerns, one file:
//
//  A. Navigation model (R6 AC1, AC6):
//     - NAV_VIEWS has exactly five entries, in order: Discordance, Live, Generalized,
//       Metabolomics, About/Methods (matched by id: discordance, live, generalized,
//       metabolomics, methods).
//     - DEFAULT_VIEW === "discordance" (Discordance is the default view).
//     - Rendering <App /> shows exactly five nav entries and the Discordance view by default.
//
//  B. Banned-identifier grep (R6 AC3):
//     - No navigable component source contains the hardcoded mock gene identifiers
//       PPARGC1A, SOD2, or COL1A1. This reads each navigable component file from disk with
//       node:fs (readFileSync) and asserts none contain the banned identifiers. Documentation
//       (DESIGN_PROVENANCE.md) and *.test.* files are intentionally NOT navigable views and are
//       out of scope for this grep.
//
// Requirements: 6.1, 6.3, 6.6.

import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, render, screen, within } from "@testing-library/react";

import { DEFAULT_VIEW, NAV_VIEWS } from "./App";

// ─────────────────────────── Mock the API client module ───────────────────────────
// Rendering <App /> mounts the default Discordance view, which calls api.discordanceCatalog()
// (and, best-effort, discordanceModel()/discordancePtmParent()). Mock the client so no child view
// makes a real network call. The nav assertions do not depend on child data, so rejecting mocks are
// fine — Discordance simply falls into its error state while the App shell and nav still render.
// ApiError stays a real class because Discordance's errText does `e instanceof ApiError`.
vi.mock("./api/client", async () => {
  const actual = await vi.importActual<typeof import("./api/client")>("./api/client");
  const reject = () => Promise.reject(new actual.ApiError(503, "unavailable in test"));
  return {
    ApiError: actual.ApiError,
    api: {
      discordanceCatalog: reject,
      discordanceModel: reject,
      discordancePtmParent: reject,
    },
  };
});

// Import AFTER vi.mock so App (and the Discordance child it renders) bind to the mocked client.
import App from "./App";

afterEach(() => {
  cleanup();
});

// ─────────────────────────── A. Navigation model ───────────────────────────

describe("App navigation model (R6 AC1, AC6)", () => {
  test("NAV_VIEWS has exactly five entries in the required order (R6 AC6)", () => {
    expect(NAV_VIEWS).toHaveLength(5);
    expect(NAV_VIEWS.map((entry) => entry.id)).toEqual([
      "discordance",
      "live",
      "generalized",
      "metabolomics",
      "methods",
    ]);
  });

  test("DEFAULT_VIEW is the Discordance view (R6 AC1/AC6)", () => {
    expect(DEFAULT_VIEW).toBe("discordance");
    // The default must also be one of the five nav entries.
    expect(NAV_VIEWS.some((entry) => entry.id === DEFAULT_VIEW)).toBe(true);
  });

  test("renders exactly five nav entries and Discordance by default (R6 AC1, AC6)", async () => {
    render(<App />);

    // Query the navigation landmark (the `.view-nav` element carries role=navigation).
    const nav = screen.getByRole("navigation");
    const navButtons = within(nav).getAllByRole("button");
    expect(navButtons).toHaveLength(5);

    // Every nav label from the model renders, in order, inside the nav.
    const navText = navButtons.map((b) => b.textContent ?? "");
    for (const entry of NAV_VIEWS) {
      expect(navText.some((t) => t.includes(entry.label))).toBe(true);
    }

    // Discordance is the default view: its nav entry is the active one on initial load.
    const discordanceLabel = NAV_VIEWS.find((e) => e.id === "discordance")!.label;
    const activeButton = navButtons.find((b) => b.className.includes("view-nav__item--active"));
    expect(activeButton).toBeTruthy();
    expect(activeButton!.textContent).toContain(discordanceLabel);

    // The Discordance view itself is mounted by default. Its API calls are mocked to reject, so it
    // resolves into its empty/error state — confirm a Discordance-specific marker is present rather
    // than another view's content. The error state carries a stable Discordance test id.
    expect(await screen.findByTestId("discordance-empty-or-error")).toBeTruthy();
  });
});

// ─────────────────────────── B. Banned-identifier grep (R6 AC3) ───────────────────────────

// Resolve paths relative to this test file (apps/web/src/App.test.tsx).
const HERE = path.dirname(fileURLToPath(import.meta.url));
const COMPONENTS = path.join(HERE, "components");

// The navigable surface: App.tsx plus every component it renders (directly or transitively) in a
// navigable view. DESIGN_PROVENANCE.md (documentation) and *.test.* files are deliberately excluded.
const NAVIGABLE_SOURCES: string[] = [
  path.join(HERE, "App.tsx"),
  path.join(COMPONENTS, "Discordance.tsx"),
  path.join(COMPONENTS, "LiveDashboard.tsx"),
  path.join(COMPONENTS, "QueryBuilder.tsx"),
  path.join(COMPONENTS, "GeneralizedQuery.tsx"),
  path.join(COMPONENTS, "MetabolomicsCaseStudy.tsx"),
  path.join(COMPONENTS, "ResultsTable.tsx"),
];

const BANNED_IDENTIFIERS = ["PPARGC1A", "SOD2", "COL1A1"] as const;

describe("No hardcoded mock gene identifiers in any navigable view (R6 AC3)", () => {
  test.each(NAVIGABLE_SOURCES)("%s contains no banned mock gene identifier", (source) => {
    const contents = readFileSync(source, "utf8");
    for (const identifier of BANNED_IDENTIFIERS) {
      expect(contents).not.toContain(identifier);
    }
  });

  test("no navigable source contains PPARGC1A, SOD2, or COL1A1 (aggregate)", () => {
    const offenders: string[] = [];
    for (const source of NAVIGABLE_SOURCES) {
      const contents = readFileSync(source, "utf8");
      for (const identifier of BANNED_IDENTIFIERS) {
        if (contents.includes(identifier)) {
          offenders.push(`${path.basename(source)} contains ${identifier}`);
        }
      }
    }
    expect(offenders).toEqual([]);
  });
});
