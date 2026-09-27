// App navigation tests + banned-identifier grep — task 8.4 (updated at the motrpac-explorer merge).
//
// Two concerns, one file:
//
//  A. Navigation model (R6 AC1, AC6; updated post-merge):
//     - The mock design views stay RETIRED (R6). The app presents live views only.
//     - After merging origin/motrpac-explorer the app has SIX live views, in order:
//       Explorer, Live results (API), Discordance, Generalized, Metabolomics, About/Methods
//       (ids: explorer, live, discordance, generalized, metabolomics, methods).
//     - DEFAULT_VIEW === "explorer" — the generalized MoTrPAC explorer is the primary landing
//       view; Live results (LiveDashboard) is the backup live view and Generalized (bundled
//       example signatures) a further backup. (Product decision recorded in DECISIONS.md.)
//     - Rendering <App /> shows exactly six nav entries and the Explorer view by default.
//
//  B. Banned-identifier grep (R6 AC3):
//     - No navigable component source contains the hardcoded mock gene identifiers
//       PPARGC1A, SOD2, or COL1A1 as fabricated result data. This now also covers the merged
//       explorer surface (MotrpacExplorer/liveCharts/ui). Documentation (DESIGN_PROVENANCE.md)
//       and *.test.* files are intentionally NOT navigable views and are out of scope.
//
// Requirements: 6.1, 6.3, 6.6.

import { readFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import { afterEach, describe, expect, test, vi } from "vitest";
import { cleanup, render, screen, within } from "@testing-library/react";

import { DEFAULT_VIEW, NAV_VIEWS } from "./App";

// ─────────────────────────── Mock the API client module ───────────────────────────
// Rendering <App /> mounts the default Explorer view, which calls api.explorerExamples() on mount.
// Mock the client so no child view makes a real network call. The nav assertions do not depend on
// child data, so rejecting/empty mocks are fine — the mounted view simply renders its idle/error
// shell while the App shell and nav still render. ApiError stays a real class because child views'
// error handling does `e instanceof ApiError`.
vi.mock("./api/client", async () => {
  const actual = await vi.importActual<typeof import("./api/client")>("./api/client");
  const reject = () => Promise.reject(new actual.ApiError(503, "unavailable in test"));
  return {
    ApiError: actual.ApiError,
    api: {
      // Explorer (default view) — resolve examples empty so the view mounts cleanly.
      explorerExamples: () => Promise.resolve({ examples: [] }),
      explorerReadFile: reject,
      explorerAnalyse: reject,
      // Discordance (also reachable) — reject so it falls into its empty/error state.
      discordanceCatalog: reject,
      discordanceModel: reject,
      discordancePtmParent: reject,
    },
  };
});

// Import AFTER vi.mock so App (and the children it renders) bind to the mocked client.
import App from "./App";

afterEach(() => {
  cleanup();
});

// ─────────────────────────── A. Navigation model ───────────────────────────

const EXPECTED_ORDER = ["explorer", "live", "discordance", "generalized", "metabolomics", "methods"];

describe("App navigation model (R6 AC1, AC6; post motrpac-explorer merge)", () => {
  test("NAV_VIEWS has exactly six entries in the required order", () => {
    expect(NAV_VIEWS).toHaveLength(6);
    expect(NAV_VIEWS.map((entry) => entry.id)).toEqual(EXPECTED_ORDER);
  });

  test("the retired mock views are not present in the nav (R6)", () => {
    const ids = new Set(NAV_VIEWS.map((e) => e.id));
    for (const retired of ["architecture", "technical", "workflow", "dashboard", "slide"]) {
      expect(ids.has(retired as never)).toBe(false);
    }
  });

  test("DEFAULT_VIEW is the Explorer view", () => {
    expect(DEFAULT_VIEW).toBe("explorer");
    // The default must also be one of the nav entries.
    expect(NAV_VIEWS.some((entry) => entry.id === DEFAULT_VIEW)).toBe(true);
  });

  test("renders exactly six nav entries and Explorer by default", async () => {
    render(<App />);

    // Query the navigation landmark (the `.view-nav` element carries role=navigation).
    const nav = screen.getByRole("navigation");
    const navButtons = within(nav).getAllByRole("button");
    expect(navButtons).toHaveLength(6);

    // Every nav label from the model renders, in order, inside the nav.
    const navText = navButtons.map((b) => b.textContent ?? "");
    for (const entry of NAV_VIEWS) {
      expect(navText.some((t) => t.includes(entry.label))).toBe(true);
    }

    // Explorer is the default view: its nav entry is the active one on initial load.
    const explorerLabel = NAV_VIEWS.find((e) => e.id === "explorer")!.label;
    const activeButton = navButtons.find((b) => b.className.includes("view-nav__item--active"));
    expect(activeButton).toBeTruthy();
    expect(activeButton!.textContent).toContain(explorerLabel);
  });
});

// ─────────────────────────── B. Banned-identifier grep (R6 AC3) ───────────────────────────

// Resolve paths relative to this test file (apps/web/src/App.test.tsx).
const HERE = path.dirname(fileURLToPath(import.meta.url));
const COMPONENTS = path.join(HERE, "components");

// The navigable surface: App.tsx plus every component it renders (directly or transitively) in a
// navigable view — now including the merged explorer surface. DESIGN_PROVENANCE.md (documentation)
// and *.test.* files are deliberately excluded.
const NAVIGABLE_SOURCES: string[] = [
  path.join(HERE, "App.tsx"),
  path.join(COMPONENTS, "MotrpacExplorer.tsx"),
  path.join(COMPONENTS, "liveCharts.tsx"),
  path.join(COMPONENTS, "ui.tsx"),
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
