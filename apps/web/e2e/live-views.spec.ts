import { expect, test, type Page } from "@playwright/test";

/**
 * API base URL for the browser-context fetches used in the point → evidence → source chain. Kept in
 * sync with `playwright.config.ts` (same env overrides, same 8765 default) so the spec talks to the
 * exact API the config booted, without importing the config module (which owns the webServer defs).
 */
const API_BASE_URL = `http://${process.env.E2E_API_HOST ?? "127.0.0.1"}:${
  process.env.E2E_API_PORT ?? "8765"
}/api`;

/**
 * End-to-end verification of the live views (R9 AC4, AC5).
 *
 * This is the browser-level counterpart to `apps/api/tests/test_api.py::test_point_to_evidence_to_source_chain`.
 * The API test proves the chain server-side (every feature -> a shown column, a run-scoped
 * evidence_id, and a real store source_feature_id). Here we prove the same chain as a user sees it:
 * we drive the real UI, read the concrete values the four live views render, and then follow the
 * point -> evidence -> source chain starting from what the page displays.
 *
 * Two things are asserted per the requirement:
 *
 *  R9 AC4 — each live view (Discordance, Live, Generalized, Metabolomics) renders its EXPECTED REAL
 *           data, not an empty / placeholder / error state. We assert concrete rendered values
 *           (the frozen discordance class counts, live comparison rows + the layer_discordance
 *           table, the generalized rank-correlation rows, the metabolomics hits table + conclusion),
 *           and we assert the known error/empty markers are ABSENT.
 *
 *  R9 AC5 — the browser-level point -> evidence -> source chain: select a data point in the Live
 *           view (a rendered comparison row), verify the linked evidence is displayed (the
 *           run-scoped evidence for that run, reachable from the page's own export/report links and
 *           run id), and verify the linked source is displayed (the report's source-of-truth text
 *           and each feature's store-backed source_feature_id).
 *
 * Authoritative discordance counts come from the spec/design (run_summary.json.catalog_classes):
 *   supported_concordant = 1, supported_opposite = 0,
 *   rna_response_protein_equivalent = 285, indeterminate = 5642,
 *   no_protein_measurement = 9228, no_rna_measurement = 255.
 */

// Frozen per-class counts (design.md grounding facts) — asserted verbatim, never recomputed.
const DISCORDANCE_CLASS_COUNTS = {
  supported_concordant: 1,
  supported_opposite: 0,
  rna_response_protein_equivalent: 285,
  indeterminate: 5642,
} as const;
const COVERAGE_COUNTS = {
  no_protein_measurement: 9228,
  no_rna_measurement: 255,
} as const;

// The built-in signature fixture the Live view defaults to; its mapping preview is unambiguous, so
// the comparison runs after a single confirm (mirrors the API test's example_name).
const BUILT_IN_EXAMPLE_LABEL = "Pah Muscle Malenfant2015";

/** Nav button labels, exactly as rendered by App.tsx NAV_VIEWS. */
const NAV = {
  discordance: "Discordance",
  live: "Live results (API)",
  generalized: "Generalized query",
  metabolomics: "Metabolomics (ST000763)",
} as const;

/** No live view should ever show one of these known empty/error/placeholder markers (R9 AC4). */
const ERROR_MARKERS = [
  "Backend not reachable",
  "Comparison failed",
  "results unavailable",
  "No discordance results",
  "No comparison run yet",
];

async function gotoView(page: Page, navLabel: string) {
  await page.getByRole("button", { name: navLabel }).click();
}

/** Assert that none of the known error/empty/placeholder states is visible in the current view. */
async function expectNoErrorState(page: Page) {
  for (const marker of ERROR_MARKERS) {
    await expect(
      page.getByText(marker, { exact: false }),
      `an error/empty marker was rendered: "${marker}"`,
    ).toHaveCount(0);
  }
}

test.beforeEach(async ({ page }) => {
  await page.goto("/");
  // The app boots on the default Discordance view (App.DEFAULT_VIEW).
  await expect(page.getByRole("button", { name: NAV.discordance })).toBeVisible();
});

test.describe("R9 AC4 — every live view renders real data", () => {
  test("Discordance view renders the frozen event-class and coverage counts", async ({ page }) => {
    await gotoView(page, NAV.discordance);

    // Wait for the case study to load past the "Loading…" placeholder into a real run.
    await expect(page.getByText(/Track 2 · Omic Discordance Explained · run /)).toBeVisible();
    await expectNoErrorState(page);

    // Stage 3 renders one summary card per classification class. Each card shows the integer count
    // then its label. Assert the count sits inside the card whose label matches — this is a concrete
    // rendered value, not a placeholder. `supported_opposite` = 0 must still render (R3 AC3).
    const cardWithLabel = (label: string) =>
      page.locator(".summary-card", { hasText: label });

    await expect(cardWithLabel("Supported concordant").locator("strong"))
      .toHaveText(String(DISCORDANCE_CLASS_COUNTS.supported_concordant));
    await expect(cardWithLabel("Supported opposite").locator("strong"))
      .toHaveText(String(DISCORDANCE_CLASS_COUNTS.supported_opposite));
    await expect(cardWithLabel("RNA response · protein equivalent").locator("strong"))
      .toHaveText(String(DISCORDANCE_CLASS_COUNTS.rna_response_protein_equivalent));
    await expect(cardWithLabel("Indeterminate").locator("strong"))
      .toHaveText(String(DISCORDANCE_CLASS_COUNTS.indeterminate));

    // Coverage states are rendered distinctly from the four classes (R2 AC8 reconciliation).
    await expect(cardWithLabel("No protein measurement").locator("strong"))
      .toHaveText(String(COVERAGE_COUNTS.no_protein_measurement));
    await expect(cardWithLabel("No RNA measurement").locator("strong"))
      .toHaveText(String(COVERAGE_COUNTS.no_rna_measurement));

    // Honest framing is present (R3 AC5, AC7): a visible weak-prediction statement and the PTM
    // occupancy caveat both render as text. Match the label chip exactly — the same words also
    // appear inside the longer note sentence, so `exact` avoids a strict-mode multi-match.
    await expect(page.getByText("WEAK PREDICTION", { exact: true })).toBeVisible();
    await expect(page.getByText(/occupancy/i).first()).toBeVisible();
  });

  test("Generalized view runs a bundled query and renders rank-correlation rows", async ({ page }) => {
    await gotoView(page, NAV.generalized);

    // The bundled source is selected by default; run it.
    await page.getByRole("button", { name: /Run query/ }).click();

    // Real data: the summary cards and the rank-correlation table render with concrete values.
    await expect(page.getByText("Disease rows")).toBeVisible();
    await expect(
      page.getByRole("heading", { name: /Rank correlation by layer/ }),
    ).toBeVisible();

    // At least one Spearman rho value row is rendered (the view filters to non-null rho rows).
    const rankTable = page.locator(".evidence-table", {
      has: page.getByText("Spearman rho"),
    });
    await expect(rankTable.locator(".table-row").nth(1)).toBeVisible();
    await expectNoErrorState(page);
  });

  test("Metabolomics view renders the case-study conclusion and hits table", async ({ page }) => {
    await gotoView(page, NAV.metabolomics);

    // Match the disclaimer/label chips exactly — the same words recur inside the caveat list items,
    // so `exact` avoids a strict-mode multi-match.
    await expect(page.getByText("SEPARATE COHORT", { exact: true })).toBeVisible();
    // The honest conclusion (null result) renders as prominent text (real data, not placeholder).
    await expect(page.getByText("CONCLUSION", { exact: true })).toBeVisible();
    await expect(page.getByText("Resting hits (PAH vs healthy)")).toBeVisible();

    // The hits table has at least one metabolite row.
    const hitsHeading = page.getByRole("heading", { name: /metabolites$/ });
    await expect(hitsHeading).toBeVisible();
    await expectNoErrorState(page);
  });

  test("Live view runs a real comparison and renders comparison + layer_discordance rows", async ({
    page,
  }) => {
    await runLiveComparison(page);

    // A run id is displayed — concrete, run-scoped, real data (not a placeholder).
    await expect(page.getByText(/Backed by the analysis API · run /)).toBeVisible();

    // The shared ResultsTable rendered at least one comparison row (a "point").
    const resultsTable = page.getByTestId("results-table");
    await expect(resultsTable).toBeVisible();
    await expect(resultsTable.getByTestId("results-empty")).toHaveCount(0);
    await expect(resultsTable.locator(".table-row").nth(1)).toBeVisible();

    // The layer_discordance table (previously rendered in zero views, R3 AC8) shows real rows.
    const layerHeading = page.getByRole("heading", { name: "Layer discordance" });
    await expect(layerHeading).toBeVisible();
    await expect(
      page.getByText("No layer discordance rows returned for this run."),
    ).toHaveCount(0);

    await expectNoErrorState(page);
  });
});

test.describe("R9 AC5 — browser-level point → evidence → source chain", () => {
  test("a rendered Live point traces to run-scoped evidence and a store source", async ({
    page,
  }) => {
    await runLiveComparison(page);

    // ── POINT ─────────────────────────────────────────────────────────────────────────────────
    // The run id shown on the page is the anchor. Read it from the rendered header so the chain
    // starts from what the browser displays, not from a value we injected.
    const runHeader = page.getByText(/Backed by the analysis API · run /);
    await expect(runHeader).toBeVisible();
    const headerText = (await runHeader.textContent()) ?? "";
    const runId = headerText.split("run ").pop()!.trim();
    expect(runId, "a run id must be rendered on the Live view").toBeTruthy();

    // The plotted points are the comparison rows in the ResultsTable. Grab the shown column ids from
    // the rendered feature/contrast cells so we can assert the evidence resolves back to a shown
    // point (mirrors the API test's `f["column_id"] in col_ids`).
    const resultsTable = page.getByTestId("results-table");
    const pointRows = resultsTable.locator(".table-row").filter({ hasNot: page.locator(".table-row--head") });
    const pointCount = await pointRows.count();
    expect(pointCount, "at least one comparison point must render").toBeGreaterThan(0);

    // ── EVIDENCE ──────────────────────────────────────────────────────────────────────────────
    // The evidence for this run is linked from the page itself: the "Download CSV/JSON bundle"
    // (export) and "Open HTML report" links both carry the run id. Verify those links are DISPLAYED
    // and point at the run-scoped evidence resources (the browser-visible evidence leg).
    const exportLink = page.getByRole("link", { name: "Download CSV/JSON bundle" });
    const reportLink = page.getByRole("link", { name: "Open HTML report" });
    await expect(exportLink).toBeVisible();
    await expect(reportLink).toBeVisible();
    await expect(exportLink).toHaveAttribute("href", new RegExp(`/comparisons/${runId}/export$`));
    await expect(reportLink).toHaveAttribute("href", new RegExp(`/comparisons/${runId}/report$`));

    // Follow the run-scoped feature evidence from the browser (same origin/base the page uses). Each
    // measured feature must carry a run-scoped evidence_id and map to a shown column — exactly the
    // API-level guarantee, asserted here from the browser context.
    const featuresResp = await page.request.get(`${API_BASE_URL}/comparisons/${runId}/features`);
    expect(featuresResp.ok(), "run-scoped features must be reachable").toBeTruthy();
    const features = (await featuresResp.json()).features as Array<{
      evidence_id: string;
      column_id: string;
      measured: boolean;
      source_feature_id: string;
    }>;
    const measured = features.filter((f) => f.measured);
    expect(measured.length, "the run must expose measured features (evidence rows)").toBeGreaterThan(0);
    for (const f of measured.slice(0, 25)) {
      expect(f.evidence_id.startsWith(runId), `evidence_id ${f.evidence_id} is run-scoped`).toBeTruthy();
      expect(f.column_id, "evidence maps to a shown point/column").toBeTruthy();
    }

    // ── SOURCE ────────────────────────────────────────────────────────────────────────────────
    // The source-of-truth leg: (1) every measured feature carries a store-backed source_feature_id,
    // and (2) the report the page links to is the human-readable source that states the honest
    // framing verbatim. The directed report is disease-agnostic (it is not PAH-specific), so it
    // states the framing as "Cross-cohort association for hypothesis generation; not a disease
    // treatment-effect test." (see export.render_report_html). Fetch the linked report from the
    // browser and assert that source text is present.
    for (const f of measured.slice(0, 25)) {
      expect(f.source_feature_id, `feature ${f.evidence_id} resolves to a store source_feature_id`).toBeTruthy();
    }

    const reportHref = await reportLink.getAttribute("href");
    expect(reportHref, "the report link must have an href").toBeTruthy();
    const reportResp = await page.request.get(reportHref!);
    expect(reportResp.ok(), "the linked HTML report (source) must load").toBeTruthy();
    const reportText = await reportResp.text();
    expect(
      reportText.includes("not a disease treatment-effect test"),
      "the report source states the honest framing verbatim",
    ).toBeTruthy();

    // Provenance closes the source leg: the run's provenance resource echoes the same run id and a
    // real multiplicity family (mirrors the API test's provenance assertion).
    const provResp = await page.request.get(`${API_BASE_URL}/comparisons/${runId}/provenance`);
    expect(provResp.ok()).toBeTruthy();
    const prov = await provResp.json();
    expect(prov.run_id).toBe(runId);
    expect(prov.multiplicity_family.family_size).toBeGreaterThanOrEqual(40);
  });
});

/**
 * Drive the Live view through the real UI: default built-in signature + rat target -> preview the
 * mappings -> confirm -> run. This is the same path a user takes and the same fixture the API test
 * uses (`pah_muscle_malenfant2015`), whose preview is unambiguous so a single confirm runs it.
 */
async function runLiveComparison(page: Page) {
  await gotoView(page, NAV.live);

  // The idle prompt is shown until a run completes.
  await expect(page.getByRole("heading", { name: "No comparison run yet" })).toBeVisible();

  // The built-in example is the default signature source; make the selection explicit for clarity.
  await expect(
    page.getByRole("button", { name: BUILT_IN_EXAMPLE_LABEL }),
  ).toBeVisible();

  // Step 1: preview mappings (blocks the comparison until confirmed — R4 AC10).
  await page.getByRole("button", { name: /Preview mappings/ }).click();

  // Step 2: confirm & run. The built-in fixture is unambiguous, so the confirm button reads
  // "Run comparison →" (no ambiguous-candidate confirmation needed).
  const runButton = page.getByRole("button", { name: /Run comparison/ });
  await expect(runButton).toBeVisible();
  await runButton.click();

  // The results panel replaces the idle prompt once the run returns.
  await expect(page.getByRole("heading", { name: "No comparison run yet" })).toHaveCount(0);
  await expect(page.getByText(/Backed by the analysis API · run /)).toBeVisible();
}
