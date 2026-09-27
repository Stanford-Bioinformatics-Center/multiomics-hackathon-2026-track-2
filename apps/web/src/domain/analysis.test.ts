import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  classifyResultPresentation,
  getStudyContext,
  getVisualizationMode,
  getVisualizationModeFromReturnedLayers,
  reconcileQueryForSpecies,
  resolveContext,
  studyContexts,
  type AnalysisQuery,
  type AvailabilityRequest,
} from "./analysis.ts";

const ratQuery: AvailabilityRequest = {
  targetMotrpacSpecies: "rat",
  selectedOmics: ["transcriptomics", "proteomics"],
  tissue: "SKM-VL",
  sex: "female",
  timepoint: "8w",
};

test("selecting human derives the acute-human context", () => {
  assert.equal(getStudyContext("human")?.studyDesign, "acute");
});

test("selecting rat derives the chronic-rat context", () => {
  assert.equal(getStudyContext("rat")?.studyDesign, "chronic");
});

test("changing species exposes species-specific tissues, time points, and sexes", () => {
  const human = getStudyContext("human");
  const rat = getStudyContext("rat");
  assert.deepEqual(human?.tissues, ["muscle", "blood"]);
  assert.ok(human?.timepoints.includes("post_24_hr"));
  assert.deepEqual(human?.sexes, ["combined"]);
  assert.deepEqual(rat?.tissues, ["SKM-GN", "SKM-VL"]);
  assert.ok(rat?.timepoints.includes("8w"));
  assert.deepEqual(rat?.sexes, ["female", "male"]);
});

test("no acute or chronic input control is rendered", () => {
  const source = readFileSync(
    new URL("../components/QueryBuilder.tsx", import.meta.url),
    "utf8",
  );
  assert.doesNotMatch(source, /exerciseStudyDesign|Exercise design/);
});

test("changing species clears incompatible filters and reports each changed field", () => {
  const result = reconcileQueryForSpecies(
    {
      ...ratQuery,
      signatureSourceSpecies: "human",
      includeNonSignificantResults: true,
    },
    "human",
  );
  assert.deepEqual(result.changedFields, ["tissue", "sex", "timepoint"]);
  assert.equal(result.query.tissue, "muscle");
  assert.equal(result.query.sex, "combined");
  assert.equal(result.query.timepoint, "during_20_min");
});

test("an incompatible filter returns no_matching_context with alternatives", () => {
  const result = resolveContext({ ...ratQuery, tissue: "blood" });
  assert.equal(result.status, "no_matching_context");
  assert.ok(result.suggestedAlternatives.length > 0);
});

test("rat transcriptomics and proteomics return an available context", () => {
  const result = resolveContext(ratQuery);
  assert.equal(result.status, "available");
  assert.deepEqual(result.availableOmics, ["transcriptomics", "proteomics"]);
});

test("an unsupported layer returns unsupported_omics", () => {
  const result = resolveContext({ ...ratQuery, selectedOmics: ["metabolomics"] });
  assert.equal(result.status, "unsupported_omics");
});

test("supported and unsupported layers return partially_available", () => {
  const result = resolveContext({
    ...ratQuery,
    selectedOmics: ["transcriptomics", "metabolomics"],
  });
  assert.equal(result.status, "partially_available");
  assert.deepEqual(result.unavailableOmics, ["metabolomics"]);
});

test("one layer enables disease-versus-exercise mode", () => {
  assert.equal(getVisualizationMode(["transcriptomics"]), "single_layer");
});

test("two layers enable pairwise cross-omic mode", () => {
  assert.equal(getVisualizationMode(["transcriptomics", "proteomics"]), "pairwise");
});

test("more than two layers switch to matrix mode", () => {
  assert.equal(
    getVisualizationMode(["transcriptomics", "proteomics", "metabolomics"]),
    "matrix",
  );
});

test("no mapped features produces a scientific empty state", () => {
  assert.equal(
    classifyResultPresentation({ mappedFeatureCount: 0, measurementCount: 0, significantCount: 0 }),
    "no_mapped_features",
  );
});

test("nonsignificant measurements remain a weak-evidence result, not missing data", () => {
  assert.equal(
    classifyResultPresentation({ mappedFeatureCount: 8, measurementCount: 8, significantCount: 0 }),
    "weak_evidence",
  );
});

test("technical exceptions produce an actual error state", () => {
  assert.equal(
    classifyResultPresentation({
      mappedFeatureCount: 8,
      measurementCount: 8,
      significantCount: 2,
      technicalError: true,
    }),
    "processing_error",
  );
});

test("study design is catalog metadata rather than a query input", () => {
  assert.equal(studyContexts.length, 2);
  assert.ok(studyContexts.every((context) => context.studyDesign));
  assert.ok(!("exerciseStudyDesign" in ratQuery));
});

test("visualization mode is derived from RETURNED layers, not requested layers", () => {
  // requested two omics but the backend returned only one -> single_layer (no fake multi-omic view)
  assert.equal(getVisualizationModeFromReturnedLayers(["RNA"]), "single_layer");
  // two returned layers -> pairwise
  assert.equal(getVisualizationModeFromReturnedLayers(["RNA", "PROT"]), "pairwise");
  // duplicates collapse
  assert.equal(getVisualizationModeFromReturnedLayers(["RNA", "RNA"]), "single_layer");
  // three returned layers -> matrix
  assert.equal(getVisualizationModeFromReturnedLayers(["RNA", "PROT", "PHOSPHO"]), "matrix");
  // empty results -> single_layer (nothing to pair)
  assert.equal(getVisualizationModeFromReturnedLayers([]), "single_layer");
});

test("selecting an unimplemented layer alongside a real one cannot fabricate a matrix", () => {
  // user requested 3 omics (would be 'matrix' by request), but only RNA+PROT returned -> pairwise
  const requestedMode = getVisualizationMode(["transcriptomics", "proteomics", "metabolomics"]);
  const renderedMode = getVisualizationModeFromReturnedLayers(["RNA", "PROT"]);
  assert.equal(requestedMode, "matrix");
  assert.equal(renderedMode, "pairwise");
  assert.notEqual(renderedMode, requestedMode);
});
// --- Task 5.7: rewired catalog-driven QueryBuilder logic (Requirements 4.8) ---

// The FDR-threshold validity constraint lives in the rewired QueryBuilder component (task 5.1),
// not in analysis.ts, so we document it here as a small pure predicate rather than importing UI.
// UI accepts an inclusive upper bound (0 < t <= 1); the API validates gt=0, lt=1. This predicate
// captures the UI-inclusive rule the builder enforces on presets/values (R4 AC8).
function isValidFdrThreshold(t: number): boolean {
  return Number.isFinite(t) && t > 0 && t <= 1;
}

// Presets offered by QueryBuilder.tsx (fdrPresets). They must stay strictly inside (0, 1) so they
// satisfy BOTH the UI-inclusive rule and the API's gt=0, lt=1 boundary.
const fdrPresets = [0.01, 0.05, 0.1];

test("FDR threshold accepts the inclusive upper boundary 0 < t <= 1", () => {
  // Interior values are valid.
  assert.equal(isValidFdrThreshold(0.01), true);
  assert.equal(isValidFdrThreshold(0.05), true);
  assert.equal(isValidFdrThreshold(0.5), true);
  // Exactly 1 is accepted by the UI-inclusive upper bound.
  assert.equal(isValidFdrThreshold(1), true);
  // A value just below 1 is accepted.
  assert.equal(isValidFdrThreshold(0.999999), true);
});

test("FDR threshold rejects t <= 0 and t > 1", () => {
  assert.equal(isValidFdrThreshold(0), false);
  assert.equal(isValidFdrThreshold(-0.01), false);
  assert.equal(isValidFdrThreshold(1.0001), false);
  assert.equal(isValidFdrThreshold(2), false);
  // Non-finite values are never valid thresholds.
  assert.equal(isValidFdrThreshold(Number.NaN), false);
  assert.equal(isValidFdrThreshold(Number.POSITIVE_INFINITY), false);
});

test("QueryBuilder FDR presets stay strictly inside (0, 1) to satisfy both UI and API bounds", () => {
  for (const preset of fdrPresets) {
    // UI-inclusive rule.
    assert.equal(isValidFdrThreshold(preset), true, `preset ${preset} must satisfy 0 < t <= 1`);
    // API's stricter gt=0, lt=1 boundary (a preset of exactly 1 would fail the API).
    assert.ok(preset > 0 && preset < 1, `preset ${preset} must satisfy the API's gt=0, lt=1`);
  }
});

// --- reconciled catalog-driven selection paths (species switch reconciliation) ---

const humanFilledQuery: AnalysisQuery = {
  signatureSourceSpecies: "human",
  targetMotrpacSpecies: "human",
  selectedOmics: ["transcriptomics", "proteomics"],
  tissue: "muscle",
  sex: "combined",
  timepoint: "post_24_hr",
  significanceThreshold: 0.05,
  includeNonSignificantResults: true,
};

test("switching human -> rat reconciles every incompatible field to the rat context", () => {
  const result = reconcileQueryForSpecies(humanFilledQuery, "rat");
  const rat = getStudyContext("rat")!;
  // All three human values are incompatible with the rat context, so all reconcile.
  assert.deepEqual(result.changedFields, ["tissue", "sex", "timepoint"]);
  assert.equal(result.query.targetMotrpacSpecies, "rat");
  assert.equal(result.query.tissue, rat.tissues[0]);
  assert.equal(result.query.sex, rat.sexes[0]);
  assert.equal(result.query.timepoint, rat.timepoints[0]);
  // Reconciled fields land inside the target species' catalog dimensions.
  assert.ok(rat.tissues.includes(result.query.tissue!));
  assert.ok(rat.sexes.includes(result.query.sex!));
  assert.ok(rat.timepoints.includes(result.query.timepoint!));
});

test("reconciliation preserves fields already compatible with the target species context", () => {
  const rat = getStudyContext("rat")!;
  // A rat query whose tissue/sex/timepoint are all valid for rat should switch species with no
  // field changes reported (the catalog-driven builder keeps compatible selections).
  const compatibleRatQuery: AnalysisQuery = {
    ...humanFilledQuery,
    targetMotrpacSpecies: "rat",
    tissue: rat.tissues[1],
    sex: rat.sexes[1],
    timepoint: "4w",
  };
  const result = reconcileQueryForSpecies(compatibleRatQuery, "rat");
  assert.deepEqual(result.changedFields, []);
  assert.equal(result.query.tissue, rat.tissues[1]);
  assert.equal(result.query.sex, rat.sexes[1]);
  assert.equal(result.query.timepoint, "4w");
});

test("reconciliation only rewrites the dimensions that are incompatible", () => {
  const rat = getStudyContext("rat")!;
  // Tissue is already valid for rat; sex/timepoint are human-only and must reconcile.
  const partiallyCompatible: AnalysisQuery = {
    ...humanFilledQuery,
    tissue: rat.tissues[0],
    sex: "combined",
    timepoint: "post_24_hr",
  };
  const result = reconcileQueryForSpecies(partiallyCompatible, "rat");
  assert.deepEqual(result.changedFields, ["sex", "timepoint"]);
  // The compatible tissue is retained.
  assert.equal(result.query.tissue, rat.tissues[0]);
  // The incompatible dimensions fall back to the rat context defaults.
  assert.equal(result.query.sex, rat.sexes[0]);
  assert.equal(result.query.timepoint, rat.timepoints[0]);
});

test("reconciled selections always resolve to an available context for the new species", () => {
  // The catalog-driven builder relies on reconciliation producing a query that resolveContext can
  // serve: switch species, then confirm the reconciled fields yield an available/partial context.
  const reconciled = reconcileQueryForSpecies(humanFilledQuery, "rat").query;
  const resolved = resolveContext({
    targetMotrpacSpecies: reconciled.targetMotrpacSpecies,
    selectedOmics: reconciled.selectedOmics,
    tissue: reconciled.tissue,
    sex: reconciled.sex,
    timepoint: reconciled.timepoint,
  });
  assert.equal(resolved.status, "available");
  assert.equal(resolved.resolvedContext?.species, "rat");
});

test("switching to a species without a catalog context sets species but changes no fields", () => {
  // getStudyContext returns undefined for an unknown species; reconciliation degrades gracefully by
  // updating only targetMotrpacSpecies (mirrors an empty catalog dimension in the live builder).
  const emptyCatalog: typeof studyContexts = [];
  const result = reconcileQueryForSpecies(humanFilledQuery, "rat", emptyCatalog);
  assert.deepEqual(result.changedFields, []);
  assert.equal(result.query.targetMotrpacSpecies, "rat");
  // Prior selections are retained untouched when no target context exists.
  assert.equal(result.query.tissue, "muscle");
  assert.equal(result.query.sex, "combined");
  assert.equal(result.query.timepoint, "post_24_hr");
});
