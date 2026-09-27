import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import test from "node:test";
import {
  classifyResultPresentation,
  getStudyContext,
  getVisualizationMode,
  reconcileQueryForSpecies,
  resolveContext,
  studyContexts,
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
