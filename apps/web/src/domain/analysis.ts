export type MotrpacSpecies = "human" | "rat";

export type SignatureSourceSpecies = "human" | "rat" | "mouse" | "other";

export type ExerciseStudyDesign = "acute" | "chronic";

export type OmicLayer =
  | "genomics"
  | "epigenomics"
  | "transcriptomics"
  | "proteomics"
  | "metabolomics";

export type MolecularEntityType =
  | "variant"
  | "gene"
  | "transcript"
  | "protein"
  | "ptm_site"
  | "metabolite";

export type CrossOmicState =
  | "concordant_up"
  | "concordant_down"
  | "rna_up_protein_down"
  | "rna_down_protein_up"
  | "rna_only"
  | "protein_only"
  | "near_zero"
  | "insufficient_data";

export interface DirectedSignature {
  id: string;
  sourceSpecies: SignatureSourceSpecies;
  featureCount: number;
}

export interface AnalysisQuery {
  signatureId?: string;
  uploadedSignature?: DirectedSignature;
  signatureSourceSpecies: SignatureSourceSpecies;
  targetMotrpacSpecies: MotrpacSpecies;
  selectedOmics: OmicLayer[];
  tissue?: string;
  sex?: "female" | "male" | "combined";
  timepoint?: string;
  significanceThreshold?: number;
  includeNonSignificantResults: boolean;
}

export interface MotrpacStudyContext {
  species: MotrpacSpecies;
  studyDesign: ExerciseStudyDesign;
  displayStudyDesign: string;
  tissues: string[];
  omics: OmicLayer[];
  sexes: Array<"female" | "male" | "combined">;
  timepoints: string[];
  datasetId: string;
  datasetVersion: string;
}

export interface LayerCapability {
  layer: OmicLayer;
  status: "implemented" | "partial" | "planned";
  entityTypes: MolecularEntityType[];
}

export interface AvailabilityRequest {
  targetMotrpacSpecies: MotrpacSpecies;
  selectedOmics: OmicLayer[];
  tissue?: string;
  sex?: "female" | "male" | "combined";
  timepoint?: string;
}

export interface AvailabilityResult {
  status:
    | "available"
    | "partially_available"
    | "unsupported_omics"
    | "no_matching_context";
  resolvedContext?: MotrpacStudyContext;
  availableOmics: OmicLayer[];
  unavailableOmics: OmicLayer[];
  message: string;
  suggestedAlternatives: Array<{
    species?: MotrpacSpecies;
    tissue?: string;
    timepoint?: string;
    omics?: OmicLayer[];
  }>;
}

export const studyContexts: MotrpacStudyContext[] = [
  {
    species: "human",
    studyDesign: "acute",
    displayStudyDesign: "Acute exercise study",
    tissues: ["muscle", "blood"],
    omics: ["transcriptomics", "proteomics"],
    sexes: ["combined"],
    timepoints: [
      "during_20_min",
      "during_40_min",
      "post_10_min",
      "post_15_30_45_min",
      "post_3.5_4_hr",
      "post_24_hr",
    ],
    datasetId: "motrpac-human-pre-suspension",
    datasetVersion: "catalog-supplied",
  },
  {
    species: "rat",
    studyDesign: "chronic",
    displayStudyDesign: "Chronic exercise-training study",
    tissues: ["SKM-GN", "SKM-VL"],
    omics: ["transcriptomics", "proteomics"],
    sexes: ["female", "male"],
    timepoints: ["1w", "2w", "4w", "8w"],
    datasetId: "motrpac-rat-training",
    datasetVersion: "catalog-supplied",
  },
];

export const layerCapabilities: LayerCapability[] = [
  { layer: "transcriptomics", status: "implemented", entityTypes: ["gene", "transcript"] },
  { layer: "proteomics", status: "implemented", entityTypes: ["protein", "ptm_site"] },
  { layer: "genomics", status: "planned", entityTypes: ["variant", "gene"] },
  { layer: "epigenomics", status: "planned", entityTypes: ["gene"] },
  { layer: "metabolomics", status: "planned", entityTypes: ["metabolite"] },
];

export function getStudyContext(species: MotrpacSpecies, catalog = studyContexts) {
  return catalog.find((context) => context.species === species);
}

export function reconcileQueryForSpecies(
  query: AnalysisQuery,
  species: MotrpacSpecies,
  catalog = studyContexts,
): { query: AnalysisQuery; changedFields: Array<"tissue" | "sex" | "timepoint"> } {
  const context = getStudyContext(species, catalog);
  if (!context) return { query: { ...query, targetMotrpacSpecies: species }, changedFields: [] };

  const changedFields: Array<"tissue" | "sex" | "timepoint"> = [];
  const tissue = context.tissues.includes(query.tissue ?? "") ? query.tissue : context.tissues[0];
  const sex = context.sexes.includes(query.sex ?? "combined") ? query.sex : context.sexes[0];
  const timepoint = context.timepoints.includes(query.timepoint ?? "")
    ? query.timepoint
    : context.timepoints[0];

  if (tissue !== query.tissue) changedFields.push("tissue");
  if (sex !== query.sex) changedFields.push("sex");
  if (timepoint !== query.timepoint) changedFields.push("timepoint");

  return {
    query: { ...query, targetMotrpacSpecies: species, tissue, sex, timepoint },
    changedFields,
  };
}

export function resolveContext(
  query: AvailabilityRequest,
  catalog = studyContexts,
  capabilities = layerCapabilities,
): AvailabilityResult {
  const speciesContext = getStudyContext(query.targetMotrpacSpecies, catalog);

  if (!speciesContext) {
    return {
      status: "no_matching_context",
      availableOmics: [],
      unavailableOmics: query.selectedOmics,
      message: "No MoTrPAC study is available for the selected species.",
      suggestedAlternatives: [],
    };
  }

  const implemented = new Set(
    capabilities.filter((item) => item.status === "implemented").map((item) => item.layer),
  );
  const availableOmics = query.selectedOmics.filter(
    (layer) => implemented.has(layer) && speciesContext.omics.includes(layer),
  );
  const unavailableOmics = query.selectedOmics.filter((layer) => !availableOmics.includes(layer));

  const filterMismatch =
    (query.tissue && !speciesContext.tissues.includes(query.tissue)) ||
    (query.sex && !speciesContext.sexes.includes(query.sex)) ||
    (query.timepoint && !speciesContext.timepoints.includes(query.timepoint));

  if (filterMismatch) {
    return {
      status: "no_matching_context",
      resolvedContext: speciesContext,
      availableOmics,
      unavailableOmics,
      message: `The ${speciesContext.displayStudyDesign.toLowerCase()} is available, but no result matches the selected tissue, sex, or time point.`,
      suggestedAlternatives: [
        {
          species: speciesContext.species,
          tissue: speciesContext.tissues[0],
          timepoint: speciesContext.timepoints[0],
          omics: speciesContext.omics,
        },
      ],
    };
  }

  if (availableOmics.length === 0) {
    return {
      status: "unsupported_omics",
      resolvedContext: speciesContext,
      availableOmics,
      unavailableOmics,
      message: "The selected omic layers are not currently connected for this species context.",
      suggestedAlternatives: [{ species: speciesContext.species, omics: speciesContext.omics }],
    };
  }

  return {
    status: unavailableOmics.length ? "partially_available" : "available",
    resolvedContext: speciesContext,
    availableOmics,
    unavailableOmics,
    message: unavailableOmics.length
      ? "Supported layers can run; planned layers remain excluded and visible."
      : "A compatible MoTrPAC context was found.",
    suggestedAlternatives: [],
  };
}

export type VisualizationMode = "single_layer" | "pairwise" | "matrix";

/**
 * Request-time hint only (how many layers the user asked for). Do NOT use this to decide what to
 * render: selecting an unimplemented layer (e.g. metabolomics) alongside a real one must not
 * fabricate a multi-omic view. Use getVisualizationModeFromReturnedLayers for rendering.
 */
export function getVisualizationMode(selectedOmics: OmicLayer[]): VisualizationMode {
  if (selectedOmics.length <= 1) return "single_layer";
  if (selectedOmics.length === 2) return "pairwise";
  return "matrix";
}

/**
 * The mode that must drive rendering: it depends on the layers the backend actually RETURNED
 * (present in the results), not on what the user selected. This prevents a fake multi-omic matrix
 * when a requested layer produced no data. Layer codes here are the engine/store codes returned by
 * the API (e.g. "RNA", "PROT", "PHOSPHO", "METAB").
 */
export function getVisualizationModeFromReturnedLayers(returnedLayers: string[]): VisualizationMode {
  const distinct = Array.from(new Set(returnedLayers));
  if (distinct.length <= 1) return "single_layer";
  if (distinct.length === 2) return "pairwise";
  return "matrix";
}

export type ResultPresentationState =
  | "results"
  | "no_mapped_features"
  | "weak_evidence"
  | "processing_error";

export function classifyResultPresentation(input: {
  mappedFeatureCount: number;
  measurementCount: number;
  significantCount: number;
  technicalError?: boolean;
}): ResultPresentationState {
  if (input.technicalError) return "processing_error";
  if (input.mappedFeatureCount === 0) return "no_mapped_features";
  if (input.measurementCount > 0 && input.significantCount === 0) return "weak_evidence";
  return "results";
}
