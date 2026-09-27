/**
 * Typed client for the Exercise Signature Explorer API (apps/api, FastAPI).
 *
 * The backend's OpenAPI schema is the transport source of truth; these types converge to it and to
 * the existing domain/analysis.ts names. mprobe remains the scientific source of truth — this client
 * performs no computation, it only calls the service that wraps the engine.
 *
 * Base URL: VITE_API_BASE_URL (default "/api", so a single container can serve web + api).
 */

const BASE = (import.meta as { env?: Record<string, string> }).env?.VITE_API_BASE_URL ?? "/api";

export interface CapabilityFlags {
  engine: boolean;
  api: boolean;
  react: boolean;
  demo_validated: boolean;
}

export interface CatalogTissue {
  tissue: string;
  tissue_label: string;
  layers: string[];
  omics: string[];
  sexes: string[];
  timepoints: string[];
}

export interface CatalogContext {
  species: string;
  dataset: string;
  study_design: string;
  study_design_label: string;
  tissues: CatalogTissue[];
}

export interface Catalog {
  contexts: CatalogContext[];
  capability_matrix: Record<string, CapabilityFlags>;
  supported_source_species: string[];
  unsupported_source_species: string[];
  store_hash: string;
}

export interface AvailabilityResponse {
  status: "available" | "partially_available" | "unsupported_omics" | "no_matching_context";
  resolved_context: { species: string; dataset: string; study_design: string; study_design_label: string } | null;
  available_omics: string[];
  unavailable_omics: string[];
  message: string;
  suggested_alternatives: Array<Record<string, unknown>>;
}

export interface MappingAuditRow {
  input_row: number;
  input_id: string;
  selected_symbol: string;
  selected_method: string;
  status: string;
  n_candidates: number;
  ambiguous: boolean;
  candidates: string;
}

export interface MappingPreview {
  counts: Record<string, number>;
  rows: MappingAuditRow[];
  requires_confirmation?: boolean;
}

export interface ColumnResult {
  column_id: string;
  comparison_label: string;
  dataset: string;
  species: string;
  tissue: string;
  layer: string;
  sex: string;
  timepoint: string;
  n_measured: number;
  n_opposed: number;
  n_same: number;
  camera_t: number | null;
  camera_p: number | null;
  camera_fdr: number | null;
  camera_bonferroni: number | null;
  verdict: string;
  significant: boolean;
}

export interface StoryRecord {
  id: string;
  panel_id: string;
  record_type: string;
  species: string | null;
  tissue: string | null;
  layer: string | null;
  exercise: string | null;
  contrast: string | null;
  timepoint: string | null;
  sex: string | null;
  molecule: string | null;
  disease_direction: string | null;
  disease_contrast?: string | null;
  disease_layer?: string | null;
  disease_unit?: string | null;
  motrpac_direction: string | null;
  motrpac_measured?: boolean;
  disease_value?: number | null;
  disease_p_value?: number | null;
  disease_q_value?: number | null;
  value: number | null;
  unit: string | null;
  n?: number | null;
  p_value: number | null;
  q_value: number | null;
  bonferroni_p?: number | null;
  correction_method?: string | null;
  correction_family?: string | null;
  label: string;
  source: string;
  metrics?: Record<string, unknown>;
}

export interface StoryPanel {
  id: string;
  title: string;
  question: string;
  observation: string;
  interpretation: string;
  hypothesis: string;
  caveat: string;
  records: StoryRecord[];
}

export interface StoryResponse {
  panels: StoryPanel[];
  availability: {
    species: string[];
    tissues: string[];
    layers: string[];
    exercises: string[];
    contrasts: string[];
    timepoints: string[];
    sexes: string[];
    contexts: Array<Pick<StoryRecord, "panel_id" | "species" | "tissue" | "layer" | "exercise" | "contrast" | "timepoint" | "sex">>;
  };
  provenance?: Record<string, string>;
}

export interface FeatureEvidence {
  evidence_id: string;
  gene: string;
  disease_direction: number;
  column_id: string;
  comparison_label: string;
  dataset: string;
  species: string;
  tissue: string;
  layer: string;
  sex: string;
  timepoint: string;
  exercise_logfc: number | null;
  exercise_stat: number | null;
  fdr_bh: number | null;
  opposed: boolean | null;
  measured: boolean;
  mapping_decision_id: string;
  source_feature_id: string;
  n_collapsed: number | null;
  aggregation_method: string;
}

export interface AnalysisResponse {
  run_id: string;
  schema_version: string;
  status: "ok" | "no_compatible_data" | "no_mapped_features";
  message: string;
  signature_name: string;
  n_input_rows: number;
  n_counted_genes: number;
  n_mapped_rows: number;
  returned_layers: string[];
  contrasts: Record<string, Record<string, string>>;
  mapping_audit: MappingAuditRow[];
  columns: ColumnResult[];
  features: FeatureEvidence[];
  layer_discordance: Array<Record<string, unknown>>;
  multiplicity_family: { method: string; threshold: number; family_size: number; n_tests: number; column_ids: string[] };
  headline: Record<string, unknown>;
  caveats: string[];
  guardrails: Array<{ safe: string; unsafe: string }>;
  provenance: Record<string, unknown>;
}

export interface GeneralizedSignature {
  id: string;
  label: string;
  tissue: string;
  reference: string;
  reference_contrast_category: string;
}

export interface GeneralizedRankRow {
  layer: string;
  timepoint: string;
  contrast_category: string;
  n_used: number;
  spearman_rho: number | null;
  status: string;
}

export interface GeneralizedQueryResult {
  analysis_type: string;
  schema_version: string;
  signature_id: string;
  signature_label: string;
  tissue: string;
  reference: string;
  counts: Record<string, number>;
  filters: Record<string, unknown>;
  thresholds: Record<string, unknown>;
  interpretation: string;
  rank_correlation: GeneralizedRankRow[];
  coverage_rows: number;
  provenance: Record<string, unknown>;
}

/**
 * Request body for `POST /generalized/query` (task 7.3): run an uploaded/pasted disease signature
 * through query_core. The backend validates `signature_csv_text` against the query_core schema and
 * returns the same normalized payload as the bundled GET route, so the response is `GeneralizedQueryResult`.
 */
export interface GeneralizedQueryUploadInput {
  signature_csv_text: string;
  tissue: string;
  reference_contrast_category: string;
  reference?: string;
}

export interface MetabHit {
  evidence_id: string;
  refmet_name: string;
  pah_direction: string;
  pah_log2_effect: number | null;
  pah_q_value: number | null;
  context_label: string;
  motrpac_blood_match_basis: string;
  blood_ee_con_logfc_20min: number | null;
  blood_peak_logfc: number | null;
  blood_max_abs_control_logfc: number | null;
  pah_effect_over_control_drift: number | null;
  exercise_vs_pah_direction: string;
  caveat: string;
}

export interface MetabCaseStudy {
  run_id: string;
  schema_version: string;
  analysis_type: string;
  case_study_id: string;
  cohort_label: string;
  title: string;
  n_hits: number;
  n_hits_matched_blood: number;
  n_background_matched_blood: number;
  hit_labels: Record<string, number>;
  null_result: Record<string, Record<string, number>>;
  median_pah_effect_over_control_drift: number | null;
  contrasts: Record<string, Record<string, string>>;
  conclusion: string;
  caveats: string[];
  hits: MetabHit[];
  provenance: Record<string, unknown>;
}

/**
 * Discordance case-study responses (Option-1 read-only adapter: `discordance_casestudy.py`).
 * These mirror the adapter's dataclasses EXACTLY (CatalogSummaryResponse, ModelResponse,
 * PtmParentResponse). The adapter serves committed demo outputs and computes no statistics.
 */
export interface DiscordanceCatalogResponse {
  run_id: string;
  schema_version: string;
  analysis_type: "discordance_catalog";
  case_study_id: string;
  cohort_label: string;
  title: string;
  tissue: string;
  contrast_category: string;
  total_events: number;
  /** The four supported classification classes, zero-filled (supported_opposite may be 0). */
  classification_counts: {
    supported_concordant: number;
    supported_opposite: number;
    rna_response_protein_equivalent: number;
    indeterminate: number;
  };
  /** Coverage states, kept separate from the classification classes. */
  coverage_counts: {
    no_protein_measurement: number;
    no_rna_measurement: number;
  };
  occupancy_caveat: string;
  limitations: string[];
  provenance: Record<string, unknown>;
}

export interface DiscordanceModelMetricRow {
  subset: "all_mapped" | "rna_responsive";
  model: "zero" | "rna_only" | "temporal";
  n_rows: number | null;
  n_genes: number | null;
  mae: number | null;
  rmse: number | null;
  r2: number | null;
}

export interface DiscordanceModelResponse {
  run_id: string;
  schema_version: string;
  analysis_type: "discordance_model";
  case_study_id: string;
  tissue: string;
  contrast_category: string;
  /** Exactly 6 rows = {all_mapped, rna_responsive} x {zero, rna_only, temporal}. */
  metrics: DiscordanceModelMetricRow[];
  /** Visible honest-framing statement: the out-of-fold R2 is near zero (weak prediction). */
  weak_prediction_note: string;
  limitations: string[];
  provenance: Record<string, unknown>;
}

export interface PtmParentSummaryRow {
  tissue: string;
  contrast_category: string;
  timepoint: string;
  raw_phosphosite_rows: number | null;
  after_mapping_checks_with_uniprot: number | null;
  excluded_missing_or_ambiguous_site_mapping: number | null;
  matched_to_single_gene_uniprot_single_feature_parent: number | null;
  matched_with_both_site_and_parent_ci: number | null;
  supported_phosphosite_rows_among_matched: number | null;
  parent_equivalent_rows_among_matched: number | null;
  supported_phosphosite_and_parent_equivalent_rows: number | null;
  candidate_genes: number | null;
}

export interface PtmParentResponse {
  run_id: string;
  schema_version: string;
  analysis_type: "ptm_parent_audit";
  /** One row per timepoint. */
  summary_rows: PtmParentSummaryRow[];
  candidate_count: number;
  occupancy_caveat: string;
  provenance: Record<string, unknown>;
}

export interface MetabConvergence {
  contrast: string;
  motrpac_package_version: string;
  source_collection: string;
  unambiguous_matched_cells: number;
  n_within_tol: number;
  max_abs_diff: number;
  identical: boolean;
  excluded_ambiguous_module_rows: number;
  excluded_ambiguous_engine_rows: number;
  note: string;
}

export interface MetabolitePreview {
  n_input_rows: number;
  n_blood_matched: number;
  unmapped_ids: string[];
  rows: Array<{
    input_row: number;
    list_id: string;
    name: string;
    pah_direction: string;
    pah_log2_effect: number | null;
    pah_q_value: number | null;
    is_reported_hit: boolean;
    blood_matched: boolean;
    muscle_matched: boolean;
    context_label: string;
    observations: Array<{
      tissue: string;
      feature_id: string;
      match_basis: string;
      contrast: string;
      timepoint: string;
      logFC: number | null;
      p_value: number | null;
      bh_q: number | null;
    }>;
    source: string;
    caveat: string;
  }>;
  rule: string;
  note: string;
}


export type ExplorerKind = "genes" | "metabolites" | "pathways";
export interface ExplorerExample { id: string; kind: ExplorerKind; label: string; source: string }
export interface ExplorerInput {
  id: number; name: string; kind: ExplorerKind; source: string | null; n_rows: number; n_mapped: number; n_names?: number;
  unmapped: string[]; directed: boolean; ranked: boolean;
}
export interface ExplorerColumn {
  id: string; species: string; tissue: string; tissue_label: string; layer: string; category: string; exercise: string;
  time: string; time_label: string; time_rank: number; sex: string; dataset: string; study_label: string; n_tested?: number; n_measured?: number;
  set_t?: number | null; set_p?: number | null; set_bh?: number | null; set_bonferroni?: number | null;
  n_opposite?: number; n_same?: number; n_up?: number; n_down?: number;
  rho?: number | null; rho_p?: number | null; rho_bh?: number | null; rho_bonferroni?: number | null; rho_n?: number;
}
/** [molecule index, column index, value (log2 FC or pathway z), raw P, BH q, Bonferroni] */
export type ExplorerValue = [number, number, number | null, number | null, number | null, number | null];
export interface ExplorerLayer {
  id: string; label: string; available: boolean; reason?: string; n_matched: number; input: number | null;
  molecules?: Array<{ name: string; dir: number; score: number | null; set?: string }>; n_shown?: number;
  columns?: ExplorerColumn[]; values?: ExplorerValue[]; value_unit?: string;
  family?: { n_tests: number | null; scope: string; methods: string[] }; molecule_rule?: string;
}
export interface ExplorerResponse { inputs: ExplorerInput[]; layers: ExplorerLayer[]; saved?: { generated_at: string; example: string; note: string } }

export interface SignatureRowInput {
  gene_symbol?: string;
  uniprot?: string;
  ensembl?: string;
  rat_symbol?: string;
  direction?: string;
  group?: string;
  weight?: string;
  source?: string;
}

export interface AnalysisRequestInput {
  signature_rows?: SignatureRowInput[];
  signature_csv_text?: string;
  example_name?: string;
  signature_name?: string;
  target_species: "rat" | "human";
  selected_omics: string[];
  tissue?: string;
  sex?: string;
  timepoint?: string;
  fdr_threshold?: number;
  include_nonsignificant?: boolean;
}

async function getJSON<T>(path: string): Promise<T> {
  const r = await fetch(`${BASE}${path}`);
  if (!r.ok) throw new ApiError(r.status, await r.text());
  return r.json() as Promise<T>;
}

async function postJSON<T>(path: string, body: unknown): Promise<T> {
  const r = await fetch(`${BASE}${path}`, {
    method: "POST",
    headers: { "Content-Type": "application/json" },
    body: JSON.stringify(body),
  });
  if (!r.ok) throw new ApiError(r.status, await r.text());
  return r.json() as Promise<T>;
}

export class ApiError extends Error {
  constructor(public status: number, public detail: string) {
    super(`API ${status}: ${detail}`);
  }
}

export const api = {
  story: () => getJSON<StoryResponse>("/story"),
  explorerExamples: () => getJSON<{ examples: ExplorerExample[] }>("/explorer/examples"),
  explorerReadFile: (body: { filename: string; content_base64: string }) =>
    postJSON<{ text: string; kind: ExplorerKind | null; n_rows: number; columns: string[] }>("/explorer/read-file", body),
  explorerAnalyse: (body: { lists: Array<{ kind: ExplorerKind; name?: string; text?: string; example?: string }> }) =>
    postJSON<ExplorerResponse>("/explorer/analyse", body),
  health: () => getJSON<{ status: string; schema_version: string; store_hash: string }>("/health"),
  catalog: () => getJSON<Catalog>("/catalog"),
  availability: (body: {
    target_species: string;
    selected_omics: string[];
    tissue?: string;
    sex?: string;
    timepoint?: string;
  }) => postJSON<AvailabilityResponse>("/availability", body),
  validate: (body: { example_name?: string; signature_csv_text?: string; signature_rows?: SignatureRowInput[] }) =>
    postJSON<MappingPreview>("/signatures/validate", body),
  mappingsPreview: (body: { example_name?: string; signature_csv_text?: string; signature_rows?: SignatureRowInput[] }) =>
    postJSON<MappingPreview>("/mappings/preview", body),
  runComparison: (body: AnalysisRequestInput) => postJSON<AnalysisResponse>("/comparisons", body),
  getComparison: (runId: string) => getJSON<AnalysisResponse>(`/comparisons/${runId}`),
  features: (runId: string) => getJSON<{ run_id: string; features: FeatureEvidence[] }>(`/comparisons/${runId}/features`),
  provenance: (runId: string) => getJSON<Record<string, unknown>>(`/comparisons/${runId}/provenance`),
  reportUrl: (runId: string) => `${BASE}/comparisons/${runId}/report`,
  exportUrl: (runId: string) => `${BASE}/comparisons/${runId}/export`,
  metabCaseStudy: () => getJSON<MetabCaseStudy>("/metabolomics/casestudy"),
  metabolitePreview: (body: { csv_text: string }) => postJSON<MetabolitePreview>("/metabolomics/preview", body),
  metabConvergence: () => getJSON<MetabConvergence>("/metabolomics/convergence"),
  metabExportUrl: () => `${BASE}/metabolomics/casestudy/export`,
  generalizedSignatures: () => getJSON<{ signatures: GeneralizedSignature[]; note: string }>("/generalized/signatures"),
  generalizedQuery: (id: string) => getJSON<GeneralizedQueryResult>(`/generalized/query/${id}`),
  generalizedQueryUpload: (body: GeneralizedQueryUploadInput) =>
    postJSON<GeneralizedQueryResult>("/generalized/query", body),
  discordanceCatalog: () => getJSON<DiscordanceCatalogResponse>("/discordance/catalog"),
  discordanceModel: () => getJSON<DiscordanceModelResponse>("/discordance/model"),
  discordancePtmParent: () => getJSON<PtmParentResponse>("/discordance/ptm-parent"),
};
