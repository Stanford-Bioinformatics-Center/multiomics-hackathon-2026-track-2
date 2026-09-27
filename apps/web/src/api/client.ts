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
  verdict: string;
  significant: boolean;
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
  metabConvergence: () => getJSON<MetabConvergence>("/metabolomics/convergence"),
  metabExportUrl: () => `${BASE}/metabolomics/casestudy/export`,
  generalizedSignatures: () => getJSON<{ signatures: GeneralizedSignature[]; note: string }>("/generalized/signatures"),
  generalizedQuery: (id: string) => getJSON<GeneralizedQueryResult>(`/generalized/query/${id}`),
};
