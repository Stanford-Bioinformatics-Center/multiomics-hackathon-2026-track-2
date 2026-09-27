import type { ColumnResult, StoryRecord, StoryResponse } from "../api/client";

export type Dimension = "exercise" | "contrast" | "timepoint" | "layer" | "species" | "sex" | "tissue";
export type Significance = "raw" | "bh" | "bonferroni";
export type LiveFilters = Record<Dimension, string> & { significance: Significance };
export type FilterContext = Pick<StoryRecord, Dimension | "motrpac_measured">;

export const ALL = "all";
export const DEFAULT_FILTERS: LiveFilters = {
  exercise: ALL,
  contrast: ALL,
  timepoint: ALL,
  significance: "bh",
  layer: ALL,
  species: ALL,
  sex: ALL,
  tissue: ALL,
};

export const DIMENSIONS: Dimension[] = [
  "exercise", "contrast", "timepoint", "layer", "species", "sex", "tissue",
];

export const PATIENT_ONLY = "patient study only";
export const TRAINED = "TRAIN-SED";
export const COMBINED_SEX = "sexes combined";
const MOTRPAC_CONTRASTS = new Set(["EE-CON", "RE-CON", "EE-EE", "RE-RE", "CON-CON"]);

export function exerciseKind(value: string | null | undefined): string {
  if (!value) return "not applicable";
  const lower = value.toLowerCase();
  if (lower.includes("control") && !lower.includes("endurance") && !lower.includes("resistance")) return "control-only drift";
  if (lower.includes("endurance")) return "endurance";
  if (lower.includes("resistance")) return "resistance";
  return "exercise type unspecified";
}

/** Map a record onto the small set of values each filter offers. The raw value stays in the table. */
export function dimensionValue(record: FilterContext, dimension: Dimension): string {
  if (dimension === "exercise") return exerciseKind(record.exercise);
  if (dimension === "contrast") {
    if (record.motrpac_measured === false) return PATIENT_ONLY;
    if (record.species === "rat") return TRAINED;
    return record.contrast && MOTRPAC_CONTRASTS.has(record.contrast) ? record.contrast : PATIENT_ONLY;
  }
  if (dimension === "sex") {
    if (!record.sex || record.sex === "all" || record.sex === "both") return COMBINED_SEX;
    return record.sex;
  }
  return record[dimension] ?? "not applicable";
}

export function matchesFilters(record: FilterContext, filters: LiveFilters, except?: Dimension): boolean {
  return DIMENSIONS.every((dimension) =>
    dimension === except ||
    filters[dimension] === ALL ||
    dimensionValue(record, dimension) === filters[dimension]
  );
}

export function filterRecords<T extends FilterContext>(records: T[], filters: LiveFilters): T[] {
  return records.filter((record) => matchesFilters(record, filters));
}

export function optionAvailable(
  contexts: FilterContext[], filters: LiveFilters, dimension: Dimension, option: string,
): boolean {
  if (option === ALL) return true;
  return contexts.some((context) =>
    dimensionValue(context, dimension) === option && matchesFilters(context, filters, dimension)
  );
}

const ORDER: Record<Dimension, string[]> = {
  exercise: ["endurance", "resistance", "control-only drift", "exercise type unspecified", "not applicable"],
  contrast: ["EE-CON", "RE-CON", "EE-EE", "CON-CON", TRAINED, PATIENT_ONLY],
  timepoint: ["rest", "rest to peak", "during_20_min", "during_40_min", "post_10_min", "post_15_30_45_min",
    "post_3.5_4_hr", "post_24_hr", "1 wk", "2 wk", "4 wk", "8 wk", "all sampled weeks"],
  layer: ["transcriptomics", "proteomics", "phosphoproteomics", "metabolomics", "methylation"],
  species: ["human", "rat"],
  sex: [COMBINED_SEX, "female", "male"],
  tissue: ["skeletal muscle", "blood", "heart", "adipose"],
};

export function optionValues(contexts: FilterContext[], dimension: Dimension): string[] {
  const observed = new Set(contexts.map((context) => dimensionValue(context, dimension)));
  const requested = dimension === "exercise" ? ["endurance", "resistance", "control-only drift"]
    : dimension === "contrast" ? ["EE-CON", "EE-EE", "CON-CON"]
    : dimension === "layer" ? ["transcriptomics", "proteomics", "phosphoproteomics", "metabolomics", "methylation"]
    : dimension === "tissue" ? ["skeletal muscle", "blood", "adipose"] : [];
  const values = Array.from(new Set([...requested, ...observed]));
  const rank = (value: string) => {
    const index = ORDER[dimension].indexOf(value);
    return index === -1 ? 999 : index;
  };
  return [ALL, ...values.sort((a, b) => rank(a) - rank(b) || a.localeCompare(b))];
}

export const OPTION_LABELS: Record<string, string> = {
  "EE-CON": "EE-CON · endurance vs control",
  "RE-CON": "RE-CON · resistance vs control",
  "EE-EE": "EE-EE · endurance post − pre",
  "CON-CON": "CON-CON · control drift",
  [TRAINED]: "Trained vs sedentary (rat)",
  [PATIENT_ONLY]: "Patient study only",
  during_20_min: "20 min during",
  during_40_min: "40 min during",
  post_10_min: "10 min after",
  post_15_30_45_min: "15–45 min after",
  "post_3.5_4_hr": "3.5–4 h after",
  post_24_hr: "24 h after",
  "rest to peak": "Rest → peak",
};

export function optionLabel(value: string): string {
  return OPTION_LABELS[value] ?? value.replace(/_/g, " ");
}

export function isSummaryRecord(record: StoryRecord): boolean {
  return record.record_type.includes("summary");
}

export function significanceValue(record: StoryRecord, method: Significance): number | null {
  const value = method === "raw" ? record.p_value : method === "bh" ? record.q_value : record.bonferroni_p ?? null;
  return typeof value === "number" && Number.isFinite(value) ? value : null;
}

/** true / false against 0.05, or null when the selected rule cannot be checked for this row. */
export function significancePass(record: StoryRecord, method: Significance): boolean | null {
  // Summary rows carry the minimum q of a tested family, not a q for their count.
  if (isSummaryRecord(record) || record.motrpac_measured === false) return null;
  const value = significanceValue(record, method);
  return value === null ? null : value < 0.05;
}

export function unavailableRuleReason(record: StoryRecord, method: Significance): string {
  if (isSummaryRecord(record)) return "Summary row: counts, not a single test.";
  if (record.motrpac_measured === false) return "No MoTrPAC measurement for this row.";
  const note = record.metrics?.bonferroni_note;
  if (method === "bonferroni" && typeof note === "string") return note;
  if (method === "bonferroni") return "Bonferroni unavailable: the committed export does not reproduce this BH family.";
  if (method === "raw") return "The committed source table has no raw P for this row.";
  return "The committed source table has no BH q for this row.";
}

function direction(value: string | null | undefined): "up" | "down" | null {
  if (!value) return null;
  const lower = value.toLowerCase();
  if (lower === "up" || lower === "positive") return "up";
  if (lower === "down" || lower === "negative") return "down";
  return null;
}

export const AGREEMENT = {
  same: "same direction as disease",
  opposite: "opposite direction to disease",
  nonsignificant: "measured, no significant change",
  unchecked: "measured, rule not checkable",
  drift: "control drift, not an exercise effect",
  missing: "not measured",
} as const;

/** Disease direction versus MoTrPAC direction for one molecule, under the selected rule. */
export function agreementClass(record: StoryRecord, method: Significance): string {
  if (record.motrpac_measured === false || record.value === null || record.molecule === null) return AGREEMENT.missing;
  const disease = direction(record.disease_direction);
  const exercise = direction(record.motrpac_direction);
  if (!disease || !exercise) return AGREEMENT.missing;
  if (record.contrast === "CON-CON") return AGREEMENT.drift;
  const significant = significancePass(record, method);
  if (significant === null) return AGREEMENT.unchecked;
  if (!significant) return AGREEMENT.nonsignificant;
  return disease === exercise ? AGREEMENT.same : AGREEMENT.opposite;
}

export function unavailableReason(
  dimension: Dimension, option: string, contexts: FilterContext[], filters: LiveFilters,
): string {
  const existsAnywhere = contexts.some((context) => dimensionValue(context, dimension) === option);
  if (!existsAnywhere) {
    if (option === "EE-EE" || option === "CON-CON") return "Only the metabolite time course carries EE-EE and CON-CON; the committed gene and protein exports are EE-CON.";
    if (option === "resistance" && filters.species === "rat") return "Rats were endurance-trained only.";
    if (option === "methylation" || option === "phosphoproteomics") return "Not part of these PAH results; use a live run for layers the engine scores.";
    if (option === "adipose") return "No PAH result here uses adipose tissue.";
    return "No result is available for this option.";
  }
  return "No result combines this option with the other filters. Set a conflicting filter to All.";
}

/** Engine comparison columns → the same filter values as story records. */
const TISSUE: Record<string, string> = { VL: "skeletal muscle", "SKM-GN": "skeletal muscle", "SKM-VL": "skeletal muscle", HEART: "heart" };
const LAYER: Record<string, string> = { RNA: "transcriptomics", PROT: "proteomics", PHOSPHO: "phosphoproteomics", METAB: "metabolomics" };
const HUMAN_TIME: Record<string, string> = { "15-45m": "post_15_30_45_min", "3.5h": "post_3.5_4_hr", "24h": "post_24_hr", "10m": "post_10_min" };

export function columnContext(column: ColumnResult): FilterContext {
  const rat = column.dataset === "rat_train";
  const resistance = column.column_id.includes("|RE_");
  return {
    species: column.species,
    tissue: TISSUE[column.tissue] ?? column.tissue,
    layer: LAYER[column.layer] ?? column.layer,
    exercise: rat ? "endurance training" : resistance ? "resistance" : "endurance",
    contrast: rat ? TRAINED : resistance ? "RE-CON" : "EE-CON",
    timepoint: rat ? column.timepoint.replace("w", " wk") : HUMAN_TIME[column.timepoint] ?? column.timepoint,
    sex: column.sex === "F" ? "female" : column.sex === "M" ? "male" : "all",
    motrpac_measured: true,
  };
}

export function columnSignificance(column: ColumnResult, method: Significance): number | null {
  const value = method === "raw" ? column.camera_p : method === "bh" ? column.camera_fdr : column.camera_bonferroni;
  return typeof value === "number" && Number.isFinite(value) ? value : null;
}

/** Set-level direction of a disease signature in one MoTrPAC comparison, in neutral words. */
export function columnDirection(column: ColumnResult, method: Significance): string {
  const value = columnSignificance(column, method);
  if (value === null || column.camera_t === null) return "not testable";
  if (value >= 0.05) return "no significant set-level shift";
  return column.camera_t > 0 ? "opposite direction to disease" : "same direction as disease";
}

export type StoryAvailability = StoryResponse["availability"];
