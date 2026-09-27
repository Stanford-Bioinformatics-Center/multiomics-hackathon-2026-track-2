import { createElement, useCallback, useEffect, useMemo, useState, type ReactNode } from "react";
import AvailabilityState from "./AvailabilityState";
import {
  api,
  ApiError,
  type AnalysisRequestInput,
  type Catalog,
  type CatalogContext,
  type MappingPreview,
  type SignatureRowInput,
} from "../api/client";
import {
  getStudyContext,
  getVisualizationModeFromReturnedLayers,
  layerCapabilities,
  reconcileQueryForSpecies,
  type AnalysisQuery,
  type MotrpacSpecies,
  type OmicLayer,
  type SignatureSourceSpecies,
} from "../domain/analysis";

// getVisualizationModeFromReturnedLayers is re-exported for downstream consumers (LiveDashboard,
// task 5.2) that render off the layers the backend actually returned rather than what was requested.
export { getVisualizationModeFromReturnedLayers };

/** Maximum accepted uploaded-signature size (bytes). Empty or larger files are rejected (R4 AC4). */
export const MAX_SIGNATURE_UPLOAD_BYTES = 1_000_000;

function Button({
  className = "",
  children,
  ...props
}: {
  className?: string;
  children: ReactNode;
  [key: string]: unknown;
}) {
  return createElement("button", { className, type: "button", ...props }, children);
}

export const initialAnalysisQuery: AnalysisQuery = {
  signatureId: "pah-skm-v1.2",
  signatureSourceSpecies: "human",
  targetMotrpacSpecies: "rat",
  selectedOmics: ["transcriptomics", "proteomics"],
  tissue: "SKM-VL",
  sex: "female",
  timepoint: "8w",
  significanceThreshold: 0.05,
  includeNonSignificantResults: true,
};

/** Exactly one signature source is active at a time (R4 AC1). */
export type SignatureSource =
  | { kind: "example"; exampleName: string }
  | { kind: "paste"; rows: SignatureRowInput[]; text: string }
  | { kind: "upload"; csvText: string; fileName: string };

/**
 * The run request the builder emits once the user has confirmed the mapping preview. Task 5.2 wires
 * this into `LiveDashboard` (which calls `api.runComparison`); the builder itself never calls
 * `/api/comparisons` directly, so the confirmation gate (R4 AC10) stays enforced in one place.
 */
export interface QueryBuilderRunContract {
  request: AnalysisRequestInput;
  preview: MappingPreview;
}

function title(value: string): string {
  // ES2020 lib target: avoid String.prototype.replaceAll (ES2021). split/join is equivalent.
  return value
    .split("_")
    .join(" ")
    .replace(/\b\w/g, (character: string) => character.toUpperCase());
}

/** Parse a pasted signature into rows. Accepts a header row plus comma/tab-delimited columns. */
function parsePastedSignature(text: string): SignatureRowInput[] {
  const lines = text
    .split(/\r?\n/)
    .map((line) => line.trim())
    .filter((line) => line.length > 0);
  if (lines.length === 0) return [];
  const split = (line: string) => line.split(/[\t,]/).map((cell) => cell.trim());
  const header = split(lines[0]).map((cell) => cell.toLowerCase());
  const hasHeader = header.some((cell) =>
    ["gene_symbol", "gene", "symbol", "uniprot", "ensembl", "rat_symbol", "direction"].includes(cell),
  );
  const bodyLines = hasHeader ? lines.slice(1) : lines;
  const cols = hasHeader ? header : ["gene_symbol", "direction"];
  const rows: SignatureRowInput[] = [];
  for (const line of bodyLines) {
    const cells = split(line);
    if (cells.length === 0 || cells.every((cell) => cell === "")) continue;
    const row: SignatureRowInput = {};
    cols.forEach((col, index) => {
      const value = cells[index];
      if (value === undefined || value === "") return;
      switch (col) {
        case "gene_symbol":
        case "gene":
        case "symbol":
          row.gene_symbol = value;
          break;
        case "uniprot":
          row.uniprot = value;
          break;
        case "ensembl":
          row.ensembl = value;
          break;
        case "rat_symbol":
          row.rat_symbol = value;
          break;
        case "direction":
          row.direction = value;
          break;
        case "group":
          row.group = value;
          break;
        case "weight":
          row.weight = value;
          break;
        case "source":
          row.source = value;
          break;
        default:
          if (!row.gene_symbol) row.gene_symbol = value;
      }
    });
    if (row.gene_symbol || row.uniprot || row.ensembl || row.rat_symbol) rows.push(row);
  }
  return rows;
}

function Field({
  label,
  children,
}: {
  label: string;
  children: ReactNode;
}) {
  return (
    <div className="query-field">
      <span className="query-field__label">{label}</span>
      <div className="query-options">{children}</div>
    </div>
  );
}

type CatalogState =
  | { kind: "loading" }
  | { kind: "error"; message: string }
  | { kind: "ready"; catalog: Catalog };

/** Distinct dimension values available for a species context in the live catalog. */
interface CatalogDimensions {
  omics: string[];
  tissues: string[];
  sexes: string[];
  timepoints: string[];
}

function dimensionsForContext(context: CatalogContext | undefined): CatalogDimensions {
  if (!context) return { omics: [], tissues: [], sexes: [], timepoints: [] };
  const omics = new Set<string>();
  const tissues = new Set<string>();
  const sexes = new Set<string>();
  const timepoints = new Set<string>();
  for (const tissue of context.tissues) {
    tissues.add(tissue.tissue);
    tissue.omics.forEach((value) => omics.add(value));
    tissue.sexes.forEach((value) => sexes.add(value));
    tissue.timepoints.forEach((value) => timepoints.add(value));
  }
  return {
    omics: [...omics],
    tissues: [...tissues],
    sexes: [...sexes],
    timepoints: [...timepoints],
  };
}

/** Map a catalog omic string ("transcriptomics"/"RNA"/…) onto a domain OmicLayer where possible. */
function omicToLayer(value: string): OmicLayer | null {
  const normalized = value.toLowerCase();
  const direct: Record<string, OmicLayer> = {
    transcriptomics: "transcriptomics",
    proteomics: "proteomics",
    genomics: "genomics",
    epigenomics: "epigenomics",
    metabolomics: "metabolomics",
    rna: "transcriptomics",
    prot: "proteomics",
    protein: "proteomics",
  };
  return direct[normalized] ?? null;
}

export default function QueryBuilder({
  onQueryChange,
  onRun,
}: {
  onQueryChange: (query: AnalysisQuery) => void;
  /** Emitted once the mapping preview is confirmed. Task 5.2 wires this to `api.runComparison`. */
  onRun?: (contract: QueryBuilderRunContract) => void;
}) {
  const [query, setQuery] = useState(initialAnalysisQuery);
  const [changeNotice, setChangeNotice] = useState<string | null>(null);

  const [catalogState, setCatalogState] = useState<CatalogState>({ kind: "loading" });

  const [signatureSource, setSignatureSource] = useState<SignatureSource>({
    kind: "example",
    exampleName: "pah_muscle_malenfant2015",
  });
  const [pasteText, setPasteText] = useState("");
  const [uploadError, setUploadError] = useState<string | null>(null);

  const [preview, setPreview] = useState<MappingPreview | null>(null);
  const [previewError, setPreviewError] = useState<string | null>(null);
  const [previewLoading, setPreviewLoading] = useState(false);
  const [confirmed, setConfirmed] = useState(false);

  const loadCatalog = useCallback(async () => {
    setCatalogState({ kind: "loading" });
    try {
      const catalog = await api.catalog();
      setCatalogState({ kind: "ready", catalog });
    } catch (error) {
      const message = error instanceof ApiError ? `${error.status}: ${error.detail}` : String(error);
      setCatalogState({ kind: "error", message });
    }
  }, []);

  useEffect(() => {
    loadCatalog().catch(() => undefined);
  }, [loadCatalog]);

  // Species-derived study design comes from the in-memory domain (rat->chronic, human->acute); this
  // stays valid and is shown read-only. The backend derives study from species alone (R11.2/R11.3),
  // so there is no acute/chronic input control here (R11.3).
  const context = getStudyContext(query.targetMotrpacSpecies)!;

  const catalog = catalogState.kind === "ready" ? catalogState.catalog : null;
  const catalogContext = useMemo<CatalogContext | undefined>(
    () => catalog?.contexts.find((entry) => entry.species === query.targetMotrpacSpecies),
    [catalog, query.targetMotrpacSpecies],
  );
  const dimensions = useMemo(() => dimensionsForContext(catalogContext), [catalogContext]);

  // Species options come from the catalog when available; the union across contexts constrained to
  // the two MoTrPAC species the domain understands.
  const speciesOptions = useMemo<MotrpacSpecies[]>(() => {
    if (!catalog) return ["human", "rat"];
    const seen = new Set<string>();
    catalog.contexts.forEach((entry) => seen.add(entry.species));
    const filtered = (["human", "rat"] as MotrpacSpecies[]).filter((species) => seen.has(species));
    return filtered.length ? filtered : ["human", "rat"];
  }, [catalog]);

  const supportedSourceSpecies = useMemo<SignatureSourceSpecies[]>(() => {
    const all: SignatureSourceSpecies[] = ["human", "rat", "mouse", "other"];
    if (!catalog || catalog.supported_source_species.length === 0) return all;
    const supported = new Set(catalog.supported_source_species.map((value) => value.toLowerCase()));
    const filtered = all.filter((species) => supported.has(species));
    return filtered.length ? filtered : all;
  }, [catalog]);

  const totalCatalogFailure = catalogState.kind === "error";
  const selectorsDisabled = totalCatalogFailure || catalogState.kind === "loading";

  // Any change to the query or the signature source invalidates a prior confirmation (R4 AC10): the
  // user must re-preview + re-confirm before a comparison can run.
  function invalidateConfirmation() {
    setPreview(null);
    setPreviewError(null);
    setConfirmed(false);
  }

  const update = useCallback(
    (next: AnalysisQuery) => {
      setQuery(next);
      onQueryChange(next);
      invalidateConfirmation();
    },
    [onQueryChange],
  );

  function changeSpecies(species: MotrpacSpecies) {
    if (species === query.targetMotrpacSpecies) return;
    const nextContext = getStudyContext(species)!;
    const previousContext = context;
    const reconciled = reconcileQueryForSpecies(query, species);
    setChangeNotice(
      `Context changed: ${previousContext.displayStudyDesign} \u2192 ${nextContext.displayStudyDesign}. Updated incompatible fields: ${reconciled.changedFields.join(", ")}.`,
    );
    update(reconciled.query);
  }

  function toggleOmic(layer: OmicLayer) {
    const exists = query.selectedOmics.includes(layer);
    const selectedOmics = exists
      ? query.selectedOmics.filter((item) => item !== layer)
      : [...query.selectedOmics, layer];
    update({ ...query, selectedOmics });
  }

  // --- signature source handling (exactly one active source) -------------------------------------

  function selectExample(exampleName: string) {
    setSignatureSource({ kind: "example", exampleName });
    setUploadError(null);
    invalidateConfirmation();
  }

  function applyPaste(text: string) {
    setPasteText(text);
    const rows = parsePastedSignature(text);
    setSignatureSource({ kind: "paste", rows, text });
    setUploadError(null);
    invalidateConfirmation();
  }

  function handleUpload(file: File | undefined) {
    if (!file) return;
    if (file.size === 0) {
      // Reject empty upload; retain the prior source unchanged (R4 AC4).
      setUploadError("The uploaded file is empty.");
      return;
    }
    if (file.size > MAX_SIGNATURE_UPLOAD_BYTES) {
      setUploadError(
        `The uploaded file is too large (max ${Math.floor(MAX_SIGNATURE_UPLOAD_BYTES / 1000)} KB).`,
      );
      return;
    }
    const reader = new FileReader();
    reader.onerror = () => setUploadError("The uploaded file could not be read.");
    reader.onload = () => {
      const csvText = typeof reader.result === "string" ? reader.result : "";
      if (csvText.trim().length === 0 || parsePastedSignature(csvText).length === 0) {
        // Unparseable / no rows: reject and retain the prior source (R4 AC4).
        setUploadError("The uploaded file could not be parsed as a signature CSV.");
        return;
      }
      setUploadError(null);
      setSignatureSource({ kind: "upload", csvText, fileName: file.name });
      invalidateConfirmation();
    };
    reader.readAsText(file);
  }

  // --- request assembly --------------------------------------------------------------------------

  // Example/paste -> signature_rows (or example_name); upload -> signature_csv_text (R4 AC2, AC3).
  function signatureRequestFields(): Pick<
    AnalysisRequestInput,
    "signature_rows" | "signature_csv_text" | "example_name"
  > {
    switch (signatureSource.kind) {
      case "example":
        return { example_name: signatureSource.exampleName };
      case "paste":
        return { signature_rows: signatureSource.rows };
      case "upload":
        return { signature_csv_text: signatureSource.csvText };
    }
  }

  function previewRequestBody() {
    const fields = signatureRequestFields();
    return {
      example_name: fields.example_name,
      signature_rows: fields.signature_rows,
      signature_csv_text: fields.signature_csv_text,
    };
  }

  function analysisRequest(): AnalysisRequestInput {
    return {
      ...signatureRequestFields(),
      target_species: query.targetMotrpacSpecies,
      selected_omics: query.selectedOmics,
      tissue: query.tissue,
      sex: query.sex,
      timepoint: query.timepoint,
      // FDR threshold is an analysis parameter recorded in the run_id via service._run_id (R4 AC13);
      // not re-implemented here. UI constrains it to 0 < t <= 1 (R4 AC8).
      fdr_threshold: query.significanceThreshold,
      // include-nonsignificant is presentation-only and already excluded from run_id (R4 AC14).
      include_nonsignificant: query.includeNonSignificantResults,
    };
  }

  // --- mapping-preview confirmation gate (R4 AC10, AC11, AC12) -----------------------------------

  async function requestPreview() {
    setPreviewLoading(true);
    setPreviewError(null);
    setConfirmed(false);
    try {
      const result = await api.mappingsPreview(previewRequestBody());
      setPreview(result);
    } catch (error) {
      // Preview failure blocks the comparison, retains selections, shows an error (R4 AC11).
      const message = error instanceof ApiError ? `${error.status}: ${error.detail}` : String(error);
      setPreviewError(message);
      setPreview(null);
    } finally {
      setPreviewLoading(false);
    }
  }

  function confirmAndRun() {
    if (!preview) return;
    setConfirmed(true);
    // The builder never calls /api/comparisons itself; it emits the confirmed contract so the gate
    // lives in one place (R4 AC10). Task 5.2's LiveDashboard performs the actual runComparison call.
    onRun?.({ request: analysisRequest(), preview });
  }

  const structurallyInvalid = query.selectedOmics.length === 0;
  const signatureReady =
    (signatureSource.kind === "example" && signatureSource.exampleName.length > 0) ||
    (signatureSource.kind === "paste" && signatureSource.rows.length > 0) ||
    (signatureSource.kind === "upload" && signatureSource.csvText.length > 0);

  // FDR threshold constrained to 0 < t <= 1 (R4 AC8). UI note: this upper bound is inclusive in the
  // UI, whereas the API validates gt=0, lt=1 — a value of exactly 1.0 is a UI-only convenience and
  // should be avoided when the backend rejects it; the offered presets stay strictly inside (0,1).
  const fdrPresets = [0.01, 0.05, 0.1];
  const ambiguousRows = (preview?.rows ?? []).filter((row) => row.ambiguous || row.n_candidates > 1);
  const requiresConfirmation = preview?.requires_confirmation ?? ambiguousRows.length > 0;

  return (
    <div className="query-builder">
      <div className="query-builder__title">
        <div>
          <span className="query-builder__eyebrow">Analysis query</span>
          <strong>User selections</strong>
        </div>
        <span className="query-builder__version">
          {catalog ? `STORE ${catalog.store_hash.slice(0, 8)}` : "LIVE API"}
        </span>
      </div>

      {totalCatalogFailure && (
        <div className="query-catalog-error" role="alert">
          The analysis catalog is unavailable, so target selectors are disabled.{" "}
          {catalogState.kind === "error" ? catalogState.message : ""}{" "}
          <Button className="query-chip" onClick={() => loadCatalog()}>
            Retry
          </Button>
        </div>
      )}

      {/* --- signature source: exactly one of example / paste / upload (R4 AC1) --- */}
      <Field label="Signature source">
        {(["example", "paste", "upload"] as SignatureSource["kind"][]).map((kind) => (
          <Button
            className={`query-chip ${signatureSource.kind === kind ? "query-chip--selected" : ""}`}
            onClick={() => {
              if (kind === "example") selectExample(signatureSource.kind === "example" ? signatureSource.exampleName : "pah_muscle_malenfant2015");
              else if (kind === "paste") applyPaste(pasteText);
              else setSignatureSource((prev) => (prev.kind === "upload" ? prev : { kind: "upload", csvText: "", fileName: "" }));
            }}
            key={kind}
          >
            {title(kind)}
          </Button>
        ))}
      </Field>

      {signatureSource.kind === "example" && (
        <Field label="Built-in example">
          {["pah_muscle_malenfant2015"].map((name) => (
            <Button
              className={`query-chip ${
                signatureSource.exampleName === name ? "query-chip--selected" : ""
              }`}
              onClick={() => selectExample(name)}
              key={name}
            >
              {title(name)}
            </Button>
          ))}
        </Field>
      )}

      {signatureSource.kind === "paste" && (
        <div className="query-field">
          <span className="query-field__label">Pasted signature (CSV/TSV)</span>
          <textarea
            className="query-paste"
            value={pasteText}
            onChange={(event) => applyPaste(event.target.value)}
            placeholder={"gene_symbol,direction\nGENE1,-1\nGENE2,1"}
            rows={5}
          />
          <small>{signatureSource.rows.length} parsed row(s)</small>
        </div>
      )}

      {signatureSource.kind === "upload" && (
        <div className="query-field">
          <span className="query-field__label">Uploaded CSV</span>
          <input
            type="file"
            accept=".csv,text/csv,text/plain"
            onChange={(event) => handleUpload(event.target.files?.[0])}
          />
          {signatureSource.fileName && !uploadError && (
            <small>Loaded {signatureSource.fileName}</small>
          )}
        </div>
      )}

      {uploadError && (
        <div className="query-upload-error" role="alert">
          {uploadError}
        </div>
      )}

      <Field label="Disease-signature source species">
        {supportedSourceSpecies.map((species) => (
          <Button
            className={`query-chip ${query.signatureSourceSpecies === species ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, signatureSourceSpecies: species })}
            key={species}
          >
            {title(species)}
          </Button>
        ))}
      </Field>

      <Field label="Target MoTrPAC species">
        {speciesOptions.map((species) => {
          const speciesContext = getStudyContext(species)!;
          return (
            <Button
              className={`species-choice ${query.targetMotrpacSpecies === species ? "species-choice--selected" : ""}`}
              onClick={() => changeSpecies(species)}
              disabled={selectorsDisabled}
              key={species}
            >
              <strong>{title(species)}</strong>
              <span>{speciesContext.displayStudyDesign}</span>
            </Button>
          );
        })}
      </Field>

      <div className="derived-study-context">
        <span>Derived study design</span>
        <strong>{context.displayStudyDesign}</strong>
        <small>Backend-derived from species \u00b7 read-only metadata \u00b7 not a user input</small>
      </div>

      {changeNotice && <div className="query-change-notice">{changeNotice}</div>}

      <Field label="Omic layers \u00b7 multi-select">
        {layerCapabilities.map((capability) => {
          // A layer is offerable when the live catalog exposes it for the current species context.
          const catalogHasLayer =
            dimensions.omics.length === 0
              ? false
              : dimensions.omics.some((value) => omicToLayer(value) === capability.layer);
          const disabled = selectorsDisabled || (catalog !== null && !catalogHasLayer);
          return (
            <Button
              className={`omic-choice ${query.selectedOmics.includes(capability.layer) ? "omic-choice--selected" : ""} ${capability.status !== "implemented" ? "omic-choice--planned" : ""}`}
              onClick={() => toggleOmic(capability.layer)}
              disabled={disabled}
              key={capability.layer}
            >
              <span>{title(capability.layer)}</span>
              <b>{capability.status === "implemented" ? "Implemented" : "Planned"}</b>
            </Button>
          );
        })}
      </Field>

      <Field label="Tissue">
        {(dimensions.tissues.length ? dimensions.tissues : context.tissues).map((tissue) => (
          <Button
            className={`query-chip ${query.tissue === tissue ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, tissue })}
            disabled={selectorsDisabled || dimensions.tissues.length === 0}
            key={tissue}
          >
            {tissue}
          </Button>
        ))}
      </Field>

      <Field label="Sex">
        {(dimensions.sexes.length ? dimensions.sexes : context.sexes).map((sex) => (
          <Button
            className={`query-chip ${query.sex === sex ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, sex: sex as AnalysisQuery["sex"] })}
            disabled={selectorsDisabled || dimensions.sexes.length === 0}
            key={sex}
          >
            {title(sex)}
          </Button>
        ))}
      </Field>

      <Field label="Time point">
        {(dimensions.timepoints.length ? dimensions.timepoints : context.timepoints).map((timepoint) => (
          <Button
            className={`query-chip ${query.timepoint === timepoint ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, timepoint })}
            disabled={selectorsDisabled || dimensions.timepoints.length === 0}
            key={timepoint}
          >
            {title(timepoint)}
          </Button>
        ))}
      </Field>

      <Field label="FDR threshold (0 < t \u2264 1)">
        {fdrPresets.map((threshold) => (
          <Button
            className={`query-chip ${query.significanceThreshold === threshold ? "query-chip--selected" : ""}`}
            onClick={() => update({ ...query, significanceThreshold: threshold })}
            key={threshold}
          >
            q \u2264 {threshold}
          </Button>
        ))}
      </Field>

      <div className="query-toggle-row">
        <div>
          <strong>Include nonsignificant measurements</strong>
          <span>Presentation-only \u00b7 never recorded in run_id</span>
        </div>
        <Button
          className={`query-toggle ${query.includeNonSignificantResults ? "query-toggle--on" : ""}`}
          onClick={() => update({ ...query, includeNonSignificantResults: !query.includeNonSignificantResults })}
          aria-pressed={query.includeNonSignificantResults}
        >
          <span />
        </Button>
      </div>

      <AvailabilityState
        result={{
          status: dimensions.omics.length === 0 && catalog !== null ? "no_matching_context" : "available",
          resolvedContext: context,
          availableOmics: query.selectedOmics,
          unavailableOmics: [],
          message: catalog
            ? "Target selectors are sourced from the live /api/catalog."
            : "Loading the live catalog\u2026",
          suggestedAlternatives: [],
        }}
      />

      {/* --- mapping-preview confirmation gate before any comparison (R4 AC10) --- */}
      <div className="mapping-gate">
        <Button
          className="run-comparison"
          disabled={structurallyInvalid || !signatureReady || previewLoading || selectorsDisabled}
          onClick={() => requestPreview()}
        >
          {structurallyInvalid
            ? "Select at least one omic layer"
            : !signatureReady
              ? "Provide a signature source"
              : previewLoading
                ? "Previewing mappings\u2026"
                : "Preview mappings \u2192"}
        </Button>

        {previewError && (
          <div className="mapping-preview-error" role="alert">
            Mapping preview unavailable \u2014 comparison blocked. {previewError}
          </div>
        )}

        {preview && (
          <div className="mapping-preview">
            <div className="mapping-preview__head">
              <strong>Confirm mappings before running the comparison</strong>
              <small>
                {Object.entries(preview.counts)
                  .map(([key, value]) => `${title(key)}: ${value}`)
                  .join(" \u00b7 ")}
              </small>
            </div>

            {ambiguousRows.length > 0 && (
              <div className="mapping-preview__ambiguous">
                <span>
                  {ambiguousRows.length} term(s) map to multiple candidates \u2014 all candidates are
                  shown and none are auto-selected. Confirm to accept the previewed mappings.
                </span>
                <ul>
                  {ambiguousRows.map((row) => (
                    <li key={`${row.input_row}-${row.input_id}`}>
                      <strong>{row.input_id}</strong>: {row.candidates}
                    </li>
                  ))}
                </ul>
              </div>
            )}

            <Button
              className={`run-comparison ${confirmed ? "run-comparison--confirmed" : ""}`}
              onClick={() => confirmAndRun()}
              disabled={confirmed}
            >
              {confirmed
                ? "Comparison submitted \u2713"
                : requiresConfirmation
                  ? "Confirm mappings & run comparison \u2192"
                  : "Run comparison \u2192"}
            </Button>
          </div>
        )}
      </div>

      <div className="run-note">
        The comparison is blocked until the previewed mappings are confirmed. Availability never
        disables submission; only a structurally invalid query, a missing signature, or an
        unavailable catalog does.
      </div>
    </div>
  );
}
