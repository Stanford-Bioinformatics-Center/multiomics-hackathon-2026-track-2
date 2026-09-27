import { useEffect, useState } from "react";
import {
  api,
  type DiscordanceCatalogResponse,
  type DiscordanceModelResponse,
  type DiscordanceModelMetricRow,
  type PtmParentResponse,
  ApiError,
} from "../api/client";

/**
 * Discordance — the app's main/default view (Track 2 "Omic Discordance Explained").
 *
 * Organized as a three-stage flow, rendered in order:
 *   Stage 1 — Query builder over the committed demo(s). Live recompute of the catalog is an
 *             explicit out-of-scope follow-on; the selectors present the committed demo outputs.
 *   Stage 2 — Data visualization: class-count chart, model-metrics comparison, PTM-parent timing
 *             (and an optional same-time RNA/protein scatter ONLY if catalog rows are delivered).
 *   Stage 3 — Discordance interpretation: all four class counts (incl. zero), the model-metrics
 *             table (R2/MAE/RMSE for zero/rna_only/temporal), the PTM-parent divergence summary,
 *             the PTM occupancy caveat, and a visible weak-prediction statement (honest framing).
 *
 * The standalone discordance modules are the source of truth; this view serves committed demo
 * outputs verbatim through the read-only adapter and computes no statistics.
 */

type Bundle = {
  catalog: DiscordanceCatalogResponse;
  model: DiscordanceModelResponse | null;
  ptm: PtmParentResponse | null;
  /** Optional per-feature catalog rows for the same-time RNA/protein scatter. Absent by default;
   *  the scatter is omitted (never faked) when these are not delivered by the adapter. */
  catalogRows: CatalogRow[];
};

type LoadState =
  | { kind: "loading" }
  | { kind: "error"; message: string }
  | { kind: "ready"; data: Bundle };

/** Shape of an optional per-feature catalog row, if the adapter ever delivers them. */
interface CatalogRow {
  gene?: string;
  rna_log2fc?: number | null;
  protein_log2fc?: number | null;
  event_class?: string;
}

const CLASS_LABELS: Array<{ key: keyof DiscordanceCatalogResponse["classification_counts"]; label: string }> = [
  { key: "supported_concordant", label: "Supported concordant" },
  { key: "supported_opposite", label: "Supported opposite" },
  { key: "rna_response_protein_equivalent", label: "RNA response · protein equivalent" },
  { key: "indeterminate", label: "Indeterminate" },
];

const MODEL_LABELS: Record<DiscordanceModelMetricRow["model"], string> = {
  zero: "Zero-change baseline",
  rna_only: "RNA-only",
  temporal: "Temporal",
};

const MODEL_ORDER: Array<DiscordanceModelMetricRow["model"]> = ["zero", "rna_only", "temporal"];

function fmt(value: number | null | undefined, digits = 3): string {
  return value === null || value === undefined || Number.isNaN(value) ? "—" : value.toFixed(digits);
}

function errText(e: unknown): string {
  return e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e);
}

/**
 * Sum the four classification classes for a delivered catalog. Used as the empty-state signal:
 * a catalog that reports zero total events (or whose four classes all sum to zero) carries no
 * discordance result to render.
 */
function classCountTotal(catalog: DiscordanceCatalogResponse): number {
  const c = catalog.classification_counts;
  return (
    c.supported_concordant +
    c.supported_opposite +
    c.rna_response_protein_equivalent +
    c.indeterminate
  );
}

/**
 * Empty/error decision for the READY path.
 *
 * A successfully delivered catalog is still treated as an EMPTY result — and therefore renders the
 * empty-state and NONE of the counts / model metrics / PTM summary — when it carries no events:
 *   - `total_events === 0`, or
 *   - all four classification classes sum to zero (an empty catalog).
 *
 * This mirrors the failure path so tasks 4.5 (Property 14) and 4.6 (component tests) can assert a
 * single exclusive branch: on success everything renders; on failure/empty nothing does.
 */
function catalogIsEmpty(catalog: DiscordanceCatalogResponse): boolean {
  return catalog.total_events === 0 || classCountTotal(catalog) === 0;
}

/** The reusable empty/error block. Rendering this block renders NONE of the counts/metrics/PTM. */
function EmptyOrError({
  heading,
  message,
  onRetry,
}: {
  heading: string;
  message: string;
  onRetry?: () => void;
}) {
  return (
    <section className="view">
      <div className="live-state live-state--error" data-testid="discordance-empty-or-error">
        <h3>{heading}</h3>
        <p>{message}</p>
        {onRetry ? (
          <button type="button" onClick={onRetry}>
            Retry
          </button>
        ) : null}
      </div>
    </section>
  );
}

/**
 * Extract the optional per-feature catalog rows if the catalog payload happens to carry them.
 * The typed response does not require them; when absent we return an empty list so Stage 2 omits
 * the RNA/protein scatter rather than fabricating points.
 */
function extractCatalogRows(catalog: DiscordanceCatalogResponse): CatalogRow[] {
  const maybe = (catalog as unknown as { rows?: unknown }).rows;
  if (!Array.isArray(maybe)) return [];
  return (maybe as CatalogRow[]).filter(
    (r) => typeof r?.rna_log2fc === "number" && typeof r?.protein_log2fc === "number",
  );
}

export default function Discordance() {
  const [state, setState] = useState<LoadState>({ kind: "loading" });
  // Stage 1 selection: the committed demo option(s). Currently the adapter serves one committed
  // demo; the selector presents it as the selectable option and is structured to accept more.
  const [selectedDemo, setSelectedDemo] = useState<string | null>(null);

  async function load() {
    setState({ kind: "loading" });
    try {
      // The catalog is the anchor resource for the view. Model and PTM are best-effort: each is a
      // separate committed output that may be independently unavailable; the view still renders the
      // stages it can. (Task 4.3 hardens the full empty/error handling.)
      const catalog = await api.discordanceCatalog();
      let model: DiscordanceModelResponse | null = null;
      let ptm: PtmParentResponse | null = null;
      try {
        model = await api.discordanceModel();
      } catch {
        model = null;
      }
      try {
        ptm = await api.discordancePtmParent();
      } catch {
        ptm = null;
      }
      const catalogRows = extractCatalogRows(catalog);
      setSelectedDemo(catalog.case_study_id);
      setState({ kind: "ready", data: { catalog, model, ptm, catalogRows } });
    } catch (e) {
      setState({ kind: "error", message: errText(e) });
    }
  }

  useEffect(() => {
    load().catch(() => undefined);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  if (state.kind === "loading") {
    return (
      <section className="view">
        <div className="live-state"><p>Loading the discordance case study from the API…</p></div>
      </section>
    );
  }

  // ─────────────────────── Single exclusive empty/error branch ───────────────────────
  // Either the catalog fetch FAILED (error state) or the READY catalog carries NO results
  // (total_events === 0, or all four classes sum to zero). In BOTH cases we render only the
  // empty/error block and NONE of: the four class counts, the coverage counts, the model
  // metrics table/chart, or the PTM-parent summary. There is exactly one return here so a test
  // can assert the exclusivity crisply (tasks 4.5 / 4.6).
  if (state.kind === "error") {
    return (
      <EmptyOrError
        heading="Discordance results unavailable"
        message={state.message}
        onRetry={load}
      />
    );
  }

  // state.kind === "ready" from here. A delivered-but-empty catalog is treated identically to a
  // failure: render only the empty-state and NONE of the counts / metrics / PTM summary.
  if (catalogIsEmpty(state.data.catalog)) {
    return (
      <EmptyOrError
        heading="No discordance results"
        message="The discordance catalog returned no events for this selection, so there are no class counts, model metrics, or PTM-parent summary to display."
        onRetry={load}
      />
    );
  }

  const { catalog, model, ptm, catalogRows } = state.data;
  const counts = catalog.classification_counts;
  const maxClassCount = Math.max(1, ...CLASS_LABELS.map(({ key }) => counts[key]));

  // Order the model metric rows deterministically by (subset, model) for the comparison chart.
  const metrics = model?.metrics ?? [];
  const maxMae = Math.max(1e-9, ...metrics.map((m) => Math.abs(m.mae ?? 0)));

  return (
    <section className="view">
      <div className="section-intro">
        <div className="section-intro__number">DIS</div>
        <div>
          <div className="eyebrow">Track 2 · Omic Discordance Explained · run {catalog.run_id}</div>
          <h2 className="section-title">{catalog.title}</h2>
          <p className="section-copy">{catalog.cohort_label}</p>
        </div>
      </div>

      <div className="science-disclaimer">
        <span>RESEARCH USE</span> Cross-cohort association for hypothesis generation. The discordance
        catalog, predictive model, and PTM-parent audit are served as committed demo outputs.
      </div>

      {/* ───────────────────────── Stage 1 · Query builder ───────────────────────── */}
      <div className="assumption-strip">
        <span className="status-mark status-mark--assumption">STAGE 1 · QUERY</span>
        <span>
          Select a committed discordance demo to explore. Live recompute of the catalog is out of
          scope for this MVP — the query builder reads the committed module outputs verbatim.
        </span>
      </div>

      <div className="query-field">
        <span className="query-field__label">Committed discordance demo</span>
        <div className="query-options">
          <button
            type="button"
            className={`query-chip ${selectedDemo === catalog.case_study_id ? "query-chip--selected" : ""}`}
            onClick={() => setSelectedDemo(catalog.case_study_id)}
          >
            {catalog.tissue} · {catalog.contrast_category}
          </button>
        </div>
      </div>

      {/* ──────────────────────── Stage 2 · Data visualization ───────────────────── */}
      <div className="assumption-strip">
        <span className="status-mark status-mark--assumption">STAGE 2 · VISUALIZE</span>
        <span>Distribution of classified events, model-metrics comparison, and PTM-parent timing.</span>
      </div>

      {/* Class-count chart from classification_counts (horizontal bars) */}
      <div className="visual-card">
        <div className="visual-header">
          <div>
            <div className="eyebrow">Event classification · {catalog.total_events} total events</div>
            <h3>Discordance classes</h3>
          </div>
        </div>
        <div style={{ display: "flex", flexDirection: "column", gap: 8, padding: "6px 18px 20px" }}>
          {CLASS_LABELS.map(({ key, label }) => {
            const value = counts[key];
            const pct = (value / maxClassCount) * 100;
            return (
              <div key={key} style={{ display: "grid", gridTemplateColumns: "220px 1fr 60px", alignItems: "center", gap: 10 }}>
                <span style={{ fontSize: 10, fontWeight: 600 }}>{label}</span>
                <span style={{ display: "block", height: 14, background: "var(--muted)", borderRadius: 4, overflow: "hidden" }}>
                  <span
                    style={{
                      display: "block",
                      height: "100%",
                      width: `${pct}%`,
                      minWidth: value > 0 ? 2 : 0,
                      background: "var(--teal)",
                    }}
                  />
                </span>
                <span style={{ fontSize: 11, fontWeight: 700, textAlign: "right" }}>{value}</span>
              </div>
            );
          })}
        </div>
      </div>

      {/* Model-metrics comparison chart from model.metrics (MAE bars per variant) */}
      {metrics.length > 0 && (
        <div className="visual-card">
          <div className="visual-header">
            <div>
              <div className="eyebrow">Out-of-fold model comparison · mean absolute error</div>
              <h3>Predictive model variants</h3>
            </div>
          </div>
          <div style={{ display: "flex", flexDirection: "column", gap: 8, padding: "6px 18px 20px" }}>
            {metrics.map((m, i) => {
              const pct = (Math.abs(m.mae ?? 0) / maxMae) * 100;
              return (
                <div key={`${m.subset}-${m.model}-${i}`} style={{ display: "grid", gridTemplateColumns: "260px 1fr 70px", alignItems: "center", gap: 10 }}>
                  <span style={{ fontSize: 10, fontWeight: 600 }}>{MODEL_LABELS[m.model]} · {m.subset}</span>
                  <span style={{ display: "block", height: 14, background: "var(--muted)", borderRadius: 4, overflow: "hidden" }}>
                    <span style={{ display: "block", height: "100%", width: `${pct}%`, background: "var(--navy)" }} />
                  </span>
                  <span style={{ fontSize: 11, fontWeight: 700, textAlign: "right" }}>{fmt(m.mae)}</span>
                </div>
              );
            })}
          </div>
        </div>
      )}

      {/* Optional same-time RNA/protein scatter — rendered ONLY if catalog rows are delivered.
          Omitted (never faked) when the adapter does not provide per-feature rows. */}
      {catalogRows.length > 0 && (
        <div className="visual-card">
          <div className="visual-header">
            <div>
              <div className="eyebrow">Same-time RNA vs protein log2FC</div>
              <h3>Per-feature discordance scatter</h3>
            </div>
          </div>
          <div className="plot-shell">
            <div className="plot-y-label">Protein log₂FC</div>
            <div className="plot">
              <span className="axis axis--x" />
              <span className="axis axis--y" />
              {catalogRows.map((r, i) => {
                // Clamp effects to a symmetric ±3 window and map to 0–100% for placement.
                const clamp = (v: number) => Math.max(-3, Math.min(3, v));
                const x = ((clamp(r.rna_log2fc ?? 0) + 3) / 6) * 100;
                const y = (1 - (clamp(r.protein_log2fc ?? 0) + 3) / 6) * 100;
                return (
                  <span
                    key={r.gene ?? i}
                    className="plot-point plot-point--neutral"
                    style={{ left: `${x}%`, top: `${y}%` }}
                  />
                );
              })}
            </div>
            <div className="plot-x-label">RNA log₂FC</div>
          </div>
        </div>
      )}

      {/* PTM-parent timing from ptm.summary_rows */}
      {ptm && ptm.summary_rows.length > 0 && (
        <div className="visual-card">
          <div className="visual-header">
            <div>
              <div className="eyebrow">PTM phosphosite vs parent protein · by timepoint</div>
              <h3>PTM-parent timing</h3>
            </div>
          </div>
          <div className="evidence-table">
            <div className="table-row table-row--head" style={{ gridTemplateColumns: "1fr repeat(4, 1fr)" }}>
              <span>Timepoint</span>
              <span>Matched (site + parent CI)</span>
              <span>Supported phosphosite</span>
              <span>Parent equivalent</span>
              <span>Supported + parent-equivalent</span>
            </div>
            {ptm.summary_rows.map((r, i) => (
              <div className="table-row" key={i} style={{ gridTemplateColumns: "1fr repeat(4, 1fr)" }}>
                <span>{r.timepoint}</span>
                <span>{r.matched_with_both_site_and_parent_ci ?? "—"}</span>
                <span>{r.supported_phosphosite_rows_among_matched ?? "—"}</span>
                <span>{r.parent_equivalent_rows_among_matched ?? "—"}</span>
                <span>{r.supported_phosphosite_and_parent_equivalent_rows ?? "—"}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* ───────────────────── Stage 3 · Discordance interpretation ──────────────── */}
      <div className="assumption-strip">
        <span className="status-mark status-mark--assumption">STAGE 3 · INTERPRET</span>
        <span>Classified event counts, model metrics, PTM-parent divergence, and honest framing.</span>
      </div>

      {/* All four classification classes, including any that are zero */}
      <div className="summary-cards" style={{ gridTemplateColumns: "repeat(4, 1fr)" }}>
        {CLASS_LABELS.map(({ key, label }) => (
          <div className="summary-card" key={key}>
            <strong>{counts[key]}</strong>
            <span>{label}</span>
            <small>events</small>
          </div>
        ))}
      </div>

      {/* Coverage states, kept distinct from the four classification classes */}
      <div className="summary-cards" style={{ gridTemplateColumns: "repeat(2, 1fr)", marginTop: 8 }}>
        <div className="summary-card">
          <strong>{catalog.coverage_counts.no_protein_measurement}</strong>
          <span>No protein measurement</span>
          <small>coverage state (not a class)</small>
        </div>
        <div className="summary-card">
          <strong>{catalog.coverage_counts.no_rna_measurement}</strong>
          <span>No RNA measurement</span>
          <small>coverage state (not a class)</small>
        </div>
      </div>

      {/* Model-metrics table: R2 / MAE / RMSE for the zero, rna_only, and temporal variants */}
      {model && (
        <div className="visual-card">
          <div className="visual-header">
            <div>
              <div className="eyebrow">Out-of-fold predictive model · {model.contrast_category}</div>
              <h3>Model metrics (R² / MAE / RMSE)</h3>
            </div>
          </div>

          {/* Honest framing: a visible weak-prediction statement alongside the metrics (R12.4). */}
          <div className="assumption-strip" style={{ margin: "0 18px 10px" }}>
            <span className="status-mark status-mark--assumption">WEAK PREDICTION</span>
            <span>{model.weak_prediction_note}</span>
          </div>

          {MODEL_ORDER.map((variant) => {
            const rows = metrics.filter((m) => m.model === variant);
            if (rows.length === 0) return null;
            return (
              <div className="evidence-table" key={variant} style={{ marginTop: 8 }}>
                <div className="table-row table-row--head" style={{ gridTemplateColumns: "1.4fr 1fr repeat(3, 0.8fr)" }}>
                  <span>{MODEL_LABELS[variant]}</span>
                  <span>Subset</span>
                  <span>R²</span>
                  <span>MAE</span>
                  <span>RMSE</span>
                </div>
                {rows.map((m, i) => (
                  <div className="table-row" key={`${variant}-${m.subset}-${i}`} style={{ gridTemplateColumns: "1.4fr 1fr repeat(3, 0.8fr)" }}>
                    <span>{m.n_rows === null ? "—" : `${m.n_rows} rows`}</span>
                    <span>{m.subset}</span>
                    <span className={(m.r2 ?? 0) > 0 ? "good-text" : "warn-text"}>{fmt(m.r2)}</span>
                    <span>{fmt(m.mae)}</span>
                    <span>{fmt(m.rmse)}</span>
                  </div>
                ))}
              </div>
            );
          })}
        </div>
      )}

      {/* PTM-parent divergence summary */}
      {ptm && (
        <div className="visual-card">
          <div className="visual-header">
            <div>
              <div className="eyebrow">PTM phosphosite-vs-parent divergence audit</div>
              <h3>PTM-parent divergence summary · {ptm.candidate_count} candidate genes</h3>
            </div>
          </div>
          <div className="evidence-table">
            <div className="table-row table-row--head" style={{ gridTemplateColumns: "1fr repeat(5, 1fr)" }}>
              <span>Timepoint</span>
              <span>Raw phosphosite rows</span>
              <span>After mapping checks</span>
              <span>Matched (single gene/parent)</span>
              <span>Supported phosphosite</span>
              <span>Parent equivalent</span>
            </div>
            {ptm.summary_rows.map((r, i) => (
              <div className="table-row" key={i} style={{ gridTemplateColumns: "1fr repeat(5, 1fr)" }}>
                <span>{r.timepoint}</span>
                <span>{r.raw_phosphosite_rows ?? "—"}</span>
                <span>{r.after_mapping_checks_with_uniprot ?? "—"}</span>
                <span>{r.matched_to_single_gene_uniprot_single_feature_parent ?? "—"}</span>
                <span>{r.supported_phosphosite_rows_among_matched ?? "—"}</span>
                <span>{r.parent_equivalent_rows_among_matched ?? "—"}</span>
              </div>
            ))}
          </div>
        </div>
      )}

      {/* PTM occupancy caveat + catalog limitations as visible text */}
      <div className="methods-grid">
        <div className="method-block method-block--warning">
          <div className="eyebrow">PTM</div>
          <h3>Occupancy caveat</h3>
          <p>{ptm?.occupancy_caveat ?? catalog.occupancy_caveat}</p>
        </div>
        <div className="method-block">
          <div className="eyebrow">Limitations</div>
          <h3>Read before interpreting</h3>
          <ul>
            {catalog.limitations.map((l, i) => <li key={i}>{l}</li>)}
            {(model?.limitations ?? []).map((l, i) => <li key={`m-${i}`}>{l}</li>)}
          </ul>
        </div>
      </div>
    </section>
  );
}
