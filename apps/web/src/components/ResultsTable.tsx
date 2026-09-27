import type { ColumnResult } from "../api/client";

/**
 * Shared, reusable results table for the Live, Generalized, Discordance, and Metabolomics views.
 *
 * Invariants enforced BY THE TYPE SYSTEM and the render (R5 AC5-AC9, R11.5, R12.2/12.3/12.5):
 *  - Every rendered value carries species + dataset + contrast. `ResultRow` REQUIRES those three
 *    fields (non-optional) so a caller can never construct a row that renders a bare value.
 *  - `camera_p` (raw, nominal, unadjusted) is rendered in the column IMMEDIATELY ADJACENT to
 *    `camera_fdr` (BH q, adjusted over the frozen 52-column multiplicity family). Each carries a
 *    visible header label stating exactly that (R5 AC5).
 *  - `species` and `dataset` are distinctly labeled columns (R5 AC6).
 *  - Cross-species rows (query feature mapped to a different-species target) show the ortholog-link
 *    relation (R5 AC9).
 *  - The only headline this component surfaces is the q = 0.0584 not-significant directed result.
 *    A 16-column q = 0.0413 value is ONLY ever rendered as a labeled sensitivity analysis, never a
 *    headline (R12.2, R12.3).
 *
 * This component computes NO statistics; it reads values already produced by the engine/service.
 */

/** The frozen canonical BH multiplicity family size for the directed path (R12.1). */
export const FROZEN_FAMILY_SIZE = 52;
/** The sensitivity-only alternative family; NEVER a headline (R12.3). */
export const SENSITIVITY_FAMILY_SIZE = 16;

/** Header-label copy, exported so views and tests can assert on the exact wording (R5 AC5). */
export const CAMERA_FDR_HEADER_NOTE = `BH q — adjusted over the frozen ${FROZEN_FAMILY_SIZE}-column multiplicity family`;
export const CAMERA_P_HEADER_NOTE = "raw p — nominal (unadjusted)";

/**
 * The cross-species ortholog-link relation for a row whose query feature maps to a different-species
 * target. Rendered only when present (R5 AC9). Non-null presence is the signal that this is a
 * cross-species result.
 */
export interface OrthologLink {
  /** The query-side feature/gene symbol (in the query species). */
  querySymbol: string;
  /** The query-side species. */
  querySpecies: string;
  /** The target-side (ortholog) feature/gene symbol. */
  targetSymbol: string;
  /** The target-side species (differs from querySpecies for a true cross-species relation). */
  targetSpecies: string;
  /** How the ortholog relation was established (e.g. "HomoloGene", "UniProt", "manual"). Optional. */
  method?: string;
}

/**
 * The row model. `species`, `dataset`, and `contrast` are REQUIRED (non-optional) by design: the type
 * system prevents a caller from rendering a value without all three (R5 AC7/AC8, R11.5). Statistic
 * fields mirror `ColumnResult` and may be null (engine may not produce a value for a given cell).
 */
export interface ResultRow {
  /** Stable key for the row. */
  id: string;
  /** REQUIRED — distinctly labeled species column (R5 AC6/AC7/AC8, R11.5). */
  species: string;
  /** REQUIRED — distinctly labeled dataset column (R5 AC6/AC7/AC8, R11.5). */
  dataset: string;
  /** REQUIRED — the contrast associated with this value (R5 AC7/AC8, R11.5). */
  contrast: string;
  /** Human-readable label for the feature/gene-set/column this row describes. */
  label: string;
  /** Raw nominal (unadjusted) p-value — the `camera_p` field of `ColumnResult`. */
  camera_p: number | null;
  /** BH q-value adjusted over the frozen 52-column family — the `camera_fdr` field. */
  camera_fdr: number | null;
  /** CAMERA verdict text (e.g. "concordant"/"opposite"), optional. */
  verdict?: string;
  /** Whether this row is significant at the analysis threshold, optional. */
  significant?: boolean;
  /**
   * When present, this row maps a query feature to a DIFFERENT-species target; the ortholog-link
   * relation is rendered (R5 AC9). Absent for same-species rows.
   */
  orthologLink?: OrthologLink;
}

/** A headline surfaced above the table. The component only ever surfaces a NON-significant headline
 *  (the q = 0.0584 directed result). A sensitivity-family value must be passed via `sensitivity`. */
export interface HeadlineResult {
  /** The BH q-value of the headline (the directed path headline is 0.0584, not significant). */
  q: number;
  /** Human-readable label, e.g. "male rat SKM-GN protein, 8 weeks". */
  label: string;
  /** Species/dataset/contrast so even the headline is never a bare value (R11.5). */
  species: string;
  dataset: string;
  contrast: string;
}

/** A labeled sensitivity analysis (e.g. the 16-column q = 0.0413). NEVER surfaced as a headline. */
export interface SensitivityResult {
  /** The alternative-family q-value (e.g. 0.0413). */
  q: number;
  /** The alternative family size (e.g. 16). */
  familySize: number;
  /** Human-readable label for the sensitivity comparison. */
  label: string;
}

export interface ResultsTableProps {
  /** The rows to render. Every row is fully qualified by construction (see `ResultRow`). */
  rows: ResultRow[];
  /** Optional caption/title for the table. */
  caption?: string;
  /**
   * Optional headline. The component labels it "not significant" whenever `q >= threshold` (default
   * 0.05), matching the q = 0.0584 directed result. Surfacing the headline never upgrades it.
   */
  headline?: HeadlineResult;
  /** Significance threshold used only to label the headline. Defaults to 0.05. */
  headlineThreshold?: number;
  /**
   * Optional sensitivity analyses. These are rendered in a clearly labeled "Sensitivity analysis"
   * strip and are NEVER promoted to the headline slot (R12.3).
   */
  sensitivity?: SensitivityResult[];
  /** Text shown when there are no rows. */
  emptyMessage?: string;
}

/** Map an API `ColumnResult` into a fully-qualified `ResultRow`. Requires the caller to pass the
 *  contrast (from the analysis response's `contrasts`/`comparison_label`) so the row is never bare. */
export function rowFromColumnResult(
  col: ColumnResult,
  opts?: { contrast?: string; orthologLink?: OrthologLink },
): ResultRow {
  return {
    id: col.column_id,
    species: col.species,
    dataset: col.dataset,
    contrast: opts?.contrast ?? col.comparison_label,
    label: col.comparison_label,
    camera_p: col.camera_p,
    camera_fdr: col.camera_fdr,
    verdict: col.verdict,
    significant: col.significant,
    orthologLink: opts?.orthologLink,
  };
}

function fmtStat(v: number | null | undefined): string {
  if (v === null || v === undefined || Number.isNaN(v)) return "—";
  // p / q values are small; exponential for tiny values, fixed otherwise.
  return v !== 0 && Math.abs(v) < 1e-3 ? v.toExponential(1) : v.toFixed(4);
}

/** Column order is intentional: species, dataset, contrast, feature, then camera_p IMMEDIATELY
 *  ADJACENT to camera_fdr, then verdict. Keep camera_p and camera_fdr next to each other (R5 AC5). */
const GRID_TEMPLATE = "0.9fr 1.2fr 1.3fr 1.6fr 0.9fr 0.9fr 0.9fr";

export default function ResultsTable({
  rows,
  caption,
  headline,
  headlineThreshold = 0.05,
  sensitivity,
  emptyMessage = "No results to display.",
}: ResultsTableProps) {
  const headlineNotSignificant = headline ? headline.q >= headlineThreshold : false;

  return (
    <div className="results-table" data-testid="results-table">
      {caption && <div className="eyebrow results-table__caption">{caption}</div>}

      {/* Headline slot: only ever the not-significant directed result (q = 0.0584). R12.2/12.3. */}
      {headline && (
        <div
          className="assumption-strip"
          data-testid="results-headline"
          data-significant={headlineNotSignificant ? "false" : "true"}
        >
          <span className="status-mark status-mark--assumption">
            {headlineNotSignificant ? "HEADLINE · NOT SIGNIFICANT" : "HEADLINE"}
          </span>
          <span>
            {headline.label}: BH q = {headline.q.toFixed(4)}{" "}
            {headlineNotSignificant ? "(not significant)" : ""} · {headline.species} ·{" "}
            {headline.dataset} · {headline.contrast}
          </span>
        </div>
      )}

      {/* Sensitivity analyses: always labeled as such, never a headline. R12.3. */}
      {sensitivity && sensitivity.length > 0 && (
        <div className="deterministic-note" data-testid="results-sensitivity">
          <span className="deterministic-note__icon">S</span>
          <span>
            <strong>Sensitivity analysis (not a headline):</strong>{" "}
            {sensitivity
              .map((s) => `${s.label} — BH q = ${s.q.toFixed(4)} over a ${s.familySize}-column family`)
              .join("; ")}
            .
          </span>
        </div>
      )}

      {rows.length === 0 ? (
        <div className="live-state" data-testid="results-empty">
          <p>{emptyMessage}</p>
        </div>
      ) : (
        <div className="evidence-table">
          <div
            className="table-row table-row--head"
            style={{ gridTemplateColumns: GRID_TEMPLATE }}
          >
            <span title="Distinctly labeled species column">Species</span>
            <span title="Distinctly labeled dataset column">Dataset</span>
            <span title="The contrast associated with this value">Contrast</span>
            <span>Feature</span>
            {/* camera_p immediately adjacent to camera_fdr, each with a visible header label. */}
            <span title={CAMERA_P_HEADER_NOTE}>
              raw p<br />
              <small>nominal (unadjusted)</small>
            </span>
            <span title={CAMERA_FDR_HEADER_NOTE}>
              BH q<br />
              <small>adj. over frozen {FROZEN_FAMILY_SIZE}-col family</small>
            </span>
            <span>Verdict</span>
          </div>

          {rows.map((row) => (
            <div className="table-row" style={{ gridTemplateColumns: GRID_TEMPLATE }} key={row.id}>
              {/* species + dataset + contrast render on EVERY row — a value is never bare. */}
              <span data-field="species">{row.species}</span>
              <span data-field="dataset">{row.dataset}</span>
              <span data-field="contrast">
                {row.contrast}
                {row.orthologLink && (
                  <>
                    <br />
                    <span
                      className="results-table__ortholog"
                      data-testid="ortholog-link"
                      title="Cross-species ortholog-link relation"
                    >
                      ortholog: {row.orthologLink.querySymbol} ({row.orthologLink.querySpecies}) →{" "}
                      {row.orthologLink.targetSymbol} ({row.orthologLink.targetSpecies})
                      {row.orthologLink.method ? ` · ${row.orthologLink.method}` : ""}
                    </span>
                  </>
                )}
              </span>
              <span data-field="label">{row.label}</span>
              {/* raw p (camera_p) — nominal/unadjusted — immediately adjacent to BH q. */}
              <span data-field="camera_p">{fmtStat(row.camera_p)}</span>
              {/* BH q (camera_fdr) — adjusted over the frozen 52-column family. */}
              <span data-field="camera_fdr">{fmtStat(row.camera_fdr)}</span>
              <span data-field="verdict">
                {row.verdict ?? "—"}
                {row.significant !== undefined && (
                  <span className={row.significant ? "good-text" : "warn-text"}>
                    {" "}
                    {row.significant ? "(sig)" : "(n.s.)"}
                  </span>
                )}
              </span>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}
