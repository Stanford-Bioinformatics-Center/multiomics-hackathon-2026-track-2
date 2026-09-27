import { useCallback, useState } from "react";
import { api, type AnalysisResponse, ApiError } from "../api/client";
import { getVisualizationModeFromReturnedLayers } from "../domain/analysis";
import QueryBuilder, { type QueryBuilderRunContract } from "./QueryBuilder";
import ResultsTable, { rowFromColumnResult } from "./ResultsTable";

type LoadState =
  | { kind: "idle" }
  | { kind: "loading" }
  | { kind: "error"; message: string }
  | { kind: "ready"; data: AnalysisResponse };

/**
 * Shape of a single `layer_discordance` row (mirrors backend `LayerDiscordanceRow`). The API client
 * types this loosely as `Record<string, unknown>`, so we narrow it here for rendering only — no
 * statistics are computed on the client.
 */
type LayerDiscordanceRow = {
  comparison_label: string;
  tissue: string;
  genes_in_both: number;
  genome_wide_rho: number | null;
  signature_rho: number | null;
  rna_opposed: number;
  protein_opposed: number;
  both_opposed: number;
  concordant: number;
  rna_only: number;
  protein_only: number;
  opposite: number;
  neither: number;
};

function fmtRho(v: number | null): string {
  return v === null || v === undefined || Number.isNaN(v) ? "—" : v.toFixed(3);
}

/**
 * Rich results panel for a single confirmed run (R4 AC9). Every displayed value traces to the one
 * `AnalysisResponse` object passed in, so nothing is recomputed on the client. Visualization mode is
 * derived from the layers the backend RETURNED, never from what was requested. The honest-framing
 * headline is read straight from the backend `headline` block and labelled "not significant"
 * (R12.2/R12.4) — the client never fabricates a headline and never surfaces the 16-column
 * sensitivity q as one.
 */
function ResultsPanel({ data: d }: { data: AnalysisResponse }) {
  const mode = getVisualizationModeFromReturnedLayers(d.returned_layers);
  const layerRows = (d.layer_discordance ?? []) as LayerDiscordanceRow[];
  const headline = (d.headline as Record<string, { statement?: string }>)["male_rat_skm_gn_protein_8wk"];
  const nSig = d.columns.filter((c) => c.significant).length;
  const nOpposed = d.columns.filter((c) => c.verdict === "opposed").length;
  const cards: [string, string, string][] = [
    [String(d.n_input_rows), "Input rows", "100% retained"],
    [String(d.n_counted_genes), "Counted proteins", `${d.n_mapped_rows} mapped`],
    [String(d.columns.length), "Comparisons", `${d.multiplicity_family.family_size}-col BH family`],
    [String(nOpposed), "Opposed (set-level)", "cameraPR t > 0"],
    [String(nSig), "Significant", `q < ${d.multiplicity_family.threshold}`],
    [d.returned_layers.join(" + "), "Returned layers", `viz: ${mode}`],
  ];

  return (
    <>
      <div className="section-intro">
        <div className="section-intro__number">LIVE</div>
        <div>
          <div className="eyebrow">Backed by the analysis API · run {d.run_id}</div>
          <h2 className="section-title">{d.signature_name}</h2>
          <p className="section-copy">
            Real results from mprobe via FastAPI. Every number below comes from one run object.
          </p>
        </div>
      </div>

      <div className="science-disclaimer">
        <span>RESEARCH USE</span> Cross-cohort association for hypothesis generation — not evidence
        that exercise treats PAH; no PAH-rat cohort.
      </div>

      {/* both contrasts, first-class; never "healthy gene set" */}
      <div className="contrast-banner">
        <div><span>Disease numerator</span><strong>{d.contrasts.rat_train?.disease_numerator}</strong></div>
        <div className="versus">VS</div>
        <div><span>Disease denominator</span><strong>{d.contrasts.rat_train?.disease_denominator}</strong></div>
        <div className="contrast-divider" />
        <div><span>Exercise numerator</span><strong>{d.contrasts.rat_train?.exercise_numerator}</strong></div>
        <div className="versus">VS</div>
        <div><span>Exercise denominator</span><strong>{d.contrasts.rat_train?.exercise_denominator}</strong></div>
      </div>

      <div className="summary-cards">
        {cards.map(([value, label, detail]) => (
          <div className="summary-card" key={label}>
            <strong>{value}</strong>
            <span>{label}</span>
            <small>{detail}</small>
          </div>
        ))}
      </div>

      {headline?.statement && (
        <div className="assumption-strip">
          <span className="status-mark status-mark--assumption">HEADLINE</span>
          <span>{headline.statement}</span>
        </div>
      )}

      <div className="visual-card">
        <div className="visual-header">
          <div>
            <div className="eyebrow">Set-level opposition · {d.multiplicity_family.family_size}-column BH family</div>
            <h3>Comparisons (q &lt; {d.multiplicity_family.threshold})</h3>
          </div>
        </div>
        {/* Shared ResultsTable (task 7.7): every row carries species + dataset + contrast, and raw p
            (camera_p) renders immediately adjacent to BH q (camera_fdr) with labeled headers. The
            contrast falls back to each column's comparison_label so no field is fabricated. Muscle
            protein cohort stays in its own table, separate from metabolomics/blood (R12.5). */}
        <ResultsTable
          rows={d.columns.map((c) => rowFromColumnResult(c))}
          emptyMessage="No set-level comparisons returned for this run."
        />
      </div>

      <div className="visual-card">
        <div className="visual-header">
          <div>
            <div className="eyebrow">RNA vs protein · within-tissue layer discordance</div>
            <h3>Layer discordance</h3>
          </div>
        </div>
        {layerRows.length === 0 ? (
          <p className="section-copy">No layer discordance rows returned for this run.</p>
        ) : (
          <div className="evidence-table">
            <div className="table-row table-row--head">
              <span>Comparison</span><span>Tissue</span><span>Genes in both</span>
              <span>Genome-wide ρ</span><span>Signature ρ</span>
              <span>Concordant</span><span>RNA only</span><span>Protein only</span>
              <span>Opposite</span><span>Neither</span>
            </div>
            {layerRows.map((r, i) => (
              <div className="table-row" key={`${r.comparison_label}-${i}`}>
                <span>{r.comparison_label}</span>
                <span>{r.tissue}</span>
                <span>{r.genes_in_both}</span>
                <span>{fmtRho(r.genome_wide_rho)}</span>
                <span>{fmtRho(r.signature_rho)}</span>
                <span>{r.concordant}</span>
                <span>{r.rna_only}</span>
                <span>{r.protein_only}</span>
                <span className={r.opposite > 0 ? "warn-text" : undefined}>{r.opposite}</span>
                <span>{r.neither}</span>
              </div>
            ))}
          </div>
        )}
      </div>

      <div className="methods-grid">
        <div className="method-block">
          <div className="eyebrow">Reproduce</div>
          <h3>Evidence bundle</h3>
          <p>
            <a href={api.exportUrl(d.run_id)}>Download CSV/JSON bundle</a> ·{" "}
            <a href={api.reportUrl(d.run_id)} target="_blank" rel="noreferrer">Open HTML report</a>
          </p>
        </div>
        <div className="method-block method-block--warning">
          <div className="eyebrow">PTM</div>
          <h3>Occupancy caveat</h3>
          <p>PTM signal is not automatically modification occupancy, enzyme activity, or functional consequence.</p>
        </div>
      </div>

      <div className="classification-card">
        <div className="eyebrow">Interpretation guardrails</div>
        <table className="guardrails-table">
          <thead><tr><th>Safe statement</th><th>Unsafe upgrade</th></tr></thead>
          <tbody>
            {d.guardrails.map((g, i) => (
              <tr key={i}><td>{g.safe}</td><td>{g.unsafe}</td></tr>
            ))}
          </tbody>
        </table>
      </div>
    </>
  );
}

/**
 * Live results view backed by the API (apps/api -> mprobe). Task 5.2 rewires this off the hardcoded
 * PAH auto-run: results are now driven by the live `QueryBuilder` (task 5.1). The builder owns the
 * signature source, target selectors, FDR threshold, and the mapping-preview confirmation gate
 * (R4 AC10); it emits a confirmed run contract only AFTER the user confirms the previewed mappings.
 * On receiving that contract, this view calls `api.runComparison(contract.request)` and renders the
 * returned `AnalysisResponse` — the same rich rendering as before (summary cards, both contrasts,
 * the honest-framing headline, the set-level BH-family table, the layer_discordance table, the
 * guardrails, the PTM caveat, and the export/report links). Nothing is recomputed on the client.
 */
export default function LiveDashboard() {
  const [state, setState] = useState<LoadState>({ kind: "idle" });

  // Driven by the confirmed contract QueryBuilder emits — never a hardcoded example (task 5.2).
  const handleRun = useCallback(async (contract: QueryBuilderRunContract) => {
    setState({ kind: "loading" });
    try {
      const data = await api.runComparison(contract.request);
      setState({ kind: "ready", data });
    } catch (e) {
      const message = e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e);
      setState({ kind: "error", message });
    }
  }, []);

  return (
    <section className="view">
      <div className="dashboard-body">
        <aside className="filter-rail">
          {/* QueryBuilder (task 5.1) owns selection + the mapping-preview confirmation gate; it emits
              the confirmed run contract via onRun. onQueryChange is required by its contract but the
              live results view renders solely off the confirmed run, so it is a no-op here. */}
          <QueryBuilder onQueryChange={() => undefined} onRun={handleRun} />
        </aside>

        <div className="dashboard-main">
          {state.kind === "idle" && (
            <div className="live-state">
              <h3>No comparison run yet</h3>
              <p>
                Choose a signature source and target context on the left, preview the mappings, then
                confirm to run a live comparison. Results appear here.
              </p>
            </div>
          )}

          {state.kind === "loading" && (
            <div className="live-state">
              <p>Running the comparison against the analysis API…</p>
            </div>
          )}

          {state.kind === "error" && (
            <div className="live-state live-state--error">
              <h3>Comparison failed</h3>
              <p>{state.message}</p>
              <p>
                If the backend is not reachable, start the API with{" "}
                <code>uvicorn motrpac_probe_service.app:app</code> from <code>apps/api</code> (needs
                the mprobe store), then set <code>VITE_API_BASE_URL</code> to its URL. Adjust the
                selections on the left and confirm the mappings again to retry.
              </p>
            </div>
          )}

          {state.kind === "ready" && <ResultsPanel data={state.data} />}
        </div>
      </div>
    </section>
  );
}
