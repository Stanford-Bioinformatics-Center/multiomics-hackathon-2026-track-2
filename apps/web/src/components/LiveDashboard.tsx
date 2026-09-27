import { useEffect, useState } from "react";
import { api, type AnalysisResponse, ApiError } from "../api/client";
import { getVisualizationModeFromReturnedLayers } from "../domain/analysis";

type LoadState =
  | { kind: "idle" }
  | { kind: "loading" }
  | { kind: "error"; message: string }
  | { kind: "ready"; data: AnalysisResponse };

/**
 * Live results view backed by the API (apps/api -> mprobe). Unlike the static design dashboard, this
 * shows REAL numbers: summary cards, the set-level table (with the frozen family q-values), the
 * headline result, both contrasts, and the guardrails — all from one backend response object, so
 * every displayed value traces to the same run. Visualization mode is derived from the layers the
 * backend RETURNED, never from what was requested.
 */
export default function LiveDashboard() {
  const [state, setState] = useState<LoadState>({ kind: "idle" });
  const [exampleName] = useState("pah_muscle_malenfant2015");

  async function run() {
    setState({ kind: "loading" });
    try {
      const data = await api.runComparison({
        example_name: exampleName,
        target_species: "rat",
        selected_omics: ["transcriptomics", "proteomics"],
        fdr_threshold: 0.05,
      });
      setState({ kind: "ready", data });
    } catch (e) {
      const message = e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e);
      setState({ kind: "error", message });
    }
  }

  useEffect(() => {
    // attempt an automatic run; if the API is not running this surfaces a clear connection state
    run().catch(() => undefined);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, []);

  if (state.kind === "idle" || state.kind === "loading") {
    return (
      <section className="view">
        <div className="live-state">
          <p>Loading live results from the analysis API…</p>
          <button type="button" onClick={run}>Run PAH 19-protein comparison</button>
        </div>
      </section>
    );
  }

  if (state.kind === "error") {
    return (
      <section className="view">
        <div className="live-state live-state--error">
          <h3>Backend not reachable</h3>
          <p>{state.message}</p>
          <p>
            Start the API with <code>uvicorn motrpac_probe_service.app:app</code> from{" "}
            <code>apps/api</code> (needs the mprobe store), then set{" "}
            <code>VITE_API_BASE_URL</code> to its URL.
          </p>
          <button type="button" onClick={run}>Retry</button>
        </div>
      </section>
    );
  }

  const d = state.data;
  const mode = getVisualizationModeFromReturnedLayers(d.returned_layers);
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
    <section className="view">
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
        <div className="evidence-table">
          <div className="table-row table-row--head">
            <span>Comparison</span><span>Opposed / measured</span><span>cameraPR t</span>
            <span>BH q (family)</span><span>Significant</span><span>Verdict</span>
          </div>
          {d.columns.map((c) => (
            <div className="table-row" key={c.column_id}>
              <span>{c.comparison_label}</span>
              <span>{c.n_opposed}/{c.n_measured}</span>
              <span>{c.camera_t === null ? "—" : c.camera_t.toFixed(2)}</span>
              <span>{c.camera_fdr === null ? "—" : c.camera_fdr.toFixed(4)}</span>
              <span className={c.significant ? "good-text" : "warn-text"}>{c.significant ? "yes" : "no"}</span>
              <span>{c.verdict}</span>
            </div>
          ))}
        </div>
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
    </section>
  );
}
