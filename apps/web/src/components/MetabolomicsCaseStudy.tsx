import { useEffect, useState } from "react";
import { api, type MetabCaseStudy, type MetabConvergence, ApiError } from "../api/client";

type LoadState =
  | { kind: "loading" }
  | { kind: "error"; message: string }
  | { kind: "ready"; data: MetabCaseStudy; convergence: MetabConvergence | null };

/**
 * Metabolomics case study (ST000763) — a SEPARATE, read-only analysis type. This is a distinct
 * human SSc-PAH plasma cohort, NOT the Malenfant muscle protein signature and NOT the rat exercise
 * data; the view never blends them. It surfaces the standalone module's committed result verbatim,
 * including its honest conclusion (sampling setting / scleroderma, not PAH; null acute response) and
 * the module<->engine convergence cross-check on the MoTrPAC side.
 */
export default function MetabolomicsCaseStudy() {
  const [state, setState] = useState<LoadState>({ kind: "loading" });

  async function load() {
    setState({ kind: "loading" });
    try {
      const data = await api.metabCaseStudy();
      let convergence: MetabConvergence | null = null;
      try {
        convergence = await api.metabConvergence();
      } catch {
        convergence = null; // convergence needs the engine store; case study still renders without it
      }
      setState({ kind: "ready", data, convergence });
    } catch (e) {
      setState({ kind: "error", message: e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e) });
    }
  }

  useEffect(() => {
    load().catch(() => undefined);
  }, []);

  if (state.kind === "loading") {
    return (
      <section className="view">
        <div className="live-state"><p>Loading the metabolomics case study from the API…</p></div>
      </section>
    );
  }
  if (state.kind === "error") {
    return (
      <section className="view">
        <div className="live-state live-state--error">
          <h3>Backend not reachable</h3>
          <p>{state.message}</p>
          <button type="button" onClick={load}>Retry</button>
        </div>
      </section>
    );
  }

  const d = state.data;
  const conv = state.convergence;
  const es = d.null_result.blood_exercise_sensitive ?? {};

  return (
    <section className="view">
      <div className="section-intro">
        <div className="section-intro__number">MET</div>
        <div>
          <div className="eyebrow">Separate analysis type · read-only case study · run {d.run_id}</div>
          <h2 className="section-title">{d.title}</h2>
          <p className="section-copy">{d.cohort_label}</p>
        </div>
      </div>

      {/* distinct-cohort banner so this is never confused with the muscle/rat PAH work */}
      <div className="science-disclaimer">
        <span>SEPARATE COHORT</span> Human SSc-PAH plasma (ST000763) — distinct from the Malenfant
        muscle protein signature and the MoTrPAC rat exercise data.
      </div>

      {/* honest conclusion, prominent */}
      <div className="assumption-strip">
        <span className="status-mark status-mark--assumption">CONCLUSION</span>
        <span>{d.conclusion}</span>
      </div>

      <div className="summary-cards">
        <div className="summary-card"><strong>{d.n_hits}</strong><span>Resting hits (PAH vs healthy)</span><small>q &lt; 0.05, |log2FC| ≥ 0.5</small></div>
        <div className="summary-card"><strong>{d.n_hits_matched_blood}</strong><span>Matched in MoTrPAC blood</span><small>of {d.n_hits}</small></div>
        <div className="summary-card summary-card--focus"><strong>{(es.hits_fraction ?? 0) * 100}%</strong><span>Exercise-sensitive hits</span><small>Fisher p {es.fisher_p}</small></div>
        <div className="summary-card"><strong>{d.median_pah_effect_over_control_drift?.toFixed(1)}×</strong><span>PAH effect vs control drift</span><small>median</small></div>
        {conv && (
          <div className="summary-card summary-card--focus">
            <strong>{conv.identical ? "✓" : "≠"}</strong>
            <span>Module ↔ engine agree</span>
            <small>{conv.unambiguous_matched_cells} cells, max diff {conv.max_abs_diff}</small>
          </div>
        )}
      </div>

      {/* contrasts, first-class; EE-EE like-for-like caveat */}
      <div className="contrast-banner">
        <div><span>Disease numerator</span><strong>{d.contrasts.disease?.numerator}</strong></div>
        <div className="versus">VS</div>
        <div><span>Disease denominator</span><strong>{d.contrasts.disease?.denominator}</strong></div>
        <div className="contrast-divider" />
        <div><span>Same-setting comparator</span><strong>{d.contrasts.disease?.same_setting_comparator}</strong></div>
      </div>
      <div className="workflow-safety">
        <div className="workflow-safety__mark">!</div>
        <div>
          <h3>EE-EE like-for-like caveat</h3>
          <p>{d.contrasts.motrpac_exercise?.ee_ee_caveat}</p>
        </div>
      </div>

      {conv && (
        <div className="deterministic-note">
          <span className="deterministic-note__icon">✓</span>
          <span>
            <strong>Convergence cross-check.</strong> {conv.note} ({conv.n_within_tol}/
            {conv.unambiguous_matched_cells} EE-CON cells identical; {conv.excluded_ambiguous_module_rows}+
            {conv.excluded_ambiguous_engine_rows} ambiguous-key rows excluded.)
          </span>
        </div>
      )}

      <div className="visual-card">
        <div className="visual-header">
          <div>
            <div className="eyebrow">Resting PAH-vs-healthy hits with MoTrPAC context</div>
            <h3>{d.n_hits} metabolites</h3>
          </div>
        </div>
        <div className="evidence-table">
          <div className="table-row table-row--head">
            <span>Metabolite (RefMet)</span><span>Dir</span><span>log2 effect</span><span>q</span>
            <span>MoTrPAC context</span><span>Effect/drift</span>
          </div>
          {d.hits.map((h) => (
            <div className="table-row" key={h.evidence_id}>
              <span>{h.refmet_name}</span>
              <span>{h.pah_direction}</span>
              <span>{h.pah_log2_effect === null ? "—" : h.pah_log2_effect.toFixed(2)}</span>
              <span>{h.pah_q_value === null ? "—" : h.pah_q_value.toExponential(1)}</span>
              <span>{h.context_label}</span>
              <span>{h.pah_effect_over_control_drift === null ? "—" : `${h.pah_effect_over_control_drift.toFixed(1)}×`}</span>
            </div>
          ))}
        </div>
      </div>

      <div className="methods-grid">
        <div className="method-block">
          <div className="eyebrow">Reproduce</div>
          <h3>Case-study bundle</h3>
          <p><a href={api.metabExportUrl()}>Download JSON bundle</a> (summary, hits, conclusion, provenance)</p>
        </div>
        <div className="method-block method-block--warning">
          <div className="eyebrow">Caveats</div>
          <h3>Read before interpreting</h3>
          <ul>
            {d.caveats.map((c, i) => <li key={i}>{c}</li>)}
          </ul>
        </div>
      </div>
    </section>
  );
}
