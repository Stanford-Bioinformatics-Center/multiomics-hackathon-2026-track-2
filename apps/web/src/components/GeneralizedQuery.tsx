import { useEffect, useState } from "react";
import {
  api, type GeneralizedSignature, type GeneralizedQueryResult, ApiError,
} from "../api/client";

/**
 * Generalized query view: run any normalized disease signature through the standalone query_core
 * engine (the generalized Track 2 backend) against the MoTrPAC reference. The bundled PAH examples
 * reproduce the legacy hardcoded results, which is the point of the generalization. The engine is
 * the source of truth; this view only displays what it returns.
 */
export default function GeneralizedQuery() {
  const [sigs, setSigs] = useState<GeneralizedSignature[] | null>(null);
  const [selected, setSelected] = useState<string>("pah_blood_rna");
  const [result, setResult] = useState<GeneralizedQueryResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [running, setRunning] = useState(false);

  useEffect(() => {
    api.generalizedSignatures()
      .then((r) => setSigs(r.signatures))
      .catch((e) => setError(e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e)));
  }, []);

  async function run(id: string) {
    setRunning(true);
    setError(null);
    setResult(null);
    try {
      setResult(await api.generalizedQuery(id));
    } catch (e) {
      setError(e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e));
    } finally {
      setRunning(false);
    }
  }

  if (error && !sigs) {
    return (
      <section className="view">
        <div className="live-state live-state--error">
          <h3>Backend not reachable</h3>
          <p>{error}</p>
        </div>
      </section>
    );
  }

  const rank = result?.rank_correlation.filter((r) => r.spearman_rho !== null) ?? [];

  return (
    <section className="view">
      <div className="section-intro">
        <div className="section-intro__number">GEN</div>
        <div>
          <div className="eyebrow">Generalized query_core engine · any signature vs MoTrPAC</div>
          <h2 className="section-title">Generalized disease-signature query</h2>
          <p className="section-copy">
            The generalized Track 2 backend. Pick a bundled signature and run it against the MoTrPAC
            reference. The PAH examples reproduce the legacy hardcoded results — that is the
            generalization working.
          </p>
        </div>
      </div>

      <div className="query-field">
        <span className="query-field__label">Bundled signature</span>
        <div className="query-options">
          {(sigs ?? []).map((s) => (
            <button
              type="button"
              key={s.id}
              className={`query-chip ${selected === s.id ? "query-chip--selected" : ""}`}
              onClick={() => setSelected(s.id)}
            >
              {s.label}
            </button>
          ))}
        </div>
      </div>
      <button className="run-comparison" type="button" disabled={running} onClick={() => run(selected)}>
        {running ? "Running query_core…" : "Run query →"}
      </button>

      {error && result === null && <div className="live-state live-state--error"><p>{error}</p></div>}

      {result && (
        <>
          <div className="summary-cards">
            <div className="summary-card"><strong>{result.counts.disease_rows}</strong><span>Disease rows</span><small>{result.signature_id}</small></div>
            <div className="summary-card"><strong>{result.counts.matched_candidate_rows}</strong><span>Matched candidates</span><small>vs {result.reference}</small></div>
            <div className="summary-card"><strong>{result.counts.disease_rows_unmatched_everywhere}</strong><span>Unmatched</span><small>fully audited</small></div>
            <div className="summary-card"><strong>{result.counts.ambiguous_candidate_rows}</strong><span>Ambiguous</span><small>never collapsed</small></div>
          </div>

          {rank.length > 0 && (
            <div className="visual-card">
              <div className="visual-header">
                <div>
                  <div className="eyebrow">Descriptive Spearman rho (disease vs MoTrPAC exercise ranking)</div>
                  <h3>Rank correlation by layer × timepoint</h3>
                </div>
              </div>
              <div className="evidence-table">
                <div className="table-row table-row--head">
                  <span>Layer</span><span>Timepoint</span><span>n used</span><span>Spearman rho</span>
                </div>
                {rank.map((r, i) => (
                  <div className="table-row" key={i}>
                    <span>{r.layer}</span><span>{r.timepoint}</span><span>{r.n_used}</span>
                    <span className={(r.spearman_rho ?? 0) >= 0 ? "good-text" : "warn-text"}>
                      {r.spearman_rho === null ? "—" : r.spearman_rho.toFixed(4)}
                    </span>
                  </div>
                ))}
              </div>
            </div>
          )}

          <div className="assumption-strip">
            <span className="status-mark status-mark--assumption">INTERPRETATION</span>
            <span>{result.interpretation}</span>
          </div>

          <div className="methods-grid">
            <div className="method-block">
              <div className="eyebrow">Provenance</div>
              <h3>Auditable</h3>
              <p>
                disease sha256 {String(result.provenance.disease_sha256).slice(0, 12)}… · reference{" "}
                {String(result.provenance.reference_sha256).slice(0, 12)}… · engine: query_core
                (source of truth)
              </p>
            </div>
          </div>
        </>
      )}
    </section>
  );
}
