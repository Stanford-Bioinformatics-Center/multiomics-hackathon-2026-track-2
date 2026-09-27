import { useEffect, useState } from "react";
import {
  api, type GeneralizedSignature, type GeneralizedQueryResult, ApiError,
} from "../api/client";

/**
 * Generalized query view: run any normalized disease signature through the standalone query_core
 * engine (the generalized Track 2 backend) against the MoTrPAC reference. The bundled PAH examples
 * reproduce the legacy hardcoded results, which is the point of the generalization. The engine is
 * the source of truth; this view only displays what it returns.
 *
 * Three signature sources are offered (R5 AC3), all rendering into the SAME result panel:
 *   - bundled: pick a committed example (GET /generalized/query/{id})
 *   - paste:   paste a query_core-format CSV (POST /generalized/query)
 *   - upload:  read a CSV file's text and submit it (POST /generalized/query)
 *
 * For paste/upload, the backend validates against the query_core schema and returns 422 naming the
 * failing requirement on non-conforming input; that error is surfaced to the analyst. Per R5 AC4, a
 * pasted signature that cannot be parsed is rejected with a visible parse-failure indication and the
 * analyst's pasted text is NOT cleared.
 */

type SourceMode = "bundled" | "paste" | "upload";

// Sensible defaults matching the bundled EE-CON examples (the c2.0 endurance-exercise-vs-control demo).
const DEFAULT_TISSUE = "muscle";
const DEFAULT_REFERENCE_CONTRAST_CATEGORY = "EE-CON";
const DEFAULT_REFERENCE = "motrpac_muscle";

/**
 * Client-side pre-check for an obviously-invalid query_core CSV (R5 AC4). This is intentionally
 * minimal: it rejects empty/whitespace-only text and text that has no data rows beyond a header.
 * The authoritative schema validation lives in the backend, which returns 422 naming the failing
 * requirement — that error is surfaced too. A pre-check keeps a doomed request from clearing state.
 */
function looksParseable(text: string): boolean {
  const lines = text
    .split(/\r?\n/)
    .map((l) => l.trim())
    .filter((l) => l.length > 0);
  // Need at least a header line and one data line.
  return lines.length >= 2;
}

export default function GeneralizedQuery() {
  const [sigs, setSigs] = useState<GeneralizedSignature[] | null>(null);
  const [selected, setSelected] = useState<string>("pah_blood_rna");
  const [result, setResult] = useState<GeneralizedQueryResult | null>(null);
  const [error, setError] = useState<string | null>(null);
  const [running, setRunning] = useState(false);

  // Source mode + upload/paste inputs.
  const [mode, setMode] = useState<SourceMode>("bundled");
  const [pastedText, setPastedText] = useState<string>("");
  const [uploadText, setUploadText] = useState<string>("");
  const [uploadName, setUploadName] = useState<string>("");
  const [tissue, setTissue] = useState<string>(DEFAULT_TISSUE);
  const [refContrast, setRefContrast] = useState<string>(DEFAULT_REFERENCE_CONTRAST_CATEGORY);
  const [reference, setReference] = useState<string>(DEFAULT_REFERENCE);

  useEffect(() => {
    api.generalizedSignatures()
      .then((r) => setSigs(r.signatures))
      .catch((e) => setError(e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e)));
  }, []);

  function errText(e: unknown): string {
    return e instanceof ApiError ? `${e.status}: ${e.detail}` : String(e);
  }

  async function runBundled(id: string) {
    setRunning(true);
    setError(null);
    setResult(null);
    try {
      setResult(await api.generalizedQuery(id));
    } catch (e) {
      setError(errText(e));
    } finally {
      setRunning(false);
    }
  }

  // Submit pasted/uploaded CSV text via POST /generalized/query. Note: on failure we do NOT clear
  // the source text (R5 AC4) — only `setResult(null)` is cleared, never the pasted/uploaded text.
  async function runUpload(signatureCsvText: string) {
    if (!looksParseable(signatureCsvText)) {
      // Client-side parse-failure pre-check: reject WITHOUT clearing the pasted/uploaded text.
      setError("Parse failure: the signature could not be read as a query_core CSV (need a header row and at least one data row). Fix the text and resubmit.");
      setResult(null);
      return;
    }
    setRunning(true);
    setError(null);
    setResult(null);
    try {
      setResult(
        await api.generalizedQueryUpload({
          signature_csv_text: signatureCsvText,
          tissue,
          reference_contrast_category: refContrast,
          reference,
        }),
      );
    } catch (e) {
      // Surface the backend error (422 names the failing schema requirement) WITHOUT clearing text.
      setError(errText(e));
    } finally {
      setRunning(false);
    }
  }

  async function onFileChange(ev: React.ChangeEvent<HTMLInputElement>) {
    const file = ev.target.files?.[0];
    if (!file) return;
    setUploadName(file.name);
    setError(null);
    try {
      const text = await file.text();
      setUploadText(text);
    } catch (e) {
      setError(`Could not read file: ${errText(e)}`);
    }
  }

  function onSubmit() {
    if (mode === "bundled") return runBundled(selected);
    if (mode === "paste") return runUpload(pastedText);
    return runUpload(uploadText);
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

  const submitDisabled =
    running ||
    (mode === "paste" && pastedText.trim().length === 0) ||
    (mode === "upload" && uploadText.trim().length === 0);

  return (
    <section className="view">
      <div className="section-intro">
        <div className="section-intro__number">GEN</div>
        <div>
          <div className="eyebrow">Generalized query_core engine · any signature vs MoTrPAC</div>
          <h2 className="section-title">Generalized disease-signature query</h2>
          <p className="section-copy">
            The generalized Track 2 backend. Supply a signature — pick a bundled example, paste a
            query_core-format CSV, or upload a CSV file — and run it against the MoTrPAC reference.
            The PAH examples reproduce the legacy hardcoded results — that is the generalization
            working.
          </p>
        </div>
      </div>

      {/* Source-mode selector: keep all three modes clear to the analyst. */}
      <div className="query-field">
        <span className="query-field__label">Signature source</span>
        <div className="query-options" role="tablist" aria-label="Signature source">
          <button
            type="button"
            role="tab"
            aria-selected={mode === "bundled"}
            className={`query-chip ${mode === "bundled" ? "query-chip--selected" : ""}`}
            onClick={() => { setMode("bundled"); setError(null); }}
          >
            Bundled example
          </button>
          <button
            type="button"
            role="tab"
            aria-selected={mode === "paste"}
            className={`query-chip ${mode === "paste" ? "query-chip--selected" : ""}`}
            onClick={() => { setMode("paste"); setError(null); }}
          >
            Paste CSV
          </button>
          <button
            type="button"
            role="tab"
            aria-selected={mode === "upload"}
            className={`query-chip ${mode === "upload" ? "query-chip--selected" : ""}`}
            onClick={() => { setMode("upload"); setError(null); }}
          >
            Upload CSV
          </button>
        </div>
      </div>

      {mode === "bundled" && (
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
      )}

      {(mode === "paste" || mode === "upload") && (
        <>
          {mode === "paste" && (
            <div className="query-field">
              <label className="query-field__label" htmlFor="gen-paste">
                Paste query_core-format signature CSV
              </label>
              <textarea
                id="gen-paste"
                className="query-paste"
                rows={8}
                spellCheck={false}
                placeholder={"feature_id,log2fc,q_value\n…"}
                value={pastedText}
                onChange={(e) => setPastedText(e.target.value)}
              />
            </div>
          )}

          {mode === "upload" && (
            <div className="query-field">
              <label className="query-field__label" htmlFor="gen-upload">
                Upload signature CSV
              </label>
              <input id="gen-upload" type="file" accept=".csv,text/csv" onChange={onFileChange} />
              {uploadName && (
                <small className="section-copy">
                  Loaded {uploadName} ({uploadText.length} chars)
                </small>
              )}
            </div>
          )}

          {/* Required query_core inputs for uploaded/pasted signatures (defaults match EE-CON demos). */}
          <div className="query-field">
            <label className="query-field__label" htmlFor="gen-tissue">Tissue</label>
            <input
              id="gen-tissue"
              className="query-paste"
              value={tissue}
              onChange={(e) => setTissue(e.target.value)}
            />
          </div>
          <div className="query-field">
            <label className="query-field__label" htmlFor="gen-refcontrast">Reference contrast category</label>
            <input
              id="gen-refcontrast"
              className="query-paste"
              value={refContrast}
              onChange={(e) => setRefContrast(e.target.value)}
            />
          </div>
          <div className="query-field">
            <label className="query-field__label" htmlFor="gen-reference">Reference</label>
            <input
              id="gen-reference"
              className="query-paste"
              value={reference}
              onChange={(e) => setReference(e.target.value)}
            />
          </div>
        </>
      )}

      <button className="run-comparison" type="button" disabled={submitDisabled} onClick={onSubmit}>
        {running ? "Running query_core…" : "Run query →"}
      </button>

      {error && result === null && (
        <div className="live-state live-state--error" role="alert">
          <p>{error}</p>
        </div>
      )}

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
