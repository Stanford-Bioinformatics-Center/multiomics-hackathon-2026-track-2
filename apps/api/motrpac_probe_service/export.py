"""Evidence bundle + a self-contained HTML report, both built from an AnalysisResponse.

Bundle contents (fixed contract):
  input_signature.csv, mapping_audit.csv, feature_evidence.csv, layer_comparison.csv,
  summary.json, analysis_parameters.json, provenance.json, methods_and_limitations.md, report.html

Every row in feature_evidence.csv carries its lineage (evidence_id, source_feature_id, n_collapsed,
aggregation_method) so a plotted point traces to an exported record and then to a store feature.
"""
from __future__ import annotations

import csv
import html
import io
import json
from dataclasses import asdict


def _csv(rows: list[dict], columns: list[str]) -> str:
    buf = io.StringIO()
    w = csv.DictWriter(buf, fieldnames=columns, extrasaction="ignore")
    w.writeheader()
    for r in rows:
        w.writerow(r)
    return buf.getvalue()


def _input_signature_csv(resp) -> str:
    rows = [{"input_row": a.input_row, "input_id": a.input_id, "selected_symbol": a.selected_symbol,
             "status": a.status} for a in resp.mapping_audit]
    return _csv(rows, ["input_row", "input_id", "selected_symbol", "status"])


def _mapping_audit_csv(resp) -> str:
    rows = [asdict(a) for a in resp.mapping_audit]
    return _csv(rows, ["input_row", "input_id", "selected_symbol", "selected_method", "status",
                       "n_candidates", "ambiguous", "candidates"])


def _feature_evidence_csv(resp) -> str:
    rows = [asdict(f) for f in resp.features]
    cols = ["evidence_id", "gene", "disease_direction", "column_id", "comparison_label", "dataset",
            "species", "tissue", "layer", "sex", "timepoint", "exercise_logfc", "exercise_stat",
            "fdr_bh", "opposed", "measured", "mapping_decision_id", "source_feature_id",
            "n_collapsed", "aggregation_method"]
    return _csv(rows, cols)


def _layer_comparison_csv(resp) -> str:
    rows = [asdict(l) for l in resp.layer_discordance]
    cols = ["comparison_label", "tissue", "genes_in_both", "genome_wide_rho", "signature_rho",
            "rna_opposed", "protein_opposed", "both_opposed", "concordant", "rna_only",
            "protein_only", "opposite", "neither"]
    return _csv(rows, cols)


def _summary_json(resp) -> str:
    return json.dumps({
        "run_id": resp.run_id, "schema_version": resp.schema_version, "status": resp.status,
        "signature_name": resp.signature_name, "n_input_rows": resp.n_input_rows,
        "n_counted_genes": resp.n_counted_genes, "n_mapped_rows": resp.n_mapped_rows,
        "returned_layers": resp.returned_layers, "headline": resp.headline,
        "multiplicity_family": asdict(resp.multiplicity_family),
        "contrasts": resp.contrasts,
        "n_columns": len(resp.columns), "n_significant_columns": sum(c.significant for c in resp.columns),
    }, indent=2)


def _analysis_parameters_json(resp) -> str:
    return json.dumps(resp.provenance.get("parameters", {}), indent=2)


def _methods_md(resp) -> str:
    g = "\n".join(f"| {row['safe']} | {row['unsafe']} |" for row in resp.guardrails)
    caveats = "\n".join(f"- {c}" for c in resp.caveats)
    fam = resp.multiplicity_family
    return (f"# Methods and limitations — {resp.signature_name}\n\n"
            f"Run id: `{resp.run_id}`  ·  schema {resp.schema_version}\n\n"
            "## What this compares\n"
            "The supplied source signature is directionally compared with responses in healthy MoTrPAC "
            "cohorts (trained vs sex-matched sedentary rats; acute human bout). This is a cross-cohort "
            "association for hypothesis generation, **not** a disease treatment-effect test. "
            "An exercise-RNA result does not establish a disease RNA baseline.\n\n"
            "## Contrasts (both first-class)\n"
            + "".join(f"- **{ds}** — disease: {c['disease_numerator']} vs {c['disease_denominator']}; "
                      f"exercise: {c['exercise_numerator']} vs {c['exercise_denominator']}\n"
                      for ds, c in resp.contrasts.items())
            + "\n## Multiple testing\n"
            f"cameraPR BH and Bonferroni are applied across the same comparison family of "
            f"{fam.family_size} columns ({fam.n_tests} non-missing tests) BEFORE any tissue/sex/timepoint "
            "filtering. Display filtering never re-runs either correction.\n\n"
            "## Interpretation guardrails\n\n"
            "| Safe statement | Unsafe upgrade |\n|---|---|\n" + g + "\n\n"
            "## Caveats\n" + caveats + "\n\n"
            "PTM signal is not automatically a measure of modification occupancy, enzyme activity, or "
            "functional consequence.\n")


def bundle_files(resp) -> dict[str, str]:
    return {
        "input_signature.csv": _input_signature_csv(resp),
        "mapping_audit.csv": _mapping_audit_csv(resp),
        "feature_evidence.csv": _feature_evidence_csv(resp),
        "layer_comparison.csv": _layer_comparison_csv(resp),
        "summary.json": _summary_json(resp),
        "analysis_parameters.json": _analysis_parameters_json(resp),
        "provenance.json": json.dumps(resp.provenance, indent=2, default=str),
        "methods_and_limitations.md": _methods_md(resp),
        "report.html": render_report_html(resp),
    }


def render_report_html(resp) -> str:
    """A minimal, deterministic HTML report (all prose from computed numbers)."""
    h = resp.headline.get("male_rat_skm_gn_protein_8wk", {})
    rows = "".join(
        f"<tr><td>{html.escape(c.comparison_label)}</td><td>{c.n_opposed}/{c.n_measured}</td>"
        f"<td>{'' if c.camera_t is None else f'{c.camera_t:+.2f}'}</td>"
        f"<td>{'' if c.camera_fdr is None else f'{c.camera_fdr:.4f}'}</td>"
        f"<td>{'' if c.camera_bonferroni is None else f'{c.camera_bonferroni:.4f}'}</td>"
        f"<td>{'yes' if c.significant else 'no'}</td><td>{html.escape(c.verdict)}</td></tr>"
        for c in resp.columns)
    guard = "".join(f"<tr><td>{html.escape(g['safe'])}</td><td>{html.escape(g['unsafe'])}</td></tr>"
                    for g in resp.guardrails)
    return (
        "<!doctype html><html><head><meta charset='utf-8'>"
        f"<title>Exercise Signature Explorer — {html.escape(resp.signature_name)}</title>"
        "<style>body{font-family:system-ui,sans-serif;margin:2rem;max-width:60rem}"
        "table{border-collapse:collapse;width:100%;font-size:13px;margin:1rem 0}"
        "th,td{border:1px solid #ddd;padding:4px 8px;text-align:left}"
        ".warn{background:#fde8e8;padding:.75rem;border-radius:6px}</style></head><body>"
        f"<h1>{html.escape(resp.signature_name)}</h1>"
        f"<p>Run <code>{resp.run_id}</code> · schema {resp.schema_version} · "
        f"{resp.n_counted_genes} counted genes of {resp.n_input_rows} rows.</p>"
        "<div class='warn'><strong>Cross-cohort association for hypothesis generation; "
        "not a disease treatment-effect test.</strong></div>"
        + (f"<h2>Headline</h2><p>{html.escape(h.get('statement',''))}</p>" if h else "")
        + f"<h2>Set-level results ({resp.multiplicity_family.family_size}-column BH family, "
          f"q&lt;{resp.multiplicity_family.threshold})</h2>"
        "<table><tr><th>Comparison</th><th>Opposed/measured</th><th>cameraPR t</th>"
        "<th>BH q (family)</th><th>Bonferroni p (family)</th><th>significant</th><th>verdict</th></tr>" + rows + "</table>"
        "<h2>Interpretation guardrails</h2>"
        "<table><tr><th>Safe statement</th><th>Unsafe upgrade</th></tr>" + guard + "</table>"
        "<p><small>All prose is templated from computed numbers; PTM signal is not modification "
        "occupancy, enzyme activity, or functional consequence.</small></p>"
        "</body></html>")
