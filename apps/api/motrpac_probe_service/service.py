"""run_analysis(request) -> AnalysisResponse: the single UI-independent entry point.

Wraps motrpac_probe.run.compute(). Performs NO statistics. Returns only JSON-safe values.

Key properties (see DECISIONS.md):
- Content-addressed input: the signature is written to a temp file whose name embeds the sha256 of
  its bytes; the service never accepts an arbitrary server file path from a client. This also fixes
  the compute()+provenance() KeyError('_path') by setting opts['_path'] to that temp file.
- Frozen multiplicity family: the engine applies BH across the full core-column set inside
  compute(); the service reports camera_fdr as-is and NEVER re-runs BH when filtering for display.
- Deterministic run_id: sha256 over signature bytes + analysis parameters + the ordered family
  column IDs + mapping/store/code versions + seed + schema version. Presentation-only options
  (include_nonsignificant, display tissue/sex/timepoint) are excluded from the run_id.
"""
from __future__ import annotations

import hashlib
import json
import math
import tempfile
from dataclasses import asdict
from pathlib import Path

from motrpac_probe import caveats as _caveats
from motrpac_probe import core, legacy, run, signature as _sig, store

from .schema import (
    AnalysisRequest, AnalysisResponse, ColumnResult, ContrastPair, FeatureEvidence,
    LayerDiscordanceRow, MappingAuditRow, MultiplicityFamily,
)

SCHEMA_VERSION = "1.1.0"
CODE_VERSION = "motrpac_probe_service/0.2.0"

# Analysis-affecting request fields (the run_id inputs). Presentation-only fields are excluded.
_ANALYSIS_PARAM_KEYS = ("target_species", "selected_omics", "fdr_threshold")

# Exercise contrast wording per dataset (derived from species; never "healthy gene set").
_CONTRASTS = {
    "human_acute": ("post-exercise", "matched pre-exercise / control"),
    "rat_train": ("trained rat", "sex-matched sedentary rat"),
}
# Comparators below come from the committed fixture descriptions in examples/README.md.
# An upload has no verified source design, so its comparator is deliberately generic.
_EXAMPLE_CONTRASTS = {
    "pah_muscle_malenfant2015": ("idiopathic PAH vastus lateralis (4 patients)",
                                  "control vastus lateralis (4 people)"),
    "pah_muscle_lower9_malenfant2015": ("idiopathic PAH vastus lateralis (4 patients)",
                                         "control vastus lateralis (4 people)"),
    "pah_blood_cheadle2012_eds": ("idiopathic PAH PBMC (30 people)", "healthy PBMC (41 people)"),
    "pah_blood_gse33463_ranked": ("idiopathic PAH PBMC (30 people)", "healthy PBMC (41 people)"),
    "type2_diabetes_muscle_mootha2003": (
        "impaired-glucose-tolerance or type-2-diabetes vastus lateralis, OXPHOS gene set",
        "normal-glucose-tolerance vastus lateralis"),
    "heart_failure_lv_hannenhalli2006": ("failing human left ventricle", "non-failing human left ventricle"),
    "aging_muscle_liu2013": ("older-adult biceps brachii", "young-adult biceps brachii"),
    "septic_shock_blood_cvijanovich2008": ("pediatric septic-shock blood, day 1", "control blood"),
    "hostrup2022_hiit_proteome": ("post-HIIT-training vastus lateralis", "pre-training vastus lateralis"),
}
_UNVERIFIED_CONTRAST = ("uploaded/source signature condition", "source comparator as supplied")

# omic layer name (React/UI) -> engine store layer code
_OMIC_TO_LAYER = {"transcriptomics": "RNA", "proteomics": "PROT", "phosphoproteomics": "PHOSPHO",
                  "metabolomics": "METAB"}


def _clean(v):
    """JSON-safe scalar: NaN/inf -> None; numpy scalar -> python; else as-is."""
    if v is None:
        return None
    if isinstance(v, float):
        return None if (math.isnan(v) or math.isinf(v)) else v
    if hasattr(v, "item"):          # numpy scalar
        v = v.item()
        return None if (isinstance(v, float) and (math.isnan(v) or math.isinf(v))) else v
    return v


def _materialize_signature(request: AnalysisRequest) -> tuple[Path, bytes, str, bool]:
    """Return (path, raw_bytes, name, is_temp). Content-addressed: temp file name embeds the sha256."""
    from motrpac_probe.paths import EXAMPLES
    if request.example_name:
        p = EXAMPLES / (request.example_name if request.example_name.endswith(".csv")
                        else request.example_name + ".csv")
        if not p.is_file():
            raise FileNotFoundError(f"unknown example: {request.example_name}")
        return p, p.read_bytes(), request.signature_name or p.stem, False
    if request.signature_csv_text is not None:
        text = request.signature_csv_text
    elif request.signature_rows:
        cols = ["gene_symbol", "uniprot", "ensembl", "rat_symbol", "direction", "group", "weight", "source"]
        lines = [",".join(cols)]
        for r in request.signature_rows:
            d = asdict(r)
            lines.append(",".join("" if d.get(c) is None else str(d[c]).replace(",", " ") for c in cols))
        text = "\n".join(lines) + "\n"
    else:
        raise ValueError("request must supply signature_rows, signature_csv_text, or example_name")
    raw = text.encode("utf-8")
    digest = hashlib.sha256(raw).hexdigest()
    tmp = Path(tempfile.gettempdir()) / f"mprobe_sig_{digest[:16]}.csv"
    tmp.write_bytes(raw)
    return tmp, raw, request.signature_name or "uploaded_signature", True


def _run_id(sig_bytes: bytes, request: AnalysisRequest, family_ids: list[str],
            mapping_version: str, store_hash: str) -> str:
    params = {k: getattr(request, k) for k in _ANALYSIS_PARAM_KEYS}
    material = json.dumps({
        "signature_sha256": hashlib.sha256(sig_bytes).hexdigest(),
        "params": params,
        "family": family_ids,                 # ordered
        "mapping_version": mapping_version,
        "store_hash": store_hash,
        "code_version": CODE_VERSION,
        "seed": run.DEFAULTS["seed"],
        "schema_version": SCHEMA_VERSION,
    }, sort_keys=True).encode("utf-8")
    return hashlib.sha256(material).hexdigest()[:24]


def _col_meta(S, cid):
    c = S.cols.loc[cid]
    return dict(dataset=str(c.dataset), species=str(c.species), tissue=str(c.tissue), layer=str(c.layer),
                sex=str(c.sex), timepoint=str(c.time), label=core.col_words(c))


def _disease_contrast(request: AnalysisRequest) -> tuple[str, str]:
    if request.example_name:
        name = Path(request.example_name).name.removesuffix(".gz").removesuffix(".csv")
        return _EXAMPLE_CONTRASTS.get(name, _UNVERIFIED_CONTRAST)
    return _UNVERIFIED_CONTRAST


def run_analysis(request: AnalysisRequest) -> AnalysisResponse:
    sig_path, sig_bytes, name, is_temp = _materialize_signature(request)
    try:
        R = run.compute(str(sig_path), name=name, cutoff=request.fdr_threshold, quiet=True)
        # Fix the KeyError('_path'): provenance() reads opts['_path']; compute() never set it.
        R["opts"]["_path"] = str(sig_path)

        S = R["S"]
        family_ids = list(R["core_cols"])
        store_hash = store.store_hash()
        mapping_version = f"idmap@{store_hash}"
        run_id = _run_id(sig_bytes, request, family_ids, mapping_version, store_hash)

        # ---- mapping audit (from the committed engine helper; never silences duplication)
        audit_res = _sig.mapping_audit(sig_path, universe=S.genes)
        audit_rows = [MappingAuditRow(
            input_row=int(r.input_row), input_id=str(r.input_id), selected_symbol=str(r.selected_symbol),
            selected_method=str(r.selected_method), status=str(r.status), n_candidates=int(r.n_candidates),
            ambiguous=bool(r.ambiguous), candidates=str(r.candidates)) for r in audit_res["audit"].itertuples()]

        # ---- set-level column results (camera_fdr already BH over the whole family)
        sc = R["scores"]
        cutoff = request.fdr_threshold
        columns = []
        for r in sc.itertuples():
            m = _col_meta(S, r.column_id)
            q = _clean(r.camera_fdr)
            t = _clean(r.camera_t)
            verdict = run.verdict(r.camera_t, r.camera_fdr, cutoff)
            columns.append(ColumnResult(
                column_id=str(r.column_id), comparison_label=m["label"], dataset=m["dataset"], species=m["species"],
                tissue=m["tissue"], layer=m["layer"], sex=m["sex"], timepoint=m["timepoint"],
                n_measured=int(r.n_measured), n_opposed=int(r.n_opposed), n_same=int(r.n_same),
                camera_t=t, camera_p=_clean(r.camera_p), camera_fdr=q,
                camera_bonferroni=_clean(r.camera_bonferroni), verdict=verdict,
                significant=bool(q is not None and q < cutoff)))

        # ---- per-feature evidence (grid_long: one row per shown gene per column), with full lineage
        # source_feature_id + n_collapsed come from the store row backing each gene in each column
        # (collapse rule: max |stat| when several assay features map to one gene).
        feature_ids = {}          # (column_id, gene) -> store feature_id
        for cid in family_ids:
            fr = S.frame(cid)
            feature_ids.update(
                {(cid, g): fid for g, fid in zip(fr.gene_symbol_human.astype(str), fr.feature_id.astype(str))})
        grid = R["grid"]
        features = []
        for r in grid.itertuples():
            m = _col_meta(S, r.column_id)
            lfc = _clean(r.logFC)
            direction = int(r.direction)
            measured = lfc is not None
            opposed = (None if not measured else bool((1 if lfc > 0 else -1) * direction < 0))
            evidence_id = f"{run_id}:{r.column_id}:{r.gene}"
            n_coll = _clean(r.n_collapsed)
            n_coll = int(n_coll) if n_coll is not None else None
            agg = ("single_feature" if (n_coll or 0) <= 1 else "max_abs_stat") if measured else ""
            features.append(FeatureEvidence(
                evidence_id=evidence_id, gene=str(r.gene), disease_direction=direction, column_id=str(r.column_id),
                comparison_label=m["label"], dataset=m["dataset"], species=m["species"], tissue=m["tissue"],
                layer=m["layer"], sex=m["sex"], timepoint=m["timepoint"], exercise_logfc=lfc,
                exercise_stat=_clean(r.stat), fdr_bh=_clean(r.fdr_bh), opposed=opposed, measured=measured,
                mapping_decision_id=f"{run_id}:map:{r.gene}",
                source_feature_id=(feature_ids.get((str(r.column_id), str(r.gene)), "") if measured else ""),
                n_collapsed=n_coll, aggregation_method=agg))

        # ---- layer discordance (RNA vs protein)
        disc = R.get("disc")
        layer_rows = []
        if disc is not None and len(disc):
            for r in disc.itertuples():
                if getattr(r, "cross", False):
                    continue
                layer_rows.append(LayerDiscordanceRow(
                    comparison_label=str(r.label), tissue=str(S.cols.loc[r.prot, "tissue"]),
                    genes_in_both=int(r.n_genes_both), genome_wide_rho=_clean(r.rho_all),
                    signature_rho=_clean(r.rho_sig), rna_opposed=int(r.rna_opposed),
                    protein_opposed=int(r.prot_opposed), both_opposed=int(r.both_opposed),
                    concordant=int(r.class_concordant), rna_only=int(r.class_RNA_only),
                    protein_only=int(r.class_PROT_only), opposite=int(r.class_opposite),
                    neither=int(r.class_ns)))

        # ---- returned layers (drives viz mode downstream): layers actually present in the columns
        returned_layers = sorted({c.layer for c in columns})

        # ---- multiplicity family
        n_tests = int(sc.camera_p.notna().sum())
        family = MultiplicityFamily(method="BH", methods=["BH", "Bonferroni"], threshold=cutoff,
                                    family_size=len(family_ids),
                                    n_tests=n_tests, column_ids=family_ids)

        # ---- headline (the frozen result), computed, not hardcoded
        headline = _headline(S, sc, cutoff)

        # ---- provenance (engine builds it; opts['_path'] now set)
        prov = run.provenance(R, command="run_analysis")
        prov["multiplicity_family"] = {"method": "BH", "methods": ["BH", "Bonferroni"],
                                       "family_size": len(family_ids), "n_tests": n_tests,
                                       "threshold": cutoff, "column_ids": family_ids}
        prov["run_id"] = run_id
        prov["schema_version"] = SCHEMA_VERSION

        disease_num, disease_den = _disease_contrast(request)
        prov["disease_contrast"] = {"numerator": disease_num, "denominator": disease_den,
                                    "source": "committed example metadata" if request.example_name else
                                              "user-supplied signature; comparator unverified"}
        contrasts = {ds: asdict(ContrastPair(
            disease_numerator=disease_num, disease_denominator=disease_den,
            exercise_numerator=ex[0], exercise_denominator=ex[1])) for ds, ex in _CONTRASTS.items()}

        status = "ok"
        message = "Analysis complete."
        if R["sig"].genes == []:
            status, message = "no_mapped_features", "No signature gene mapped to a measured MoTrPAC feature."
        elif not columns:
            status, message = "no_compatible_data", "No compatible MoTrPAC comparison for this request."

        caveats = list(_caveats.fixed_caveats())
        caveats[0] = ("Direction agreement compares the supplied source signature with a MoTrPAC group "
                      "contrast measured in different people or animals. It is not a treatment-effect test.")

        return AnalysisResponse(
            run_id=run_id, schema_version=SCHEMA_VERSION, status=status, message=message,
            signature_name=str(R["sig"].name), n_input_rows=audit_res["counts"]["n_input_rows"],
            n_counted_genes=len(R["sig"].genes), n_mapped_rows=audit_res["counts"]["n_mapped"],
            returned_layers=returned_layers, contrasts=contrasts, mapping_audit=audit_rows,
            columns=columns, features=features, layer_discordance=layer_rows, multiplicity_family=family,
            headline=headline, caveats=caveats,
            guardrails=[{"safe": s, "unsafe": u} for s, u in _caveats.GUARDRAILS], provenance=prov)
    finally:
        if is_temp:
            try:
                sig_path.unlink()
            except OSError:
                pass


def _headline(S, sc, cutoff) -> dict:
    """Frozen headline: rat gastrocnemius protein, male, 8 wk, under the family BH. Computed."""
    out = {"family_note": "cameraPR BH is applied across the whole comparison family before any "
                          "tissue/sex/timepoint filtering; display filtering never re-runs BH."}
    s = sc.set_index("column_id")
    for cid in s.index:
        c = S.cols.loc[cid]
        if (c.dataset == "rat_train" and c.tissue == "SKM-GN" and c.layer == "PROT"
                and c.sex == "M" and str(c.time).startswith("8")):
            r = s.loc[cid]
            q = _clean(r.camera_fdr)
            q_text = "not estimable" if q is None else f"{q:.3g}"
            out["male_rat_skm_gn_protein_8wk"] = {
                "column_id": cid, "camera_t": _clean(r.camera_t), "camera_p": _clean(r.camera_p),
                "camera_fdr_family": q, "n_opposed": int(r.n_opposed), "n_measured": int(r.n_measured),
                "significant_at_threshold": bool(q is not None and q < cutoff),
                "interpretation": run.verdict(r.camera_t, r.camera_fdr, cutoff),
                "statement": (f"{int(r.n_opposed)}/{int(r.n_measured)} measured signature members have an "
                              "opposite-direction MoTrPAC estimate; "
                              f"set-level BH q={q_text} across {len(s)} columns "
                              f"({'significant' if q is not None and q < cutoff else 'not significant'} "
                              f"at q<{cutoff}). Cross-cohort comparison, not a treatment-effect test."),
            }
    return out
