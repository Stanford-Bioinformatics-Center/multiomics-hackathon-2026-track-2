"""Metabolomics case study (ST000763), Option-1 adapter: serve the standalone module's COMMITTED
outputs, normalized to the service's evidence/provenance schema. No statistics, no R at request time.

This is a distinct analysis TYPE from the directed-signature runs (run_analysis): it is a read-only,
descriptive plasma case study of a SEPARATE cohort (SSc-PAH plasma, Metabolomics Workbench ST000763),
not the Malenfant muscle protein signature and not the rat exercise data. Its honest conclusion - the
PAH-vs-healthy plasma differences are confounded by sampling setting / scleroderma, and the acute
exercise response is null - is carried verbatim.

Source of truth: the standalone module under `MoTrPAC Hackathon/Metabolomics/` (ADR-0011 framing:
the module is the case study; the engine METAB layer is the reusable scorer; they are complementary).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

_REPO = Path(__file__).resolve().parents[3]
MODULE_DIR = _REPO / "MoTrPAC Hackathon" / "Metabolomics"
OUTPUT_DIR = MODULE_DIR / "output"
HITS_CSV = OUTPUT_DIR / "06_context_table_hits.csv"
SUMMARY_JSON = OUTPUT_DIR / "06_context_summary.json"
MANIFEST_JSON = MODULE_DIR / "data" / "source_manifest.json"

SCHEMA_VERSION = "1.0.0"
CASE_STUDY_ID = "st000763_rest_pah_vs_healthy"

CONTRASTS = {
    "disease": {
        "numerator": "SSc-PAH plasma at rest (catheterization setting)",
        "denominator": "healthy plasma at rest (non-invasive setting)",
        "same_setting_comparator": "normal-pressure SSc at rest (same catheterization setting) - 0 hits",
    },
    "motrpac_exercise": {
        "EE-CON": "(exercise post - pre) - (control post - pre): the exercise effect",
        "EE-EE": "exercise post - pre only: like-for-like with patient peak - rest",
        "CON-CON": "control post - pre only: drift without exercise",
        "ee_ee_caveat": "Under the like-for-like EE-EE contrast, PAH is no more MoTrPAC-like than "
                        "healthy; PAH's apparent extra likeness under EE-CON comes from subtracting "
                        "control-side drift.",
    },
}

CONCLUSION = (
    "In ST000763, no metabolite's acute exercise response differs between SSc-PAH and healthy people "
    "or other SSc groups. At rest PAH differs strongly from healthy, but normal-pressure SSc sampled "
    "in the same catheterization setting shows almost the same differences (PAH vs that group: 0 hits). "
    "MoTrPAC shows these metabolites are largely stable in healthy people during and without exercise, "
    "and the resting differences are several times larger than normal physiological drift. Sampling "
    "setting and/or scleroderma is the leading explanation, NOT PAH. These are cross-study descriptions "
    "from inferred pairs and public labels, not proof of mechanism, and not evidence exercise treats PAH."
)

CAVEATS = [
    "Separate cohort: human SSc-PAH plasma (ST000763), distinct from the Malenfant vastus-lateralis "
    "protein signature and the MoTrPAC rat exercise data - do not blend them into one 'PAH result'.",
    "PAH-vs-healthy always mixes three things: PAH, scleroderma, and sampling setting (all PAH samples "
    "are catheterization-setting; all healthy are non-invasive).",
    "Rest/peak pairs are inferred from group/sex/setting order (no participant IDs): 108 plausible, "
    "unverified pairs.",
    "The acute exercise response is null (0 metabolites at q<0.05 for PAH vs healthy change).",
    "PTM signal is not automatically modification occupancy, enzyme activity, or functional consequence.",
]


@dataclass
class MetabEvidenceRow:
    evidence_id: str
    refmet_name: str
    pah_direction: str
    pah_log2_effect: Optional[float]
    pah_q_value: Optional[float]
    context_label: str
    motrpac_blood_match_basis: str
    blood_ee_con_logfc_20min: Optional[float]
    blood_peak_logfc: Optional[float]
    blood_max_abs_control_logfc: Optional[float]
    pah_effect_over_control_drift: Optional[float]
    exercise_vs_pah_direction: str
    caveat: str


@dataclass
class MetabCaseStudyResponse:
    run_id: str
    schema_version: str
    analysis_type: str
    case_study_id: str
    cohort_label: str
    title: str
    n_hits: int
    n_hits_matched_blood: int
    n_background_matched_blood: int
    hit_labels: dict
    null_result: dict
    median_pah_effect_over_control_drift: Optional[float]
    contrasts: dict
    conclusion: str
    caveats: list
    hits: list
    provenance: dict


def _num(v) -> Optional[float]:
    if v is None or v == "":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(f) or math.isinf(f)) else f


def _sha16(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def available() -> bool:
    return HITS_CSV.is_file() and SUMMARY_JSON.is_file()


def build_case_study() -> MetabCaseStudyResponse:
    if not available():
        raise FileNotFoundError("metabolomics module outputs not found; run the Metabolomics pipeline "
                                "(MoTrPAC Hackathon/Metabolomics/run_pipeline.py)")
    summary = json.loads(SUMMARY_JSON.read_text())
    manifest = json.loads(MANIFEST_JSON.read_text()) if MANIFEST_JSON.is_file() else {}
    lst = summary.get("lists", {}).get(CASE_STUDY_ID, {})

    with HITS_CSV.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    content = b"".join(p.read_bytes() for p in (HITS_CSV, SUMMARY_JSON) if p.is_file())
    run_id = hashlib.sha256(
        content + summary.get("motrpac_package_version", "").encode() + SCHEMA_VERSION.encode()
    ).hexdigest()[:24]

    hits = []
    for r in rows:
        evidence_id = f"{run_id}:{CASE_STUDY_ID}:{r.get('refmet_name','')}"
        hits.append(MetabEvidenceRow(
            evidence_id=evidence_id,
            refmet_name=r.get("refmet_name", ""),
            pah_direction=r.get("pah_direction", ""),
            pah_log2_effect=_num(r.get("pah_log2_effect")),
            pah_q_value=_num(r.get("pah_q_value")),
            context_label=r.get("context_label", ""),
            motrpac_blood_match_basis=r.get("motrpac_blood_match_basis", ""),
            blood_ee_con_logfc_20min=_num(r.get("blood_ee_con_logFC_during_20_min")),
            blood_peak_logfc=_num(r.get("blood_peak_logFC")),
            blood_max_abs_control_logfc=_num(r.get("blood_max_abs_control_logFC")),
            pah_effect_over_control_drift=_num(r.get("pah_effect_over_control_drift")),
            exercise_vs_pah_direction=r.get("exercise_vs_pah_direction", ""),
            caveat=r.get("caveat", ""),
        ))

    provenance = {
        "analysis_type": "metabolomics_case_study",
        "case_study_id": CASE_STUDY_ID,
        "run_id": run_id,
        "schema_version": SCHEMA_VERSION,
        "source_of_truth": "standalone module MoTrPAC Hackathon/Metabolomics/ (Option-1 adapter: "
                           "committed outputs served read-only; not recomputed in the service)",
        "motrpac_package_version": summary.get("motrpac_package_version"),
        "label_rule": summary.get("label_rule"),
        "output_files": {
            "06_context_table_hits.csv": _sha16(HITS_CSV),
            "06_context_summary.json": _sha16(SUMMARY_JSON),
        },
        "source_datasets": manifest,
    }

    return MetabCaseStudyResponse(
        run_id=run_id, schema_version=SCHEMA_VERSION, analysis_type="metabolomics_case_study",
        case_study_id=CASE_STUDY_ID, cohort_label="Human SSc-PAH plasma (Metabolomics Workbench ST000763)",
        title="PAH plasma metabolites in MoTrPAC context (ST000763)",
        n_hits=int(lst.get("n_hits", len(rows))),
        n_hits_matched_blood=int(lst.get("n_hits_matched_blood", 0)),
        n_background_matched_blood=int(lst.get("n_background_matched_blood", 0)),
        hit_labels=lst.get("hit_labels", {}),
        null_result={"blood_exercise_sensitive": lst.get("blood_exercise_sensitive", {}),
                     "blood_control_drift": lst.get("blood_control_drift", {})},
        median_pah_effect_over_control_drift=_num(lst.get("median_pah_effect_over_control_drift_hits")),
        contrasts=CONTRASTS, conclusion=CONCLUSION, caveats=CAVEATS,
        hits=[asdict(h) for h in hits], provenance=provenance)
