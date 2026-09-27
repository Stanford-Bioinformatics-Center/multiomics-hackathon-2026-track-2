"""Discordance case study, Option-1 adapter: serve the standalone discordance module's COMMITTED
outputs, normalized to the service's JSON-safe schema. No statistics computed at request time.

This mirrors `metab_casestudy.py`: a read-only adapter that reads committed demo outputs produced by
the standalone discordance modules (`build_discordance.py`, `audit_ptm_parent.py`) under
`MoTrPAC Hackathon/generalized/discordance/`. Those modules remain the SOLE source of truth for every
statistic; this adapter only reads, shapes (JSON-safe), and returns. It computes NO statistics.

Track 2 deliverables served here:
  - Concordant/Discordant Event Catalog  (catalog.csv, run_summary.json)      -> build_catalog_summary()
  - Predictive Model of Discordance      (model_metrics.csv, model_predictions.csv) -> build_model()
  - PTM phosphosite-vs-parent audit      (ptm_parent_summary.csv, ptm_parent_candidates.csv) -> build_ptm_parent()

Standing caveat carried verbatim: PTM signal is not automatically a measure of modification occupancy,
enzyme activity, or functional consequence.

Path coupling: `Path(__file__).resolve().parents[3]` climbs
  discordance_casestudy.py -> motrpac_probe_service -> api -> apps -> repo root.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

_REPO = Path(__file__).resolve().parents[3]
DISCORDANCE_DIR = _REPO / "MoTrPAC Hackathon" / "generalized" / "discordance"
DEMO_DIR = DISCORDANCE_DIR / "demo_muscle_ee"
PTM_DIR = DISCORDANCE_DIR / "demo_ptm_parent_ee"

CATALOG_CSV = DEMO_DIR / "catalog.csv"
MODEL_METRICS_CSV = DEMO_DIR / "model_metrics.csv"
MODEL_PREDICTIONS_CSV = DEMO_DIR / "model_predictions.csv"
RUN_SUMMARY_JSON = DEMO_DIR / "run_summary.json"
PTM_SUMMARY_CSV = PTM_DIR / "ptm_parent_summary.csv"
PTM_CANDIDATES_CSV = PTM_DIR / "ptm_parent_candidates.csv"

SCHEMA_VERSION = "1.0.0"
CASE_STUDY_ID = "muscle_ee_discordance"

# The four supported classification classes (R2 AC6). `supported_opposite` is a valid class
# even when zero rows carry it in a given demo; it is always surfaced (zero-filled).
CLASSIFICATION_CLASSES = (
    "supported_concordant",
    "supported_opposite",
    "rna_response_protein_equivalent",
    "indeterminate",
)

# The two coverage states, kept separate from the classification classes.
COVERAGE_STATES = (
    "no_protein_measurement",
    "no_rna_measurement",
)

# The full event-class domain the catalog gate validates against (R2 AC8). Any catalog row whose
# `event_class` is missing, malformed, or not one of these makes the catalog unavailable.
KNOWN_EVENT_CLASSES = frozenset(CLASSIFICATION_CLASSES) | frozenset(COVERAGE_STATES)

# The fixed model grid: {all_mapped, rna_responsive} x {zero, rna_only, temporal} = 6 rows (Property 3).
MODEL_SUBSETS = ("all_mapped", "rna_responsive")
MODEL_VARIANTS = ("zero", "rna_only", "temporal")

WEAK_PREDICTION_NOTE = (
    "The predictive model of discordance is a WEAK prediction reported honestly as a result, not a "
    "success claim: out-of-fold R2 is near zero, so early RNA and PTM signals explain little of the "
    "same-time protein change beyond a zero-change baseline. This is an association result, not a "
    "mechanistic or causal claim."
)

OCCUPANCY_CAVEAT = (
    "PTM signal is not automatically a measure of modification occupancy, enzyme activity, or "
    "functional consequence. Site-minus-parent protein log2 fold-change is descriptive only."
)


def _num(v) -> Optional[float]:
    """Numeric-normalization contract: NaN/Infinity/unparseable/empty -> None."""
    if v is None or v == "":
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(f) or math.isinf(f)) else f


def _int(v) -> Optional[int]:
    """Integer-normalization: parse via float first (handles '5586.0'); NaN/inf/unparseable -> None."""
    f = _num(v)
    return None if f is None else int(f)


def _clean(s) -> str:
    """String-normalization: None -> '' so payloads never carry a non-JSON scalar."""
    return "" if s is None else str(s)


def _sha16(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# Availability guards
# --------------------------------------------------------------------------- #
def catalog_available() -> bool:
    """True iff the catalog + run summary exist AND every catalog row's `event_class` is a member
    of KNOWN_EVENT_CLASSES. A single missing/malformed/out-of-set label makes the catalog
    unavailable (R2 AC8, Property 4). Reads module-level constants at call time so tests can
    repoint DEMO_DIR/CATALOG_CSV/RUN_SUMMARY_JSON via monkeypatch."""
    if not (CATALOG_CSV.is_file() and RUN_SUMMARY_JSON.is_file()):
        return False
    try:
        with CATALOG_CSV.open(newline="") as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames is None or "event_class" not in reader.fieldnames:
                return False
            for row in reader:
                label = row.get("event_class")
                if label is None or label not in KNOWN_EVENT_CLASSES:
                    return False
    except (OSError, csv.Error):
        return False
    return True


def model_available() -> bool:
    return MODEL_METRICS_CSV.is_file()


def ptm_available() -> bool:
    return PTM_SUMMARY_CSV.is_file()


# --------------------------------------------------------------------------- #
# Dataclasses (JSON-safe by construction)
# --------------------------------------------------------------------------- #
@dataclass
class CatalogSummaryResponse:
    run_id: str
    schema_version: str
    analysis_type: str
    case_study_id: str
    cohort_label: str
    title: str
    tissue: str
    contrast_category: str
    total_events: int
    classification_counts: dict
    coverage_counts: dict
    occupancy_caveat: str
    limitations: list
    provenance: dict


@dataclass
class ModelMetricRow:
    subset: str
    model: str
    n_rows: Optional[int]
    n_genes: Optional[int]
    mae: Optional[float]
    rmse: Optional[float]
    r2: Optional[float]


@dataclass
class ModelResponse:
    run_id: str
    schema_version: str
    analysis_type: str
    case_study_id: str
    tissue: str
    contrast_category: str
    metrics: list
    weak_prediction_note: str
    limitations: list
    provenance: dict


@dataclass
class PtmParentSummaryRow:
    tissue: str
    contrast_category: str
    timepoint: str
    raw_phosphosite_rows: Optional[int]
    after_mapping_checks_with_uniprot: Optional[int]
    excluded_missing_or_ambiguous_site_mapping: Optional[int]
    matched_to_single_gene_uniprot_single_feature_parent: Optional[int]
    matched_with_both_site_and_parent_ci: Optional[int]
    supported_phosphosite_rows_among_matched: Optional[int]
    parent_equivalent_rows_among_matched: Optional[int]
    supported_phosphosite_and_parent_equivalent_rows: Optional[int]
    candidate_genes: Optional[int]


@dataclass
class PtmParentResponse:
    run_id: str
    schema_version: str
    analysis_type: str
    case_study_id: str
    summary_rows: list
    candidate_count: int
    occupancy_caveat: str
    provenance: dict


# --------------------------------------------------------------------------- #
# Builders — read committed outputs, compute NO statistics
# --------------------------------------------------------------------------- #
def build_catalog_summary() -> CatalogSummaryResponse:
    if not catalog_available():
        raise FileNotFoundError(
            "discordance catalog outputs not found or event-class domain invalid; run the discordance "
            "pipeline (MoTrPAC Hackathon/generalized/discordance/build_discordance.py)"
        )
    summary = json.loads(RUN_SUMMARY_JSON.read_text())
    catalog_classes = summary.get("catalog_classes", {}) or {}

    # classification_counts: the four supported classes, zero-filled.
    classification_counts = {cls: int(catalog_classes.get(cls, 0) or 0) for cls in CLASSIFICATION_CLASSES}
    # coverage_counts: the two coverage states, kept separate.
    coverage_counts = {st: int(catalog_classes.get(st, 0) or 0) for st in COVERAGE_STATES}
    total_events = sum(classification_counts.values()) + sum(coverage_counts.values())

    run_id = _run_id(CATALOG_CSV, RUN_SUMMARY_JSON)
    provenance = _provenance(
        "discordance_catalog",
        {"catalog.csv": _sha16(CATALOG_CSV), "run_summary.json": _sha16(RUN_SUMMARY_JSON)},
    )

    return CatalogSummaryResponse(
        run_id=run_id,
        schema_version=SCHEMA_VERSION,
        analysis_type="discordance_catalog",
        case_study_id=CASE_STUDY_ID,
        cohort_label="MoTrPAC muscle, endurance-exercise-vs-control (EE-CON)",
        title="Concordant/discordant event catalog (muscle EE-CON)",
        tissue=_clean(summary.get("tissue")),
        contrast_category=_clean(summary.get("contrast_category")),
        total_events=total_events,
        classification_counts=classification_counts,
        coverage_counts=coverage_counts,
        occupancy_caveat=OCCUPANCY_CAVEAT,
        limitations=[_clean(x) for x in (summary.get("limitations") or [])],
        provenance=provenance,
    )


def build_model() -> ModelResponse:
    if not model_available():
        raise FileNotFoundError(
            "discordance model outputs not found; run the discordance pipeline "
            "(MoTrPAC Hackathon/generalized/discordance/build_discordance.py)"
        )
    summary = json.loads(RUN_SUMMARY_JSON.read_text()) if RUN_SUMMARY_JSON.is_file() else {}

    with MODEL_METRICS_CSV.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    metrics = [
        ModelMetricRow(
            subset=_clean(r.get("subset")),
            model=_clean(r.get("model")),
            n_rows=_int(r.get("n_rows")),
            n_genes=_int(r.get("n_genes")),
            mae=_num(r.get("mae")),
            rmse=_num(r.get("rmse")),
            r2=_num(r.get("r2")),
        )
        for r in rows
    ]

    files = {"model_metrics.csv": _sha16(MODEL_METRICS_CSV)}
    if MODEL_PREDICTIONS_CSV.is_file():
        files["model_predictions.csv"] = _sha16(MODEL_PREDICTIONS_CSV)
    run_id = _run_id(MODEL_METRICS_CSV, RUN_SUMMARY_JSON)

    return ModelResponse(
        run_id=run_id,
        schema_version=SCHEMA_VERSION,
        analysis_type="discordance_model",
        case_study_id=CASE_STUDY_ID,
        tissue=_clean(summary.get("tissue")),
        contrast_category=_clean(summary.get("contrast_category")),
        metrics=[asdict(m) for m in metrics],
        weak_prediction_note=WEAK_PREDICTION_NOTE,
        limitations=[_clean(x) for x in (summary.get("limitations") or [])],
        provenance=_provenance("discordance_model", files),
    )


def build_ptm_parent() -> PtmParentResponse:
    if not ptm_available():
        raise FileNotFoundError(
            "PTM phosphosite-vs-parent audit outputs not found; run the audit "
            "(MoTrPAC Hackathon/generalized/discordance/audit_ptm_parent.py)"
        )
    with PTM_SUMMARY_CSV.open(newline="") as fh:
        rows = list(csv.DictReader(fh))

    summary_rows = [
        PtmParentSummaryRow(
            tissue=_clean(r.get("tissue")),
            contrast_category=_clean(r.get("contrast_category")),
            timepoint=_clean(r.get("timepoint")),
            raw_phosphosite_rows=_int(r.get("raw_phosphosite_rows")),
            after_mapping_checks_with_uniprot=_int(r.get("after_mapping_checks_with_uniprot")),
            excluded_missing_or_ambiguous_site_mapping=_int(
                r.get("excluded_missing_or_ambiguous_site_mapping")
            ),
            matched_to_single_gene_uniprot_single_feature_parent=_int(
                r.get("matched_to_single_gene_uniprot_single_feature_parent")
            ),
            matched_with_both_site_and_parent_ci=_int(r.get("matched_with_both_site_and_parent_ci")),
            supported_phosphosite_rows_among_matched=_int(
                r.get("supported_phosphosite_rows_among_matched")
            ),
            parent_equivalent_rows_among_matched=_int(r.get("parent_equivalent_rows_among_matched")),
            supported_phosphosite_and_parent_equivalent_rows=_int(
                r.get("supported_phosphosite_and_parent_equivalent_rows")
            ),
            candidate_genes=_int(r.get("candidate_genes")),
        )
        for r in rows
    ]

    candidate_count = 0
    files = {"ptm_parent_summary.csv": _sha16(PTM_SUMMARY_CSV)}
    if PTM_CANDIDATES_CSV.is_file():
        with PTM_CANDIDATES_CSV.open(newline="") as fh:
            candidate_count = sum(1 for _ in csv.DictReader(fh))
        files["ptm_parent_candidates.csv"] = _sha16(PTM_CANDIDATES_CSV)

    run_id = _run_id(PTM_SUMMARY_CSV, PTM_CANDIDATES_CSV if PTM_CANDIDATES_CSV.is_file() else None)

    return PtmParentResponse(
        run_id=run_id,
        schema_version=SCHEMA_VERSION,
        analysis_type="ptm_parent_audit",
        case_study_id=CASE_STUDY_ID,
        summary_rows=[asdict(r) for r in summary_rows],
        candidate_count=candidate_count,
        occupancy_caveat=OCCUPANCY_CAVEAT,
        provenance=_provenance("ptm_parent_audit", files),
    )


# --------------------------------------------------------------------------- #
# Shared provenance / run_id helpers
# --------------------------------------------------------------------------- #
def _run_id(*paths: Optional[Path]) -> str:
    content = b"".join(p.read_bytes() for p in paths if p is not None and p.is_file())
    return hashlib.sha256(content + SCHEMA_VERSION.encode()).hexdigest()[:24]


def _provenance(analysis_type: str, output_files: dict) -> dict:
    return {
        "analysis_type": analysis_type,
        "case_study_id": CASE_STUDY_ID,
        "schema_version": SCHEMA_VERSION,
        "source_of_truth": (
            "standalone discordance modules MoTrPAC Hackathon/generalized/discordance/ "
            "(Option-1 adapter: committed outputs served read-only; not recomputed in the service)"
        ),
        "output_files": output_files,
    }
