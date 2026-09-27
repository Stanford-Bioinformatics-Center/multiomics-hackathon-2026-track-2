"""Generalized query adapter: expose the standalone `query_core` engine (the generalized Track 2
backend, MoTrPAC Hackathon/generalized/) through the service, as a THIRD analysis type.

query_core takes any normalized disease-signature CSV and matches it to MoTrPAC contrast summaries.
It is CLI/importable-only in the repo; this adapter runs it in-process into a temp dir and normalizes
its outputs (summary.json + rank_correlation.csv) to a JSON-safe response. No statistics here - the
engine computes everything. The generalized engine (query_core) is the source of truth, exactly as
mprobe is for the directed-signature path and the standalone module is for the metabolomics case study.
"""
from __future__ import annotations

import csv
import importlib
import json
import math
import sys
import tempfile
from pathlib import Path

_REPO = Path(__file__).resolve().parents[3]
GENERALIZED_DIR = _REPO / "MoTrPAC Hackathon" / "generalized"
SIGNATURES_DIR = GENERALIZED_DIR / "examples" / "signatures"
REFERENCE_DEFAULT = GENERALIZED_DIR / "query_core" / "motrpac_reference.csv.gz"
REFERENCE_FULL = GENERALIZED_DIR / "query_core" / "motrpac_reference_full.csv.gz"

SCHEMA_VERSION = "1.0.0"

# Bundled example signatures reconstructed from the legacy hardcoded scripts. Each is the same input
# the corresponding hardcoded pipeline used, so query_core reproduces the legacy result.
BUNDLED_SIGNATURES = {
    "pah_blood_rna": {
        "file": "pah_blood_rna.csv.gz", "tissue": "blood", "reference": "default",
        "label": "PAH blood PBMC RNA (GSE33463) - reproduces the legacy 6 Spearman values",
        "reference_contrast_category": "EE-CON"},
    "pah_muscle_protein": {
        "file": "pah_muscle_protein.csv.gz", "tissue": "muscle", "reference": "full",
        "label": "PAH muscle 9-protein (Malenfant Table 2)",
        "reference_contrast_category": "EE-CON"},
    "cheadle2012_biocarta_ipah": {
        "file": "cheadle2012_biocarta_ipah.csv.gz", "tissue": "blood", "reference": "default",
        "label": "PAH blood BioCarta pathways (Cheadle 2012)",
        "reference_contrast_category": "EE-CON"},
    "st000763_rest_pah_vs_healthy": {
        "file": "st000763_rest_pah_vs_healthy.csv.gz", "tissue": "blood", "reference": "default",
        "label": "ST000763 plasma metabolites, PAH vs healthy (confounded by setting)",
        "reference_contrast_category": "EE-CON"},
    "st000763_rest_pah_vs_normal_pressure": {
        "file": "st000763_rest_pah_vs_normal_pressure.csv.gz", "tissue": "blood", "reference": "default",
        "label": "ST000763 plasma metabolites, PAH vs normal-pressure SSc (same setting)",
        "reference_contrast_category": "EE-CON"},
}


def available() -> bool:
    return (SIGNATURES_DIR.is_dir() and REFERENCE_DEFAULT.is_file()
            and (GENERALIZED_DIR / "query_core" / "engine.py").is_file())


def _run_query(**kwargs):
    """Import and call query_core.engine.run_query (adds the generalized dir to sys.path once)."""
    p = str(GENERALIZED_DIR)
    if p not in sys.path:
        sys.path.insert(0, p)
    engine = importlib.import_module("query_core.engine")
    return engine.run_query(**kwargs)


def _clean(v):
    if v in (None, ""):
        return None
    try:
        f = float(v)
    except (TypeError, ValueError):
        return v
    return None if (math.isnan(f) or math.isinf(f)) else f


def list_signatures() -> dict:
    """The bundled example signatures the UI can offer (read-only picklist)."""
    return {"signatures": [{"id": k, **{kk: vv for kk, vv in v.items() if kk != "file"}}
                           for k, v in BUNDLED_SIGNATURES.items()],
            "note": "Each bundled signature is the same input the corresponding legacy hardcoded "
                    "pipeline used; query_core reproduces the legacy result from it."}


def run_bundled(signature_id: str) -> dict:
    if not available():
        raise FileNotFoundError("generalized query_core module not present")
    if signature_id not in BUNDLED_SIGNATURES:
        raise KeyError(signature_id)
    spec = BUNDLED_SIGNATURES[signature_id]
    disease = SIGNATURES_DIR / spec["file"]
    reference = REFERENCE_FULL if spec["reference"] == "full" else REFERENCE_DEFAULT
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "q"
        summary = _run_query(
            disease_path=disease, reference_path=reference, out_dir=out,
            tissue=spec["tissue"], reference_contrast_category=spec["reference_contrast_category"],
        )
        rank = _read_rank_correlation(out / "rank_correlation.csv")
        coverage = _read_counts(out / "coverage.csv")
    # normalize summary (already plain dict/JSON-ish from the engine) + attach rank + provenance
    return {
        "analysis_type": "generalized_query",
        "schema_version": SCHEMA_VERSION,
        "signature_id": signature_id,
        "signature_label": spec["label"],
        "tissue": spec["tissue"],
        "reference": "motrpac_reference_full" if spec["reference"] == "full" else "motrpac_reference",
        "counts": summary.get("counts", {}),
        "filters": summary.get("filters", {}),
        "thresholds": summary.get("thresholds", {}),
        "interpretation": summary.get("interpretation", ""),
        "rank_correlation": rank,
        "coverage_rows": coverage,
        "provenance": {
            "engine": "query_core (generalized Track 2 backend), source of truth",
            "disease_input": summary.get("disease_input"),
            "disease_sha256": summary.get("disease_sha256"),
            "reference_input": summary.get("reference_input"),
            "reference_sha256": summary.get("reference_sha256"),
            "reference_releases": summary.get("reference_releases"),
        },
    }


def _read_rank_correlation(path: Path) -> list:
    if not path.is_file():
        return []
    rows = []
    with open(path, encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            if not r.get("spearman_rho"):
                continue
            rows.append({
                "layer": r.get("reference_layer", ""),
                "timepoint": r.get("reference_timepoint", ""),
                "contrast_category": r.get("reference_contrast_category", ""),
                "n_used": int(r["n_used"]) if r.get("n_used") else 0,
                "spearman_rho": _clean(r.get("spearman_rho")),
                "status": r.get("status", ""),
            })
    return rows


def _read_counts(path: Path) -> int:
    if not path.is_file():
        return 0
    with open(path, encoding="utf-8-sig", newline="") as fh:
        return sum(1 for _ in csv.DictReader(fh))
