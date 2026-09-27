"""FastAPI app: the transport source of truth (OpenAPI + Pydantic). No statistics here.

Endpoints wrap the service (`run_analysis`, `build_catalog`, `resolve_availability`,
`signature.mapping_audit`). The OpenAPI schema this app publishes is what the TypeScript client is
generated from. "No compatible data" / "zero mapped features" are valid scientific outcomes returned
with HTTP 200; only malformed input is 4xx.

Runs (deterministic run_id) are cached in-process so /comparisons/{run_id}/... can serve features,
provenance, report, and the export bundle without recomputation.
"""
from __future__ import annotations

import io
import os
import json
import zipfile
from dataclasses import asdict
from typing import Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse, StreamingResponse
from pydantic import BaseModel, Field, model_validator

from . import catalog as _catalog
from . import convergence as _convergence
from . import generalized_query as _genq
from . import metab_casestudy as _metab
from . import metab_upload as _metab_upload
from . import story as _story
from . import explorer as _explorer
from . import export as _export
from .schema import AnalysisRequest as _AnalysisRequest
from .schema import SignatureRow as _SignatureRow
from .service import SCHEMA_VERSION, run_analysis

app = FastAPI(title="Exercise Signature Explorer API", version="0.1.0",
              description="Thin transport over the motrpac_probe engine. mprobe is the sole "
                          "scientific engine; this layer performs no statistics.")

# The web dev server (Vite) runs on a different origin than the API, so the browser sends a CORS
# preflight (OPTIONS) before each request. Without CORS middleware those preflights 405 and every
# fetch fails with "Failed to fetch". Allow local dev origins (and any localhost/127.0.0.1 port);
# MPROBE_CORS_ORIGINS overrides for other deployments.
_default_origins = [
    "http://localhost:8443", "http://127.0.0.1:8443",
    "http://localhost:5173", "http://127.0.0.1:5173",
    "http://localhost:3000", "http://127.0.0.1:3000",
]
_env_origins = os.environ.get("MPROBE_CORS_ORIGINS", "")
_allow_origins = [o.strip() for o in _env_origins.split(",") if o.strip()] or _default_origins
app.add_middleware(
    CORSMiddleware,
    allow_origins=_allow_origins,
    allow_origin_regex=r"http://(localhost|127\.0\.0\.1):\d+",
    allow_methods=["*"],
    allow_headers=["*"],
)

# in-process run cache: run_id -> AnalysisResponse (as a dataclass instance)
_RUNS: dict[str, object] = {}


# ---- Pydantic transport models (OpenAPI source of truth) ----------------------------------------

class SignatureRowIn(BaseModel):
    gene_symbol: Optional[str] = None
    uniprot: Optional[str] = None
    ensembl: Optional[str] = None
    rat_symbol: Optional[str] = None
    direction: Optional[str] = None
    group: Optional[str] = None
    weight: Optional[str] = None
    source: Optional[str] = None


class AnalysisRequestIn(BaseModel):
    signature_rows: Optional[list[SignatureRowIn]] = None
    signature_csv_text: Optional[str] = None
    example_name: Optional[str] = None
    signature_name: Optional[str] = None
    target_species: str = Field("rat", pattern="^(rat|human)$")
    selected_omics: list[str] = Field(default_factory=lambda: ["transcriptomics", "proteomics"])
    tissue: Optional[str] = None
    sex: Optional[str] = None
    timepoint: Optional[str] = None
    fdr_threshold: float = Field(0.05, gt=0.0, lt=1.0)
    include_nonsignificant: bool = True

    @model_validator(mode="after")
    def _one_signature_source(self):
        provided = [self.signature_rows is not None, bool(self.signature_csv_text), bool(self.example_name)]
        if sum(provided) != 1:
            raise ValueError("provide exactly one of: signature_rows, signature_csv_text, example_name")
        return self

    def to_dataclass(self) -> _AnalysisRequest:
        rows = None
        if self.signature_rows is not None:
            rows = [_SignatureRow(**r.model_dump()) for r in self.signature_rows]
        return _AnalysisRequest(
            signature_rows=rows, signature_csv_text=self.signature_csv_text, example_name=self.example_name,
            signature_name=self.signature_name, target_species=self.target_species,
            selected_omics=self.selected_omics, tissue=self.tissue, sex=self.sex, timepoint=self.timepoint,
            fdr_threshold=self.fdr_threshold, include_nonsignificant=self.include_nonsignificant)


class AvailabilityRequestIn(BaseModel):
    target_species: str
    selected_omics: list[str]
    tissue: Optional[str] = None
    sex: Optional[str] = None
    timepoint: Optional[str] = None
    contrast_category: Optional[str] = None
    exercise_type: Optional[str] = None


class ValidateRequestIn(BaseModel):
    signature_rows: Optional[list[SignatureRowIn]] = None
    signature_csv_text: Optional[str] = None
    example_name: Optional[str] = None

    @model_validator(mode="after")
    def _one(self):
        if sum([self.signature_rows is not None, bool(self.signature_csv_text), bool(self.example_name)]) != 1:
            raise ValueError("provide exactly one of: signature_rows, signature_csv_text, example_name")
        return self


class MetabolitePreviewRequestIn(BaseModel):
    csv_text: str = Field(min_length=1, max_length=2_000_000)


# ---- helpers ------------------------------------------------------------------------------------

def _run_and_cache(req: AnalysisRequestIn):
    resp = run_analysis(req.to_dataclass())
    _RUNS[resp.run_id] = resp
    return resp


def _mapping_audit_payload(vr: ValidateRequestIn) -> dict:
    """Run only the mapping audit (validate / mappings preview) without a full analysis."""
    import tempfile, hashlib
    from pathlib import Path
    from motrpac_probe import core, signature as _sig
    from motrpac_probe.paths import EXAMPLES

    if vr.example_name:
        p = EXAMPLES / (vr.example_name if vr.example_name.endswith(".csv") else vr.example_name + ".csv")
        if not p.is_file():
            raise HTTPException(status_code=422, detail=f"unknown example: {vr.example_name}")
    else:
        if vr.signature_csv_text is not None:
            text = vr.signature_csv_text
        else:
            cols = ["gene_symbol", "uniprot", "ensembl", "rat_symbol", "direction", "group", "weight", "source"]
            lines = [",".join(cols)]
            for r in vr.signature_rows:
                d = r.model_dump()
                lines.append(",".join("" if d.get(c) is None else str(d[c]).replace(",", " ") for c in cols))
            text = "\n".join(lines) + "\n"
        raw = text.encode("utf-8")
        p = Path(tempfile.gettempdir()) / f"mprobe_val_{hashlib.sha256(raw).hexdigest()[:16]}.csv"
        p.write_bytes(raw)
    try:
        res = _sig.mapping_audit(p, universe=core.get_store().genes)
    except Exception as e:  # malformed CSV (missing id/direction columns, etc.)
        raise HTTPException(status_code=422, detail=f"invalid signature: {e}")
    return {"counts": res["counts"], "rows": res["audit"].to_dict("records")}


# ---- endpoints ----------------------------------------------------------------------------------

@app.get("/api/health")
def health():
    from motrpac_probe import store
    return {"status": "ok", "schema_version": SCHEMA_VERSION, "store_hash": store.store_hash()}


@app.get("/api/catalog")
def get_catalog():
    return _catalog.build_catalog()


@app.get("/api/story", response_model=_story.StoryResponse)
def get_story():
    """Read-only narrative panels, each record traced to a committed source output."""
    try:
        return _story.build_story()
    except FileNotFoundError as error:
        raise HTTPException(status_code=404, detail=str(error))


class ExplorerListIn(BaseModel):
    kind: str
    name: Optional[str] = None
    text: Optional[str] = Field(default=None, max_length=5_000_000)
    example: Optional[str] = None


class ExplorerRequestIn(BaseModel):
    lists: list[ExplorerListIn] = Field(min_length=1, max_length=6)


@app.get("/api/explorer/examples")
def explorer_examples():
    """Neutral example inputs; each carries only a source citation."""
    return {"examples": _explorer.list_examples()}


class ExplorerFileIn(BaseModel):
    filename: str = Field(min_length=1, max_length=300)
    content_base64: str = Field(min_length=1, max_length=14_000_000)


@app.post("/api/explorer/read-file")
def explorer_read_file(req: ExplorerFileIn):
    """Uploaded spreadsheet or text file -> CSV text plus a guessed list kind."""
    import base64
    try:
        return _explorer.read_upload(req.filename, base64.b64decode(req.content_base64))
    except (ValueError, UnicodeDecodeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except Exception as exc:  # unreadable workbook
        raise HTTPException(status_code=422, detail=f"Could not read {req.filename}: {exc}") from exc


@app.post("/api/explorer/analyse")
def explorer_analyse(req: ExplorerRequestIn):
    """Any gene/protein, metabolite or pathway list -> every matching MoTrPAC comparison by layer."""
    try:
        return _explorer.analyse([item.model_dump() for item in req.lists])
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc


@app.post("/api/availability")
def post_availability(req: AvailabilityRequestIn):
    # a scientific empty state is a 200, not an error
    return _catalog.resolve_availability(req.target_species, req.selected_omics, req.tissue, req.sex, req.timepoint,
                                         contrast_category=req.contrast_category, exercise_type=req.exercise_type)


@app.post("/api/signatures/validate")
def post_validate(req: ValidateRequestIn):
    return _mapping_audit_payload(req)


@app.post("/api/mappings/preview")
def post_mappings_preview(req: ValidateRequestIn):
    payload = _mapping_audit_payload(req)
    payload["requires_confirmation"] = any(r["ambiguous"] for r in payload["rows"])
    return payload


@app.post("/api/comparisons")
def post_comparisons(req: AnalysisRequestIn):
    try:
        resp = _run_and_cache(req)
    except FileNotFoundError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    return asdict(resp)   # includes status ("ok" | "no_compatible_data" | "no_mapped_features")


def _get_run(run_id: str):
    resp = _RUNS.get(run_id)
    if resp is None:
        raise HTTPException(status_code=404, detail=f"unknown run_id: {run_id}")
    return resp


@app.get("/api/comparisons/{run_id}")
def get_comparison(run_id: str):
    return asdict(_get_run(run_id))


@app.get("/api/comparisons/{run_id}/features")
def get_features(run_id: str):
    return {"run_id": run_id, "features": [asdict(f) for f in _get_run(run_id).features]}


@app.get("/api/comparisons/{run_id}/provenance")
def get_provenance(run_id: str):
    return _get_run(run_id).provenance


@app.get("/api/comparisons/{run_id}/report", response_class=HTMLResponse)
def get_report(run_id: str):
    return HTMLResponse(_export.render_report_html(_get_run(run_id)))


@app.get("/api/comparisons/{run_id}/export")
def get_export(run_id: str):
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for name, content in _export.bundle_files(_get_run(run_id)).items():
            z.writestr(name, content)
    buf.seek(0)
    return StreamingResponse(buf, media_type="application/zip",
                             headers={"Content-Disposition": f'attachment; filename="evidence_bundle_{run_id}.zip"'})


# ---- metabolomics case study (ST000763): a SEPARATE, read-only analysis type -----------------------

@app.post("/api/metabolomics/preview")
def preview_metabolite_list(request: MetabolitePreviewRequestIn):
    """Match a documented user list to precomputed MoTrPAC metabolite contrasts."""
    try:
        return _metab_upload.preview_metabolite_list(request.csv_text)
    except ValueError as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc
    except FileNotFoundError as exc:
        raise HTTPException(status_code=503, detail="MoTrPAC metabolite export is unavailable") from exc

@app.get("/api/metabolomics/casestudy")
def get_metab_casestudy():
    if not _metab.available():
        raise HTTPException(status_code=404,
                            detail="metabolomics module outputs not present (run the Metabolomics pipeline)")
    return asdict(_metab.build_case_study())


@app.get("/api/metabolomics/casestudy/export")
def get_metab_casestudy_export():
    if not _metab.available():
        raise HTTPException(status_code=404, detail="metabolomics module outputs not present")
    resp = _metab.build_case_study()
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("summary.json", json.dumps({
            "run_id": resp.run_id, "cohort": resp.cohort_label, "n_hits": resp.n_hits,
            "n_hits_matched_blood": resp.n_hits_matched_blood, "hit_labels": resp.hit_labels,
            "null_result": resp.null_result, "contrasts": resp.contrasts}, indent=2))
        z.writestr("hits.json", json.dumps(resp.hits, indent=2))
        z.writestr("conclusion.txt", resp.conclusion + "\n\n" + "\n".join(f"- {c}" for c in resp.caveats))
        z.writestr("provenance.json", json.dumps(resp.provenance, indent=2, default=str))
    buf.seek(0)
    return StreamingResponse(buf, media_type="application/zip",
                             headers={"Content-Disposition":
                                      f'attachment; filename="metabolomics_casestudy_{resp.run_id}.zip"'})


@app.get("/api/metabolomics/convergence")
def get_metab_convergence():
    """Integrity cross-check: module MoTrPAC export vs engine METAB store on shared EE-CON cells."""
    try:
        return _convergence.cross_check()
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))

# ---- generalized query (query_core): THIRD analysis type - run any signature vs MoTrPAC ------------

@app.get("/api/generalized/signatures")
def get_generalized_signatures():
    """Bundled example signatures the generalized query_core engine can run (read-only picklist)."""
    if not _genq.available():
        raise HTTPException(status_code=404, detail="generalized query_core module not present")
    return _genq.list_signatures()


@app.get("/api/generalized/query/{signature_id}")
def run_generalized_query(signature_id: str):
    """Run a bundled signature through query_core against the MoTrPAC reference. The engine is the
    source of truth; this reproduces the legacy hardcoded result for the PAH examples."""
    if not _genq.available():
        raise HTTPException(status_code=404, detail="generalized query_core module not present")
    try:
        return _genq.run_bundled(signature_id)
    except KeyError:
        raise HTTPException(status_code=404, detail=f"unknown signature: {signature_id}")
