"""Catalog & availability derived from the store's columns (never from hardcoded constants).

The React app shipped an in-memory catalog that (a) incorrectly listed proteomics for rat SKM-VL
(the store has SKM-VL RNA only; SKM-GN has RNA+PROT+PHOSPHO) and (b) offered mouse/other source
species the mapper cannot handle. This module reads `store.load_columns()` so the catalog is always
the ground truth, and reports capability status per layer (engine / api / react / demo_validated).
"""
from __future__ import annotations

from motrpac_probe import core

# omic layer (React/UI name) <-> engine store layer code
OMIC_TO_LAYER = {"transcriptomics": "RNA", "proteomics": "PROT", "phosphoproteomics": "PHOSPHO",
                 "metabolomics": "METAB", "plasma_protein_olink": "PROT_OLINK"}
LAYER_TO_OMIC = {v: k for k, v in OMIC_TO_LAYER.items()}

# study design is DERIVED from species (no acute/chronic input).
SPECIES_DESIGN = {"human": ("acute", "Acute exercise study"),
                  "rat": ("chronic", "Chronic exercise-training study")}
DATASET_FOR_SPECIES = {"human": "human_acute", "rat": "rat_train"}

# Source species the mapper actually supports (human symbols/UniProt/Ensembl; rat symbols via RGD).
SUPPORTED_SOURCE_SPECIES = ("human", "rat")
UNSUPPORTED_SOURCE_SPECIES = ("mouse", "other")

# Capability matrix: engine (implemented in mprobe scoring), api (exposed by this service),
# react (surfaced in the UI), demo_validated (has a validated committed demo). Kept explicit so a
# layer that is engine-supported but not UI-exposed is not mislabeled globally "planned".
CAPABILITY_MATRIX = {
    "transcriptomics":     {"engine": True,  "api": True,  "react": True,  "demo_validated": True},
    "proteomics":          {"engine": True,  "api": True,  "react": True,  "demo_validated": True},
    "phosphoproteomics":   {"engine": True,  "api": True,  "react": False, "demo_validated": False},
    "metabolomics":        {"engine": True,  "api": False, "react": False, "demo_validated": False},
    "genomics":            {"engine": False, "api": False, "react": False, "demo_validated": False},
    "epigenomics":         {"engine": False, "api": False, "react": False, "demo_validated": False},
}

# Analysis types the app exposes. The directed-signature path (run_analysis) scores a gene/protein
# signature against the MoTrPAC family. The metabolomics case study is a SEPARATE, read-only analysis
# type: the standalone ST000763 module's committed outputs served via an Option-1 adapter. This is why
# the metabolomics *engine live-scoring* row above stays api/react False while the case study is
# exposed - they are different analyses, not the same one toggled on.
ANALYSIS_TYPES = {
    "directed_signature": {
        "label": "Directed disease signature vs exercise",
        "cohort": "Malenfant vastus-lateralis PAH proteins (muscle) vs MoTrPAC exercise",
        "engine": "motrpac_probe (run_analysis)",
        "api": True, "react": True,
    },
    "metabolomics_case_study": {
        "label": "Metabolomics case study (ST000763)",
        "cohort": "Human SSc-PAH plasma (Metabolomics Workbench ST000763) - SEPARATE cohort",
        "engine": "standalone MoTrPAC Hackathon/Metabolomics module (read-only adapter)",
        "api": True, "react": True,
    },
    "generalized_query": {
        "label": "Generalized query (any signature vs MoTrPAC)",
        "cohort": "Any normalized disease signature (bundled PAH examples reproduce the legacy results)",
        "engine": "standalone MoTrPAC Hackathon/generalized query_core (source of truth)",
        "api": True, "react": True,
    },
    "discordance_catalog": {
        "label": "Discordance catalog + predictive model",
        "cohort": "MoTrPAC muscle EE-CON: RNA vs total protein vs PTM within tissue (committed demo)",
        "engine": "standalone MoTrPAC Hackathon/generalized/discordance module (read-only adapter)",
        "api": True, "react": True,
        "note": "Serves committed demo outputs; live recompute is a documented follow-on, not in this MVP.",
    },
}


def _exercise_columns(S):
    c = S.cols
    return c[c.kind.isin(core.EXERCISE_KINDS)]


def build_catalog(S=None) -> dict:
    """Full catalog from the store: per species, the tissues and, per tissue, the layers/omics,
    sexes, and timepoints that ACTUALLY exist in exercise-kind columns."""
    S = S or core.get_store()
    ex = _exercise_columns(S)
    contexts = []
    for species, ds in DATASET_FOR_SPECIES.items():
        design_id, design_label = SPECIES_DESIGN[species]
        d = ex[ex.dataset == ds]
        tissues = []
        for tissue in sorted(d.tissue.unique()):
            td = d[d.tissue == tissue]
            layers = sorted(td.layer.unique())
            omics = [LAYER_TO_OMIC.get(l, l) for l in layers]
            tissues.append({
                "tissue": tissue,
                "tissue_label": core.TISSUE_WORDS.get(tissue, tissue),
                "layers": layers,
                "omics": omics,
                "sexes": sorted(td.sex.unique()) or ["combined"],
                "timepoints": sorted(td.time.unique()),
            })
        contexts.append({
            "species": species, "dataset": ds, "study_design": design_id,
            "study_design_label": design_label, "tissues": tissues,
        })
    return {
        "contexts": contexts,
        "capability_matrix": CAPABILITY_MATRIX,
        "supported_source_species": list(SUPPORTED_SOURCE_SPECIES),
        "unsupported_source_species": list(UNSUPPORTED_SOURCE_SPECIES),
        "store_hash": __import__("motrpac_probe.store", fromlist=["store_hash"]).store_hash(),
        "analysis_types": ANALYSIS_TYPES,
    }


def resolve_availability(target_species: str, selected_omics: list[str], tissue: str | None = None,
                         sex: str | None = None, timepoint: str | None = None, S=None) -> dict:
    """Availability for a specific request, computed from the store. 'no compatible data' and
    'unsupported omics' are valid scientific outcomes (the API returns them with HTTP 200)."""
    S = S or core.get_store()
    if target_species not in DATASET_FOR_SPECIES:
        return {"status": "no_matching_context", "resolved_context": None, "available_omics": [],
                "unavailable_omics": selected_omics,
                "message": f"No MoTrPAC study for species '{target_species}'.", "suggested_alternatives": []}
    ds = DATASET_FOR_SPECIES[target_species]
    ex = _exercise_columns(S)
    d = ex[ex.dataset == ds]
    design_id, design_label = SPECIES_DESIGN[target_species]

    # engine-implemented layers only (from the capability matrix)
    engine_layers = {OMIC_TO_LAYER[o] for o, cap in CAPABILITY_MATRIX.items()
                     if cap["engine"] and o in OMIC_TO_LAYER}
    requested_layers = {o: OMIC_TO_LAYER.get(o) for o in selected_omics}

    scope = d if tissue is None else d[d.tissue == tissue]
    present_layers = set(scope.layer.unique())

    available_omics, unavailable_omics = [], []
    for omic, layer in requested_layers.items():
        if layer is not None and layer in engine_layers and layer in present_layers:
            available_omics.append(omic)
        else:
            unavailable_omics.append(omic)

    # filter mismatch (tissue/sex/timepoint not present)
    mismatch = ((tissue is not None and tissue not in set(d.tissue)) or
                (sex is not None and sex not in set(scope.sex)) or
                (timepoint is not None and timepoint not in set(scope.time)))
    resolved = {"species": target_species, "dataset": ds, "study_design": design_id,
                "study_design_label": design_label}
    if mismatch:
        return {"status": "no_matching_context", "resolved_context": resolved,
                "available_omics": available_omics, "unavailable_omics": unavailable_omics,
                "message": f"The {design_label.lower()} exists, but no result matches the selected "
                           "tissue, sex, or timepoint.",
                "suggested_alternatives": [{"tissue": t} for t in sorted(d.tissue.unique())][:6]}
    if not available_omics:
        return {"status": "unsupported_omics", "resolved_context": resolved,
                "available_omics": [], "unavailable_omics": unavailable_omics,
                "message": "The selected omic layers are not available for this species/tissue.",
                "suggested_alternatives": [{"omics": sorted({LAYER_TO_OMIC.get(l, l)
                                                             for l in present_layers})}]}
    return {"status": "partially_available" if unavailable_omics else "available",
            "resolved_context": resolved, "available_omics": available_omics,
            "unavailable_omics": unavailable_omics,
            "message": ("Supported layers can run; unavailable layers are excluded and shown."
                        if unavailable_omics else "A compatible MoTrPAC context was found."),
            "suggested_alternatives": []}
