"""Catalog & availability derived from store.columns; capability matrix; source-species gating."""
import pytest

from motrpac_probe_service import (
    build_catalog, resolve_availability, CAPABILITY_MATRIX,
    SUPPORTED_SOURCE_SPECIES, UNSUPPORTED_SOURCE_SPECIES,
)

pytestmark = pytest.mark.skipif(
    __import__("motrpac_probe.paths", fromlist=["CONTRASTS"]).CONTRASTS.exists() is False,
    reason="no mprobe store: run `mprobe store fetch`")


def _tissue(catalog, species, tissue):
    ctx = next(c for c in catalog["contexts"] if c["species"] == species)
    return next((t for t in ctx["tissues"] if t["tissue"] == tissue), None)


def test_catalog_built_from_store_species_design_derived():
    cat = build_catalog()
    species = {c["species"]: c for c in cat["contexts"]}
    assert species["human"]["study_design"] == "acute"
    assert species["rat"]["study_design"] == "chronic"


def test_skmgn_has_protein_and_ptm_skmvl_lacks_them():
    """The METAB index adds a separate layer, but SKM-VL still lacks protein/PTM."""
    cat = build_catalog()
    skmgn = _tissue(cat, "rat", "SKM-GN")
    skmvl = _tissue(cat, "rat", "SKM-VL")
    assert set(skmgn["layers"]) == {"RNA", "PROT", "PHOSPHO", "METAB"}
    assert set(skmvl["layers"]) == {"RNA", "METAB"}
    assert "proteomics" not in skmvl["omics"]
    assert "proteomics" in skmgn["omics"] and "transcriptomics" in skmgn["omics"]


def test_availability_proteomics_unavailable_for_skmvl():
    res = resolve_availability("rat", ["transcriptomics", "proteomics"], tissue="SKM-VL")
    assert "proteomics" in res["unavailable_omics"]
    assert "transcriptomics" in res["available_omics"]
    assert res["status"] == "partially_available"


def test_availability_full_for_skmgn():
    res = resolve_availability("rat", ["transcriptomics", "proteomics"], tissue="SKM-GN")
    assert res["status"] == "available"
    assert set(res["available_omics"]) == {"transcriptomics", "proteomics"}


def test_catalog_exposes_exact_store_combinations():
    cat = build_catalog()
    combos = cat["combinations"]
    assert any(c["species"] == "human" and c["tissue"] == "BLOOD" and c["layer"] == "RNA"
               and c["contrast_category"] == "EE-CON" and c["timepoint"] == "20m-during" for c in combos)
    assert not any(c["species"] == "human" and c["tissue"] == "BLOOD"
                   and c["exercise_type"] == "resistance" and c["timepoint"] == "20m-during" for c in combos)
    assert any(c["species"] == "rat" and c["contrast_category"] == "TRAIN-SED"
               and c["exercise_type"] == "endurance" for c in combos)
    assert any(c["layer"] == "METAB" for c in combos)


def test_availability_filters_axes_jointly():
    yes = resolve_availability("human", ["transcriptomics"], tissue="BLOOD", timepoint="20m-during",
                               contrast_category="EE-CON", exercise_type="endurance")
    assert yes["status"] == "available"
    no = resolve_availability("human", ["transcriptomics"], tissue="BLOOD", timepoint="20m-during",
                              contrast_category="RE-CON", exercise_type="resistance")
    assert no["status"] == "no_matching_context" and no["available_omics"] == []
    rat_no = resolve_availability("rat", ["proteomics"], tissue="SKM-GN", sex="female", timepoint="8w",
                                  exercise_type="resistance")
    assert rat_no["status"] == "no_matching_context"


def test_measured_metabolomics_distinguished_from_live_api_support():
    res = resolve_availability("human", ["metabolomics"], tissue="BLOOD", timepoint="20m-during",
                               contrast_category="EE-CON")
    assert res["status"] == "unsupported_omics"
    assert "Measured in MoTrPAC" in res["unavailable_reasons"]["metabolomics"]


def test_availability_empty_state_is_not_an_error():
    # a tissue that does not exist -> scientific empty state (not an exception)
    res = resolve_availability("rat", ["transcriptomics"], tissue="NON_EXISTENT_TISSUE")
    assert res["status"] == "no_matching_context"
    assert res["suggested_alternatives"]


def test_capability_matrix_layers():
    # engine supports RNA/PROT/PHOSPHO/METAB; genomics/epigenomics unsupported
    assert CAPABILITY_MATRIX["transcriptomics"]["engine"] and CAPABILITY_MATRIX["transcriptomics"]["react"]
    assert CAPABILITY_MATRIX["phosphoproteomics"]["engine"] and not CAPABILITY_MATRIX["phosphoproteomics"]["react"]
    assert CAPABILITY_MATRIX["metabolomics"]["engine"] and not CAPABILITY_MATRIX["metabolomics"]["api"]
    assert not CAPABILITY_MATRIX["genomics"]["engine"]
    assert not CAPABILITY_MATRIX["epigenomics"]["engine"]


def test_source_species_gating():
    assert set(SUPPORTED_SOURCE_SPECIES) == {"human", "rat"}
    assert set(UNSUPPORTED_SOURCE_SPECIES) == {"mouse", "other"}
    cat = build_catalog()
    assert cat["supported_source_species"] == list(SUPPORTED_SOURCE_SPECIES)


def test_unknown_species_no_matching_context():
    res = resolve_availability("mouse", ["transcriptomics"])
    assert res["status"] == "no_matching_context"
