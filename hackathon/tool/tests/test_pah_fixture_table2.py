"""Validate the PAH muscle signature fixtures against Malenfant 2015 Table 2.

These tests read only committed CSVs (no store, no local hackathon data), so they run
on a fresh clone. They enforce the Gate-2 contract:

* The full-19 fixture (`pah_muscle_malenfant2015.csv`) has 25 rows = 19 counted proteins
  (from Table 2) + 6 phenotype/context rows, with the directions and ratios of Table 2.
* The companion provenance table
  (`pah_muscle_malenfant2015.provenance.csv`) losslessly records paper vs canonical
  identifiers, the ratio, disease_log2fc = log2(ratio), the reported (nominal) p-value
  text, and a normalization/correction note, and reconciles row-for-row with the fixture.
* The lower-9 regression fixture (`pah_muscle_lower9_malenfant2015.csv`) is UNCHANGED:
  exactly the nine downregulated Table-2 proteins with the paper symbols ATP5L / ATP5B.

Rationale (DECISIONS.md ADR-0005): enrich the existing fixture via a companion table;
never fork a competing 9+10 source of truth; preserve paper provenance for integrity.
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
FULL19 = EXAMPLES / "pah_muscle_malenfant2015.csv"
PROV = EXAMPLES / "pah_muscle_malenfant2015.provenance.csv"
LOWER9 = EXAMPLES / "pah_muscle_lower9_malenfant2015.csv"

# Malenfant 2015 J Mol Med, Table 2 (iTRAQ vastus lateralis, 4 IPAH vs 4 controls).
# canonical_gene_symbol -> (pah_to_control_ratio, direction, reported_p_value_text).
# Nine lower proteins (ratio < 1, direction -1) and ten higher proteins (ratio > 1, +1).
TABLE2 = {
    # lower in PAH
    "NDUFA9": (0.71, -1, "0.004"),
    "UQCRC2": (0.74, -1, "0.006"),
    "UQCRC1": (0.74, -1, "0.001"),
    "ATP5MG": (0.70, -1, "0.01"),     # paper symbol ATP5L
    "ATP5F1B": (0.74, -1, "<0.0001"),  # paper symbol ATP5B
    "IDH2": (0.74, -1, "0.05"),
    "OGDH": (0.77, -1, "0.002"),
    "SLC25A4": (0.65, -1, "<0.0001"),
    "ECH1": (0.57, -1, "0.02"),
    # higher in PAH
    "GLO1": (1.39, 1, "0.03"),
    "FBP2": (1.33, 1, "0.0006"),
    "PEBP1": (1.29, 1, "0.007"),
    "ATP2A1": (1.33, 1, "0.0007"),
    "MYH1": (2.42, 1, "<0.0001"),
    "MYH2": (1.26, 1, "<0.0001"),
    "MYH7": (1.26, 1, "<0.0001"),
    "MYLPF": (1.58, 1, "0.03"),
    "APOA1": (1.29, 1, "<0.0001"),
    "PDLIM3": (1.58, 1, "<0.0001"),
}

# The six phenotype/context rows carried in the full-19 fixture (not counted in the stats).
CONTEXT_ROWS = {"CS": -1, "TFAM": -1, "PDHA1": -1, "PDHB": -1, "LDHA": 1, "LDHB": 1}

# Identifier corrections that MUST be recorded in the provenance table.
# canonical_gene_symbol -> (paper_uniprot, canonical_uniprot, must-mention substring in note)
IDENTIFIER_CORRECTIONS = {
    "GLO1": ("Q04760-2", "Q04760", "Q04760"),
    "ATP2A1": ("O14983-2", "O14983", "O14983"),
    "MYH1": ("Q9UKX2", "P12882", "P12882"),
    "MYH2": ("Q9UKX2", "Q9UKX2", None),
}


def test_full19_fixture_shape_and_table2():
    df = pd.read_csv(FULL19)
    assert len(df) == 25, "full-19 fixture must be 19 counted + 6 context = 25 rows"
    counted = df[df.group != "phenotype"]
    context = df[df.group == "phenotype"]
    assert len(counted) == 19
    assert len(context) == 6
    # every counted protein is a Table-2 protein with the right direction and ratio (stored in `weight`)
    assert set(counted.gene_symbol) == set(TABLE2)
    for r in counted.itertuples():
        ratio, direction, _ = TABLE2[r.gene_symbol]
        assert r.direction == direction, f"{r.gene_symbol} direction"
        assert math.isclose(float(r.weight), ratio, abs_tol=1e-9), f"{r.gene_symbol} ratio"
    # context rows match
    assert dict(zip(context.gene_symbol, context.direction)) == CONTEXT_ROWS


def test_provenance_table_reconciles_with_fixture():
    df = pd.read_csv(FULL19)
    prov = pd.read_csv(PROV)
    assert len(prov) == 25
    merged = df.rename(columns={"gene_symbol": "canonical_gene_symbol"}).merge(
        prov, on="canonical_gene_symbol", how="outer", indicator=True)
    assert (merged._merge == "both").all(), "every fixture row must have exactly one provenance row"
    # directions agree
    assert (merged.direction_x == merged.direction_y).all()
    # canonical UniProt agrees with the fixture's `uniprot` wherever the fixture has one
    has_u = merged[merged.uniprot.notna() & (merged.uniprot.astype(str) != "")]
    assert (has_u.uniprot == has_u.canonical_uniprot_accession).all()
    # the fixture ratio (`weight`) equals the provenance ratio for the 19 counted proteins
    has_w = merged[merged.weight.notna()]
    assert len(has_w) == 19
    assert np.allclose(has_w.weight.astype(float), has_w.pah_to_control_ratio.astype(float), atol=1e-9)


def test_provenance_log2fc_and_direction_signs():
    prov = pd.read_csv(PROV)
    counted = prov[prov.pah_to_control_ratio.notna()].copy()
    assert len(counted) == 19
    expected = np.log2(counted.pah_to_control_ratio.astype(float))
    assert np.allclose(counted.disease_log2fc.astype(float), expected, atol=5e-4), \
        "disease_log2fc must equal log2(pah_to_control_ratio)"
    # sign of the log2 fold change must equal the stated direction
    assert (np.sign(counted.disease_log2fc.astype(float)).astype(int) == counted.direction.astype(int)).all()


def test_provenance_reported_pvalues_match_table2_text():
    prov = pd.read_csv(PROV, dtype=str, keep_default_na=False)
    counted = prov[prov.pah_to_control_ratio != ""]
    got = dict(zip(counted.canonical_gene_symbol, counted.reported_p_value_text))
    for gene, (_, _, ptext) in TABLE2.items():
        assert got[gene] == ptext, f"{gene}: p-value text {got[gene]!r} != Table 2 {ptext!r}"
    # censored values are preserved literally as text, not coerced to a number
    assert got["ATP5F1B"] == "<0.0001"
    assert got["MYH1"] == "<0.0001"


def test_provenance_records_identifier_corrections():
    prov = pd.read_csv(PROV, dtype=str, keep_default_na=False).set_index("canonical_gene_symbol")
    for gene, (paper_u, canon_u, note_substr) in IDENTIFIER_CORRECTIONS.items():
        row = prov.loc[gene]
        assert row.paper_uniprot_accession == paper_u, f"{gene}: paper accession preserved"
        assert row.canonical_uniprot_accession == canon_u, f"{gene}: canonical accession"
        if note_substr is not None:
            assert note_substr in row.normalization_note, f"{gene}: correction note must mention {note_substr}"
    # the MYH1 <-> MYH2 collision correction must be explicit
    assert "MYH2" in prov.loc["MYH1"].normalization_note
    assert "CORRECTION" in prov.loc["MYH1"].normalization_note.upper()


def test_lower9_regression_fixture_unchanged():
    """The committed lower-9 regression fixture must remain exactly the nine downregulated
    Table-2 proteins with the ORIGINAL paper symbols ATP5L / ATP5B (not the updated symbols)."""
    df = pd.read_csv(LOWER9)
    assert len(df) == 9
    assert (df.direction == -1).all()
    assert set(df.gene_symbol) == {
        "NDUFA9", "UQCRC2", "UQCRC1", "ATP5L", "ATP5B", "IDH2", "OGDH", "SLC25A4", "ECH1"}
    # paper symbols, not the HGNC-updated ones
    assert "ATP5L" in set(df.gene_symbol) and "ATP5MG" not in set(df.gene_symbol)
    assert "ATP5B" in set(df.gene_symbol) and "ATP5F1B" not in set(df.gene_symbol)


def test_no_conflicting_or_duplicate_canonical_symbols():
    prov = pd.read_csv(PROV)
    assert prov.canonical_gene_symbol.is_unique, "no duplicate canonical symbols in the provenance table"
    # the paper's MYH1/MYH2 collision (both Q9UKX2) must be resolved to distinct canonical accessions
    myo = prov.set_index("canonical_gene_symbol")
    assert myo.loc["MYH1"].canonical_uniprot_accession != myo.loc["MYH2"].canonical_uniprot_accession


def test_sensitivity_views_partition_and_membership():
    """Predefined sensitivity views must be self-consistent with the fixture (planned analysis,
    not post-hoc gene removal): lower9 and upper10 partition full19; fibre is a subset of full19."""
    import json
    views = {v["id"]: v["genes"] for v in json.loads(
        (EXAMPLES / "pah_muscle_sensitivity_views.json").read_text())["views"]}
    full19 = set(views["full19"])
    lower9 = set(views["lower9"])
    upper10 = set(views["upper10"])
    fibre = set(views["fibre_myosin"])
    counted = set(pd.read_csv(FULL19).query("group != 'phenotype'").gene_symbol)

    assert full19 == counted, "full19 view must equal the 19 counted proteins"
    assert len(views["full19"]) == 19 and len(views["lower9"]) == 9 and len(views["upper10"]) == 10
    assert lower9.isdisjoint(upper10), "lower9 and upper10 must not overlap"
    assert lower9 | upper10 == full19, "lower9 and upper10 must partition full19"
    assert fibre.issubset(full19), "fibre subgroup must be a subset of full19"
    # every gene in every view exists in the fixture's counted set
    for gid, genes in views.items():
        assert set(genes).issubset(counted), f"view {gid} has genes not in the fixture"
