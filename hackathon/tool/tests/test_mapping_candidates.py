"""Mapping-confirmation gate: candidate exposure and audit row-count reconciliation.

The scientific analysis resolves each identifier to one human symbol (first-wins). The
application's upload workflow must additionally surface EVERY candidate before a
selection is locked in, because 176 Ensembl ids in the released ID map resolve to more
than one human symbol. These tests enforce that ambiguity is surfaced (not silenced) and
that the audit conserves input-row counts (never drop_duplicates() to hide duplication).
"""
import pandas as pd
import pytest

from motrpac_probe import signature


def _write(tmp_path, rows, columns, name="sig.csv"):
    p = tmp_path / name
    pd.DataFrame(rows, columns=columns).to_csv(p, index=False)
    return p


def _first_ambiguous_ensembl(S):
    """An Ensembl id that maps to >1 human symbol, both measured in the store if possible."""
    from motrpac_probe import store
    m = store.load_idmap()
    ens = m[m.id_type == "ensembl"]
    counts = ens.groupby("id").gene_symbol_human.nunique()
    for eid in counts[counts > 1].index:
        yield eid, sorted(ens[ens.id == eid].gene_symbol_human.unique())


def test_candidates_expose_ambiguous_ensembl(S):
    eid, syms = next(_first_ambiguous_ensembl(S))
    cands = signature.candidates_for({"ensembl": eid}, S.genes)
    got = sorted({c["candidate_symbol"] for c in cands})
    # every ID-map candidate for this Ensembl id is exposed (not just the first)
    assert set(syms).issubset(set(got)), f"{eid}: expected {syms}, exposed {got}"
    assert len({c["candidate_symbol"] for c in cands}) > 1


def test_audit_flags_ambiguous_row(S, tmp_path):
    eid, syms = next(_first_ambiguous_ensembl(S))
    # gene_symbol column holds NDUFA9; a dedicated ensembl column holds the ambiguous id
    p = _write(tmp_path, [("", eid, 1), ("NDUFA9", "", -1)], ["gene_symbol", "ensembl", "direction"])
    res = signature.mapping_audit(p, universe=S.genes)
    audit, counts = res["audit"], res["counts"]
    assert counts["n_audit_rows"] == counts["n_input_rows"] == 2  # row conservation
    amb = audit[audit.input_id == eid].iloc[0]
    assert amb.ambiguous and amb.n_candidates > 1
    # the unambiguous NDUFA9 row is not flagged ambiguous
    clean = audit[audit.input_id == "NDUFA9"].iloc[0]
    assert not clean.ambiguous and clean.selected_symbol == "NDUFA9"


def test_audit_conserves_rows_with_duplicates(S, tmp_path):
    # duplicates and conflicts must remain visible as separate audit rows, never collapsed away
    p = _write(tmp_path, [("NDUFA9", -1), ("ndufa9", -1), ("UQCRC1", 1), ("UQCRC1", -1), ("NOTAGENE", 1)],
               ["gene_symbol", "direction"])
    res = signature.mapping_audit(p, universe=S.genes)
    audit, counts = res["audit"], res["counts"]
    assert counts["n_input_rows"] == counts["n_audit_rows"] == 5
    assert counts["n_unmapped"] >= 1  # NOTAGENE
    # status classes partition the input rows
    assert counts["n_mapped"] + counts["n_unmapped"] + counts["n_dropped"] == counts["n_input_rows"]


def test_audit_on_full19_fixture(S):
    from motrpac_probe.paths import EXAMPLES
    res = signature.mapping_audit(EXAMPLES / "pah_muscle_malenfant2015.csv", universe=S.genes)
    counts = res["counts"]
    assert counts["n_input_rows"] == 25
    assert counts["n_audit_rows"] == 25
    # the 19 counted proteins plus context rows that map; none should be silently dropped as duplicates
    assert counts["n_mapped"] >= 19
