"""Signature CSV reading and gene mapping (signature.load / signature.map_genes)."""
import numpy as np
import pandas as pd
import pytest

from motrpac_probe import signature


def _load(tmp_path, rows, columns=None, name="sig.csv", **kw):
    df = pd.DataFrame(rows, columns=columns)
    p = tmp_path / name
    df.to_csv(p, index=False)
    return signature.load(p, **kw)


def _status(sig):
    t = sig.table
    return {r.key: (r.gene, r.map_method, r.status) for r in t.assign(key=t.iloc[:, 0].astype(str)).itertuples()}


def test_human_symbols_and_case(S, tmp_path):
    sig = _load(tmp_path, [("NDUFA9", -1), ("uqcrc2", -1), ("Atp5f1b", -1), ("MYH7", 1)],
                ["gene_symbol", "direction"])
    assert sig.genes == ["NDUFA9", "UQCRC2", "ATP5F1B", "MYH7"]
    assert list(sig.dirs) == [-1, -1, -1, 1]
    st = _status(sig)
    assert st["NDUFA9"][1] == "human symbol"
    assert st["uqcrc2"][1] == "human symbol (case-normalised)"
    assert st["Atp5f1b"][1] == "human symbol (case-normalised)"


def test_rat_symbols(S, tmp_path):
    # Ndufa9 upper-cases to the human symbol; Zfp316 needs the RGD ortholog table (human ZNF316)
    sig = _load(tmp_path, [("Ndufa9", -1), ("Zfp316", 1), ("zfp316x_not_a_gene", 1)], ["gene_symbol", "direction"])
    st = _status(sig)
    assert st["Ndufa9"][0] == "NDUFA9"
    assert st["Zfp316"][:2] == ("ZNF316", "rat symbol -> human ortholog (RGD)")
    assert st["zfp316x_not_a_gene"][2].startswith("unmapped")
    # a dedicated rat_symbol column
    sig = _load(tmp_path, [("", "Zfp316", 1), ("", "Ndufa9", -1)], ["gene_symbol", "rat_symbol", "direction"],
                name="rat.csv")
    assert sig.genes == ["ZNF316", "NDUFA9"]
    assert set(sig.table.map_method) == {"rat symbol -> human ortholog (RGD)"}


def test_uniprot(S, tmp_path):
    sig = _load(tmp_path, [("P04217", "", -1), ("Q16795-1", "", -1), ("", "P04217", 1), ("XYZ", "Q16795", -1)],
                ["gene_symbol", "uniprot", "direction"])
    t = sig.table
    assert t.gene.tolist() == ["A1BG", "NDUFA9", "A1BG", "NDUFA9"]
    assert t.map_method.tolist()[:3] == ["UniProt in gene_symbol", "UniProt in gene_symbol",
                                          "UniProt (MoTrPAC feature map)"]
    assert t.map_method.iloc[3] == "UniProt (MoTrPAC feature map)"
    # A1BG listed -1 then +1 -> conflicting, dropped; NDUFA9 listed twice at -1 -> kept once
    assert sig.genes == ["NDUFA9"]
    assert t.status.tolist() == ["dropped: gene listed with conflicting directions", "mapped",
                                 "dropped: gene listed with conflicting directions",
                                 "dropped: duplicate of an earlier row"]


def test_uniprot_column_only(S, tmp_path):
    sig = _load(tmp_path, [("P04217", -1)], ["uniprot", "direction"])
    assert sig.genes == ["A1BG"]
    assert sig.table.map_method.iloc[0] == "UniProt (MoTrPAC feature map)"


def test_ensembl(S, tmp_path):
    sig = _load(tmp_path, [("ENSG00000121410", "", 1), ("ENSG00000139180.12", "", -1), ("", "ENSG00000183978.5", 1),
                           ("ENSG00000000000", "", 1)], ["gene_symbol", "ensembl", "direction"])
    t = sig.table
    assert t.gene.iloc[0] == "A1BG" and t.map_method.iloc[0] == "Ensembl id in gene_symbol"
    assert t.gene.iloc[1] == "NDUFA9" and t.map_method.iloc[1] == "Ensembl id in gene_symbol"
    assert t.map_method.iloc[2] == "Ensembl" and t.status.iloc[2] == "mapped"
    assert t.status.iloc[3].startswith("unmapped")
    assert "A1BG" in sig.genes and "NDUFA9" in sig.genes


def test_unknown_genes_reported(S, tmp_path):
    sig = _load(tmp_path, [("NOTAGENE123", 1), ("NDUFA9", -1), ("ZZZ9999", -1)], ["gene_symbol", "direction"])
    t = sig.table
    assert sig.genes == ["NDUFA9"]
    assert len(t) == 3  # every input row is reported
    un = t[t.status != "mapped"]
    assert un.gene_symbol.tolist() == ["NOTAGENE123", "ZZZ9999"]
    assert un.status.str.startswith("unmapped").all()
    assert un.gene.isna().all()


def test_duplicates(S, tmp_path):
    sig = _load(tmp_path, [("NDUFA9", -1), ("ndufa9", -1), ("Q16795", -1), ("UQCRC1", 1), ("UQCRC1", -1),
                           ("MYH7", 1)], ["gene_symbol", "direction"])
    t = sig.table
    assert sig.genes == ["NDUFA9", "MYH7"]
    assert t.status.tolist() == ["mapped", "dropped: duplicate of an earlier row", "dropped: duplicate of an earlier row",
                                 "dropped: gene listed with conflicting directions",
                                 "dropped: gene listed with conflicting directions", "mapped"]


@pytest.mark.parametrize("word, d", [("up", 1), ("UP", 1), (" Up ", 1), ("+", 1), ("1", 1), ("2.5", 1),
                                     ("down", -1), ("Down", -1), ("dn", -1), ("-", -1), ("-1", -1), ("-0.3", -1)])
def test_direction_words(S, tmp_path, word, d):
    sig = _load(tmp_path, [("NDUFA9", word)], ["gene_symbol", "direction"])
    assert list(sig.dirs) == [d]


@pytest.mark.parametrize("word", ["sideways", "0", ""])
def test_bad_direction_dropped(S, tmp_path, word):
    sig = _load(tmp_path, [("NDUFA9", word), ("MYH7", "up")], ["gene_symbol", "direction"])
    assert sig.genes == ["MYH7"]
    assert sig.table.status.iloc[0] == "dropped: direction not +1/-1"


def test_missing_direction_column(S, tmp_path):
    p = tmp_path / "nodir.csv"
    pd.DataFrame({"gene_symbol": ["NDUFA9"], "weight": [1.0]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="direction"):
        signature.load(p)


def test_missing_id_column(S, tmp_path):
    p = tmp_path / "noid.csv"
    pd.DataFrame({"refmet_name": ["Citric acid"], "direction": [-1]}).to_csv(p, index=False)
    with pytest.raises(ValueError, match="gene_symbol"):
        signature.load(p)


def test_header_case_and_context_groups(S, tmp_path):
    sig = _load(tmp_path, [("NDUFA9", -1, "OXPHOS"), ("CS", -1, "phenotype"), ("LDHA", 1, "phenotype")],
                ["Gene_Symbol", " Direction ", "GROUP"], name="MySig.csv")
    assert sig.name == "MySig"
    assert sig.genes == ["NDUFA9"]
    assert sig.shown == ["NDUFA9", "CS", "LDHA"]
    assert list(sig.shown_dirs) == [-1, -1, 1]
    assert sig.groups == {"NDUFA9": "OXPHOS", "CS": "phenotype", "LDHA": "phenotype"}
    all_counted = signature.load(tmp_path / "MySig.csv", context_groups=())
    assert all_counted.genes == ["NDUFA9", "CS", "LDHA"]


def test_empty_group_becomes_signature(S, tmp_path):
    sig = _load(tmp_path, [("NDUFA9", -1)], ["gene_symbol", "direction"])
    assert sig.groups == {"NDUFA9": "signature"}
    assert isinstance(sig.dirs, np.ndarray)
