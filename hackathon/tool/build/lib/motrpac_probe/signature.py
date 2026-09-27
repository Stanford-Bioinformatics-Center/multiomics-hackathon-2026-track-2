"""Read a signature CSV and map it to the human gene symbols used in the store.

Input columns (case-insensitive): gene_symbol and/or uniprot / ensembl / rat_symbol, direction (+1 up in disease,
-1 down; 'up'/'down' accepted), optional group, weight, source. Mapping order per row, first hit wins:
  1. gene_symbol exactly as given, then upper-cased, against the store's human symbols;
  2. a gene_symbol value that looks like an Ensembl (ENSG/ENSRNOG) or UniProt accession is mapped as such;
  3. rat symbol -> human ortholog (RGD table used by scripts/03_join.py), case-insensitive;
  4. the uniprot column (MoTrPAC human proteomics feature map; isoform suffix dropped if needed);
  5. the ensembl column; 6. the rat_symbol column.
Duplicated genes keep the first row if directions agree and are dropped if they conflict. Every input row is
reported with its mapping status.
"""
import re

import numpy as np
import pandas as pd

from . import store

UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$")
ENSEMBL_RE = re.compile(r"^ENS(RNO)?G\d{11}(\.\d+)?$")
DIR_WORDS = {"up": 1, "+": 1, "down": -1, "-": -1, "dn": -1}


def _dir(v):
    if isinstance(v, str):
        s = v.strip().lower()
        if s in DIR_WORDS:
            return DIR_WORDS[s]
        try:
            v = float(s)
        except ValueError:
            return 0
    try:
        v = float(v)
    except (TypeError, ValueError):
        return 0
    return int(np.sign(v)) if np.isfinite(v) else 0


def _maps():
    m = store.load_idmap()
    uni = m[m.id_type == "uniprot"].drop_duplicates("id").set_index("id").gene_symbol_human
    ens = m[m.id_type == "ensembl"].drop_duplicates("id").set_index("id").gene_symbol_human
    rat = m[m.id_type == "rat_symbol"]
    rat = rat.assign(k=rat.id.str.lower()).drop_duplicates("k").set_index("k").gene_symbol_human
    return uni, ens, rat


def read_table(path):
    t = pd.read_csv(path, dtype=str, keep_default_na=False, na_values=[""])
    t.columns = [c.strip().lower() for c in t.columns]
    return t


def map_genes(t, universe):
    """Returns the input table with gene (human symbol or NaN), map_method, status columns."""
    uni, ens, rat = _maps()
    uset = set(universe)

    def one(r):
        def val(c):
            v = r.get(c)
            return v.strip() if isinstance(v, str) and v.strip() else None
        g = val("gene_symbol")
        if g:
            if g in uset:
                return g, "human symbol"
            if g.upper() in uset:
                return g.upper(), "human symbol (case-normalised)"
            if ENSEMBL_RE.match(g) and g.split(".")[0] in ens.index:
                return ens[g.split(".")[0]], "Ensembl id in gene_symbol"
            if UNIPROT_RE.match(g):
                for k in (g, g.split("-")[0]):
                    if k in uni.index:
                        return uni[k], "UniProt in gene_symbol"
            if g.lower() in rat.index:
                return rat[g.lower()], "rat symbol -> human ortholog (RGD)"
        u = val("uniprot")
        if u:
            for k in (u, u.split("-")[0]):
                if k in uni.index:
                    return uni[k], "UniProt (MoTrPAC feature map)"
        e = val("ensembl")
        if e and e.split(".")[0] in ens.index:
            return ens[e.split(".")[0]], "Ensembl"
        rs = val("rat_symbol")
        if rs and rs.lower() in rat.index:
            return rat[rs.lower()], "rat symbol -> human ortholog (RGD)"
        return None, None

    got = [one(r) for r in t.to_dict("records")]
    t = t.copy()
    t["gene"] = [g for g, _ in got]
    t["map_method"] = [m for _, m in got]
    t["direction_parsed"] = [_dir(v) for v in t.get("direction", pd.Series([None] * len(t)))]
    t["status"] = "mapped"
    t.loc[t.gene.isna(), "status"] = "unmapped: not found in MoTrPAC gene/protein annotation"
    t.loc[t.gene.notna() & ~t.gene.isin(uset), "status"] = "unmapped: mapped symbol not measured in any MoTrPAC column"
    t.loc[t.direction_parsed == 0, "status"] = "dropped: direction not +1/-1"
    ok = t.status == "mapped"
    dup = t[ok].groupby("gene").direction_parsed.agg(["nunique", "size"])
    conflict = set(dup.index[dup["nunique"] > 1])
    t.loc[ok & t.gene.isin(conflict), "status"] = "dropped: gene listed with conflicting directions"
    ok = t.status == "mapped"
    first = ~t[ok].duplicated("gene")
    t.loc[first.index[~first], "status"] = "dropped: duplicate of an earlier row"
    return t


def load(path, name=None, context_groups=("phenotype",), universe=None):
    from .core import Sig, get_store
    t = read_table(path)
    idcols = {"gene_symbol", "uniprot", "ensembl", "rat_symbol"}
    if not idcols & set(t.columns):
        raise ValueError(f"{path}: needs a gene_symbol, uniprot, ensembl or rat_symbol column")
    if "direction" not in t.columns:
        raise ValueError(f"{path}: needs a direction column (+1 up in disease, -1 down)")
    if "group" not in t.columns:
        t["group"] = ""
    t["group"] = t.group.fillna("")
    universe = get_store().genes if universe is None else universe
    t = map_genes(t, universe)
    ok = t[t.status == "mapped"]
    ctx = ok.group.isin(context_groups)
    counted = ok[~ctx]
    shown = pd.concat([counted, ok[ctx]])
    if "weight" in ok.columns:
        w = pd.to_numeric(ok.weight, errors="coerce")
        t["weight_num"] = pd.to_numeric(t.weight, errors="coerce")
    return Sig(name=name or str(path).rsplit("/", 1)[-1].rsplit(".", 1)[0], table=t,
               genes=counted.gene.tolist(), dirs=counted.direction_parsed.to_numpy(int),
               shown=shown.gene.tolist(), shown_dirs=shown.direction_parsed.to_numpy(int),
               groups=dict(zip(shown.gene, shown.group.replace("", "signature"))))
