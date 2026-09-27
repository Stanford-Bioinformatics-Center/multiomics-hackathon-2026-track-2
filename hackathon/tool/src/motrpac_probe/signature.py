"""Read a signature CSV (symbols, UniProt, Ensembl or rat symbols; +1 / -1 directions) and map it to the
store's human gene symbols, reporting every row that could not be used.

Method details: docs/METHODS.md#signaturepy
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


# --------------------------------------------------------------------------------------------------
# Candidate exposure for the mandatory mapping-confirmation gate.
#
# `map_genes()` / `load()` resolve each identifier to ONE human symbol (first-wins on the ID map),
# which is what the scientific analysis and its golden outputs use. The application's upload workflow,
# however, must show the researcher EVERY candidate before a selection is locked in, because 176
# Ensembl ids in the released ID map resolve to more than one human symbol. The functions below expose
# all candidates and never silence duplication with drop_duplicates(); they do not change scoring.

def _all_maps():
    """Full ID map WITHOUT drop_duplicates: id_type -> {id: [human_symbol, ...]} (all candidates)."""
    m = store.load_idmap()
    out = {}
    for id_type, sub in m.groupby("id_type"):
        key = sub.id.str.lower() if id_type == "rat_symbol" else sub.id
        d = {}
        for k, g in zip(key, sub.gene_symbol_human):
            d.setdefault(k, [])
            if g not in d[k]:
                d[k].append(g)
        out[id_type] = d
    return out


def candidates_for(row, universe, all_maps=None):
    """All candidate human symbols for one input row, in resolution-priority order, each tagged with
    the id namespace and route it came from. Returns a list of dicts:
    {candidate_symbol, id_type, source_column, source_value, measured (bool)}.
    `row` is a dict-like with lower-cased keys (gene_symbol/uniprot/ensembl/rat_symbol)."""
    all_maps = all_maps or _all_maps()
    uset = set(universe)
    uni, ens, rat = all_maps.get("uniprot", {}), all_maps.get("ensembl", {}), all_maps.get("rat_symbol", {})

    def val(c):
        v = row.get(c)
        return v.strip() if isinstance(v, str) and v.strip() else None

    seen, cands = set(), []

    def add(sym, id_type, col, value):
        keyc = (sym, id_type, col)
        if keyc in seen:
            return
        seen.add(keyc)
        cands.append({"candidate_symbol": sym, "id_type": id_type, "source_column": col,
                      "source_value": value, "measured": sym in uset})

    g = val("gene_symbol")
    if g:
        if g in uset:
            add(g, "human symbol", "gene_symbol", g)
        elif g.upper() in uset:
            add(g.upper(), "human symbol (case-normalised)", "gene_symbol", g)
        if ENSEMBL_RE.match(g):
            for s in ens.get(g.split(".")[0], []):
                add(s, "Ensembl id in gene_symbol", "gene_symbol", g)
        if UNIPROT_RE.match(g):
            for k in (g, g.split("-")[0]):
                for s in uni.get(k, []):
                    add(s, "UniProt in gene_symbol", "gene_symbol", g)
        for s in rat.get(g.lower(), []):
            add(s, "rat symbol -> human ortholog (RGD)", "gene_symbol", g)
    u = val("uniprot")
    if u:
        for k in (u, u.split("-")[0]):
            for s in uni.get(k, []):
                add(s, "UniProt (MoTrPAC feature map)", "uniprot", u)
    e = val("ensembl")
    if e:
        for s in ens.get(e.split(".")[0], []):
            add(s, "Ensembl", "ensembl", e)
    rs = val("rat_symbol")
    if rs:
        for s in rat.get(rs.lower(), []):
            add(s, "rat symbol -> human ortholog (RGD)", "rat_symbol", rs)
    return cands


def mapping_audit(path, universe=None):
    """A per-input-row mapping audit for the confirmation gate. Reports EVERY input row with:
    the selected (first-wins) mapping used by the analysis, the number of candidates, whether the
    mapping is ambiguous (>1 distinct candidate symbol), and the full candidate list. Row counts are
    conserved (one audit row per input row); duplication is never silenced. Returns a dict with the
    audit DataFrame and reconciliation counts."""
    from .core import get_store
    t = read_table(path)
    if "direction" not in t.columns:
        t["direction"] = ""
    universe = get_store().genes if universe is None else universe
    resolved = map_genes(t, universe)   # the selection the analysis actually uses (unchanged)
    all_maps = _all_maps()
    records = []
    for i, row in enumerate(t.to_dict("records")):
        cands = candidates_for(row, universe, all_maps)
        measured = [c for c in cands if c["measured"]]
        distinct = sorted({c["candidate_symbol"] for c in measured}) or sorted({c["candidate_symbol"] for c in cands})
        sel = resolved.iloc[i]
        records.append({
            "input_row": i,
            "input_id": next((str(row[c]) for c in ("gene_symbol", "uniprot", "ensembl", "rat_symbol")
                              if c in row and isinstance(row[c], str) and row[c].strip()), ""),
            "selected_symbol": sel.gene if pd.notna(sel.gene) else "",
            "selected_method": sel.map_method if pd.notna(sel.map_method) else "",
            "status": sel.status,
            "n_candidates": len(distinct),
            "ambiguous": len(distinct) > 1,
            "candidates": ";".join(distinct),
        })
    audit = pd.DataFrame(records)
    counts = {
        "n_input_rows": len(t),
        "n_audit_rows": len(audit),
        "n_mapped": int((resolved.status == "mapped").sum()),
        "n_ambiguous": int(audit.ambiguous.sum()),
        "n_unmapped": int(resolved.status.str.startswith("unmapped").sum()),
        "n_dropped": int(resolved.status.str.startswith("dropped").sum()),
    }
    # row-count reconciliation: one audit row per input row, and the status classes partition the rows
    assert counts["n_audit_rows"] == counts["n_input_rows"], "mapping audit must conserve input rows"
    return {"audit": audit, "counts": counts}


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
