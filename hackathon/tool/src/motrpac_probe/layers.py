"""Which layer answered? RNA vs protein per gene at matched timepoints, next to the detection-power flags
(baseline mRNA, proteomics missingness, GTEx) that decide whether a layer could have answered.

Method details: docs/METHODS.md#layerspy
"""
import numpy as np
import pandas as pd

from . import legacy
from .core import EXERCISE_KINDS, TISSUE_WORDS, col_words

DETECT_TISSUES = [("human_acute", "VL"), ("rat_train", "SKM-GN"), ("rat_train", "SKM-VL"), ("rat_train", "HEART")]


def _tissue_rows(S, ds, tissue, layer):
    c = S.cols
    cids = c[(c.dataset == ds) & (c.tissue == tissue) & (c.layer == layer) & c.kind.isin(EXERCISE_KINDS)].column_id
    if not len(cids):
        return None, 0
    pos = np.concatenate([S.col(x) for x in cids])
    return pd.DataFrame({"gene": S.genes[S.gcode[pos]], "base": S.base[pos], "miss": S.miss[pos],
                         "ncoll": S.ncoll[pos]}), len(cids)


def _tertile(v):
    out = pd.Series(np.nan, index=v.index)
    ok = v.notna()
    if ok.sum() >= 3:
        out[ok] = pd.qcut(v[ok].rank(method="first"), 3, labels=[1, 2, 3]).astype(float)
    return out.map({1.0: "low", 2.0: "mid", 3.0: "high"})


def detection_flags(S, genes):
    """One row per gene x core tissue with the detection-power flags."""
    from .store import GTEX_TISSUE
    out = []
    for ds, tissue in DETECT_TISSUES:
        r, nr = _tissue_rows(S, ds, tissue, "RNA")
        p, npr = _tissue_rows(S, ds, tissue, "PROT")
        rec = pd.DataFrame({"gene": genes})
        rec.insert(0, "tissue", f"{'Human' if ds == 'human_acute' else 'Rat'} {TISSUE_WORDS[tissue]}")
        if r is not None:
            rb = r.groupby("gene").base.median()
            rc = r.groupby("gene").size()
            rec["rna_columns"] = rec.gene.map(rc).fillna(0).astype(int).astype(str) + f"/{nr}"
            rec["rna_baseline"] = rec.gene.map(rb)
            rec["rna_baseline_tertile"] = rec.gene.map(_tertile(rb))
            rec["n_collapsed_rna"] = rec.gene.map(r.groupby("gene").ncoll.max())
        else:
            rec["rna_columns"] = "0/0"
        if p is not None:
            pc = p.groupby("gene").size()
            rec["prot_columns"] = rec.gene.map(pc).fillna(0).astype(int).astype(str) + f"/{npr}"
            rec["prot_n_missing"] = rec.gene.map(p.groupby("gene").miss.median())
            rec["n_collapsed_prot"] = rec.gene.map(p.groupby("gene").ncoll.max())
        else:
            rec["prot_columns"] = "0/0"
        gt = GTEX_TISSUE[tissue]
        tpm = S.ann[f"gtex_tpm|{gt}"]
        rec["gtex_tpm"] = rec.gene.map(tpm)
        if r is not None:
            universe = np.log2(tpm.reindex(r.gene.unique()) + 1)
            rec["gtex_tertile"] = rec.gene.map(_tertile(universe))
        rec["baseline_unit"] = "log2 CPM (AveExpr)" if ds == "human_acute" else "log2(sedentary count + 1)"
        out.append(rec)
    return pd.concat(out, ignore_index=True)


def coverage(S, sig, cids):
    """Gene x tissue-layer block: number of columns in which the gene is measured."""
    c = S.cols.loc[list(cids)]
    blocks = c.groupby(["dataset", "tissue", "layer"], sort=False).column_id.apply(list)
    rows = []
    codes = dict(zip(sig.shown, S.codes(sig.shown)))
    for gene in sig.shown:
        rec = {"gene": gene, "group": sig.groups.get(gene, "")}
        for (ds, t, l), ids in blocks.items():
            n = sum(int((S.gcode[S.col(x)] == codes[gene]).any()) for x in ids) if codes[gene] >= 0 else 0
            rec[(ds, t, l)] = f"{n}/{len(ids)}"
        rows.append(rec)
    return pd.DataFrame(rows), list(blocks.index)


# --------------------------------------------------------------------------------------------- E
def layer_pairs(S, kinds=EXERCISE_KINDS, include_cross_muscle=True):
    """(RNA column, PROT column) pairs in the same dataset x tissue x contrast; plus rat SKM-VL RNA vs SKM-GN PROT
    (scripts/06_discordance.py PAIRS) flagged as cross-muscle."""
    c = S.cols[S.cols.kind.isin(kinds)]
    pr = c[c.layer == "PROT"]
    out = []
    for p in pr.itertuples():
        r = c[(c.dataset == p.dataset) & (c.tissue == p.tissue) & (c.layer == "RNA") & (c.contrast == p.contrast)]
        if len(r):
            out.append(dict(rna=r.column_id.iloc[0], prot=p.column_id, cross=False))
        if include_cross_muscle and p.tissue == "SKM-GN":
            r = c[(c.tissue == "SKM-VL") & (c.layer == "RNA") & (c.contrast == p.contrast)]
            if len(r):
                out.append(dict(rna=r.column_id.iloc[0], prot=p.column_id, cross=True))
    return pd.DataFrame(out)


def _merged(S, rna, prot):
    a, b = S.frame(rna), S.frame(prot)
    return a[["gene_symbol_human", "logFC", "stat", "fdr_bh"]].merge(
        b[["gene_symbol_human", "logFC", "stat", "fdr_bh"]], on="gene_symbol_human", suffixes=("_rna", "_prot")
    ).dropna(subset=["stat_rna", "stat_prot"])


def classify(m, cutoff):
    """scripts/08_discordance_model.py label rule."""
    rs, ps = m.fdr_bh_rna < cutoff, m.fdr_bh_prot < cutoff
    same = np.sign(m.stat_rna) == np.sign(m.stat_prot)
    return np.select([rs & ps & same, rs & ps & ~same, rs & ~ps, ~rs & ps],
                     ["concordant", "opposite", "RNA_only", "PROT_only"], "ns")


def discordance(S, sig, pairs, cutoff=0.05):
    """Per pair: genome-wide and signature rho (scripts/06 rho_ci), opposed counts per layer, class counts.
    Also the per-gene table for signature genes measured in both layers."""
    dirs = dict(zip(sig.genes, sig.dirs))
    summ, genes = [], []
    for pr in pairs.itertuples():
        m = _merged(S, pr.rna, pr.prot)
        n_all, r_all, p_all, lo_all, hi_all = legacy.rho_ci(m.stat_rna, m.stat_prot)
        ms = m[m.gene_symbol_human.isin(dirs)].copy()
        ms["direction"] = ms.gene_symbol_human.map(dirs)
        ms["rna_opposed"] = np.sign(ms.logFC_rna) * ms.direction < 0
        ms["prot_opposed"] = np.sign(ms.logFC_prot) * ms.direction < 0
        ms["layer_class"] = classify(ms, cutoff)
        n_s, r_s, p_s, lo_s, hi_s = legacy.rho_ci(ms.stat_rna, ms.stat_prot)
        cr, cp = S.cols.loc[pr.rna], S.cols.loc[pr.prot]
        label = col_words(cp).replace(" protein", "")
        if pr.cross:
            label = label.replace("gastrocnemius", "vastus lateralis RNA vs gastrocnemius protein")
        rec = dict(rna=pr.rna, prot=pr.prot, cross=pr.cross, label=label, dataset=cp.dataset, tissue=cp.tissue,
                   contrast=cp.contrast, early=bool(cp.early_human_muscle),
                   n_genes_both=n_all, rho_all=r_all, rho_all_lo=lo_all, rho_all_hi=hi_all,
                   n_sig_both=len(ms), rna_opposed=int(ms.rna_opposed.sum()), prot_opposed=int(ms.prot_opposed.sum()),
                   both_opposed=int((ms.rna_opposed & ms.prot_opposed).sum()),
                   rho_sig=r_s, rho_sig_lo=lo_s, rho_sig_hi=hi_s,
                   frac_same_sign_sig=float((np.sign(ms.logFC_rna) == np.sign(ms.logFC_prot)).mean()) if len(ms) else np.nan)
        for k in ["concordant", "RNA_only", "PROT_only", "opposite", "ns"]:
            rec[f"class_{k}"] = int((ms.layer_class == k).sum())
        allc = classify(m, cutoff)
        for k in ["concordant", "RNA_only", "PROT_only", "opposite"]:
            rec[f"genome_{k}"] = int((allc == k).sum())
        summ.append(rec)
        genes.append(ms.assign(pair=label, rna=pr.rna, prot=pr.prot, cross=pr.cross))
    G = pd.concat(genes, ignore_index=True) if genes else pd.DataFrame()
    return pd.DataFrame(summ), G


def sentence(r):
    """Templated E sentence (only numbers from the table)."""
    if not r.n_sig_both:
        return f"At {r.label}, no signature gene is measured in both layers."
    rho = "not estimable" if not np.isfinite(r.rho_all) else f"ρ = {r.rho_all:.2f}"
    return (f"At {r.label}, RNA opposed {r.rna_opposed}/{r.n_sig_both} and protein {r.prot_opposed}/{r.n_sig_both}; "
            f"genome-wide RNA–protein agreement here is {rho}.")
