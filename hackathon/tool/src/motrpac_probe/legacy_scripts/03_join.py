#!/usr/bin/env python3
"""03_join.py — one long gene-level table of exercise effects, human + rat, RNA + PROT.

Run from hackathon/:  python3 scripts/03_join.py
Input : data/raw/*_v1.csv written by 01_load_human.R / 02_load_rat.R
Output: data/join_table_v2.csv  (write-once; refuses to overwrite — bump OUT and note it in LOG.md)
        data/rat_rna_excluded_v1.csv (the rat RNA rows removed by the v2 filter, for transparency)

Versions:
  v1  as published, every feature kept.
  v2  rat RNA feature x contrast rows are dropped BEFORE collapsing when either group's mean normalized
      count is < 1 (min of comparison/reference_average_intensity) or |logFC| > 5. DESeq2 gives these
      on/off features (all-zero in one group; recurrent across tissues, e.g. Myh9, Mast4, Nat8b,
      Aldh18a1 — probably genotype, not training) logFC up to +/-60 and |z| ~ 10, so they won the
      max-|stat| collapse over the real gene. Human and rat PROT rows are unchanged.

One row = dataset x tissue x layer x contrast x human gene symbol.
  dataset   human_acute | rat_train
  tissue    human: VL (vastus lateralis biopsy); rat: SKM-GN, SKM-VL, HEART
  layer     RNA | PROT
  contrast  human: {EE,RE}_vs_CON_{15to45min,3.5to4h,24h}  (delta-delta vs control, see FACTSHEET)
            rat:   {F,M}_{1w,2w,4w,8w}                      (trained vs sex-matched sedentary)
  gene_symbol_human  human: HUMAN_FEATURE_TO_GENE; rat: FEATURE_TO_GENE -> RAT_TO_HUMAN_GENE (RGD)
  feature_id  the feature kept after collapsing (max |stat| among features mapping to the gene)
  logFC, stat, p, fdr  from that feature, as published
  stat_type   human: moderated t (limma/dream 't'); rat RNA: 'zscore' (DESeq2 Wald); rat PROT: 't' (limma)
  fdr         published adjusted p. human: BH within tissue x assay x contrast.
              rat: BY across ALL tissues within the assay (very conservative; see FACTSHEET)
  fdr_bh      BH recomputed by us within dataset x tissue x layer x contrast over all measured
              features (before collapsing) — comparable between human and rat
  training_fdr  rat only: selection_fdr, IHW-adjusted overall training test p (DESeq2 LRT for RNA, limma F-test for PROT; per sex, Fisher-combined across sexes; one per feature)
  n_collapsed number of features that mapped to this gene in this slice
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

RAW = "data/raw"
OUT = "data/join_table_v2.csv"
EXCL = "data/rat_rna_excluded_v1.csv"

if os.path.exists(OUT):
    print(f"{OUT} exists and is immutable; nothing to do. To regenerate, bump OUT to the next _vN and log it.")
    sys.exit(0)


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    out[ok] = false_discovery_control(p[ok], method="bh")
    return out


# ---------------------------------------------------------------- human
H_CONTRASTS = {
    "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)": "EE_vs_CON_15to45min",
    "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)": "EE_vs_CON_3.5to4h",
    "Endur.post_24_hr - Control.post_24_hr (delta-delta)": "EE_vs_CON_24h",
    "Resist.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)": "RE_vs_CON_15to45min",
    "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)": "RE_vs_CON_3.5to4h",
    "Resist.post_24_hr - Control.post_24_hr (delta-delta)": "RE_vs_CON_24h",
}
hf2g = pd.read_csv(f"{RAW}/human_feature_to_gene_muscle_v1.csv", low_memory=False)
hf2g = hf2g.dropna(subset=["gene_symbol"]).drop_duplicates(["assay", "feature_id"])

human = []
for layer, fn, assay in [("RNA", "rna", "transcript-rna-seq"), ("PROT", "prot", "prot-pr")]:
    d = pd.read_csv(f"{RAW}/human_muscle_{fn}_DA_v1.csv",
                    usecols=["contrast_short", "feature_id", "logFC", "t", "p_value", "adj_p_value"])
    d = d[d.contrast_short.isin(H_CONTRASTS)].copy()
    d["contrast"] = d.contrast_short.map(H_CONTRASTS)
    d["fdr_bh"] = d.groupby("contrast").p_value.transform(bh)
    sym = hf2g[hf2g.assay == assay].set_index("feature_id").gene_symbol
    d["gene_symbol_human"] = d.feature_id.map(sym)
    d = d.rename(columns={"t": "stat", "p_value": "p", "adj_p_value": "fdr"})
    d = d.assign(dataset="human_acute", tissue="VL", layer=layer, stat_type="t_moderated",
                 training_fdr=np.nan)
    human.append(d.drop(columns="contrast_short"))
human = pd.concat(human, ignore_index=True)

# ---------------------------------------------------------------- rat
rf2g = pd.read_csv(f"{RAW}/rat_feature_to_gene_v1.csv", low_memory=False)
rf2g = rf2g.dropna(subset=["gene_symbol"])[["feature_ID", "gene_symbol"]].drop_duplicates()
orth = pd.read_csv(f"{RAW}/rat_to_human_gene_v1.csv", low_memory=False)
orth = orth.dropna(subset=["HUMAN_ORTHOLOG_SYMBOL"])[["RAT_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL"]].drop_duplicates()
fmap = rf2g.merge(orth, left_on="gene_symbol", right_on="RAT_SYMBOL")[["feature_ID", "HUMAN_ORTHOLOG_SYMBOL"]]
fmap = fmap.drop_duplicates()

rat, excluded = [], []
for tissue, tcode, layer, assay, statcol in [
    ("SKM-GN", "skmgn", "RNA", "trnscrpt", "zscore"),
    ("SKM-GN", "skmgn", "PROT", "prot", "tscore"),
    ("SKM-VL", "skmvl", "RNA", "trnscrpt", "zscore"),
    ("HEART", "heart", "RNA", "trnscrpt", "zscore"),
    ("HEART", "heart", "PROT", "prot", "tscore"),
]:
    d = pd.read_csv(f"{RAW}/rat_{tcode}_{assay}_DA_v1.csv",
                    usecols=["feature_ID", "sex", "comparison_group", "logFC", statcol, "p_value",
                             "adj_p_value", "selection_fdr", "comparison_average_intensity",
                             "reference_average_intensity"])
    d["contrast"] = d.sex.str[0].str.upper() + "_" + d.comparison_group
    d["fdr_bh"] = d.groupby("contrast").p_value.transform(bh)  # BH over all tested features, as published
    if layer == "RNA":
        low = d[["comparison_average_intensity", "reference_average_intensity"]].min(axis=1)
        bad = (low < 1) | (d.logFC.abs() > 5)
        excluded.append(d[bad].assign(tissue=tissue))
        d = d[~bad]
    d = d.drop(columns=["comparison_average_intensity", "reference_average_intensity"])
    d = d.merge(fmap, on="feature_ID", how="left")  # a feature may map to >1 human gene -> one row each
    d = d.rename(columns={"feature_ID": "feature_id", statcol: "stat", "p_value": "p",
                          "adj_p_value": "fdr", "selection_fdr": "training_fdr",
                          "HUMAN_ORTHOLOG_SYMBOL": "gene_symbol_human"})
    d = d.assign(dataset="rat_train", tissue=tissue, layer=layer,
                 stat_type="zscore" if statcol == "zscore" else "t")
    rat.append(d.drop(columns=["sex", "comparison_group"]))
rat = pd.concat(rat, ignore_index=True)
excluded = pd.concat(excluded, ignore_index=True)
print(f"rat RNA rows excluded (low count in one group or |logFC|>5): {len(excluded)}")
print(excluded.groupby(["tissue", "contrast"]).size().unstack().to_string())
if not os.path.exists(EXCL):
    excluded.to_csv(EXCL, index=False)
    os.chmod(EXCL, 0o444)

# ---------------------------------------------------------------- collapse to gene
allx = pd.concat([human, rat], ignore_index=True)
unmapped = allx[allx.gene_symbol_human.isna()]
print("rows without a human gene symbol (dropped):")
print(unmapped.groupby(["dataset", "tissue", "layer"]).feature_id.nunique().to_string())
allx = allx.dropna(subset=["gene_symbol_human", "stat"])

key = ["dataset", "tissue", "layer", "contrast", "gene_symbol_human"]
allx["abs_stat"] = allx.stat.abs()
allx["n_collapsed"] = allx.groupby(key).feature_id.transform("nunique")
allx = allx.sort_values("abs_stat", ascending=False).drop_duplicates(key)

cols = ["dataset", "tissue", "layer", "contrast", "gene_symbol_human", "feature_id", "logFC", "stat",
        "p", "fdr", "n_collapsed", "stat_type", "fdr_bh", "training_fdr"]
allx = allx[cols].sort_values(["dataset", "tissue", "layer", "contrast", "gene_symbol_human"])
allx.to_csv(OUT, index=False)
os.chmod(OUT, 0o444)

summ = allx.groupby(["dataset", "tissue", "layer", "contrast"]).agg(
    n_genes=("gene_symbol_human", "size"), n_fdr05=("fdr", lambda s: int((s < 0.05).sum())),
    n_fdrbh05=("fdr_bh", lambda s: int((s < 0.05).sum())), multi_feature_genes=("n_collapsed", lambda s: int((s > 1).sum())))
print(summ.to_string())
print(f"\nwrote {OUT}: {len(allx)} rows")
print("READY: data/join_table_v2.csv, {} rows, tissues={}, contrasts={}".format(
    len(allx), ",".join(sorted(allx.tissue.unique())), ",".join(allx.contrast.unique())))
