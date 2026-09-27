#!/usr/bin/env python3
"""10_everywhere.py — the PAH-signature grid in EVERY available MoTrPAC comparison ("everywhere" grid).

Run from hackathon/ after scripts/10_everywhere.R:  python3 scripts/10_everywhere.py
Inputs (read-only):
  data/raw/rat_alltissue_{trnscrpt,prot,phospho}_DA_v1.csv, rat_alltissue_feature_to_gene_v1.csv,
  data/raw/rat_to_human_gene_v1.csv                       rat, 19 tissues (10_everywhere.R)
  data/raw/human_muscle_{rna,prot,phospho}_DA_v1.csv       human muscle (01_load_human.R)
  everywhere/raw/human_{adipose,blood}_*_DA_v1.csv.gz, human_feature_to_gene_all_v1.csv.gz
  data/join_table_v2.csv                                   cross-check only
Outputs (everywhere/):
  pah_everywhere_long.csv            one row per signature gene x comparison column (NA rows kept)
  pah_everywhere_columns.csv         per column: opposed / same / measured, sign test
  pah_everywhere_ranking.csv/.md     tissue x layer x time ranked by opposition and by same-direction fraction
  pah_everywhere_blocks.csv          per tissue x layer block: opposed / measured (pooled cells)
  pah_everywhere.png, captions.md, crosscheck_join_v2.md

Gene-level collapsing reuses 03_join.py exactly:
  fdr_bh = BH over ALL tested features within dataset x tissue x layer x contrast, before collapsing;
  rat RNA rows with min(comparison, reference mean normalized count) < 1 or |logFC| > 5 dropped (03_join v2 filter);
  rat feature -> FEATURE_TO_GENE symbol -> RAT_TO_HUMAN_GENE (RGD) ortholog; human feature -> HUMAN_FEATURE_TO_GENE;
  one row per gene = the feature with max |stat| (rat PHOSPHO: the site with max |stat|).
Signature and agreement exactly as 05_pah_grid.py: agreement = sign(exercise logFC) x PAH direction
(-1 = exercise moves the gene AGAINST PAH); colour = agreement x min(|stat|, 4); dot = fdr_bh < 0.05;
opposed / measured over the 19 Malenfant proteomic genes (the 6 enzyme-phenotype genes are shown, not counted).
"""
import os
import re

import numpy as np
import pandas as pd
from scipy.stats import binomtest, false_discovery_control

RAW, EV = "data/raw", "everywhere"
os.makedirs(EV, exist_ok=True)

# ------------------------------------------------------------------ signature (copied from 05_pah_grid.py)
SIGNATURE = [
    ("NDUFA9", 0.71, -1, "OXPHOS-down"), ("UQCRC2", 0.74, -1, "OXPHOS-down"),
    ("UQCRC1", 0.74, -1, "OXPHOS-down"), ("ATP5MG", 0.70, -1, "OXPHOS-down"),
    ("ATP5F1B", 0.74, -1, "OXPHOS-down"),
    ("IDH2", 0.74, -1, "TCA/transport-down"), ("OGDH", 0.77, -1, "TCA/transport-down"),
    ("SLC25A4", 0.65, -1, "TCA/transport-down"), ("ECH1", 0.57, -1, "TCA/transport-down"),
    ("GLO1", 1.39, 1, "glycolytic-up"), ("FBP2", 1.33, 1, "glycolytic-up"),
    ("PEBP1", 1.29, 1, "contractile/other-up"), ("ATP2A1", 1.33, 1, "contractile/other-up"),
    ("MYH1", 2.42, 1, "contractile/other-up"), ("MYH2", 1.26, 1, "contractile/other-up"),
    ("MYH7", 1.26, 1, "contractile/other-up"), ("MYLPF", 1.58, 1, "contractile/other-up"),
    ("APOA1", 1.29, 1, "contractile/other-up"), ("PDLIM3", 1.58, 1, "contractile/other-up"),
    ("CS", None, -1, "phenotype"), ("TFAM", None, -1, "phenotype"),
    ("PDHA1", None, -1, "phenotype"), ("PDHB", None, -1, "phenotype"),
    ("LDHA", None, 1, "phenotype"), ("LDHB", None, 1, "phenotype"),
]
sig = pd.DataFrame(SIGNATURE, columns=["gene", "pah_ratio", "pah_dir", "group"])
GENES = set(sig.gene)


def bh(p):
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    out[ok] = false_discovery_control(p[ok], method="bh")
    return out


def collapse(d, key):
    """max |stat| feature per gene (03_join.py); n_collapsed = distinct features mapping to the gene."""
    d = d.copy()
    d["abs_stat"] = d.stat.abs()
    d["n_collapsed"] = d.groupby(key).feature_id.transform("nunique")
    return d.sort_values("abs_stat", ascending=False).drop_duplicates(key).drop(columns="abs_stat")


# ------------------------------------------------------------------ rat, every tissue
rf2g = pd.read_csv(f"{RAW}/rat_alltissue_feature_to_gene_v1.csv", usecols=["feature_ID", "gene_symbol"],
                   low_memory=False).dropna().drop_duplicates()
orth = pd.read_csv(f"{RAW}/rat_to_human_gene_v1.csv", usecols=["RAT_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL"],
                   low_memory=False).dropna().drop_duplicates()
fmap = rf2g.merge(orth, left_on="gene_symbol", right_on="RAT_SYMBOL")[["feature_ID", "HUMAN_ORTHOLOG_SYMBOL"]]
fmap = fmap.drop_duplicates()
fmap_sig = fmap[fmap.HUMAN_ORTHOLOG_SYMBOL.isin(GENES)]

rat = []
for assay, layer, statcol in [("trnscrpt", "RNA", "zscore"), ("prot", "PROT", "tscore"), ("phospho", "PHOSPHO", "tscore")]:
    d = pd.read_csv(f"{RAW}/rat_alltissue_{assay}_DA_v1.csv", low_memory=False)
    d["contrast"] = d.sex.str[0].str.upper() + "_" + d.comparison_group
    d["fdr_bh"] = d.groupby(["tissue", "contrast"]).p_value.transform(bh)  # over all tested features
    if layer == "RNA":
        low = d[["comparison_average_intensity", "reference_average_intensity"]].min(axis=1)
        d = d[~((low < 1) | (d.logFC.abs() > 5))]
    d = d[d.feature_ID.isin(fmap_sig.feature_ID)].merge(fmap_sig, on="feature_ID")
    d = d.rename(columns={"feature_ID": "feature_id", statcol: "stat", "p_value": "p", "adj_p_value": "fdr",
                          "HUMAN_ORTHOLOG_SYMBOL": "gene"})
    d = d.assign(dataset="rat_train", layer=layer)[["dataset", "tissue", "layer", "contrast", "gene", "feature_id",
                                                   "logFC", "stat", "p", "fdr", "fdr_bh"]]
    rat.append(collapse(d.dropna(subset=["stat"]), ["dataset", "tissue", "layer", "contrast", "gene"]))
    print(f"rat {layer}: {d.tissue.nunique()} tissues")
rat = pd.concat(rat, ignore_index=True)

# ------------------------------------------------------------------ human, every tissue x ome x contrast
hf2g = pd.read_csv(f"{EV}/raw/human_feature_to_gene_all_v1.csv.gz", low_memory=False)
hf2g = hf2g.dropna(subset=["gene_symbol"]).drop_duplicates(["assay", "feature_id"])
HSRC = [("muscle", "RNA", "transcript-rna-seq", f"{RAW}/human_muscle_rna_DA_v1.csv"),
        ("muscle", "PROT", "prot-pr", f"{RAW}/human_muscle_prot_DA_v1.csv"),
        ("muscle", "PHOSPHO", "prot-ph", f"{RAW}/human_muscle_phospho_DA_v1.csv"),
        ("adipose", "RNA", "transcript-rna-seq", f"{EV}/raw/human_adipose_rna_DA_v1.csv.gz"),
        ("adipose", "PROT", "prot-pr", f"{EV}/raw/human_adipose_prot_DA_v1.csv.gz"),
        ("adipose", "PHOSPHO", "prot-ph", f"{EV}/raw/human_adipose_phospho_DA_v1.csv.gz"),
        ("blood", "RNA", "transcript-rna-seq", f"{EV}/raw/human_blood_rna_DA_v1.csv.gz"),
        ("blood", "PROT (Olink)", "prot-ol", f"{EV}/raw/human_blood_olink_DA_v1.csv.gz")]
human, hcon = [], []
for tissue, layer, assay, f in HSRC:
    d = pd.read_csv(f, usecols=["contrast_short", "contrast_category", "Timepoint", "feature_id", "logFC", "t",
                                "p_value", "adj_p_value"], low_memory=False)
    d["fdr_bh"] = d.groupby("contrast_short").p_value.transform(bh)
    hcon.append(d[["contrast_short", "contrast_category", "Timepoint"]].drop_duplicates().assign(tissue=tissue, layer=layer))
    sym = hf2g[hf2g.assay == assay].set_index("feature_id").gene_symbol
    d["gene"] = d.feature_id.map(sym)
    d = d[d.gene.isin(GENES)].rename(columns={"t": "stat", "p_value": "p", "adj_p_value": "fdr",
                                              "contrast_short": "contrast"})
    d = d.assign(dataset="human_acute", tissue=tissue, layer=layer)[
        ["dataset", "tissue", "layer", "contrast", "gene", "feature_id", "logFC", "stat", "p", "fdr", "fdr_bh"]]
    human.append(collapse(d.dropna(subset=["stat"]), ["dataset", "tissue", "layer", "contrast", "gene"]))
human = pd.concat(human, ignore_index=True)
hcon = pd.concat(hcon, ignore_index=True)

# ------------------------------------------------------------------ comparison columns (order = figure order)
TP = {"during_20_min": "20m-during", "during_40_min": "40m-during", "post_10_min": "10m", "post_15_30_45_min": "15-45m",
      "post_3.5_4_hr": "3.5h", "post_24_hr": "24h", "pre_exercise": "pre"}
TP_ORDER = list(TP)
GRP = {"Endur": "EE", "Resist": "RE", "Control": "CON"}


def parse_human(c):
    """-> (short label, kind, sort key). kind: exercise vs control | exercise post-pre | reference (non-exercise)."""
    lhs, rhs = [s.strip() for s in c.replace(" (delta-delta)", "").split(" - ")]
    g1, t1 = lhs.split(".", 1)
    g2, t2 = rhs.split(".", 1)
    g1, g2 = GRP[g1], GRP[g2]
    if t1 == "pre_exercise" and t2 == "pre_exercise":
        return f"{g1}−{g2} pre", "reference (non-exercise)", (2, 2, g1 + g2, 0)
    if g2 == "CON" and g1 != "CON":
        return f"{g1}−CON {TP[t1]}", "exercise vs control", (0, 0, g1, TP_ORDER.index(t1))
    if g1 == g2 and g1 != "CON":
        return f"{g1} post−pre {TP[t1]}", "exercise post−pre", (1, 0, g1, TP_ORDER.index(t1))
    if g1 == g2 == "CON":
        return f"CON post−pre {TP[t1]}", "reference (non-exercise)", (2, 0, "", TP_ORDER.index(t1))
    return f"{g1}−{g2} {TP[t1]}", "reference (non-exercise)", (2, 1, g1 + g2, TP_ORDER.index(t1))


RAT_ORDER = ["SKM-GN", "SKM-VL", "HEART", "WAT-SC", "LIVER", "KIDNEY", "LUNG", "CORTEX", "ADRNL", "BAT", "BLOOD",
             "COLON", "HIPPOC", "HYPOTH", "OVARY", "SMLINT", "SPLEEN", "TESTES", "VENACV"]
cols = []
for tissue in ["muscle", "adipose", "blood"]:
    for layer in hcon[hcon.tissue == tissue].layer.unique():
        cc = hcon[(hcon.tissue == tissue) & (hcon.layer == layer)]
        parsed = sorted([(parse_human(c), c, tp) for c, tp in zip(cc.contrast_short, cc.Timepoint)], key=lambda x: x[0][2])
        for (lab, kind, _), c, tp in parsed:
            cols.append(dict(dataset="human_acute", tissue=tissue, layer=layer, contrast=c, label=lab, kind=kind,
                             time=lab, block=f"human {tissue}"))
for tissue in RAT_ORDER:
    for layer in ["RNA", "PROT", "PHOSPHO"]:
        cc = sorted(rat[(rat.tissue == tissue) & (rat.layer == layer)].contrast.unique(),
                    key=lambda c: (c[0], int(c[2:-1])))
        for c in cc:
            cols.append(dict(dataset="rat_train", tissue=tissue, layer=layer, contrast=c, label=c.replace("_", " "),
                             kind="training vs sedentary", time=c[2:], block=f"rat {tissue}"))
cols = pd.DataFrame(cols)
cols["col"] = range(len(cols))
print(f"{len(cols)} comparison columns: " + ", ".join(f"{k}={v}" for k, v in cols.groupby("dataset").size().items()))

allx = pd.concat([rat, human], ignore_index=True)
long = sig.merge(cols, how="cross").merge(allx, on=["dataset", "tissue", "layer", "contrast", "gene"], how="left")
long["agreement"] = np.sign(long.logFC) * long.pah_dir
long = long.sort_values(["col", "gene"])
long.to_csv(f"{EV}/pah_everywhere_long.csv", index=False)

# ------------------------------------------------------------------ cross-check against join_table_v2 (03_join.py)
j = pd.read_csv("data/join_table_v2.csv", low_memory=False)
j = j[j.gene_symbol_human.isin(GENES)].rename(columns={"gene_symbol_human": "gene"})
j["tissue"] = j.tissue.replace({"VL": "muscle"})
cmp = long.merge(j[["dataset", "tissue", "layer", "contrast", "gene", "stat", "fdr_bh"]],
                 on=["dataset", "tissue", "layer", "contrast", "gene"], suffixes=("", "_join"))
cmp = cmp.dropna(subset=["stat", "stat_join"], how="all")
ok_stat = np.isclose(cmp.stat, cmp.stat_join, equal_nan=True)
ok_fdr = np.isclose(cmp.fdr_bh, cmp.fdr_bh_join, equal_nan=True)
xc = (f"Cross-check vs data/join_table_v2.csv (03_join.py) for the 25 signature genes, over the comparisons both "
      f"tables contain: {len(cmp)} gene x comparison cells; stat identical in {ok_stat.sum()}, fdr_bh identical in "
      f"{ok_fdr.sum()}.\n")
bad = cmp[~(ok_stat & ok_fdr)]
if len(bad):
    xc += "\nMismatches:\n\n" + bad[["dataset", "tissue", "layer", "contrast", "gene", "stat", "stat_join", "fdr_bh",
                                     "fdr_bh_join"]].to_markdown(index=False, floatfmt=".4g") + "\n"
open(f"{EV}/crosscheck_join_v2.md", "w").write(xc)
print(xc)


# ------------------------------------------------------------------ summaries
def counts(d):
    a = d[d.group != "phenotype"].agreement.dropna()
    a = a[a != 0]
    n, k = len(a), int((a < 0).sum())
    return n, k, n - k, (binomtest(k, n, 0.5).pvalue if n else np.nan)


rows = []
for ci, d in long.groupby("col"):
    n, k, s, p = counts(d)
    r = cols.iloc[ci].to_dict()
    rows.append({**r, "n_measured_19": n, "n_opposed_19": k, "n_same_19": s, "sign_p_19": p,
                 "n_fdrbh05": int((d.fdr_bh < 0.05).sum())})
colsum = pd.DataFrame(rows)
colsum.to_csv(f"{EV}/pah_everywhere_columns.csv", index=False)

# ranking unit = tissue x layer x time: rat pools the two sexes of a week; human = one contrast
rows = []
for (ds, t, l, tm, kind), d in long.groupby(["dataset", "tissue", "layer", "time", "kind"], sort=False):
    n, k, s, p = counts(d)
    sig19 = d[(d.group != "phenotype") & (d.fdr_bh < 0.05)]
    rows.append(dict(dataset=ds, tissue=t, layer=l, time=tm, kind=kind, n_columns=d.col.nunique(),
                     n_measured=n, n_opposed=k, n_same=s,
                     frac_opposed=k / n if n else np.nan, frac_same=s / n if n else np.nan, sign_p=p,
                     n_fdrbh05=len(sig19), n_fdrbh05_opposed=int((sig19.agreement < 0).sum()),
                     n_fdrbh05_same=int((sig19.agreement > 0).sum())))
rank = pd.DataFrame(rows)
ok = rank.sign_p.notna()
rank.loc[ok, "sign_q_bh"] = false_discovery_control(rank.loc[ok, "sign_p"], method="bh")
rank.to_csv(f"{EV}/pah_everywhere_ranking.csv", index=False)

blocks = []
for (b, l), d in long.groupby(["block", "layer"], sort=False):
    n, k, s, p = counts(d)
    blocks.append(dict(block=b, layer=l, n_columns=d.col.nunique(), cells_measured=n, cells_opposed=k,
                       frac_opposed=k / n if n else np.nan))
blocks = pd.DataFrame(blocks)
blocks.to_csv(f"{EV}/pah_everywhere_blocks.csv", index=False)

MIN_N = 10
view = rank[rank.n_measured >= MIN_N].copy()
view["unit"] = (view.dataset.map({"rat_train": "rat", "human_acute": "human"}) + " " + view.tissue + " " +
                view.layer + " " + view.time)
fmt = lambda df, c: df[["unit", "kind", "n_measured", "n_opposed", "n_same", c, "sign_p", "sign_q_bh",
                        "n_fdrbh05_opposed", "n_fdrbh05_same"]].to_markdown(index=False, floatfmt=".3g")
by_opp = view.sort_values(["frac_opposed", "n_measured", "sign_p"], ascending=[False, False, True])
by_same = view.sort_values(["frac_same", "n_measured", "sign_p"], ascending=[False, False, True])
view["majority"] = np.where(view.frac_opposed > 0.5, "opposed", "same")
n_units, n_units_all = len(view), len(rank)
qo = int(((view.sign_q_bh < 0.05) & (view.frac_opposed > 0.5)).sum())
qs = int(((view.sign_q_bh < 0.05) & (view.frac_same > 0.5)).sum())
po = int(((view.sign_p < 0.05) & (view.frac_opposed > 0.5)).sum())
ps = int(((view.sign_p < 0.05) & (view.frac_same > 0.5)).sum())
md = [
    "# PAH signature in every available comparison — ranking",
    "",
    f"Unit = tissue × layer × time (rat: the two sexes of one training week pooled, so up to 2 × 19 cells; human: one "
    f"contrast, up to 19 cells). Counts use the 19 Malenfant 2015 proteomic genes; a cell is *opposed* when "
    f"sign(exercise logFC) × PAH direction = −1 and *same* when +1. sign_p = two-sided binomial test of opposed vs 0.5; "
    f"sign_q_bh = BH over all {n_units_all} units. Genes (and the two sexes) are not independent, so p is optimistic. "
    f"Units with fewer than {MIN_N} measured cells are excluded from the ranking ({n_units_all - n_units} units, e.g. "
    f"blood Olink measures few signature genes); they remain in pah_everywhere_ranking.csv.",
    "",
    f"**Multiplicity.** {n_units} ranked units. Nominal sign_p < 0.05: {po} opposed-majority and {ps} "
    f"same-majority units (about {0.05 * n_units:.0f} expected by chance alone if units were independent). After BH: "
    f"{qo} opposed-majority and {qs} same-majority units at q < 0.05. "
    "Median opposition fraction by comparison kind: " + "; ".join(
        f"{k} {g.frac_opposed.median():.2f} ({len(g)} units)" for k, g in view.groupby("kind")) + ". "
    "BH-significant units by kind: " + "; ".join(
        f"{d}-majority / {k}: {n}" for (d, k), n in view[view.sign_q_bh < 0.05].groupby(["majority", "kind"]).size().items()) + ". "
    "Reference (non-exercise) comparisons — control time course, baseline group differences, EE−RE — are included as "
    "a negative control: they show how often the signature looks 'opposed' when no exercise contrast is involved.",
    "",
    "## Ranked by opposition fraction (exercise moves the gene AGAINST PAH) — top 30",
    "",
    fmt(by_opp.head(30), "frac_opposed"),
    "",
    "## Ranked by same-direction fraction (exercise moves the gene the SAME way as PAH) — top 30",
    "",
    fmt(by_same.head(30), "frac_same"),
    "",
    "## Per block (tissue × layer, all columns pooled)",
    "",
    blocks.to_markdown(index=False, floatfmt=".3g"),
    "",
]
open(f"{EV}/pah_everywhere_ranking.md", "w").write("\n".join(md))
print("\n".join(md[:8]))

# ------------------------------------------------------------------ figure
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle

INK, INK2, MUTED = "#0b0b0b", "#52514e", "#c3c2b7"
KIND_COL = {"exercise vs control": "#0b0b0b", "exercise post−pre": "#8a8985", "training vs sedentary": "#0b0b0b",
            "reference (non-exercise)": "#e4e3df"}
cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a61e1e"])
cmap.set_bad("white")
genes = sig.gene.tolist()
grp = sig.group.tolist()
gidx = {g: i for i, g in enumerate(genes)}
V = np.full((len(genes), len(cols)), np.nan)
Dt = np.zeros_like(V, dtype=bool)
for r in long.itertuples():
    if pd.notna(r.stat):
        V[gidx[r.gene], r.col] = r.agreement * min(abs(r.stat), 4)
        Dt[gidx[r.gene], r.col] = r.fdr_bh < 0.05
blk_stats = {(b.block, b.layer): b for b in blocks.itertuples()}

# pack whole tissue blocks into strips of <= MAXC columns (plus 1-column gaps between blocks)
MAXC = 100
strips, cur, width = [], [], 0
for b in cols.block.unique():
    nb = int((cols.block == b).sum())
    if cur and width + nb + 1 > MAXC:
        strips.append(cur); cur, width = [], 0
    cur.append(b); width += nb + 1
strips.append(cur)

CELL = 0.125  # inches
strip_w = [sum(int((cols.block == b).sum()) + 1 for b in s) for s in strips]
W = max(strip_w) * CELL + 2.6
H_STRIP = (len(genes) + 1.2) * CELL + 1.75
fig = plt.figure(figsize=(W, H_STRIP * len(strips) + 1.6))
top = 1 - 0.95 / fig.get_figheight()
h_ax = (len(genes) + 1.2) * CELL / fig.get_figheight()
for si, s in enumerate(strips):
    x0 = 1.9 / W
    y_top = top - si * H_STRIP / fig.get_figheight() - 0.75 / fig.get_figheight()
    ax = fig.add_axes([x0, y_top - h_ax, strip_w[si] * CELL / W, h_ax])
    # assemble the strip matrix with a blank spacer column after every block; row -1 = contrast-kind band
    xs, xcols = [], []
    pos = 0
    for b in s:
        for ci in cols.index[cols.block == b]:
            xs.append(pos); xcols.append(ci); pos += 1
        pos += 1
    ncol = pos
    M = np.full((len(genes), ncol), np.nan)
    Dm = np.zeros_like(M, dtype=bool)
    for x, ci in zip(xs, xcols):
        M[:, x] = V[:, ci]; Dm[:, x] = Dt[:, ci]
    ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-4, vmax=4, aspect="equal",
              extent=(-.5, ncol - .5, len(genes) - .5, -.5))
    for x, ci in zip(xs, xcols):
        for y in np.where(np.isnan(V[:, ci]))[0]:
            ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False, hatch="//////", edgecolor=MUTED, lw=0))
        ax.add_patch(Rectangle((x - .5, -1.6), 1, 0.9, color=KIND_COL[cols.kind[ci]], lw=0, clip_on=False))
    yy, xx = np.where(Dm)
    ax.scatter(xx, yy, s=4.5, color=INK, zorder=3, lw=0)
    ax.set_xlim(-.5, ncol - .5); ax.set_ylim(len(genes) - .5, -1.7)
    ax.set_xticks(np.arange(-.5, ncol), minor=True); ax.set_yticks(np.arange(-.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", lw=0.5); ax.tick_params(which="minor", length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_yticks(range(len(genes)), [f"{g} ({'↓' if d < 0 else '↑'})" for g, d in zip(genes, sig.pah_dir)], fontsize=6)
    ax.set_xticks(xs, [cols.label[ci] for ci in xcols], rotation=90, fontsize=5)
    for tl, ci in zip(ax.get_xticklabels(), xcols):
        if cols.kind[ci] == "reference (non-exercise)":
            tl.set_color("#8a8985")
    ax.tick_params(axis="both", length=0, pad=1.5)
    # block (tissue) and layer headers with pooled opposed / measured
    for b in s:
        bx = [x for x, ci in zip(xs, xcols) if cols.block[ci] == b]
        ax.vlines(bx[0] - .5, -.5, len(genes) - .5, color=INK, lw=1.0)
        ax.vlines(bx[-1] + .5, -.5, len(genes) - .5, color=INK, lw=1.0)
        for i in range(1, len(genes)):
            if grp[i] != grp[i - 1]:
                ax.hlines(i - .5, bx[0] - .5, bx[-1] + .5, color=INK, lw=0.8)
        ax.text((bx[0] + bx[-1]) / 2, -5.3, b, ha="center", va="bottom", fontsize=8, fontweight="bold")
        for l in cols[cols.block == b].layer.unique():
            lx = [x for x, ci in zip(xs, xcols) if cols.block[ci] == b and cols.layer[ci] == l]
            st = blk_stats[(b, l)]
            fr = "" if not st.cells_measured else f" {100 * st.frac_opposed:.0f}%"
            ax.text((lx[0] + lx[-1]) / 2, -3.0, f"{l}\n{st.cells_opposed}/{st.cells_measured}{fr}", ha="center",
                    va="bottom", fontsize=5.5, color=INK2, linespacing=1.0)
            if lx[0] != bx[0]:
                ax.vlines(lx[0] - .5, -.5, len(genes) - .5, color=INK, lw=0.5)

# legend + colourbar + caption
cax = fig.add_axes([1.9 / W, 1 - 0.55 / fig.get_figheight(), 2.6 / W, 0.1 / fig.get_figheight()])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(-4, 4), cmap=cmap), cax=cax, orientation="horizontal")
cb.ax.tick_params(labelsize=6)
cb.set_label("agreement × |stat| (capped at 4):  ← exercise OPPOSES PAH  |  same as PAH →", fontsize=6.5)
handles = [plt.Line2D([], [], marker="o", ls="", color=INK, ms=3, label="BH FDR < 0.05 (fdr_bh)"),
           Patch(facecolor="white", edgecolor=MUTED, hatch="//////", label="not measured"),
           Patch(color=KIND_COL["exercise vs control"], label="band: exercise vs control / rat trained vs sedentary"),
           Patch(color=KIND_COL["exercise post−pre"], label="band: human exercise post − pre (within group)"),
           Patch(color=KIND_COL["reference (non-exercise)"], label="band: non-exercise reference (control time course, "
                                                                    "baseline, EE−RE); grey labels")]
fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(5.0 / W, 1 - 0.25 / fig.get_figheight()), ncol=3,
           fontsize=6.5, frameon=False)
fig.text(1.9 / W, 1 - 0.12 / fig.get_figheight(),
         "PAH skeletal-muscle signature (Malenfant 2015) vs the exercise response in every MoTrPAC tissue, layer and "
         "comparison", fontsize=10, fontweight="bold", va="top")
cap = ("Signal searched in every available comparison. Rows: the 19 Malenfant 2015 PAH-muscle proteins + 6 enzyme-"
       "phenotype genes (↓/↑ = direction in PAH). Columns: every gene-level comparison available — human acute exercise "
       f"(muscle, adipose, blood; all {int((cols.dataset == 'human_acute').sum())} tissue × layer × contrast columns, "
       f"including non-exercise reference contrasts) and rat endurance training (all 19 tissues, RNA / PROT / PHOSPHO, "
       f"sex × week; {int((cols.dataset == 'rat_train').sum())} columns). Colour = sign(exercise logFC) × PAH direction × "
       "min(|stat|, 4); dot = BH FDR < 0.05 recomputed within tissue × layer × contrast. Header numbers: cells opposed / "
       "measured among the 19 proteomic genes, pooled over the block's columns. Group-level contrasts across genes, not "
       "per-animal or per-person predictions; with this many comparisons some blocks will look opposed or aligned by "
       "chance — see pah_everywhere_ranking.md.")
import textwrap
fig.text(1.9 / W, 0.12 / fig.get_figheight(), "\n".join(textwrap.wrap(cap, 230)), fontsize=6.5, va="bottom", color=INK2)
fig.savefig(f"{EV}/pah_everywhere.png", dpi=200)
plt.close(fig)
open(f"{EV}/captions.md", "w").write("# Captions\n\n- `pah_everywhere.png` — " + cap + "\n")
print(f"wrote {EV}/pah_everywhere.png ({len(strips)} strips, {len(cols)} columns)")
