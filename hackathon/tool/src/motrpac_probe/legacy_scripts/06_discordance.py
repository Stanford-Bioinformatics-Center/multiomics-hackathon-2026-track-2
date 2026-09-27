#!/usr/bin/env python3
"""06_discordance.py — how well does the RNA response predict the PROT response, by time and gene subset?

Run from hackathon/:  python3 scripts/06_discordance.py
Input : data/join_table_v2.csv, data/raw/genesets_v1.csv
Output: data/discordance_spearman_v1.csv (write-once), fig/tables/discordance_spearman.md, FACTSHEET.md (re-rendered),
        fig/discordance_human.png, fig/discordance_rat.png, captions in fig/captions.md

For each dataset x contrast: Spearman rho between the gene-level RNA stat and PROT stat over
genes measured in both layers, for (a) all genes, (b) KEGG OXPHOS (hsa00190), (c) the 19 Malenfant
PAH proteins. 95% CI from Fisher z with the Bonett-Wright SE for Spearman, sqrt((1 + rho^2/2)/(n-3)).
Pairs: human VL RNA x VL PROT; rat SKM-GN RNA x SKM-GN PROT; rat SKM-VL RNA x SKM-GN PROT (cross-muscle,
same animals); rat HEART RNA x HEART PROT.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(__file__))
JOIN = "data/join_table_v2.csv"
OUT = "data/discordance_spearman_v1.csv"
TAB = "fig/tables"
PAH19 = ["NDUFA9", "UQCRC2", "UQCRC1", "ATP5MG", "ATP5F1B", "IDH2", "OGDH", "SLC25A4", "ECH1", "GLO1", "FBP2",
         "PEBP1", "ATP2A1", "MYH1", "MYH2", "MYH7", "MYLPF", "APOA1", "PDLIM3"]

j = pd.read_csv(JOIN)
gs = pd.read_csv("data/raw/genesets_v1.csv")
SUBSETS = {"all genes": None,
           "OXPHOS (KEGG)": set(gs.gene_symbol[gs.gs_name == "KEGG_OXIDATIVE_PHOSPHORYLATION"]),
           "PAH 19": set(PAH19)}
PAIRS = [("human_acute", "VL", "VL"), ("rat_train", "SKM-GN", "SKM-GN"), ("rat_train", "SKM-VL", "SKM-GN"),
         ("rat_train", "HEART", "HEART")]


def rho_ci(x, y):
    n = len(x)
    if n < 5:
        return n, np.nan, np.nan, np.nan, np.nan
    r, p = spearmanr(x, y)
    se = np.sqrt((1 + r ** 2 / 2) / (n - 3))
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    return n, r, p, np.tanh(z - 1.96 * se), np.tanh(z + 1.96 * se)


rows = []
for ds, tr, tp in PAIRS:
    for c in j[(j.dataset == ds) & (j.tissue == tp)].contrast.unique():
        r = j[(j.dataset == ds) & (j.tissue == tr) & (j.layer == "RNA") & (j.contrast == c)]
        p = j[(j.dataset == ds) & (j.tissue == tp) & (j.layer == "PROT") & (j.contrast == c)]
        m = r.merge(p, on="gene_symbol_human", suffixes=("_rna", "_prot")).dropna(subset=["stat_rna", "stat_prot"])
        for sname, genes in SUBSETS.items():
            mm = m if genes is None else m[m.gene_symbol_human.isin(genes)]
            n, rho, pv, lo, hi = rho_ci(mm.stat_rna, mm.stat_prot)
            rows.append(dict(dataset=ds, RNA_tissue=tr, PROT_tissue=tp, contrast=c, subset=sname, n_genes=n,
                             spearman=rho, p=pv, ci_low=lo, ci_high=hi,
                             frac_same_sign=float((np.sign(mm.logFC_rna) == np.sign(mm.logFC_prot)).mean())))
res = pd.DataFrame(rows)
if not os.path.exists(OUT):
    res.to_csv(OUT, index=False)
    os.chmod(OUT, 0o444)
    print(f"wrote {OUT} ({len(res)} rows)")
wide = res.assign(cell=res.apply(lambda r: f"{r.spearman:+.2f} [{r.ci_low:+.2f},{r.ci_high:+.2f}] n={r.n_genes}", axis=1)) \
    .pivot_table(index=["dataset", "RNA_tissue", "PROT_tissue", "contrast"], columns="subset", values="cell",
                 aggfunc="first", sort=False)[list(SUBSETS)].reset_index()
s = wide.to_markdown(index=False)
open(f"{TAB}/discordance_spearman.md", "w").write(s + "\n")
print(s)

# ---------------------------------------------------------------- figures
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
INK2, GRID = "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.spines.top": False,
                     "axes.spines.right": False, "figure.dpi": 150, "savefig.bbox": "tight", "legend.frameon": False})
SUB_C = dict(zip(SUBSETS, C))


def caption(fn, text):
    caps = {}
    if os.path.exists("fig/captions.md"):
        for line in open("fig/captions.md"):
            if line.startswith("- `"):
                caps[line.split("`")[1]] = line.rstrip("\n")
    caps[fn] = f"- `{fn}` — {text}"
    with open("fig/captions.md", "w") as f:
        f.write("# Figure captions\n\n" + "\n".join(caps[k] for k in sorted(caps)) + "\n")


# human: x = contrast, series = gene subset
hc = ["EE_vs_CON_15to45min", "EE_vs_CON_3.5to4h", "EE_vs_CON_24h", "RE_vs_CON_15to45min", "RE_vs_CON_3.5to4h",
      "RE_vs_CON_24h"]
fig, ax = plt.subplots(figsize=(8, 4))
for k, (sname, col) in enumerate(SUB_C.items()):
    d = res[(res.dataset == "human_acute") & (res.subset == sname)].set_index("contrast").loc[hc]
    x = np.arange(len(hc)) + (k - 1) * 0.22
    ax.errorbar(x, d.spearman, yerr=[d.spearman - d.ci_low, d.ci_high - d.spearman], fmt="o", ms=6, color=col,
                elinewidth=1.2, capsize=0, label=f"{sname} (n≈{int(d.n_genes.median())})")
ax.axhline(0, color=INK2, lw=0.8); ax.axvline(2.5, color="#0b0b0b", lw=0.8)
ax.set_xticks(range(len(hc)), [c.replace("_vs_CON_", " vs CON\n") for c in hc], fontsize=8)
ax.set_ylabel("Spearman ρ (RNA stat vs PROT stat)"); ax.set_ylim(-1, 1)
ax.set_title("Human vastus lateralis, one exercise bout: RNA–protein agreement", loc="left", fontsize=10)
ax.legend(loc="lower left", fontsize=8, ncol=3); ax.grid(axis="x", visible=False)
fig.tight_layout(); fig.savefig("fig/discordance_human.png"); plt.close(fig)
caption("discordance_human.png", "Human VL: Spearman rho between gene-level RNA and PROT test statistics (moderated t) "
        "per delta-delta contrast, for all genes measured in both layers, KEGG OXPHOS genes, and the 19 PAH-signature "
        "proteins; bar = 95% CI (Fisher z, Bonett-Wright SE). Point = one contrast x gene subset.")

# rat: x = weeks, rows = sex, cols = subset, series = tissue pair
fig, axes = plt.subplots(2, 3, figsize=(11, 6.5), sharex=True, sharey=True)
pairs_r = [p for p in PAIRS if p[0] == "rat_train"]
for ri, sx in enumerate("FM"):
    for ci, sname in enumerate(SUBSETS):
        ax = axes[ri, ci]
        for k, (ds, tr, tp) in enumerate(pairs_r):
            d = res[(res.RNA_tissue == tr) & (res.PROT_tissue == tp) & (res.subset == sname) &
                    (res.contrast.str.startswith(sx))].copy()
            d["wk"] = d.contrast.str[2:-1].astype(int)
            d = d.sort_values("wk")
            off = (k - 1) * 0.18
            lab = f"RNA {tr} × PROT {tp}"
            ax.plot(d.wk + off, d.spearman, lw=2, color=C[k], ls="-" if tr == tp else "--", marker="o", ms=5,
                    label=lab)
            ax.vlines(d.wk + off, d.ci_low, d.ci_high, color=C[k], lw=1, alpha=0.7)
        ax.axhline(0, color=INK2, lw=0.8); ax.set_ylim(-1, 1); ax.set_xticks([1, 2, 4, 8])
        if ri == 0:
            ax.set_title(sname, loc="left", fontsize=10)
        if ci == 0:
            ax.set_ylabel(f"{'Females' if sx == 'F' else 'Males'}\nSpearman ρ (RNA vs PROT)")
        if ri == 1:
            ax.set_xlabel("weeks of training")
axes[0, 0].legend(loc="lower left", fontsize=7.5)
fig.tight_layout(); fig.savefig("fig/discordance_rat.png"); plt.close(fig)
caption("discordance_rat.png", "Rat training: Spearman rho between gene-level RNA (DESeq2 z) and PROT (limma t) "
        "statistics vs weeks of training, by sex (rows) and gene subset (columns: all genes, KEGG OXPHOS, PAH 19); "
        "line = tissue pairing (dashed = SKM-VL RNA vs SKM-GN PROT, different muscles of the same animals); "
        "vertical bar = 95% CI (Fisher z, Bonett-Wright SE).")

# ---------------------------------------------------------------- re-render FACTSHEET.md
# FACTSHEET.md = fig/tables/FACTSHEET_template.md with every {{name}} replaced by fig/tables/name.md,
# so it always reflects the latest run of 04/05/06. Edit prose in the template, not in FACTSHEET.md.
import re
tmpl = open(f"{TAB}/FACTSHEET_template.md").read()
open("FACTSHEET.md", "w").write(re.sub(r"\{\{(\w+)\}\}", lambda m: open(f"{TAB}/{m.group(1)}.md").read().strip(), tmpl))
print("rendered FACTSHEET.md")


# ======================================================================= NARRATIVE.md + deck_extracts/
# Every number in NARRATIVE.md comes from a registry REG via {{n:key}}; tables via {{t:name}} (their numeric
# cells are registered too). Rendering fails on an unknown key, so nothing uncited can slip in.
# Output: NARRATIVE.md, deck_extracts/{*.png, deck_numbers.csv, captions.md}
import shutil
from scipy.stats import binomtest

REG = {}


def reg(key, value, source, figure="", note=""):
    if isinstance(value, (np.floating, float)):
        value = float(value)
    elif isinstance(value, (np.integer,)):
        value = int(value)
    REG[key] = dict(key=key, value=value, source_script=source, figure=figure, note=note)
    return value


for fn in ["explore_numbers.csv", "robust_numbers.csv", "crosscheck_numbers.csv", "hyperemia_numbers.csv"]:
    if os.path.exists(f"{TAB}/{fn}"):
        for r in pd.read_csv(f"{TAB}/{fn}").fillna("").itertuples():
            reg(r.key, r.value, r.source_script, r.figure, r.note)

S06, S05, S04, S03 = ("scripts/06_discordance.py deck", "scripts/05_pah_grid.py", "scripts/04_explore.py", "scripts/03_join.py")
TABLES = {}


def table(name, df, source, figure="", fmt="{:.3g}"):
    """Register every numeric cell, keep a markdown rendering."""
    out = df.copy()
    for ri, row in df.iterrows():
        for col, v in row.items():
            if isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool) and pd.notna(v):
                reg(f"{name}[{ri}|{col}]", v, source, figure)
                out.at[ri, col] = str(int(v)) if float(v).is_integer() else fmt.format(v)
    TABLES[name] = out.to_markdown(index=False)


# ---- P0 design
hc = pd.read_csv("data/raw/human_muscle_sample_counts_v1.csv")
tp = ["pre_exercise", "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
hc["layer"] = hc.assay.map({"transcript-rna-seq": "RNA", "prot-pr": "PROT", "prot-ph": "PHOSPHO"})
hc["group"] = hc.randomGroupCode.str.replace("ADU", "")
t_n = hc.pivot_table(index=["layer", "group"], columns="Timepoint", values="n_median", aggfunc="first")[tp]
t_n.columns = ["pre", "15-45 min", "3.5-4 h", "24 h"]
table("human_n", t_n.reset_index().sort_values(["layer", "group"], key=lambda s: s.map({"RNA": 0, "PROT": 1, "PHOSPHO": 2}).fillna(s.map({"Control": 0, "Endur": 1, "Resist": 2}))).reset_index(drop=True), S04, "")
for layer, fn in [("RNA", "rna"), ("PROT", "prot"), ("PHOSPHO", "phospho")]:
    x = pd.read_csv(f"data/raw/human_muscle_{fn}_DA_v1.csv", usecols=["contrast_short", "feature_id", "adj_p_value"])
    reg(f"hum_nfeat_{layer}", x.feature_id.nunique(), "scripts/01_load_human.R")
    reg("hum_n_contrasts", x.contrast_short.nunique(), "scripts/01_load_human.R")
    if layer != "PHOSPHO":
        for c, short in {v: k for k, v in {
                "Endur.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)": "EE_vs_CON_15to45min",
                "Endur.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)": "EE_vs_CON_3.5to4h",
                "Endur.post_24_hr - Control.post_24_hr (delta-delta)": "EE_vs_CON_24h",
                "Resist.post_15_30_45_min - Control.post_15_30_45_min (delta-delta)": "RE_vs_CON_15to45min",
                "Resist.post_3.5_4_hr - Control.post_3.5_4_hr (delta-delta)": "RE_vs_CON_3.5to4h",
                "Resist.post_24_hr - Control.post_24_hr (delta-delta)": "RE_vs_CON_24h"}.items()}.items():
            s_ = x[x.contrast_short == short]
            reg(f"hum_da_{layer}_{c}", int((s_.adj_p_value < 0.05).sum()), S04, "", "published BH<0.05, feature level")
            reg(f"hum_dapct_{layer}_{c}", round(100 * (s_.adj_p_value < 0.05).mean(), 2), S04)
rc = pd.read_csv("data/raw/rat_sample_counts_v1.csv")
reg("rat_n_min", int(rc.n_animals.min()), "scripts/02_load_rat.R"); reg("rat_n_max", int(rc.n_animals.max()), "scripts/02_load_rat.R")
hcol = ["EE_vs_CON_15to45min", "EE_vs_CON_3.5to4h", "EE_vs_CON_24h", "RE_vs_CON_15to45min", "RE_vs_CON_3.5to4h", "RE_vs_CON_24h"]
t_da = pd.DataFrame([{"layer": l, **{c.replace("_vs_CON_", " "): f"{REG[f'hum_da_{l}_{c}']['value']} ({REG[f'hum_dapct_{l}_{c}']['value']}%)" for c in hcol}}
                     for l in ("RNA", "PROT")])
TABLES["human_da"] = t_da.to_markdown(index=False)

# ---- P0 RNA-PROT Spearman (human)
sp = res[res.dataset == "human_acute"]
for rr in sp.itertuples():
    reg(f"rho_hum_{rr.contrast}_{rr.subset.split(' ')[0]}", round(rr.spearman, 2), S06, "discordance_human.png")
    reg(f"rholo_hum_{rr.contrast}_{rr.subset.split(' ')[0]}", round(rr.ci_low, 2), S06, "discordance_human.png")
    reg(f"rhohi_hum_{rr.contrast}_{rr.subset.split(' ')[0]}", round(rr.ci_high, 2), S06, "discordance_human.png")
for rr in res[(res.dataset == "rat_train") & (res.RNA_tissue == "SKM-GN")].itertuples():
    reg(f"rho_rat_{rr.contrast}_{rr.subset.split(' ')[0]}", round(rr.spearman, 2), S06, "discordance_rat.png")
for rr in res[(res.dataset == "rat_train") & (res.RNA_tissue == "SKM-VL")].itertuples():
    reg(f"rho_ratVL_{rr.contrast}_{rr.subset.split(' ')[0]}", round(rr.spearman, 2), S06, "discordance_rat.png")

# ---- P0 artifacts + FDR choice
exc = pd.read_csv("data/rat_rna_excluded_v1.csv")
reg("onoff_rows", len(exc), S03); reg("onoff_features", exc.feature_ID.nunique(), S03)
j1 = pd.read_csv("data/join_table_v1.csv")
r1_ = j1[(j1.dataset == "rat_train") & (j1.layer == "RNA")]
ex = r1_.loc[r1_.logFC.abs().idxmax()]
reg("onoff_example_logFC", round(float(ex.logFC), 1), S03, "", f"{ex.gene_symbol_human} {ex.tissue} {ex.contrast} in join v1")
ONOFF_EX = f"{ex.gene_symbol_human} ({ex.tissue}, {ex.contrast})"
k = ["dataset", "tissue", "layer", "contrast", "gene_symbol_human"]
mm = j1.merge(j, on=k, how="left", suffixes=("_1", "_2"), indicator=True)
reg("v2_rows_removed", int((mm._merge == "left_only").sum()), S03)
reg("v2_rows_feature_changed", int(((mm._merge == "both") & (mm.feature_id_1 != mm.feature_id_2)).sum()), S03)
reg("onoff_max_abs_z", round(float(exc.zscore.abs().max()), 1), S03)
rows = []
for t_ in ("SKM-GN", "SKM-VL", "HEART"):
    for l_ in ("RNA", "PROT"):
        d = j[(j.dataset == "rat_train") & (j.tissue == t_) & (j.layer == l_)]
        if len(d):
            rows.append(dict(tissue=t_, layer=l_, published_BY_total=int((d.fdr < 0.05).sum()),
                             BH_total=int((d.fdr_bh < 0.05).sum()),
                             BY_8w=int((d[d.contrast.str.endswith("8w")].fdr < 0.05).sum()),
                             BH_8w=int((d[d.contrast.str.endswith("8w")].fdr_bh < 0.05).sum())))
table("fdr_choice", pd.DataFrame(rows), S06, "rat_nDA_over_time.png")

# ---- P1 DOWN-9 per gene (human)
D9 = ["NDUFA9", "UQCRC2", "UQCRC1", "ATP5MG", "ATP5F1B", "IDH2", "OGDH", "SLC25A4", "ECH1"]
hj = j[j.dataset == "human_acute"]
for layer in ("RNA", "PROT"):
    rows = []
    for g in D9:
        row = {"gene": g}
        for c in hcol:
            d = hj[(hj.layer == layer) & (hj.contrast == c) & (hj.gene_symbol_human == g)]
            row[c.replace("_vs_CON_", " ")] = f"{d.logFC.iloc[0]:+.2f} ({d.fdr.iloc[0]:.2g})" if len(d) else "NA"
            if len(d):
                reg(f"d9_{layer}_{g}_{c}_logFC", round(float(d.logFC.iloc[0]), 2), S06, "pah_grid_v2.png")
                reg(f"d9_{layer}_{g}_{c}_fdr", float(f"{d.fdr.iloc[0]:.2g}"), S06, "pah_grid_v2.png")
        rows.append(row)
    TABLES[f"down9_{layer}"] = pd.DataFrame(rows).to_markdown(index=False)
pgrid = pd.read_csv(f"{TAB}/pah_pathway_grid.csv")
for rr in pgrid.itertuples():
    reg(f"pw_{rr.gene_set}_{rr.tissue}_{rr.layer}_{rr.contrast}", round(rr.camera_t, 1), S05, "pah_pathway_grid.png")
    reg(f"pwq_{rr.gene_set}_{rr.tissue}_{rr.layer}_{rr.contrast}", float(f"{rr.fdr_bh_row:.2g}"), S05, "pah_pathway_grid.png")
cam_rows = []
for layer in ("RNA", "PROT"):
    for c in hcol:
        tag = f"{layer}_{c.replace('_vs_CON_', '_')}"
        cam_rows.append({"layer": layer, "contrast": c.replace("_vs_CON_", " "),
                         "MoTrPAC CAMERA_RESULTS (GOBP OXPHOS)": REG[f"XC_pkgcam_{tag}"]["value"],
                         "ours: same GOBP set, join stat": REG[f"XC_joingobp_{tag}"]["value"],
                         "ours: KEGG hsa00190, join stat (B3)": REG[f"XC_joinkegg_{tag}"]["value"]})
TABLES["camera_side_by_side"] = pd.DataFrame(cam_rows).to_markdown(index=False)

# ---- P2 the UP-10
U10 = ["GLO1", "FBP2", "PEBP1", "ATP2A1", "MYH1", "MYH2", "MYH7", "MYLPF", "APOA1", "PDLIM3"]
rows = []
for c in hcol[1:3] + hcol[4:6]:
    for layer in ("RNA", "PROT"):
        d = hj[(hj.layer == layer) & (hj.contrast == c) & hj.gene_symbol_human.isin(U10)]
        opp = int((d.logFC < 0).sum()); n = len(d)
        p = binomtest(opp, n, 0.5).pvalue
        rows.append(dict(contrast=c.replace("_vs_CON_", " "), layer=layer, n_measured=n, n_opposed=opp,
                         n_fdrbh05=int((d.fdr_bh < 0.05).sum()), sign_test_p=round(p, 3)))
table("up10", pd.DataFrame(rows), S06, "pah_grid_v2.png")
for layer in ("RNA", "PROT"):
    d = hj[(hj.layer == layer) & (hj.contrast == "EE_vs_CON_24h") & hj.gene_symbol_human.isin(D9)]
    reg(f"down9_opp_EE24_{layer}", int((d.logFC > 0).sum()), S06, "pah_grid_v2.png")

# ---- P3 rat
rows = []
for sx in "FM":
    for w in ("1w", "2w", "4w", "8w"):
        c = f"{sx}_{w}"
        row = {"contrast": c}
        for t_, l_ in [("SKM-GN", "RNA"), ("SKM-GN", "PROT"), ("SKM-VL", "RNA")]:
            v = pgrid[(pgrid.gene_set == "KEGG_OXIDATIVE_PHOSPHORYLATION") & (pgrid.tissue == t_) & (pgrid.layer == l_) & (pgrid.contrast == c)].camera_t.iloc[0]
            row[f"OXPHOS t {t_} {l_}"] = round(float(v), 1)
        row["RNA-PROT rho (SKM-GN, all genes)"] = REG[f"rho_rat_{c}_all"]["value"]
        d = j[(j.dataset == "rat_train") & (j.tissue == "SKM-GN") & (j.layer == "PROT") & (j.contrast == c)]
        row["PROT genes BH<0.05"] = int((d.fdr_bh < 0.05).sum())
        rows.append(row)
table("rat_time", pd.DataFrame(rows), S06, "pah_pathway_grid.png / discordance_rat.png")


# ---- literature inputs + extra aggregates
MAL = "Malenfant 2015 J Mol Med (inputs hard-coded in scripts/05_pah_grid.py)"
for key, v in [("mal_n_pah", 4), ("mal_n_ctrl", 4), ("mal_n_proteins", 231), ("mal_n_down", 9), ("mal_n_up", 10),
               ("mal_down_min", 0.57), ("mal_down_max", 0.77), ("mal_up_min", 1.26), ("mal_up_max", 2.42), ("mal_n_pheno", 6)]:
    reg(key, v, MAL)
ox = pgrid[pgrid.gene_set == "KEGG_OXIDATIVE_PHOSPHORYLATION"]
for (t_, l_), d in ox[ox.dataset == "rat_train"].groupby(["tissue", "layer"]):
    reg(f"rat_ox_upsig_{t_}_{l_}", int(((d.camera_t > 0) & (d.p < 0.05)).sum()), S05, "pah_pathway_grid.png", "columns with OXPHOS t>0, p<0.05")
    reg(f"rat_ox_ncol_{t_}_{l_}", len(d), S05, "pah_pathway_grid.png")
pw1 = pd.read_csv("data/pah_pathway_v1.csv")
kg = pw1[pw1.gene_set == "KEGG_GLYCOLYSIS_GLUCONEOGENESIS"].set_index(["dataset", "tissue", "layer", "contrast"])
hg = pw1[pw1.gene_set == "HALLMARK_GLYCOLYSIS"].set_index(["dataset", "tissue", "layer", "contrast"])
ks = kg[kg.p < 0.05].index
reg("glyc_kegg_nsig", len(ks), S05, "pah_pathway.png", "KEGG glycolysis columns with p<0.05")
reg("glyc_kegg_ncol", len(kg), S05, "pah_pathway.png")
reg("glyc_hall_replicated", int((hg.loc[ks].p < 0.05).sum()), S05, "pah_pathway.png", "of those, Hallmark p<0.05")
reg("glyc_hall_minp_of_ks", round(float(hg.loc[ks].p.min()), 2), S05, "pah_pathway.png")
reg("glyc_kegg_nopp", int(((kg.p < 0.05) & (kg.camera_t < 0)).sum()), S05, "pah_pathway.png", "significant and opposes PAH (down)")

TM = "teammate's independent scripts (numbers reported to us; not recomputed by them here)"
for key, v in [("tm_prot_matched", 9), ("tm_prot_sig", 0), ("tm_prot_tests", 27), ("tm_rna_pos_15", 0), ("tm_rna_pos_35", 7),
               ("tm_rna_pos_24", 9), ("tm_rna_q_24", 7), ("tm_camera_set_n", 139)]:
    reg(key, v, TM)
reg("rat_hours_after_last_bout", 48, "MoTrPAC rat training design (Nature 2024 paper; not in package data)")
reg("R2a_n_within", int(REG["R2a_n_MYH7"]["value"] // 3), S05, "robust_composition.png", "columns per tissue x layer")

# ---- P4 lag (other session; values transcribed from lag/RESULTS.md, not recomputed here)
LAGSRC = "lag/RESULTS.md (scripts/07_lag.py, other session)"
for key, v, note in [("lag_dR2_raw_pct", 0.34, "SKM-GN pooled, 2w RNA -> 8w PROT beyond 2w PROT, held-out, raw"),
                     ("lag_dR2_rint_pct", 0.78, "same, rank-INT"), ("lag_heart_rint_pct", 0.15, "HEART rank-INT held-out"),
                     ("lag_same_time_raw", 0.022, "8w RNA held-out dR2 raw"), ("lag_2w_raw", 0.0034, "2w RNA held-out dR2 raw"),
                     ("lag_4w_cond_raw", 0.0021, "4w RNA beyond 8w RNA, held-out raw"), ("lag_nperm", 1000, "permutations")]:
    reg(key, v, LAGSRC, "lag/fig/*.png", note)

# ---- render
pd.DataFrame(REG.values()).to_csv(f"{TAB}/all_numbers.csv", index=False)  # full registry (cited or not)
tmpl = open(f"{TAB}/NARRATIVE_template.md").read()
used = set()


def fmt(v):
    if isinstance(v, str):
        return v
    if float(v).is_integer():
        return f"{int(v):,}" if abs(v) >= 10000 else str(int(v))
    return f"{v:g}"


def sub(m):
    kind, key = m.group(1), m.group(2)
    if kind == "n":
        if key not in REG:
            raise KeyError(f"NARRATIVE cites unknown number key: {key}")
        used.add(key)
        return fmt(REG[key]["value"])
    if kind == "t":
        used.update(k for k in REG if k.startswith(f"{key}["))
        return TABLES[key]
    if kind == "s":
        return {"ONOFF_EX": ONOFF_EX}[key]


text = re.sub(r"\{\{(n|t|s):([^}]+)\}\}", sub, tmpl)
open("NARRATIVE.md", "w").write(text)
os.makedirs("deck_extracts", exist_ok=True)
nums = pd.DataFrame([REG[k] for k in sorted(used)])
nums.to_csv("deck_extracts/deck_numbers.csv", index=False)
print(f"rendered NARRATIVE.md: {len(used)} registered numbers cited")

# ---- deck figures (copied byte-for-byte: native size/aspect) + one-sentence captions
DECK = {
    "logfc_distributions.png": ("none (descriptive)", "Each box is the spread of exercise logFC for one layer and contrast; human protein boxes are about half the height of RNA boxes."),
    "robust_compression.png": ("paired Wilcoxon |PROT| vs |RNA|", "Points below the dashed line are genes whose protein effect is smaller than their RNA effect; bars compare the largest human effects per layer."),
    "rna_vs_prot_scatter.png": ("Spearman rho", "Each point is a gene measured in both layers; a cloud along the diagonal would mean RNA and protein change together."),
    "rna_prot_spearman_by_contrast.png": ("Spearman rho", "Bar height is how well RNA and protein changes agree across all genes in one contrast (0 = no agreement)."),
    "rat_nDA_over_time.png": ("BH vs published BY FDR < 0.05", "Lines climbing to the right mean more genes change the longer rats train; compare the two rows to see what the FDR choice does."),
    "pah_grid.png": ("two-sided binomial sign test per column", "Blue cells are genes that exercise moves against their PAH direction, red cells the same way as PAH; dots mark BH FDR < 0.05."),
    "pah_grid_v2.png": ("two-sided binomial sign test per column", "Same as pah_grid with human 3.5 h added: read the human block (left) against the rat blocks for layer x time."),
    "pah_pathway.png": ("limma cameraPR (Python port)", "Points above zero mean the gene set moved up with exercise relative to all other genes; for OXPHOS, up = opposes PAH."),
    "pah_pathway_grid.png": ("limma cameraPR (Python port), BH within row", "Red = pathway up with exercise, blue = down; compare the colour to the PAH arrow in each row label."),
    "discordance_human.png": ("Spearman rho with Fisher-z 95% CI", "Higher points mean RNA and protein changes agree more for that gene subset and human contrast."),
    "discordance_rat.png": ("Spearman rho with Fisher-z 95% CI", "Lines rising with weeks mean RNA and protein agree more as training accumulates."),
    "robust_specificity.png": ("empirical null (1,000 random sets)", "If the coloured line sits in the middle of the grey histogram, the PAH genes are no more 'opposed' than a random comparable set."),
    "robust_composition.png": ("Spearman rho; cameraPR", "Top: does OXPHOS induction track fibre-type markers across rat columns; bottom: do whole pathways move together in RNA and protein at each human time."),
    "robust_hyperemia.png": ("cameraPR with full vs muscle-intrinsic background", "If the dashed OXPHOS line sits closer to zero than the solid one, part of the drop is dilution of muscle by non-muscle material."),
}
for fn in DECK:
    shutil.copy2(f"fig/{fn}", f"deck_extracts/{fn}")
with open("deck_extracts/captions.md", "w") as f:
    f.write("# Deck figure captions\n\nFigures are copied unchanged from `fig/` (native size and aspect ratio). "
            "Full technical captions: `fig/captions.md`.\n\n")
    for fn, (test, how) in DECK.items():
        f.write(f"- **{fn}** — Test: {test}. How to read: {how}\n")
print(f"deck_extracts: {len(DECK)} figures, deck_numbers.csv ({len(nums)} rows), captions.md")
