#!/usr/bin/env python3
"""05_pah_grid.py — does the exercise response oppose the PAH skeletal-muscle signature?

Run from hackathon/:  python3 scripts/05_pah_grid.py [lookup] [grid] [pathway]   (default: all)
Input : data/join_table_v2.csv, data/raw/genesets_v1.csv
Outputs (data/ files are write-once; existing files are left untouched):
  lookup : data/pah_human_{prot,rna}_lookup_v1.csv, fig/tables/pah_human_{prot,rna}_lookup.md
  grid   : data/pah_grid_long_v1.csv, fig/pah_grid.png, fig/tables/pah_grid_summary.md
  pathway: data/pah_pathway_v1.csv, fig/pah_pathway.png, fig/tables/pah_pathway.md
Known NA: rat ATP2A1 RNA (Atp2a1 has no Ensembl gene in the rat RNA annotation), APOA1 RNA in human and
rat muscle (not expressed; the protein is plasma-derived).

PAH signature: Malenfant et al. 2015 J Mol Med, iTRAQ vastus lateralis, 4 IPAH vs 4 controls.
Value = PAH/control protein ratio. "phenotype" genes come from the same paper's enzyme
activity / blot results (CS, TFAM, PDH activity down; LDH activity up), not from its proteomics;
they are mapped to the genes encoding those enzymes and kept as a separate row group.
agreement = sign(exercise logFC) * PAH direction  ->  -1 = exercise moves the gene AGAINST PAH.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import binomtest

JOIN = "data/join_table_v2.csv"
TAB = "fig/tables"
os.makedirs(TAB, exist_ok=True)
parts = set(sys.argv[1:]) or {"lookup", "grid", "pathway", "grid_v2", "pathway_grid", "robust", "crosscheck", "hyperemia"}

# (gene, PAH/control ratio or None, PAH direction, row group). Current HGNC symbols;
# the paper's ATP5L = ATP5MG and ATP5B = ATP5F1B.
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


def write_once(df, path):
    if os.path.exists(path):
        print(f"EXISTS, not overwriting: {path}")
        return
    df.to_csv(path, index=False)
    os.chmod(path, 0o444)
    print(f"wrote {path} ({len(df)} rows)")


def caption(fn, text):
    """Insert/replace one line in fig/captions.md (same format as 04_explore.py)."""
    caps = {}
    if os.path.exists("fig/captions.md"):
        for line in open("fig/captions.md"):
            if line.startswith("- `"):
                caps[line.split("`")[1]] = line.rstrip("\n")
    caps[fn] = f"- `{fn}` — {text}"
    with open("fig/captions.md", "w") as f:
        f.write("# Figure captions\n\n" + "\n".join(caps[k] for k in sorted(caps)) + "\n")


def sign_test(agree):
    """agree: Series of +1/-1/NaN. Returns n measured, n opposed, two-sided binomial p (H0 0.5)."""
    a = agree.dropna()
    a = a[a != 0]
    n, k = len(a), int((a < 0).sum())
    return n, k, (binomtest(k, n, 0.5).pvalue if n else np.nan)


j = pd.read_csv(JOIN)

# ======================================================================= quick human lookup
if "lookup" in parts:
    cols = ["EE_vs_CON_3.5to4h", "EE_vs_CON_24h", "RE_vs_CON_3.5to4h", "RE_vs_CON_24h"]
    for layer in ["PROT", "RNA"]:
        h = j[(j.dataset == "human_acute") & (j.layer == layer) & j.contrast.isin(cols)]
        long = sig.merge(pd.DataFrame({"contrast": cols}), how="cross").merge(
            h[["gene_symbol_human", "contrast", "feature_id", "logFC", "stat", "fdr", "fdr_bh"]],
            left_on=["gene", "contrast"], right_on=["gene_symbol_human", "contrast"], how="left"
        ).drop(columns="gene_symbol_human")
        long["agreement"] = np.sign(long.logFC) * long.pah_dir
        long.insert(0, "layer", layer)
        write_once(long, f"data/pah_human_{layer.lower()}_lookup_v1.csv")

        # markdown: one row per gene, per contrast "logFC (fdr_bh) agreement"
        def cell(r):
            if pd.isna(r.logFC):
                return "NA"
            star = "*" if r.fdr_bh < 0.05 else ""
            return f"{r.logFC:+.2f} ({r.fdr_bh:.2g}){star} {'opp' if r.agreement < 0 else 'same'}"
        long["cell"] = long.apply(cell, axis=1)
        wide = long.pivot_table(index=["group", "gene", "pah_dir"], columns="contrast", values="cell",
                                aggfunc="first", sort=False)[cols].reset_index()
        wide["pah_dir"] = wide.pah_dir.map({1: "up", -1: "down"})
        summ = {"group": "**summary (19 proteomic)**", "gene": "", "pah_dir": ""}
        summ2 = {"group": "**summary (all 25)**", "gene": "", "pah_dir": ""}
        for c in cols:
            s = long[long.contrast == c]
            n, k, p = sign_test(s[s.group != "phenotype"].agreement)
            summ[c] = f"n={n}, opposed={k}, p={p:.2g}"
            n, k, p = sign_test(s.agreement)
            summ2[c] = f"n={n}, opposed={k}, p={p:.2g}"
        wide = pd.concat([wide, pd.DataFrame([summ, summ2])], ignore_index=True)
        s = wide.to_markdown(index=False)
        open(f"{TAB}/pah_human_{layer.lower()}_lookup.md", "w").write(s + "\n")
        print(f"\n### human {layer}\n{s}")


# ======================================================================= B2 grid
WK = ["1w", "2w", "4w", "8w"]
GRID_COLS = ([("human_acute", "VL", "RNA", "EE_vs_CON_24h"), ("human_acute", "VL", "PROT", "EE_vs_CON_24h"),
              ("human_acute", "VL", "RNA", "RE_vs_CON_24h"), ("human_acute", "VL", "PROT", "RE_vs_CON_24h")]
             + [("rat_train", "SKM-GN", "RNA", f"{s}_{w}") for s in "FM" for w in WK]
             + [("rat_train", "SKM-GN", "PROT", f"{s}_{w}") for s in "FM" for w in WK]
             + [("rat_train", "SKM-VL", "RNA", f"{s}_{w}") for s in "FM" for w in WK])


def col_label(ds, t, l, c):
    if ds == "human_acute":
        return f"{c[:2]} {c.split('_')[-1].replace('3.5to4h', '3.5h')}"
    return f"{c[0]} {c[2:]}"


def draw_grid(cols, blocks, png, summary_name, long_path=None, human_times="24 h", gap=None):
    """Gene x contrast PAH-agreement heatmap + per-column sign-test summary. Returns the summary frame."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    keys = pd.DataFrame(cols, columns=["dataset", "tissue", "layer", "contrast"])
    keys["col"] = range(len(keys))
    long = sig.merge(keys, how="cross").merge(
        j.rename(columns={"gene_symbol_human": "gene"}), on=["dataset", "tissue", "layer", "contrast", "gene"],
        how="left")
    long["agreement"] = np.sign(long.logFC) * long.pah_dir
    if long_path:
        write_once(long.drop(columns="col"), long_path)

    # per-column summary
    rows = []
    for (ds, t, l, c, ci), d in long.groupby(["dataset", "tissue", "layer", "contrast", "col"]):
        n19, k19, p19 = sign_test(d[d.group != "phenotype"].agreement)
        n25, k25, p25 = sign_test(d.agreement)
        rows.append(dict(col=ci, dataset=ds, tissue=t, layer=l, contrast=c, n_measured_19=n19, n_opposed_19=k19,
                         frac_opposed_19=k19 / n19 if n19 else np.nan, sign_p_19=p19,
                         n_measured_25=n25, n_opposed_25=k25, frac_opposed_25=k25 / n25 if n25 else np.nan,
                         sign_p_25=p25, n_fdrbh05=int((d.fdr_bh < 0.05).sum())))
    summ = pd.DataFrame(rows).sort_values("col").drop(columns="col")
    s = summ.to_markdown(index=False, floatfmt=".3g")
    open(f"{TAB}/{summary_name}.md", "w").write(s + "\n")
    print(s)

    # figure: rows = genes (grouped), cols = cols; colour = agreement x min(|stat|, 4)
    genes = sig.gene.tolist()
    M = np.full((len(genes), len(keys)), np.nan)
    D = np.zeros_like(M, dtype=bool)
    for r in long.itertuples():
        i = genes.index(r.gene)
        if pd.notna(r.stat):
            M[i, r.col] = r.agreement * min(abs(r.stat), 4)
            D[i, r.col] = r.fdr_bh < 0.05
    # optional blank spacer column before index `gap` (separates the human block from rat blocks)
    X = (lambda ci: ci + (1 if gap is not None and ci >= gap else 0))
    ncol = len(keys) + (1 if gap is not None else 0)
    if gap is not None:
        M = np.insert(M, gap, np.nan, axis=1); D = np.insert(D, gap, False, axis=1)
    cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a61e1e"])
    cmap.set_bad("white")
    fig, ax = plt.subplots(figsize=(13 + 0.33 * (len(cols) - 28), 10.5))
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-4, vmax=4, aspect="equal")
    yy, xx = np.where(np.isnan(M))
    for y, x in zip(yy, xx):
        if gap is not None and x == gap:
            continue
        ax.add_patch(plt.Rectangle((x - .5, y - .5), 1, 1, fill=False, hatch="////", edgecolor="#c3c2b7", lw=0))
        ax.text(x, y, "NA", ha="center", va="center", fontsize=5.5, color="#52514e")
    yy, xx = np.where(D)
    ax.scatter(xx, yy, s=14, color="#0b0b0b", zorder=3, label="BH FDR < 0.05 (fdr_bh)")
    # grid lines between cells, heavier between row groups and column blocks
    ax.set_xticks(np.arange(-.5, ncol), minor=True); ax.set_yticks(np.arange(-.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", lw=1.2); ax.tick_params(which="minor", length=0)
    grp = sig.group.tolist()
    for i in range(1, len(genes)):
        if grp[i] != grp[i - 1]:
            ax.axhline(i - .5, color="#0b0b0b", lw=1.2)
    for a, b, lab in blocks:
        if a:
            ax.axvline(X(a) - .5, color="#0b0b0b", lw=1.2)
            if gap is not None and a == gap:
                ax.axvline(X(a) - 1.5, color="#0b0b0b", lw=1.2)
        ax.text((X(a) + X(b - 1)) / 2, -1.6, lab, ha="center", va="bottom", fontsize=9)
    xl = [col_label(*k) + (f"\n{k[2]}" if k[0] == "human_acute" else "") for k in cols]
    ax.set_xticks([X(i) for i in range(len(keys))], xl, fontsize=7.5, rotation=90)
    ax.set_yticks(range(len(genes)),
                  [f"{g} ({'↓' if d < 0 else '↑'} in PAH)" for g, d in zip(genes, sig.pah_dir)], fontsize=8)
    # row-group labels on the right
    for gname in dict.fromkeys(grp):
        idx = [i for i, g in enumerate(grp) if g == gname]
        ax.text(ncol - .3, np.mean(idx), gname, va="center", fontsize=8, color="#52514e")
    # per-column summary under the grid: opposed / measured (19 proteomic genes)
    for ci, r in enumerate(summ.itertuples()):
        star = "*" if r.sign_p_19 < 0.05 else ""
        ax.text(X(ci), len(genes) + 2.3, f"{r.n_opposed_19}/{r.n_measured_19}{star}", ha="center", va="top",
                fontsize=6.5, rotation=90)
    ax.text(-1, len(genes) + 2.3, "opposed / measured\n(19 proteomic genes;\n* sign test p<0.05)", ha="right",
            va="top", fontsize=7)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.13)
    cb.set_label("agreement × |stat| (capped at 4)\n← exercise OPPOSES PAH    same as PAH →", fontsize=8)
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.08), fontsize=8, frameon=False)
    ax.set_title("Does the exercise response oppose the PAH skeletal-muscle signature (Malenfant 2015)?",
                 loc="left", fontsize=11, pad=62)
    fig.savefig(f"fig/{png}", dpi=150, bbox_inches="tight"); plt.close(fig)
    caption(png, "Rows: Malenfant 2015 PAH-muscle proteins (19) + 6 enzyme-phenotype genes, "
                f"grouped; columns: exercise contrast (human VL delta-delta vs control at {human_times}; rat trained vs "
                "sedentary by sex x week). Cell colour = sign(exercise logFC) x PAH direction x min(|stat|,4): blue = "
                "exercise moves the gene against PAH, red = same way as PAH; dot = BH FDR<0.05 (fdr_bh); hatched NA = "
                "not measured. Bottom numbers: genes opposed / measured among the 19, * = two-sided binomial sign "
                "test p<0.05 (H0 0.5), which can mean significantly opposed (high count) OR significantly same-as-PAH "
                "(low count, e.g. human PROT 2/19); genes are not independent, so p is optimistic.")


    return summ


if "grid" in parts:
    draw_grid(GRID_COLS, [(0, 4, "Human VL\nacute bout"), (4, 12, "Rat SKM-GN RNA\n(training weeks)"),
                          (12, 20, "Rat SKM-GN PROT\n(training weeks)"), (20, 28, "Rat SKM-VL RNA\n(training weeks)")],
              "pah_grid.png", "pah_grid_summary", long_path="data/pah_grid_long_v1.csv")

# v2: human 3.5 h columns added (v1 figure/table kept). No data/ output (data/ is read-only for this step).
GRID_COLS_V2 = ([("human_acute", "VL", l, f"{g}_vs_CON_{t}") for g in ("EE", "RE") for t in ("3.5to4h", "24h")
                 for l in ("RNA", "PROT")] + [c for c in GRID_COLS if c[0] == "rat_train"])
if "grid_v2" in parts:
    draw_grid(GRID_COLS_V2, [(0, 8, "Human VL\nacute bout"), (8, 16, "Rat SKM-GN RNA\n(training weeks)"),
                             (16, 24, "Rat SKM-GN PROT\n(training weeks)"), (24, 32, "Rat SKM-VL RNA\n(training weeks)")],
              "pah_grid_v2.png", "pah_grid_v2_summary", human_times="3.5-4 h and 24 h", gap=8)


# ======================================================================= B3 pathway level
def camera_pr(stat, index, inter_gene_cor=0.01):
    """Python port of limma::cameraPR.default (use.ranks=FALSE, fixed inter.gene.cor).
    stat: 1-D array over all genes; index: boolean mask of set members. Returns (t, df, two-sided p)."""
    from scipy.stats import t as tdist
    stat = np.asarray(stat, float)
    G, m = len(stat), int(index.sum())
    mean_stat, var_stat = stat.mean(), stat.var(ddof=1)
    vif = 1 + (m - 1) * inter_gene_cor
    delta = G / (G - m) * (stat[index].mean() - mean_stat)
    var_pooled = ((G - 1) * var_stat - delta ** 2 * m * (G - m) / G) / (G - 2)
    tstat = delta / np.sqrt(var_pooled * (vif / m + 1 / (G - m)))
    df = G - 2
    return tstat, df, 2 * tdist.sf(abs(tstat), df)


PATH_COLS = ([("human_acute", "VL", l, f"{g}_vs_CON_{t}") for g in ("EE", "RE") for t in ("15to45min", "3.5to4h", "24h")
              for l in ("RNA", "PROT")]
             + [c for c in GRID_COLS if c[0] == "rat_train"])
PAH_SETS = {"KEGG_OXIDATIVE_PHOSPHORYLATION": -1, "KEGG_GLYCOLYSIS_GLUCONEOGENESIS": 1,
            "HALLMARK_OXIDATIVE_PHOSPHORYLATION": -1, "HALLMARK_GLYCOLYSIS": 1}

if "pathway" in parts:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    gs = pd.read_csv("data/raw/genesets_v1.csv")
    rows = []
    for ds, t, l, c in PATH_COLS:
        d = j[(j.dataset == ds) & (j.tissue == t) & (j.layer == l) & (j.contrast == c)].dropna(subset=["stat"])
        for setname, pah_dir in PAH_SETS.items():
            idx = d.gene_symbol_human.isin(set(gs.gene_symbol[gs.gs_name == setname])).to_numpy()
            tt, df, p = camera_pr(d.stat.to_numpy(), idx)
            rows.append(dict(dataset=ds, tissue=t, layer=l, contrast=c, gene_set=setname, pah_direction=pah_dir,
                             n_genes_in_set=int(idx.sum()), n_genes_total=len(d), camera_t=tt, p=p,
                             direction="Up" if tt > 0 else "Down",
                             agreement=int(np.sign(tt) * pah_dir),
                             frac_set_up=float((d.stat[idx] > 0).mean())))
    pw = pd.DataFrame(rows)
    pw["fdr_bh_within_set"] = pw.groupby("gene_set").p.transform(
        lambda x: __import__("scipy.stats", fromlist=["x"]).false_discovery_control(x, method="bh"))
    write_once(pw, "data/pah_pathway_v1.csv")
    view = pw[pw.gene_set.str.startswith("KEGG")].copy()
    view["cell"] = view.apply(lambda r: f"{r.direction} t={r.camera_t:+.1f} p={r.p:.2g} "
                              f"{'OPPOSES' if r.agreement < 0 else 'same as'} PAH", axis=1)
    wide = view.pivot_table(index=["dataset", "tissue", "layer", "contrast"], columns="gene_set", values="cell",
                            aggfunc="first", sort=False).reset_index()
    s = wide.to_markdown(index=False)
    open(f"{TAB}/pah_pathway.md", "w").write(s + "\n")
    print(s)

    # figure: camera t per column for the two KEGG sets (and Hallmark as open markers)
    labels = [f"{'hum' if ds == 'human_acute' else t} {l} {c.replace('_vs_CON_', ' ')}" for ds, t, l, c in PATH_COLS]
    x = np.arange(len(PATH_COLS))
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    for ax, (kegg, hall, title) in zip(axes, [
            ("KEGG_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "OXPHOS (down in PAH)"),
            ("KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "HALLMARK_GLYCOLYSIS", "Glycolysis (up in PAH)")]):
        LC = {"RNA": "#2a78d6", "PROT": "#eb6834"}
        for setname, marker, filled, lab in [(kegg, "o", True, "KEGG"), (hall, "D", False, "Hallmark")]:
            d = pw[pw.gene_set == setname].set_index(["dataset", "tissue", "layer", "contrast"]).loc[PATH_COLS]
            for layer, col in LC.items():
                sel = np.array([k[2] == layer for k in PATH_COLS])
                ax.scatter(x[sel], d.camera_t[sel], marker=marker, s=34, facecolor=col if filled else "white",
                           edgecolor=col, lw=1.2, zorder=3, label=f"{layer}, {lab}")
        ax.axhline(0, color="#52514e", lw=0.8)
        for thr in (1.96, -1.96):
            ax.axhline(thr, color="#c3c2b7", lw=0.8, ls="--")
        pah = PAH_SETS[kegg]
        ax.text(1.005, 0.97 if pah < 0 else 0.03, "exercise OPPOSES PAH", transform=ax.transAxes, va="top" if pah < 0 else "bottom",
                fontsize=8, color="#104281")
        ax.text(1.005, 0.03 if pah < 0 else 0.97, "same as PAH", transform=ax.transAxes, va="bottom" if pah < 0 else "top",
                fontsize=8, color="#a61e1e")
        ax.set_ylabel("cameraPR t (set vs rest)"); ax.set_title(title, loc="left", fontsize=10)
        for b in (12, 20, 28):
            ax.axvline(b - .5, color="#0b0b0b", lw=0.8)
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.0), fontsize=7.5, ncol=4)
        ax.grid(axis="x", visible=False)
    axes[1].set_xticks(x, labels, rotation=90, fontsize=7)
    fig.tight_layout(); fig.savefig("fig/pah_pathway.png"); plt.close(fig)
    caption("pah_pathway.png", "Pathway-level test per contrast: limma cameraPR (Python port, inter-gene correlation "
            "0.01, validated to match limma 3.62.1 cameraPR p-values exactly) of the gene-level stat column (join_table_v2) for OXPHOS and "
            "glycolysis sets; point = one dataset x layer x contrast (blue = RNA, orange = PROT); filled circle = KEGG hsa00190/hsa00010, open diamond = "
            "MSigDB Hallmark; dashed = |t| 1.96 (~p 0.05). Positive t = set shifted up vs all other measured genes. "
            "OXPHOS up / glycolysis down = exercise opposes the PAH muscle phenotype.")


# ======================================================================= extra gene sets (cache in fig/tables)
# data/ is read-only for this step, so the additional msigdbr sets are cached next to the other derived tables.
EXTRA_GS = f"{TAB}/genesets_extra.csv"
RSCRIPT = os.path.expanduser("~/miniconda3/envs/motrpac/bin/Rscript")


def extra_genesets():
    if not os.path.exists(EXTRA_GS):
        import subprocess
        r = f'''suppressMessages({{library(msigdbr); library(data.table)}})
        a <- as.data.table(msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:KEGG_LEGACY"))[gs_name == "KEGG_CITRATE_CYCLE_TCA_CYCLE"]
        b <- as.data.table(msigdbr(species="Homo sapiens", collection="H"))
        d <- as.data.table(msigdbr(species="Homo sapiens", collection="C5", subcollection="GO:BP"))[gs_name == "GOBP_MITOCHONDRIAL_TRANSLATION"]
        e <- as.data.table(msigdbr(species="Homo sapiens", collection="C2", subcollection="CP:REACTOME"))[gs_name == "REACTOME_MITOCHONDRIAL_BIOGENESIS"]
        x <- unique(rbind(a, b, d, e, fill=TRUE)[, .(gs_name, gene_symbol)])
        x[, msigdbr_version := as.character(packageVersion("msigdbr"))]
        fwrite(x, "{EXTRA_GS}")'''
        subprocess.run([RSCRIPT if os.path.exists(RSCRIPT) else "Rscript", "-e", r], check=True)
    g = pd.concat([pd.read_csv("data/raw/genesets_v1.csv"), pd.read_csv(EXTRA_GS)]).drop_duplicates(["gs_name", "gene_symbol"])
    sets = {k: set(v) for k, v in g.groupby("gs_name").gene_symbol}
    sets["MYH_SET"] = {"MYH1", "MYH2", "MYH4", "MYH7"}
    return sets


def col_frame(ds, t, l, c):
    return j[(j.dataset == ds) & (j.tissue == t) & (j.layer == l) & (j.contrast == c)].dropna(subset=["stat"])


# ======================================================================= pathway heatmap
PW_ROWS = [("KEGG_OXIDATIVE_PHOSPHORYLATION", "OXPHOS (KEGG)", "↓"), ("KEGG_CITRATE_CYCLE_TCA_CYCLE", "TCA cycle (KEGG)", "↓"),
           ("KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "Glycolysis (KEGG)", "↑"),
           ("MYH_SET", "Contractile: MYH1/2/4/7", "↑"),
           ("GOBP_MITOCHONDRIAL_TRANSLATION", "Mito translation (GO:BP)", "↓*"),
           ("REACTOME_MITOCHONDRIAL_BIOGENESIS", "Mito biogenesis (Reactome)", "↓*")]
if "pathway_grid" in parts:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from scipy.stats import false_discovery_control
    sets = extra_genesets()
    rows = []
    for ci, key in enumerate(GRID_COLS_V2):
        d = col_frame(*key)
        for ri, (sname, lab, pdir) in enumerate(PW_ROWS):
            idx = d.gene_symbol_human.isin(sets[sname]).to_numpy()
            tt, df, p = camera_pr(d.stat.to_numpy(), idx) if idx.sum() >= 2 else (np.nan, np.nan, np.nan)
            rows.append(dict(zip(["dataset", "tissue", "layer", "contrast"], key), gene_set=sname, row=ri, col=ci,
                             n_in_set=int(idx.sum()), camera_t=tt, p=p))
    pg = pd.DataFrame(rows)
    pg["fdr_bh_row"] = pg.groupby("gene_set").p.transform(lambda x: false_discovery_control(x.fillna(1), method="bh"))
    pg.to_csv(f"{TAB}/pah_pathway_grid.csv", index=False)
    M = pg.pivot(index="row", columns="col", values="camera_t").to_numpy()
    D = pg.pivot(index="row", columns="col", values="fdr_bh_row").to_numpy() < 0.05
    cmap = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#a61e1e"])
    fig, ax = plt.subplots(figsize=(14, 4.6))
    im = ax.imshow(np.clip(M, -8, 8), cmap=cmap, vmin=-8, vmax=8, aspect="equal")
    yy, xx = np.where(D)
    ax.scatter(xx, yy, s=14, color="#0b0b0b", zorder=3, label="BH FDR < 0.05 (within row)")
    ax.set_xticks(np.arange(-.5, M.shape[1]), minor=True); ax.set_yticks(np.arange(-.5, M.shape[0]), minor=True)
    ax.grid(which="minor", color="white", lw=1.2); ax.tick_params(which="minor", length=0); ax.grid(which="major", visible=False)
    for a, b, lab in [(0, 8, "Human VL acute bout"), (8, 16, "Rat SKM-GN RNA"), (16, 24, "Rat SKM-GN PROT"), (24, 32, "Rat SKM-VL RNA")]:
        if a:
            ax.axvline(a - .5, color="#0b0b0b", lw=1.2)
        ax.text((a + b - 1) / 2, -0.9, lab, ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(GRID_COLS_V2)),
                  [col_label(*k) + (f"\n{k[2]}" if k[0] == "human_acute" else "") for k in GRID_COLS_V2], rotation=90, fontsize=7.5)
    ax.set_yticks(range(len(PW_ROWS)), [f"{lab} ({pdir} in PAH)" for _, lab, pdir in PW_ROWS], fontsize=8.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.02)
    cb.set_label("cameraPR t (set vs all other genes),\ncapped at ±8; red = up with exercise", fontsize=8)
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 1.12), fontsize=8, frameon=False)
    ax.set_title("Pathway-level exercise response (cameraPR on the gene-level stat)", loc="left", fontsize=11, pad=30)
    fig.savefig("fig/pah_pathway_grid.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    caption("pah_pathway_grid.png", "Pathway heatmap: cell = cameraPR t (Python port of limma, inter-gene cor 0.01) of "
            "the gene-level stat for a gene set vs all other measured genes in one contrast (join_table_v2), capped at "
            "±8; red = set up with exercise, blue = down; dot = BH FDR<0.05 within the row (32 columns). PAH direction "
            "in the row label (↓* = implied by the paper's TFAM/CS decrease, not measured as a set). Exercise opposes "
            "PAH where the colour is opposite to the arrow.")
    print(pg.pivot(index="gene_set", columns="col", values="camera_t").round(1).to_string())


# ======================================================================= robustness R1-R3
# Outputs: fig/robust_specificity.png, fig/robust_composition.png, fig/tables/robust_*.md|csv,
#          fig/tables/robust_numbers.csv (every number the narrative cites from this block).
MITO9 = ["NDUFA9", "UQCRC2", "UQCRC1", "ATP5MG", "ATP5F1B", "IDH2", "OGDH", "SLC25A4", "ECH1"]
PAH19 = sig[sig.group != "phenotype"].gene.tolist()
NUMS = []
_PART = ["robust"]  # which part produced a registered number (provenance in deck_numbers.csv)


def num(key, value, figure="", note=""):
    NUMS.append(dict(key=key, value=value, source_script=f"scripts/05_pah_grid.py {_PART[0]}", figure=figure, note=note))
    return value


if "robust" in parts:
    _PART[0] = "robust"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import spearmanr
    sets = extra_genesets()
    rng = np.random.default_rng(20260926)
    NBOOT = 1000
    C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]

    # ---------------- R1(a): PAH-9 mito genes vs random 9-gene OXPHOS+TCA sets
    pool_sets = sets["KEGG_OXIDATIVE_PHOSPHORYLATION"] | sets["KEGG_CITRATE_CYCLE_TCA_CYCLE"]

    def both_layers(ds, rna_t, prot_t):
        r = set(j[(j.dataset == ds) & (j.tissue == rna_t) & (j.layer == "RNA")].gene_symbol_human)
        p = set(j[(j.dataset == ds) & (j.tissue == prot_t) & (j.layer == "PROT")].gene_symbol_human)
        return r & p
    pools = {"human_acute": sorted(pool_sets & both_layers("human_acute", "VL", "VL")),
             "rat_train": sorted(pool_sets & both_layers("rat_train", "SKM-GN", "SKM-GN"))}
    draws = {k: [rng.choice(v, 9, replace=False) for _ in range(NBOOT)] for k, v in pools.items()}
    r1 = []
    for key in GRID_COLS_V2:
        d = col_frame(*key)
        st = d.stat.to_numpy(); g = d.gene_symbol_human.to_numpy()
        lfc = dict(zip(g, d.logFC))

        def score(genes):
            idx = np.isin(g, genes)
            vals = [lfc[x] for x in genes if x in lfc]
            return (np.mean(np.array(vals) > 0) if vals else np.nan,  # all 9 are "down in PAH": opposed = up
                    camera_pr(st, idx)[0] if idx.sum() >= 2 else np.nan)
        obs_f, obs_t = score(MITO9)
        null = np.array([score(s) for s in draws[key[0]]])
        r1.append(dict(zip(["dataset", "tissue", "layer", "contrast"], key), pool_size=len(pools[key[0]]),
                       obs_frac_opposed=obs_f, null_frac_median=np.nanmedian(null[:, 0]),
                       pct_frac=100 * np.nanmean(null[:, 0] <= obs_f),
                       obs_camera_t=obs_t, null_t_median=np.nanmedian(null[:, 1]),
                       null_t_q05=np.nanpercentile(null[:, 1], 5), null_t_q95=np.nanpercentile(null[:, 1], 95),
                       pct_t=100 * np.nanmean(null[:, 1] <= obs_t), _null=null))
    r1 = pd.DataFrame(r1)
    r1.drop(columns="_null").to_csv(f"{TAB}/robust_r1a.csv", index=False)
    open(f"{TAB}/robust_r1a.md", "w").write(r1.drop(columns="_null").to_markdown(index=False, floatfmt=".3g") + "\n")
    for ds, t, l, c, tag in [("human_acute", "VL", "RNA", "EE_vs_CON_24h", "hum_EE24_RNA"),
                             ("rat_train", "SKM-GN", "PROT", "F_8w", "rat_GN_PROT_F8"),
                             ("rat_train", "SKM-GN", "PROT", "M_8w", "rat_GN_PROT_M8")]:
        r = r1[(r1.dataset == ds) & (r1.tissue == t) & (r1.layer == l) & (r1.contrast == c)].iloc[0]
        for f in ["obs_frac_opposed", "null_frac_median", "pct_frac", "obs_camera_t", "null_t_median", "pct_t"]:
            num(f"R1a_{tag}_{f}", round(float(r[f]), 3), "robust_specificity.png")
    num("R1a_pool_human", len(pools["human_acute"]), "robust_specificity.png", "KEGG OXPHOS+TCA genes measured in both layers")
    num("R1a_pool_rat", len(pools["rat_train"]), "robust_specificity.png")
    num("R1a_mito9_in_pool", int(len(set(MITO9) & pool_sets)), "", "PAH-9 genes that are in KEGG OXPHOS/TCA")

    # ---------------- R1(b): 19 random genes matched on abundance decile, human EE 24h RNA
    raw = pd.read_csv("data/raw/human_muscle_rna_DA_v1.csv", usecols=["contrast_short", "feature_id", "AveExpr"])
    raw = raw[raw.contrast_short == "Endur.post_24_hr - Control.post_24_hr (delta-delta)"]
    d = col_frame("human_acute", "VL", "RNA", "EE_vs_CON_24h").merge(raw[["feature_id", "AveExpr"]], on="feature_id")
    d["decile"] = pd.qcut(d.AveExpr, 10, labels=False)
    pahd = sig.set_index("gene").pah_dir
    meas = d[d.gene_symbol_human.isin(PAH19)]
    obs = int((np.sign(meas.logFC) * meas.gene_symbol_human.map(pahd) < 0).sum())
    bydec = {k: v for k, v in d[~d.gene_symbol_human.isin(PAH19)].groupby("decile")}
    nullb = []
    for _ in range(NBOOT):
        k = 0
        for gname, dec in zip(meas.gene_symbol_human, meas.decile):
            pick = bydec[dec].iloc[rng.integers(len(bydec[dec]))]
            k += int(np.sign(pick.logFC) * pahd[gname] < 0)
        nullb.append(k)
    nullb = np.array(nullb)
    p_emp = (1 + (nullb >= obs).sum()) / (NBOOT + 1)
    num("R1b_obs_opposed", obs, "robust_specificity.png"); num("R1b_n_measured", len(meas), "robust_specificity.png")
    num("R1b_null_median", float(np.median(nullb)), "robust_specificity.png")
    num("R1b_null_q95", float(np.percentile(nullb, 95)), "robust_specificity.png")
    num("R1b_p_empirical", round(float(p_emp), 4), "robust_specificity.png", "P(null >= observed), 1000 decile-matched sets")
    num("R1b_decile_of_pah_genes_median", float(meas.decile.median()), "", "PAH genes are highly expressed")

    fig, axes = plt.subplots(1, 4, figsize=(13, 3.4))
    for ax, (ds, t, l, c, ttl) in zip(axes[:3], [("human_acute", "VL", "RNA", "EE_vs_CON_24h", "Human VL RNA, EE 24 h"),
                                                  ("rat_train", "SKM-GN", "PROT", "F_8w", "Rat SKM-GN PROT, F 8 wk"),
                                                  ("rat_train", "SKM-GN", "PROT", "M_8w", "Rat SKM-GN PROT, M 8 wk")]):
        r = r1[(r1.dataset == ds) & (r1.tissue == t) & (r1.layer == l) & (r1.contrast == c)].iloc[0]
        ax.hist(r._null[:, 1], bins=30, color="#b9b8b3", label="1,000 random 9-gene\nOXPHOS/TCA sets")
        ax.axvline(r.obs_camera_t, color=C[0], lw=2, label=f"PAH mito-9 (pct {r.pct_t:.0f})")
        ax.set_title(ttl, loc="left", fontsize=9); ax.set_xlabel("cameraPR t (set vs rest)"); ax.grid(axis="x", visible=False)
        ax.legend(fontsize=7, loc="upper left")
    axes[0].set_ylabel("n random sets")
    ax = axes[3]
    ax.hist(nullb, bins=np.arange(nullb.min() - .5, max(nullb.max(), obs) + 1.5), color="#b9b8b3",
            label="1,000 random 19-gene sets\n(abundance-decile matched)")
    ax.axvline(obs, color=C[1], lw=2, label=f"PAH-19: {obs}/{len(meas)} (p = {p_emp:.3f})")
    ax.set_title("Human VL RNA, EE 24 h", loc="left", fontsize=9); ax.set_xlabel("n genes opposed to PAH direction")
    ax.legend(fontsize=7, loc="upper left"); ax.grid(axis="x", visible=False)
    fig.tight_layout(); fig.savefig("fig/robust_specificity.png"); plt.close(fig)
    caption("robust_specificity.png", "Specificity nulls. Panels 1-3: cameraPR t of the 9 PAH mitochondrial genes (blue "
            "line) vs 1,000 random 9-gene sets drawn from KEGG OXPHOS+TCA genes measured in both layers (grey); the "
            "percentile says whether the PAH genes are more 'opposed' than any mitochondrial set. Panel 4: number of the "
            "18 measured PAH-19 genes moving against their PAH direction (orange) vs 1,000 abundance-decile-matched random "
            "19-gene sets given the same direction vector; p = one-sided empirical P(null >= observed).")

    # ---------------- R2(a): fibre-type markers vs OXPHOS across rat columns
    ox = sets["KEGG_OXIDATIVE_PHOSPHORYLATION"]
    rat_cols = [k for k in GRID_COLS_V2 if k[0] == "rat_train"]
    rows = []
    for key in rat_cols:
        d = col_frame(*key)
        tt = camera_pr(d.stat.to_numpy(), d.gene_symbol_human.isin(ox).to_numpy())[0]
        st = d.set_index("gene_symbol_human").stat
        rows.append(dict(zip(["dataset", "tissue", "layer", "contrast"], key), oxphos_t=tt,
                         **{m: st.get(m, np.nan) for m in ["MYH7", "MYH2", "MYH1", "MYH4"]}))
    r2a = pd.DataFrame(rows)
    r2a.to_csv(f"{TAB}/robust_r2a.csv", index=False)
    r2a_s = []
    for m in ["MYH7", "MYH2", "MYH1", "MYH4"]:
        ok = r2a[m].notna()
        rho, p = spearmanr(r2a.oxphos_t[ok], r2a[m][ok])
        r2a_s.append(dict(marker=m, fibre={"MYH7": "type I", "MYH2": "IIa", "MYH1": "IIx", "MYH4": "IIb"}[m],
                          n_columns=int(ok.sum()), spearman=rho, p=p))
        num(f"R2a_rho_{m}", round(float(rho), 3), "robust_composition.png"); num(f"R2a_p_{m}", round(float(p), 4), "robust_composition.png")
        num(f"R2a_n_{m}", int(ok.sum()), "robust_composition.png")
        for (t_, l_), g_ in r2a.groupby(["tissue", "layer"]):  # within one tissue x layer (n = 8 columns)
            rw = spearmanr(g_.oxphos_t, g_[m])[0]
            r2a_s[-1][f"rho_{t_}_{l_}"] = rw
            num(f"R2a_rho_{m}_{t_}_{l_}", round(float(rw), 3), "robust_composition.png", "within tissue x layer, n=8")
    r2a_s = pd.DataFrame(r2a_s)
    open(f"{TAB}/robust_r2a.md", "w").write(r2a_s.to_markdown(index=False, floatfmt=".3g") + "\n")

    # ---------------- R2(b): human Hallmark 50, RNA vs PROT co-movement by time
    hall = {k: v for k, v in sets.items() if k.startswith("HALLMARK_")}
    num("R2b_n_hallmark_sets", len(hall))
    rows = []
    for g_ in ("EE", "RE"):
        for tp in ("15to45min", "3.5to4h", "24h"):
            c = f"{g_}_vs_CON_{tp}"
            dr, dp = col_frame("human_acute", "VL", "RNA", c), col_frame("human_acute", "VL", "PROT", c)
            for hn, genes in hall.items():
                tr = camera_pr(dr.stat.to_numpy(), dr.gene_symbol_human.isin(genes).to_numpy())[0]
                tp_ = camera_pr(dp.stat.to_numpy(), dp.gene_symbol_human.isin(genes).to_numpy())[0]
                rows.append(dict(contrast=c, gene_set=hn, t_rna=tr, t_prot=tp_))
    r2b = pd.DataFrame(rows)
    r2b.to_csv(f"{TAB}/robust_r2b.csv", index=False)
    summ = []
    for c, d in r2b.groupby("contrast", sort=False):
        both = (d.t_rna.abs() > 2) & (d.t_prot.abs() > 2)
        same = both & (np.sign(d.t_rna) == np.sign(d.t_prot))
        rho = spearmanr(d.t_rna, d.t_prot)[0]
        summ.append(dict(contrast=c, n_sets=len(d), n_rna_abs_t_gt2=int((d.t_rna.abs() > 2).sum()),
                         n_prot_abs_t_gt2=int((d.t_prot.abs() > 2).sum()), n_both_same_sign=int(same.sum()),
                         n_both_opposite_sign=int((both & ~same).sum()), spearman_t_rna_vs_prot=rho,
                         n_same_sign_any=int((np.sign(d.t_rna) == np.sign(d.t_prot)).sum())))
        tag = c.replace("_vs_CON_", "_")
        num(f"R2b_{tag}_both_same", int(same.sum()), "robust_composition.png")
        num(f"R2b_{tag}_rho", round(float(rho), 3), "robust_composition.png")
        num(f"R2b_{tag}_prot_gt2", int((d.t_prot.abs() > 2).sum()), "robust_composition.png")
        num(f"R2b_{tag}_rna_gt2", int((d.t_rna.abs() > 2).sum()), "robust_composition.png")
    r2b_s = pd.DataFrame(summ)
    open(f"{TAB}/robust_r2b.md", "w").write(r2b_s.to_markdown(index=False, floatfmt=".3g") + "\n")
    print(r2a_s.to_string()); print(r2b_s.to_string())

    fig, axes = plt.subplots(2, 4, figsize=(13, 6.4))
    grp_c = {("SKM-GN", "RNA"): C[0], ("SKM-GN", "PROT"): C[1], ("SKM-VL", "RNA"): C[2]}
    for ax, row in zip(axes[0], r2a_s.itertuples()):
        for (t, l), col in grp_c.items():
            s_ = r2a[(r2a.tissue == t) & (r2a.layer == l)]
            ax.scatter(s_[row.marker], s_.oxphos_t, s=30, color=col, edgecolor="white", lw=0.5, label=f"rat {t} {l}")
        ax.set_title(f"{row.marker} ({row.fibre}): ρ = {row.spearman:.2f}, n = {row.n_columns}", loc="left", fontsize=9)
        ax.set_xlabel(f"{row.marker} gene stat (z RNA / t PROT)"); ax.axvline(0, color="#52514e", lw=0.6)
    axes[0, 0].set_ylabel("OXPHOS cameraPR t"); axes[0, 0].legend(fontsize=7, loc="upper left")
    for ax, c in zip(axes[1], ["EE_vs_CON_15to45min", "EE_vs_CON_24h", "RE_vs_CON_15to45min", "RE_vs_CON_24h"]):
        d = r2b[r2b.contrast == c]; s_ = r2b_s[r2b_s.contrast == c].iloc[0]
        ax.scatter(d.t_rna, d.t_prot, s=18, color=C[3], edgecolor="#52514e", lw=0.3)
        ax.axhline(0, color="#52514e", lw=0.6); ax.axvline(0, color="#52514e", lw=0.6)
        for v in (2, -2):
            ax.axhline(v, color="#c3c2b7", lw=0.6, ls="--"); ax.axvline(v, color="#c3c2b7", lw=0.6, ls="--")
        ax.set_title(f"Human {c.replace('_vs_CON_', ' ')}: ρ = {s_.spearman_t_rna_vs_prot:.2f}\n"
                     f"{s_.n_both_same_sign}/50 sets |t|>2 same sign in both", loc="left", fontsize=9)
        ax.set_xlabel("RNA cameraPR t (Hallmark set)")
    axes[1, 0].set_ylabel("PROT cameraPR t (Hallmark set)")
    fig.tight_layout(); fig.savefig("fig/robust_composition.png"); plt.close(fig)
    caption("robust_composition.png", "Composition checks. Top (R2a): across the 24 rat columns (point = one tissue x "
            "layer x sex x week), KEGG OXPHOS cameraPR t vs the gene stat of each myosin heavy-chain fibre-type marker; "
            "Spearman rho in the title (does 'OXPHOS up' track a fibre-type shift?). Bottom (R2b): human VL, cameraPR t "
            "of all 50 MSigDB Hallmark sets in RNA (x) vs PROT (y), at 15-45 min and 24 h after EE/RE (point = one "
            "Hallmark set); dashed = |t| 2. Many sets co-moving at 15-45 min but not at 24 h indicates a global/technical "
            "shift rather than gene-specific regulation.")

    # ---------------- R3: correlation among the 9 mito genes across the 28 grid columns
    mat = pd.DataFrame({f"{k[1]}|{k[2]}|{k[3]}": col_frame(*k).set_index("gene_symbol_human").stat.reindex(MITO9)
                        for k in GRID_COLS})
    corr = mat.T.corr(method="spearman", min_periods=10)
    lam = np.clip(np.linalg.eigvalsh(corr.to_numpy()), 0, None)
    M = len(lam)
    meff_nyholt = 1 + (M - 1) * (1 - np.var(lam, ddof=1) / M)
    meff_liji = float(np.sum((lam >= 1) + (lam - np.floor(lam))))
    off = corr.to_numpy()[np.triu_indices(M, 1)]
    num("R3_mean_offdiag_r", round(float(off.mean()), 3), "", "mean pairwise Spearman r, 9 mito genes x 28 columns")
    num("R3_min_offdiag_r", round(float(off.min()), 3)); num("R3_max_offdiag_r", round(float(off.max()), 3))
    num("R3_meff_nyholt", round(float(meff_nyholt), 2), "", "effective number of independent genes out of 9")
    num("R3_meff_liji", round(meff_liji, 2), "", "Li-Ji effective number out of 9")
    num("R3_top_eigen_frac", round(float(lam.max() / lam.sum()), 3), "", "variance share of first eigenvector")
    open(f"{TAB}/robust_r3.md", "w").write(corr.round(2).reset_index().rename(columns={"index": "gene"}).to_markdown(index=False) + "\n")
    print(corr.round(2).to_string()); print("Meff Nyholt", meff_nyholt, "LiJi", meff_liji)

    pd.DataFrame(NUMS).to_csv(f"{TAB}/robust_numbers.csv", index=False)
    print(pd.DataFrame(NUMS)[["key", "value"]].to_string())


# ======================================================================= package gene sets (MOLECULAR_SIGNATURES)
PKG_GS = f"{TAB}/genesets_pkg.csv"


def pkg_genesets():
    """Sets shipped inside MotrpacHumanPreSuspensionAnalysis (offline), incl. MitoCarta 3.0 MitoPathways."""
    if not os.path.exists(PKG_GS):
        import subprocess
        r = f'''m <- MotrpacHumanPreSuspensionAnalysis::MOLECULAR_SIGNATURES
        s <- list(GOBP_OXIDATIVE_PHOSPHORYLATION = m$GOBP[["GOBP_OXIDATIVE_PHOSPHORYLATION"]],
                  GOBP_COMPLEMENT_ACTIVATION = m$GOBP[["GOBP_COMPLEMENT_ACTIVATION"]],
                  GOCC_CONTRACTILE_FIBER = m$GOCC[["GOCC_CONTRACTILE_FIBER"]],
                  GOCC_SARCOPLASMIC_RETICULUM = m$GOCC[["GOCC_SARCOPLASMIC_RETICULUM"]],
                  MITOCARTA_ALL = unique(unlist(m$MITOCARTA)))
        x <- data.frame(gs_name = rep(names(s), lengths(s)), gene_symbol = unlist(s, use.names = FALSE))
        x$source <- paste0("MotrpacHumanPreSuspensionAnalysis ", packageVersion("MotrpacHumanPreSuspensionAnalysis"))
        write.csv(x, "{PKG_GS}", row.names = FALSE)'''
        subprocess.run([RSCRIPT if os.path.exists(RSCRIPT) else "Rscript", "-e", r], check=True)
    g = pd.read_csv(PKG_GS)
    return {k: set(v) for k, v in g.groupby("gs_name").gene_symbol}


H_TIMES = [("15to45min", "15-45 min"), ("3.5to4h", "3.5-4 h"), ("24h", "24 h")]
H_SHORT = {f"{g}_vs_CON_{t}": f"{'Endur' if g == 'EE' else 'Resist'}.{p} - Control.{p} (delta-delta)"
           for g in ("EE", "RE") for t, p in zip(["15to45min", "3.5to4h", "24h"], ["post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"])}
DOWN9 = MITO9

# ======================================================================= cross-check vs teammate + package CAMERA
if "crosscheck" in parts:
    _PART[0] = "crosscheck"
    pk = pkg_genesets()
    rows = []
    for g_ in ("EE", "RE"):
        for t, tl in H_TIMES:
            c = f"{g_}_vs_CON_{t}"
            p = col_frame("human_acute", "VL", "PROT", c); r = col_frame("human_acute", "VL", "RNA", c)
            p9, r9 = p[p.gene_symbol_human.isin(DOWN9)], r[r.gene_symbol_human.isin(DOWN9)]
            rows.append(dict(contrast=c, prot_n_matched=len(p9), prot_uniprot=";".join(p9.sort_values("gene_symbol_human").feature_id),
                             prot_n_multi_isoform=int((p9.n_collapsed > 1).sum()), prot_n_fdr05=int((p9.fdr < 0.05).sum()),
                             prot_n_positive=int((p9.logFC > 0).sum()),
                             rna_n_matched=len(r9), rna_n_positive=int((r9.logFC > 0).sum()), rna_n_fdr05=int((r9.fdr < 0.05).sum()),
                             rna_n_positive_fdr05=int(((r9.logFC > 0) & (r9.fdr < 0.05)).sum())))
    cc = pd.DataFrame(rows)
    open(f"{TAB}/crosscheck_down9.md", "w").write(cc.drop(columns="prot_uniprot").to_markdown(index=False) + "\n")
    print(cc.drop(columns="prot_uniprot").to_string()); print(cc.prot_uniprot.iloc[0])
    for rr in cc.itertuples():
        tag = rr.contrast.replace("_vs_CON_", "_")
        for f in ["prot_n_matched", "prot_n_fdr05", "prot_n_positive", "rna_n_positive", "rna_n_fdr05", "rna_n_positive_fdr05"]:
            num(f"XC_{tag}_{f}", int(getattr(rr, f)), "", "Malenfant DOWN-9, published BH adj p")

    # OXPHOS CAMERA: package precomputed vs our port (package-style input) vs our B3 input
    cam = pd.read_csv("data/raw/human_muscle_camera_v1.csv", low_memory=False)
    f2g = pd.read_csv("data/raw/human_feature_to_gene_muscle_v1.csv", low_memory=False)[["assay", "feature_id", "gene_symbol"]]
    gobp = pk["GOBP_OXIDATIVE_PHOSPHORYLATION"]
    kegg = extra_genesets()["KEGG_OXIDATIVE_PHOSPHORYLATION"]
    rows = []
    for layer, fn, assay in [("RNA", "rna", "transcript-rna-seq"), ("PROT", "prot", "prot-pr")]:
        raw = pd.read_csv(f"data/raw/human_muscle_{fn}_DA_v1.csv", usecols=["contrast_short", "feature_id", "z.std"])
        raw = raw.merge(f2g[f2g.assay == assay].drop(columns="assay").drop_duplicates(), on="feature_id", how="left")
        raw["new_id"] = raw.gene_symbol.fillna(raw.feature_id)
        for c, short in H_SHORT.items():
            x = raw[raw.contrast_short == short].copy()
            x = x.reindex(x["z.std"].abs().sort_values(ascending=False).index).drop_duplicates("new_id")
            t_pkgstyle = camera_pr(x["z.std"].to_numpy(), x.new_id.isin(gobp).to_numpy())[0]
            pre = cam[(cam.assay == assay) & (cam.contrast_short == short) & (cam.set == "GOBP_OXIDATIVE_PHOSPHORYLATION")]
            d = col_frame("human_acute", "VL", layer, c)
            rows.append(dict(layer=layer, contrast=c,
                             pkg_precomputed_t=float(pre.t.iloc[0]) if len(pre) else np.nan,
                             pkg_set_size=int(pre.set_size.iloc[0]) if len(pre) else np.nan,
                             our_port_pkg_input_t=t_pkgstyle, n_gobp_in_input=int(x.new_id.isin(gobp).sum()),
                             our_join_gobp_t=camera_pr(d.stat.to_numpy(), d.gene_symbol_human.isin(gobp).to_numpy())[0],
                             our_join_kegg_t=camera_pr(d.stat.to_numpy(), d.gene_symbol_human.isin(kegg).to_numpy())[0]))
    ox = pd.DataFrame(rows)
    ox["abs_diff_port_vs_pkg"] = (ox.our_port_pkg_input_t - ox.pkg_precomputed_t).abs()
    open(f"{TAB}/crosscheck_oxphos_camera.md", "w").write(ox.to_markdown(index=False, floatfmt=".3g") + "\n")
    print(ox.round(3).to_string())
    num("XC_camera_max_absdiff_port_vs_pkg", float(ox.abs_diff_port_vs_pkg.max()), "", "our cameraPR port on package-style input vs CAMERA_RESULTS")
    for rr in ox.itertuples():
        tag = f"{rr.layer}_{rr.contrast.replace('_vs_CON_', '_')}"
        num(f"XC_pkgcam_{tag}", round(rr.pkg_precomputed_t, 2)); num(f"XC_joinkegg_{tag}", round(rr.our_join_kegg_t, 2))
        num(f"XC_joingobp_{tag}", round(rr.our_join_gobp_t, 2))
    pd.DataFrame(NUMS).to_csv(f"{TAB}/crosscheck_numbers.csv", index=False)

# ======================================================================= hyperaemia / blood-contamination diagnostic
BLOOD = ["ALB", "HBB", "HBA1", "HBA2", "FGB", "FGA", "SERPINA1", "APOA1"]
if "hyperemia" in parts:
    _PART[0] = "hyperemia"
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    pk = pkg_genesets(); ex = extra_genesets()
    hcols = [(f"{g}_vs_CON_{t}", g, tl) for g in ("EE", "RE") for t, tl in H_TIMES]
    # (a) marker genes
    rows = []
    for c, g_, tl in hcols:
        for layer in ("RNA", "PROT"):
            d = col_frame("human_acute", "VL", layer, c).set_index("gene_symbol_human")
            for gene in BLOOD:
                rows.append(dict(contrast=c, layer=layer, gene=gene, logFC=d.logFC.get(gene, np.nan),
                                 stat=d.stat.get(gene, np.nan), fdr_bh=d.fdr_bh.get(gene, np.nan)))
    bm = pd.DataFrame(rows)
    bm.to_csv(f"{TAB}/hyperemia_markers.csv", index=False)
    # (b) blood sets + OXPHOS full vs muscle-intrinsic universe
    intrinsic = pk["MITOCARTA_ALL"] | pk["GOCC_CONTRACTILE_FIBER"] | pk["GOCC_SARCOPLASMIC_RETICULUM"]
    bsets = {"HALLMARK_COAGULATION": ex["HALLMARK_COAGULATION"], "HALLMARK_COMPLEMENT": ex["HALLMARK_COMPLEMENT"],
             "GOBP_COMPLEMENT_ACTIVATION": pk["GOBP_COMPLEMENT_ACTIVATION"]}
    kegg = ex["KEGG_OXIDATIVE_PHOSPHORYLATION"]
    rows = []
    for c, g_, tl in hcols:
        for layer in ("RNA", "PROT"):
            d = col_frame("human_acute", "VL", layer, c)
            rec = dict(contrast=c, layer=layer)
            for bn, bg in bsets.items():
                rec[f"t_{bn}"] = camera_pr(d.stat.to_numpy(), d.gene_symbol_human.isin(bg).to_numpy())[0]
            rec["t_OXPHOS_full"] = camera_pr(d.stat.to_numpy(), d.gene_symbol_human.isin(kegg).to_numpy())[0]
            di = d[d.gene_symbol_human.isin(intrinsic)]
            rec["n_universe_intrinsic"] = len(di); rec["n_universe_full"] = len(d)
            rec["t_OXPHOS_intrinsic"] = camera_pr(di.stat.to_numpy(), di.gene_symbol_human.isin(kegg).to_numpy())[0]
            rec["median_stat_intrinsic"] = float(di.stat.median()); rec["median_stat_nonintrinsic"] = float(d[~d.gene_symbol_human.isin(intrinsic)].stat.median())
            rows.append(rec)
    hs = pd.DataFrame(rows)
    open(f"{TAB}/hyperemia_sets.md", "w").write(hs.to_markdown(index=False, floatfmt=".3g") + "\n")
    print(hs.round(2).to_string())
    for rr in hs.itertuples():
        tag = f"{rr.layer}_{rr.contrast.replace('_vs_CON_', '_')}"
        for f in ["t_HALLMARK_COAGULATION", "t_HALLMARK_COMPLEMENT", "t_GOBP_COMPLEMENT_ACTIVATION", "t_OXPHOS_full", "t_OXPHOS_intrinsic",
                  "median_stat_intrinsic", "median_stat_nonintrinsic"]:
            num(f"HY_{tag}_{f}", round(float(getattr(rr, f)), 2), "robust_hyperemia.png")
    num("HY_n_universe_intrinsic_RNA", int(hs[hs.layer == "RNA"].n_universe_intrinsic.iloc[0]), "robust_hyperemia.png")
    num("HY_n_universe_intrinsic_PROT", int(hs[hs.layer == "PROT"].n_universe_intrinsic.iloc[0]), "robust_hyperemia.png")
    for (c, layer), d in bm.groupby(["contrast", "layer"]):
        tag = f"{layer}_{c.replace('_vs_CON_', '_')}"
        num(f"HY_{tag}_markers_up", int((d.logFC > 0).sum()), "robust_hyperemia.png", "of blood markers measured")
        num(f"HY_{tag}_markers_measured", int(d.logFC.notna().sum()), "robust_hyperemia.png")
        num(f"HY_{tag}_markers_up_fdr05", int(((d.logFC > 0) & (d.fdr_bh < 0.05)).sum()), "robust_hyperemia.png")
    print(bm.pivot_table(index=["layer", "gene"], columns="contrast", values="stat").round(1).to_string())

    # figure
    C = ["#2a78d6", "#eb6834"]
    xl = [f"{g_} {tl}" for c, g_, tl in hcols]
    def segplot(ax, y, label=None, **kw):
        y = np.asarray(y, float)
        ax.plot(range(3), y[:3], label=label, **kw); ax.plot(range(3, 6), y[3:], **kw)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axes[0]
    for k, layer in enumerate(("RNA", "PROT")):
        d = bm[bm.layer == layer]
        for gi, gene in enumerate(BLOOD):
            s_ = d[d.gene == gene].set_index("contrast").stat.reindex([h[0] for h in hcols])
            ax.scatter(np.arange(6) + (k - .5) * 0.3, s_, s=16, color=C[k], alpha=0.8, lw=0,
                       label=f"{layer} ({d[d.stat.notna()].gene.nunique()} of 8 genes measured)" if gi == 0 else None)
        med = d.groupby("contrast").stat.median().reindex([h[0] for h in hcols])
        for seg in (slice(0, 3), slice(3, 6)):  # EE and RE drawn as separate segments
            ax.plot((np.arange(6) + (k - .5) * 0.3)[seg], med.to_numpy()[seg], color=C[k], lw=2, marker="_", ms=14)
    ax.axhline(0, color="#52514e", lw=0.8); ax.axvline(2.5, color="#0b0b0b", lw=0.8)
    ax.set_xticks(range(6), xl, rotation=45, ha="right", fontsize=8); ax.set_ylabel("gene stat (moderated t)")
    ax.set_title("A. Blood/plasma marker genes\n(ALB, HBB, HBA1/2, FGA/B, SERPINA1, APOA1); line = median", loc="left", fontsize=9)
    ax.legend(fontsize=7.5, loc="upper right")
    ax = axes[1]
    for k, layer in enumerate(("RNA", "PROT")):
        d = hs[hs.layer == layer].set_index("contrast").reindex([h[0] for h in hcols])
        segplot(ax, d.t_HALLMARK_COAGULATION, color=C[k], lw=2, marker="o", label=f"{layer}: Hallmark coagulation")
        segplot(ax, d.t_HALLMARK_COMPLEMENT, color=C[k], lw=2, marker="s", ls="--", label=f"{layer}: Hallmark complement")
    ax.axhline(0, color="#52514e", lw=0.8); ax.axvline(2.5, color="#0b0b0b", lw=0.8)
    for v in (2, -2): ax.axhline(v, color="#c3c2b7", lw=0.8, ls=":")
    ax.set_xticks(range(6), xl, rotation=45, ha="right", fontsize=8); ax.set_ylabel("cameraPR t (set vs rest)")
    ax.set_title("B. Blood-derived gene sets", loc="left", fontsize=9); ax.legend(fontsize=7)
    ax = axes[2]
    for k, layer in enumerate(("RNA", "PROT")):
        d = hs[hs.layer == layer].set_index("contrast").reindex([h[0] for h in hcols])
        segplot(ax, d.t_OXPHOS_full, color=C[k], lw=2, marker="o", label=f"{layer}: vs all genes")
        segplot(ax, d.t_OXPHOS_intrinsic, color=C[k], lw=2, marker="D", ls="--", mfc="white",
                label=f"{layer}: vs muscle-intrinsic genes only")
    ax.axhline(0, color="#52514e", lw=0.8); ax.axvline(2.5, color="#0b0b0b", lw=0.8)
    for v in (2, -2): ax.axhline(v, color="#c3c2b7", lw=0.8, ls=":")
    ax.set_xticks(range(6), xl, rotation=45, ha="right", fontsize=8); ax.set_ylabel("OXPHOS cameraPR t")
    ax.set_title("C. OXPHOS drop with a muscle-intrinsic background\n(MitoCarta3 + contractile fibre + SR)", loc="left", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig("fig/robust_hyperemia.png"); plt.close(fig)
    caption("robust_hyperemia.png", "Hyperaemia/blood-contamination diagnostic, human VL (delta-delta vs control). A: moderated "
            "t of 8 blood/plasma marker genes per contrast (point = gene; line = median; blue RNA, orange PROT). B: cameraPR "
            "t of Hallmark coagulation and complement sets. C: KEGG OXPHOS cameraPR t against all measured genes (solid) vs "
            "against a muscle-intrinsic background only (MitoCarta 3.0 MitoPathways + GO:CC contractile fibre + sarcoplasmic "
            "reticulum; dashed). A drop that vanishes with the intrinsic background = dilution of muscle by non-muscle material, "
            "not OXPHOS-specific regulation.")
    pd.DataFrame(NUMS).to_csv(f"{TAB}/hyperemia_numbers.csv", index=False)
