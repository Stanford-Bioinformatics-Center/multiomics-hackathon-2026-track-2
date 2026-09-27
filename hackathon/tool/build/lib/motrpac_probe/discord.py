"""`mprobe discord`: omic-discordance report for one species x tissue (generalises scripts 04/06/08/11 + lag/).

Panels (order fixed): a timescale, b agreement by timepoint, c agreement by gene property, d which layer
responded and the detection-power decomposition, e protein-without-RNA calls, f does early RNA predict later
protein, g sex without thresholds, h caveats. Every number shown is also written to tables/*.csv.
"""
import datetime as dt
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import minimize  # noqa: E402
from scipy.stats import mannwhitneyu, spearmanr  # noqa: E402

from . import __version__, core, layers, legacy, store  # noqa: E402
from .caveats import GUARDRAILS, fixed_caveats  # noqa: E402
from .core import EXERCISE_KINDS, LAYER_WORDS, TISSUE_WORDS, col_words  # noqa: E402
from .paths import OUT  # noqa: E402

LAYER_ORDER = ["RNA", "PROT", "PHOSPHO", "METAB"]


def _cols(S, species, tissue):
    ds = "human_acute" if species == "human" else "rat_train"
    c = S.cols
    return c[(c.dataset == ds) & (c.tissue == tissue) & c.kind.isin(EXERCISE_KINDS)].sort_values("col_order")


def _series_key(c):
    """(series, x position, x label): human series = EE/RE, x = time; rat series = F/M, x = week."""
    if c.dataset == "rat_train":
        return c.sex, int(c.time[:-1]), f"{c.time[:-1]} wk"
    order = {"15to45min": 0, "3.5to4h": 1, "24h": 2, "10min": -1, "during20min": -3, "during40min": -2}
    t = c.contrast.split("_")[-1]
    return c.group, order.get(t, 9), core.H_TIME_WORDS.get(t, t)


# ------------------------------------------------------------------------------------------------ a
def timescale(S, cols, cutoff):
    rows = []
    for c in cols.itertuples():
        pos = S.col(c.column_id)
        n = len(pos)
        k = int((S.fdr[pos] < cutoff).sum())
        s, x, xl = _series_key(c)
        rows.append(dict(column_id=c.column_id, layer=c.layer, series=s, x=x, time=xl, n_measured=n, n_significant=k,
                         fraction=k / n if n else np.nan))
    return pd.DataFrame(rows)


def fig_timescale(ts, species):
    from .render import LAYER_COLORS, INK
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    for (layer, s), d in ts.groupby(["layer", "series"], sort=False):
        d = d.sort_values("x")
        y = d.fraction.replace(0, np.nan)
        ax.plot(d.x, y, marker="o", color=LAYER_COLORS.get(layer, INK),
                ls="-" if s in ("F", "EE") else "--", label=f"{LAYER_WORDS.get(layer, layer)}, {s}")
        zero = d[d.fraction == 0]
        ax.scatter(zero.x, np.full(len(zero), 1e-4), marker="v", color=LAYER_COLORS.get(layer, INK), s=24)
    ax.set_yscale("log")
    ax.set_ylim(5e-5, 1)
    xs = ts.drop_duplicates("x").sort_values("x")
    ax.set_xticks(xs.x, xs.time)
    ax.set_ylabel("fraction of measured features at BH FDR < cutoff")
    ax.set_xlabel("time after bout" if species == "human" else "weeks of training")
    ax.legend(fontsize=9, ncol=2, frameon=False)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ b
def pairs_for(S, cols):
    c = cols
    out = []
    for p in c[c.layer == "PROT"].itertuples():
        r = c[(c.layer == "RNA") & (c.contrast == p.contrast)]
        if len(r):
            out.append((r.column_id.iloc[0], p.column_id))
    return out


def agreement(S, pairs, sig_genes, seed=0):
    rng = np.random.default_rng(seed)
    rows, merged = [], {}
    for rna, prot in pairs:
        m = layers._merged(S, rna, prot)
        merged[(rna, prot)] = m
        c = S.cols.loc[prot]
        s, x, xl = _series_key(c)
        for subset, mm in [("all genes", m)] + ([("signature", m[m.gene_symbol_human.isin(sig_genes)])] if sig_genes else []):
            if len(mm) < 5:
                rows.append(dict(rna=rna, prot=prot, series=s, x=x, time=xl, subset=subset, n=len(mm)))
                continue
            bs = legacy.boot_rho(mm.stat_rna, mm.stat_prot, rng)
            rows.append(dict(rna=rna, prot=prot, series=s, x=x, time=xl, subset=subset, n=len(mm),
                             rho=legacy.rank_rho(mm.stat_rna, mm.stat_prot), lo=np.percentile(bs, 2.5),
                             hi=np.percentile(bs, 97.5)))
    return pd.DataFrame(rows), merged


def fig_agreement(ag):
    from .render import ACCENT, MUTED, INK
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for (sub, s), d in ag.dropna(subset=["rho"]).groupby(["subset", "series"], sort=False):
        d = d.sort_values("x")
        off = (-0.08 if s in ("F", "EE") else 0.08) + (0.02 if sub == "signature" else 0)
        ax.errorbar(d.x + off, d.rho, yerr=[(d.rho - d.lo).clip(lower=0), (d.hi - d.rho).clip(lower=0)], marker="o", ms=5, capsize=2,
                    color=ACCENT if sub == "signature" else MUTED, ls="-" if s in ("F", "EE") else "--",
                    label=f"{sub}, {s}")
    ax.axhline(0, color=INK, lw=0.8)
    xs = ag.drop_duplicates("x").sort_values("x")
    ax.set_xticks(xs.x, xs.time)
    ax.set_ylabel("Spearman ρ, RNA vs protein stat")
    ax.legend(fontsize=9, frameon=False)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ c
def gene_properties(S, tissue, genes):
    ann = S.ann.reindex(genes)
    from .store import GTEX_TISSUE
    gt = GTEX_TISSUE.get(tissue)
    p = pd.DataFrame(index=pd.Index(genes, name="gene"))
    p["mitocarta3"] = ann.mitocarta3.fillna(False).astype(int).to_numpy()
    p["stable_complex"] = ann.stable_complex.fillna(False).astype(int).to_numpy()
    p["secreted"] = ann.secreted_proxy.fillna(False).astype(int).to_numpy()
    p["gtex_log2_tpm"] = np.log2(ann[f"gtex_tpm|{gt}"].to_numpy(float) + 1) if gt else np.nan
    # phosphosites: number of phospho features collapsed into the gene in this tissue (0 if in proteome only)
    c = S.cols
    ph = c[(c.tissue == tissue) & (c.layer == "PHOSPHO")].column_id
    if len(ph):
        f = S.frame(ph.iloc[0])
        nsite = f.set_index(f.gene_symbol_human.astype(str)).n_collapsed
        p["phosphosites"] = pd.Series(genes, index=genes).map(nsite).fillna(0).to_numpy()
    else:
        p["phosphosites"] = np.nan
    return p


STRATA = {  # scripts/08 STRATA + baseline-count tertile
    "mitochondrial (MitoCarta3.0)": ("mitocarta3", [0, 1], ["no", "yes"], 0),
    "stable-complex subunit (GO:CC)": ("stable_complex", [0, 1], ["no", "yes"], 0),
    "secreted/extracellular (GO:CC proxy)": ("secreted", [0, 1], ["no", "yes"], 0),
    "GTEx abundance tertile (matched tissue)": ("gtex_tertile", [1.0, 2.0, 3.0], ["low", "mid", "high"], 1.0),
    "baseline RNA tertile (this dataset)": ("base_tertile", [1.0, 2.0, 3.0], ["low", "mid", "high"], 1.0),
    "phosphosites detected": ("phospho_bin", ["0", "1-3", ">3"], ["0", "1-3", ">3"], "0"),
}


def _tert(v):
    out = pd.Series(np.nan, index=v.index)
    ok = v.notna()
    if ok.sum() >= 30:
        out[ok] = pd.qcut(v[ok].rank(method="first"), 3, labels=[1, 2, 3]).astype(float)
    return out


def stratified(S, tissue, merged, pairs_sel, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for rna, prot in pairs_sel:
        m = merged[(rna, prot)].copy()
        props = gene_properties(S, tissue, m.gene_symbol_human.astype(str).tolist())
        m = pd.concat([m.reset_index(drop=True), props.reset_index(drop=True)], axis=1)
        base = S.frame(rna).set_index(S.frame(rna).gene_symbol_human.astype(str)).baseline_expr
        m["base"] = m.gene_symbol_human.astype(str).map(base).to_numpy(float)
        m["gtex_tertile"] = _tert(m.gtex_log2_tpm)
        m["base_tertile"] = _tert(m.base)
        m["phospho_bin"] = pd.cut(m.phosphosites, [-0.5, 0.5, 3.5, np.inf], labels=["0", "1-3", ">3"]).astype(object)
        lab = col_words(S.cols.loc[prot]).split(", ", 1)[1]
        bs_all = legacy.boot_rho(m.stat_rna, m.stat_prot, rng)
        rows.append(dict(comparison=lab, stratum="all genes", level="all", n=len(m),
                         rho=legacy.rank_rho(m.stat_rna, m.stat_prot), lo=np.percentile(bs_all, 2.5),
                         hi=np.percentile(bs_all, 97.5)))
        for sname, (col, levels, labs, ref) in STRATA.items():
            boots, recs = {}, {}
            for lev, lb in zip(levels, labs):
                s = m[m[col] == lev]
                if len(s) < 10:
                    recs[lev] = dict(comparison=lab, stratum=sname, level=lb, n=len(s))
                    continue
                boots[lev] = legacy.boot_rho(s.stat_rna, s.stat_prot, rng)
                recs[lev] = dict(comparison=lab, stratum=sname, level=lb, n=len(s),
                                 rho=legacy.rank_rho(s.stat_rna, s.stat_prot), lo=np.percentile(boots[lev], 2.5),
                                 hi=np.percentile(boots[lev], 97.5))
            for lev in levels:
                if lev != ref and lev in boots and ref in boots:
                    d = boots[lev] - boots[ref]
                    recs[lev].update(delta_vs_ref=recs[lev]["rho"] - recs[ref]["rho"], delta_lo=np.percentile(d, 2.5),
                                     delta_hi=np.percentile(d, 97.5))
                rows.append(recs[lev])
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ d
def _logit_fit(X, y, lam=1.0):
    w_pos = 0.5 / max(y.mean(), 1e-9)
    w_neg = 0.5 / max(1 - y.mean(), 1e-9)
    sw = np.where(y == 1, w_pos, w_neg)

    def f(b):
        z = X @ b[1:] + b[0]
        ll = sw * (y * -np.logaddexp(0, -z) + (1 - y) * -np.logaddexp(0, z))
        return -ll.sum() + lam * 0.5 * (b[1:] ** 2).sum()
    return minimize(f, np.zeros(X.shape[1] + 1), method="L-BFGS-B").x


def _prep(D, cols, tr):
    X = D[cols].to_numpy(float)
    med = np.nanmedian(X[tr], axis=0)
    ind = np.isnan(X)
    X = np.where(ind, med, X)
    mu, sd = X[tr].mean(0), X[tr].std(0)
    sd[sd == 0] = 1
    X = (X - mu) / sd
    keep = ind[tr].any(0)
    return np.hstack([X, ind[:, keep].astype(float)])


def auc(y, p):
    pos, neg = p[y == 1], p[y == 0]
    if not len(pos) or not len(neg):
        return np.nan
    return mannwhitneyu(pos, neg).statistic / (len(pos) * len(neg))


def cv_auc(D, cols, groups, folds=5, seed=0):
    """L2 logistic regression (class-balanced, median-imputed + missing indicators, standardised; scripts/08
    LogReg_L2), 5-fold CV with folds by gene; returns out-of-fold AUC."""
    rng = np.random.default_rng(seed)
    ug = np.unique(groups)
    fold_of = dict(zip(ug, rng.permutation(len(ug)) % folds))
    f = np.array([fold_of[g] for g in groups])
    y = D.y.to_numpy()
    oof = np.full(len(D), np.nan)
    for k in range(folds):
        tr, te = f != k, f == k
        if len(np.unique(y[tr])) < 2 or not te.any():
            continue
        X = _prep(D, cols, tr)
        b = _logit_fit(X[tr], y[tr])
        oof[te] = X[te] @ b[1:] + b[0]
    ok = np.isfinite(oof)
    return auc(y[ok], oof[ok])


def which_layer(S, tissue, species, merged, cutoff):
    rows, D = [], []
    for (rna, prot), m in merged.items():
        m = m.copy()
        m["label"] = layers.classify(m, cutoff)
        c = S.cols.loc[prot]
        s, x, xl = _series_key(c)
        rec = dict(comparison=col_words(c).split(", ", 1)[1], series=s, x=x, n_both=len(m))
        for k in ["concordant", "RNA_only", "PROT_only", "opposite", "ns"]:
            rec[k] = int((m.label == k).sum())
        rows.append(rec)
        fr = S.frame(rna).set_index(S.frame(rna).gene_symbol_human.astype(str))
        fp = S.frame(prot).set_index(S.frame(prot).gene_symbol_human.astype(str))
        g = m.gene_symbol_human.astype(str)
        m["rna_baseline"] = g.map(fr.baseline_expr).to_numpy(float)
        m["prot_missing"] = g.map(fp.prot_n_missing).to_numpy(float)
        m["n_collapsed_rna"] = g.map(fr.n_collapsed).to_numpy(float)
        m["n_collapsed_prot"] = g.map(fp.n_collapsed).to_numpy(float)
        m["comparison"] = rec["comparison"]
        D.append(m)
    counts = pd.DataFrame(rows)
    D = pd.concat(D, ignore_index=True)
    props = gene_properties(S, tissue, D.gene_symbol_human.astype(str).tolist())
    D = pd.concat([D, props.reset_index(drop=True)], axis=1)
    det = ["rna_baseline", "prot_missing"] if species == "rat" else ["rna_baseline", "gtex_log2_tpm"]
    full = det + ["gtex_log2_tpm", "mitocarta3", "stable_complex", "secreted", "n_collapsed_rna", "n_collapsed_prot",
                  "phosphosites"]
    full = list(dict.fromkeys(c for c in full if D[c].notna().any()))
    det = [c for c in det if D[c].notna().any()]
    T = D[D.label.isin(["RNA_only", "PROT_only"])].copy()
    T["y"] = (T.label == "RNA_only").astype(int)
    model = dict(n=len(T), n_rna_only=int(T.y.sum()), n_prot_only=int((1 - T.y).sum()), detection_features=det,
                 full_features=full)
    if min(model["n_rna_only"], model["n_prot_only"]) >= legacy.MIN_CLASS:
        grp = T.gene_symbol_human.astype(str).to_numpy()
        model["auc_detection"] = cv_auc(T, det, grp)
        model["auc_full"] = cv_auc(T, full, grp)
        rng = np.random.default_rng(1)
        Tp = T.copy()
        Tp["y"] = rng.permutation(Tp.y.to_numpy())
        model["auc_label_shuffle"] = cv_auc(Tp, full, grp)
    else:
        model["note"] = f"a class has < {legacy.MIN_CLASS} gene x comparison rows; model not fitted (scripts/08 rule)"
    # descriptive: RNA-only share by detection strata
    desc = []
    Tb = T.copy()
    Tb["baseline tertile"] = _tert(Tb.rna_baseline).map({1.0: "low", 2.0: "mid", 3.0: "high"})
    for lev, d in Tb.groupby("baseline tertile"):
        desc.append(dict(property="RNA baseline tertile", level=lev, n=len(d), frac_rna_only=d.y.mean()))
    if Tb.prot_missing.notna().any():
        for lev, d in Tb.groupby(Tb.prot_missing > 0):
            desc.append(dict(property="proteomics missing values > 0", level="yes" if lev else "no", n=len(d),
                             frac_rna_only=d.y.mean()))
    return counts, model, pd.DataFrame(desc)


# ------------------------------------------------------------------------------------------------ g
def sex_panel(S, cols, cutoff):
    if not (cols.sex == "F").any() or not (cols.sex == "M").any():
        return pd.DataFrame(), None
    rows, pts = [], {}
    for layer in [l for l in LAYER_ORDER if l in set(cols.layer)]:
        for wk in ["1w", "2w", "4w", "8w"]:
            f = cols[(cols.layer == layer) & (cols.contrast == f"F_{wk}")]
            m = cols[(cols.layer == layer) & (cols.contrast == f"M_{wk}")]
            if not len(f) or not len(m):
                continue
            a, b = S.frame(f.column_id.iloc[0]), S.frame(m.column_id.iloc[0])
            x = a[["gene_symbol_human", "logFC", "fdr_bh"]].merge(b[["gene_symbol_human", "logFC", "fdr_bh"]],
                                                                    on="gene_symbol_human", suffixes=("_F", "_M")).dropna()
            if len(x) < 10:
                continue
            rho = spearmanr(x.logFC_F, x.logFC_M)[0]
            slope = np.polyfit(x.logFC_F, x.logFC_M, 1)[0]
            sf, sm = x.fdr_bh_F < cutoff, x.fdr_bh_M < cutoff
            fonly, monly = sf & ~sm, sm & ~sf
            same_f = float((np.sign(x.logFC_F[fonly]) == np.sign(x.logFC_M[fonly])).mean()) if fonly.any() else np.nan
            same_m = float((np.sign(x.logFC_F[monly]) == np.sign(x.logFC_M[monly])).mean()) if monly.any() else np.nan
            rows.append(dict(layer=LAYER_WORDS.get(layer, layer), week=wk, n_genes=len(x), spearman_F_vs_M=rho,
                             slope_M_on_F=slope, sig_both=int((sf & sm).sum()), sig_F_only=int(fonly.sum()),
                             sig_M_only=int(monly.sum()), F_only_same_sign_in_M=same_f, M_only_same_sign_in_F=same_m))
            if wk == "8w":
                pts[layer] = x
    return pd.DataFrame(rows), pts


def fig_sex(pts):
    from .render import ACCENT, INK
    n = len(pts)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.0))
    axes = np.atleast_1d(axes)
    for ax, (layer, x) in zip(axes, pts.items()):
        ax.scatter(x.logFC_F, x.logFC_M, s=3, color="#9a9a9a", lw=0, alpha=0.5)
        s = (x.fdr_bh_F < 0.05) | (x.fdr_bh_M < 0.05)
        ax.scatter(x.logFC_F[s], x.logFC_M[s], s=5, color=ACCENT, lw=0)
        lim = np.nanpercentile(np.abs(np.r_[x.logFC_F, x.logFC_M]), 99.5)
        ax.plot([-lim, lim], [-lim, lim], color=INK, lw=0.6, ls="--")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_xlabel("female logFC, 8 wk")
        ax.set_ylabel("male logFC, 8 wk")
        ax.set_title(f"{LAYER_WORDS.get(layer, layer)}: ρ = {spearmanr(x.logFC_F, x.logFC_M)[0]:.2f}", loc="left",
                     fontsize=11)
    fig.tight_layout()
    return fig


# ------------------------------------------------------------------------------------------------ run
def compute(species, tissue, signature=None, cutoff=0.05, layers_=("RNA", "PROT", "PHOSPHO", "METAB"), quiet=False,
            universe="all"):
    from .run import Timer
    T = Timer(quiet=quiet)
    S = core.get_store()
    try:
        from . import metab
        S = metab.metab_store()
    except Exception as e:  # metabolomics optional
        print(f"[mprobe] METAB layer unavailable: {e}", file=sys.stderr)
    T("store load")
    cols = _cols(S, species, tissue)
    cols = cols[cols.layer.isin(layers_)]
    if not len(cols):
        raise SystemExit(f"no exercise comparisons for {species} {tissue}")
    sig = None
    if signature:
        from . import signature as sigmod
        sig = sigmod.load(signature)
    R = dict(species=species, tissue=tissue, cutoff=cutoff, cols=cols, sig=sig, S=S, universe=universe)
    R["timescale"] = timescale(S, cols, cutoff)
    T("a timescale")
    pairs = pairs_for(S, cols)
    R["agreement"], merged = agreement(S, pairs, sig.genes if sig else [])
    T("b agreement")
    late = [p for p in pairs if S.cols.loc[p[1], "contrast"] in ("F_4w", "F_8w", "M_4w", "M_8w", "EE_vs_CON_24h",
                                                                  "RE_vs_CON_24h", "EE_vs_CON_3.5to4h", "RE_vs_CON_3.5to4h")]
    R["stratified"] = stratified(S, tissue, merged, late or pairs)
    T("c stratified rho")
    R["classes"], R["model"], R["model_desc"] = which_layer(S, tissue, species, merged, cutoff)
    T("d which layer")
    R["extra"] = {}
    try:
        from . import discord_extra as dx
        if species == "rat":
            R["extra"]["e"] = dx.protein_without_rna(S, tissue, cutoff=cutoff)
            T("e protein without RNA")
            R["extra"]["f"] = dx.lag_panel(S, tissue)
            T("f lag")
    except Exception as e:
        R["extra_error"] = f"{type(e).__name__}: {e}"
    R["sex"], R["sex_pts"] = sex_panel(S, cols, cutoff)
    T("g sex")
    try:
        from . import metab
        R["pathways"] = metab.pathway_panel(S, tissue=tissue, species=species)
        T("3-layer pathways")
    except Exception as e:
        R["pathways_error"] = f"{type(e).__name__}: {e}"
    R["timer"] = T
    return R


def write(R, outdir, command=""):
    from . import render
    outdir = Path(outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    S, sp, tis, cut = R["S"], R["species"], R["tissue"], R["cutoff"]
    fp = render.fmt_p
    written = set()

    def T(key, title, caption, frame, formats=None, max_rows=None):
        if key not in written:
            frame.to_csv(outdir / "tables" / f"{key}.csv", index=False)
            written.add(key)
        return dict(kind="table", title=title, caption=caption, csv=f"tables/{key}.csv",
                    html=render.df_to_html(frame, table_id=key, formats=formats or {}, max_rows=max_rows))

    def F(key, fig, title, how):
        return dict(kind="figure", id=f"fig-{key}", title=title, how_to_read=how, png=f"figures/{key}.png",
                    svg=render.save_fig(fig, outdir / "figures" / f"{key}.png"))

    tw = f"{'Human' if sp == 'human' else 'Rat'} {TISSUE_WORDS.get(tis, tis)}"
    secs = []
    with render.mpl_style():
        # a
        ts = R["timescale"]
        tv = ts[["layer", "series", "time", "n_measured", "n_significant", "fraction"]].rename(columns={
            "layer": "Layer", "series": "Group / sex", "time": "Time", "n_measured": "Features measured",
            "n_significant": f"Features at BH FDR < {cut}", "fraction": "Fraction"})
        tv["Layer"] = tv.Layer.map(lambda l: LAYER_WORDS.get(l, l))
        byl = ts.groupby("layer").fraction.max()
        intro = ["Largest fraction significant per layer: " + "; ".join(
            f"{LAYER_WORDS.get(l, l)} {v:.3g}" for l, v in byl.items()) + "."]
        secs.append(dict(id="timescale", title="a. Timescale: how much of each layer changes, and when", intro=intro, items=[
            F("timescale", fig_timescale(ts, sp), f"{tw}: fraction of features changed per layer and time",
              f"Fraction of measured features with BH FDR < {cut} (fdr_bh, within comparison), log scale; solid = "
              f"{'female' if sp == 'rat' else 'EE'}, dashed = {'male' if sp == 'rat' else 'RE'}; ▼ at the bottom = "
              "none significant. Expected shape: RNA responds in hours, protein in weeks."),
            T("timescale_table", "Features changed per layer and time", "", tv, formats={"Fraction": "{:.4f}"})]))
        # b
        ag = R["agreement"]
        av = ag[["subset", "series", "time", "n", "rho", "lo", "hi"]].rename(columns={
            "subset": "Genes", "series": "Group / sex", "time": "Time", "n": "n genes", "rho": "Spearman ρ",
            "lo": "95% CI low", "hi": "95% CI high"})
        secs.append(dict(id="agreement", title="b. RNA–protein agreement by timepoint", intro=[], items=[
            F("agreement", fig_agreement(ag), f"{tw}: RNA vs protein agreement per timepoint",
              "Spearman ρ of gene-level RNA and protein statistics over genes measured in both layers; bars = 95% "
              f"bootstrap CI over genes (B = {legacy.B}, scripts/08); blue = signature genes only, if given."),
            T("agreement_table", "Agreement per timepoint", "", av,
              formats={"Spearman ρ": "{:.3f}", "95% CI low": "{:.3f}", "95% CI high": "{:.3f}"})]))
        # c
        st = R["stratified"].rename(columns={"comparison": "Comparison", "stratum": "Gene property", "level": "Level",
                                             "n": "n genes", "rho": "ρ", "lo": "CI low", "hi": "CI high",
                                             "delta_vs_ref": "Δρ vs reference", "delta_lo": "Δ CI low",
                                             "delta_hi": "Δ CI high"})
        f3 = {c: "{:.3f}" for c in ["ρ", "CI low", "CI high", "Δρ vs reference", "Δ CI low", "Δ CI high"]}
        secs.append(dict(id="properties", title="c. Agreement by gene property", intro=[
            "Stratified Spearman ρ between RNA and protein statistics (scripts/08): within each level of a gene "
            "property, with 95% bootstrap CI over genes, and the difference from the reference level (first level) "
            "with a paired bootstrap CI. Shown for the later timepoints."], items=[
            T("stratified_rho", "Stratified RNA–protein agreement", "Levels with < 10 genes are not estimated.", st,
              formats=f3, max_rows=200)]))
        # d
        cl = R["classes"].rename(columns={"comparison": "Comparison", "n_both": "Genes in both layers",
                                          "concordant": "Concordant", "RNA_only": "RNA only", "PROT_only": "Protein only",
                                          "opposite": "Opposite", "ns": "Neither", "series": "Group / sex"}).drop(columns="x")
        md = R["model"]
        mrows = [("RNA-only vs protein-only rows", md["n"]), ("RNA only", md["n_rna_only"]),
                 ("Protein only", md["n_prot_only"]),
                 ("Detection-only features", ", ".join(md["detection_features"])),
                 ("Full-model features", ", ".join(md["full_features"]))]
        intro = ["Much of 'which layer responded' is which assay could see the gene."]
        if "auc_full" in md:
            mrows += [("AUC, detection-only (CV)", f"{md['auc_detection']:.3f}"), ("AUC, full (CV)", f"{md['auc_full']:.3f}"),
                      ("AUC, full with shuffled labels", f"{md['auc_label_shuffle']:.3f}")]
            share = (md["auc_detection"] - 0.5) / (md["auc_full"] - 0.5) if md["auc_full"] > 0.5 else np.nan
            intro.append(f"Predicting RNA-only vs protein-only genes: detection power alone (baseline mRNA and "
                         f"{'proteomics missingness' if sp == 'rat' else 'GTEx abundance'}) gives out-of-fold AUC "
                         f"{md['auc_detection']:.2f}; adding gene properties gives {md['auc_full']:.2f} (shuffled "
                         f"labels {md['auc_label_shuffle']:.2f}). Detection alone carries "
                         f"{100 * share:.0f}% of the full model's discrimination above chance." if np.isfinite(share) else "")
        else:
            mrows.append(("Model", md.get("note", "")))
        mt = pd.DataFrame(mrows, columns=["Quantity", "Value"])
        mt["Value"] = mt.Value.astype(str)
        dd = R["model_desc"].rename(columns={"property": "Property", "level": "Level", "n": "n rows",
                                             "frac_rna_only": "Fraction RNA-only"})
        secs.append(dict(id="which", title="d. Which layer responded, and why", intro=intro, items=[
            T("layer_classes", "Layer classes per comparison",
              f"BH FDR < {cut} in each layer: concordant = both, same sign; opposite = both, opposite sign "
              "(scripts/08 rule).", cl),
            T("detection_model", "Detection-power decomposition",
              "L2 logistic regression (class-balanced; median imputation with missing indicators; standardised; "
              "5-fold CV with folds by gene; scripts/08 LogReg_L2), pooled over comparisons.", mt),
            T("detection_strata", "RNA-only share by detection stratum", "Rows = RNA-only or protein-only gene × "
              "comparison calls.", dd, formats={"Fraction RNA-only": "{:.2f}"})]))
        # e, f
        for key, sid, title in [("e", "protonly", "e. Protein-without-RNA calls"),
                                ("f", "lag", "f. Does early RNA predict later protein?")]:
            X = R["extra"].get(key)
            if X is None:
                msg = ("Rat only (needs the training time course)." if sp == "human"
                       else R.get("extra_error", "not available"))
                secs.append(dict(id=sid, title=title, intro=[msg], items=[]))
                continue
            items = []
            for fk, fg in X.get("figures", {}).items():
                ttl, how = X.get("captions", {}).get(fk, (fk, ""))
                items.append(F(f"{key}_{fk}", fg, ttl, how))
            for tk, tb in X.get("tables", {}).items():
                ttl = X.get("captions", {}).get(tk, (tk.replace("_", " "), ""))
                items.append(T(f"{key}_{tk}", ttl[0], ttl[1], tb, max_rows=60))
            secs.append(dict(id=sid, title=title, intro=list(X.get("notes", [])), items=items))
        # g
        sx = R["sex"]
        if len(sx):
            sv = sx.rename(columns={"layer": "Layer", "week": "Week", "n_genes": "Genes", "spearman_F_vs_M": "Spearman ρ F vs M",
                                    "slope_M_on_F": "Slope M on F", "sig_both": "Significant in both",
                                    "sig_F_only": "Female only", "sig_M_only": "Male only",
                                    "F_only_same_sign_in_M": "Female-only: same sign in male",
                                    "M_only_same_sign_in_F": "Male-only: same sign in female"})
            f2 = {c: "{:.2f}" for c in ["Spearman ρ F vs M", "Slope M on F", "Female-only: same sign in male",
                                        "Male-only: same sign in female"]}
            items = [T("sex_table", "Female vs male, all genes", f"Significant = BH FDR < {cut}.", sv, formats=f2)]
            if R["sex_pts"]:
                items.insert(0, F("sex", fig_sex(R["sex_pts"]), f"{tw}: female vs male training response at 8 weeks",
                                  "Point = gene (blue = BH FDR < 0.05 in either sex); dashed = identity. Spearman ρ "
                                  "over all genes in the title."))
            s8 = sx[sx.week == "8w"]
            intro = [f"At 8 weeks, {r.layer}: {r.sig_F_only} female-only and {r.sig_M_only} male-only hits; "
                     f"{100 * r.F_only_same_sign_in_M:.0f}% of female-only and {100 * r.M_only_same_sign_in_F:.0f}% of "
                     f"male-only hits have the same sign in the other sex (all-gene ρ = {r.spearman_F_vs_M:.2f})."
                     for r in s8.itertuples() if np.isfinite(r.F_only_same_sign_in_M) and np.isfinite(r.M_only_same_sign_in_F)]
            secs.append(dict(id="sex", title="g. Sex without thresholds", intro=intro, items=items))
        else:
            secs.append(dict(id="sex", title="g. Sex without thresholds", intro=[
                "The human acute study reports exercise-vs-control contrasts pooled over sexes; no sex contrast is "
                "available in the public summary statistics."], items=[]))
        # 3-layer pathways
        if "pathways" in R and R["pathways"] is not None and len(R["pathways"]):
            from . import metab
            pw = R["pathways"]
            items = [F("pathways", metab.pathway_heatmap(pw, tis), f"{tw}: pathway response in RNA, protein and "
                       "metabolites", metab.PATHWAY_CAPTION),
                     T("pathways_table", "Pathway cameraPR t per layer and time", "", pw, max_rows=200)]
            secs.append(dict(id="pathways", title="Metabolic pathways across three layers", intro=[
                "Metabolite pools are not flux. The pathway map is many-to-many (tool/store/pathway_map.csv); genes "
                "and metabolites are never matched one-to-one."], items=items))
        # h
        cav = fixed_caveats() + ["All panels are group contrasts across genes, not per-animal or per-person.",
                                 "Layer classes depend on the FDR cutoff and on detection power (panel d)."]
        secs.append(dict(id="caveats", title="h. Caveats", intro=cav, items=[
            T("guardrails", "Safe statements and the unsafe upgrades to avoid", "From the team's planning workbook, "
              "generalised.", pd.DataFrame(GUARDRAILS, columns=["Safe statement", "Unsafe upgrade"]))]))
    R["timer"]("figures + tables")
    prov = {"tool": f"motrpac_probe {__version__}", "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
            "command": command, "repo_git_sha": store.git_sha(), "store": store.store_hash(),
            "legacy_scripts_sha256_16": legacy.script_hashes(),
            "parameters": {"species": sp, "tissue": tis, "cutoff": cut},
            "timings_s": {k: round(v, 2) for k, v in R["timer"].rows}}
    (outdir / "provenance.json").write_text(json.dumps(prov, indent=2))
    first = [f"{tw}: {len(R['cols'])} exercise comparisons across layers "
             f"{', '.join(LAYER_WORDS.get(l, l) for l in R['cols'].layer.unique())}.",
             "Much of 'which layer responded' is which assay could see the gene (panel d)."]
    ag = R["agreement"]
    agl = ag[ag.subset == "all genes"].dropna(subset=["rho"])
    head = dict(title=f"{tw}: how discordant are the layers, and why?",
                how_to_read="Panels a–b answer 'when does each layer move and how well do they agree'; panel d "
                            "answers 'how much of the disagreement is detection power'.",
                table_html="", sentences=[
                    ("Genome-wide RNA–protein ρ ranges " f"{agl.rho.min():.2f} to {agl.rho.max():.2f} across "
                     f"{len(agl)} timepoints.") if len(agl) else "No timepoint has both RNA and protein."])
    ctx = dict(title=f"MoTrPAC discordance: {tw}", subtitle="Which omic layer responds to exercise, when, and why "
               "they disagree", meta=[("Species", sp), ("Tissue", tis), ("Cutoff", str(cut)),
                                      ("Store", store.store_hash()), ("Command", command)],
               headline=head, first_screen=first, sections=secs, glossary=_glossary(),
               provenance_json=json.dumps(prov, indent=2))
    render.render_report(ctx, outdir / "report.html")
    R["timer"]("render")
    return outdir / "report.html"


def _glossary():
    from .run import GLOSSARY
    return GLOSSARY + [("Concordant / RNA only / protein only / opposite", "per gene and comparison: significant "
                        "(BH FDR < cutoff) in both layers with the same sign / only RNA / only protein / both with "
                        "opposite signs.")]


def main_discord(species, tissue, signature=None, cutoff=0.05, layers_=None, outdir=None, command="", quiet=False,
                 universe="all"):
    R = compute(species, tissue, signature=signature, cutoff=cutoff,
                layers_=layers_ or ("RNA", "PROT", "PHOSPHO", "METAB"), quiet=quiet, universe=universe)
    out = Path(outdir) if outdir else OUT / f"discord_{tissue}"
    return write(R, out, command=command), R
