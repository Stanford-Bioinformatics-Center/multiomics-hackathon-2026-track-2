"""`mprobe run`: one signature -> tool/out/NAME/{report.html, figures/, tables/, provenance.json}.

All prose in the report is templated from numbers computed here; nothing is free text generated per run.
"""
import datetime as dt
import hashlib
import html
import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import __version__, core, everywhere, figures, interactive, layers, legacy, nulls, ranked, signature, store
from .caveats import GUARDRAILS, fixed_caveats
from .core import CORE_LAYERS, CORE_TISSUES, EXERCISE_KINDS, col_words
from .paths import OUT

DEFAULTS = dict(cutoff=0.05, cap=4.0, nboot=1000, seed=20260926, nperm=2000, exact_symbols=False, context_groups=("phenotype",), min_n=10,
                tissues=CORE_TISSUES, universe="all", pool_sets=(), pool_groups=(), sections=None)
SECTIONS = ["timecourse", "coverage", "grid", "camera", "nulls", "layers", "pathways", "everywhere", "caveats", "sensitivity"]
METAB_EXTRA_TISSUES = ("PLASMA", "BLOOD")   # metabolite signatures: also rat plasma and human blood
PATHWAY_PANELS = [("SKM-GN", "rat"), ("VL", "human")]


def verdict(t, q, cutoff):
    if not np.isfinite(t) or not np.isfinite(q) or q >= cutoff:
        return "no set-level shift"
    return "opposed" if t > 0 else "same direction"


def intrinsic_universe(S):
    g = S.gs
    return g["MITOCARTA_ALL"] | g["GOCC_CONTRACTILE_FIBER"] | g["GOCC_SARCOPLASMIC_RETICULUM"]


def _hash_file(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()[:16]


class Timer:
    def __init__(self, quiet=False):
        self.t0, self.last, self.rows, self.quiet = time.time(), time.time(), [], quiet

    def __call__(self, what):
        now = time.time()
        self.rows.append((what, now - self.last))
        if not self.quiet:
            print(f"[mprobe] {what:<32s} {now - self.last:6.1f} s", file=sys.stderr)
        self.last = now


def compute(sig_path, name=None, **kw):
    """All analyses for one signature. Returns a results dict (no files written)."""
    o = {**DEFAULTS, **{k: v for k, v in kw.items() if v is not None}}
    T = Timer(quiet=o.get("quiet", False))
    try:
        from . import metab
        metab_mode = metab.is_metab_signature(sig_path)
    except ImportError:
        metab, metab_mode = None, False
    if metab_mode:  # metabolite signature: METAB columns of the metab-extended store
        S = metab.metab_store()
        T("store load")
        sig = metab.load_metab_signature(sig_path, name=name)
        T("signature mapping")
        core_cols = metab.metab_columns(S, tissues=list(o["tissues"]) + list(METAB_EXTRA_TISSUES))
    else:
        S = core.get_store()
        T("store load")
        rinfo = None
        if ranked.is_ranked(sig_path):  # full disease ranking (gene + t): team blood pipeline, generalised
            rk, sig, rinfo = ranked.load_ranked(sig_path, name=name, cutoff=o["cutoff"],
                                                exact_symbols=o.get("exact_symbols", False))
        else:
            sig = signature.load(sig_path, name=name, context_groups=tuple(o["context_groups"]))
        T("signature mapping")
        core_cols = S.select(tissues=o["tissues"], layers=CORE_LAYERS, kinds=EXERCISE_KINDS).column_id.tolist()
    cutoff, cap = o["cutoff"], o["cap"]
    R = dict(opts=o, sig=sig, core_cols=core_cols, S=S, metab_mode=metab_mode,
             ranked_mode=not metab_mode and rinfo is not None)
    if R["ranked_mode"]:
        R["rinfo"] = rinfo
        c = S.cols
        hx = (c.dataset == "human_acute") & (c.kind == "exercise vs control")
        focus = c[(hx & c.layer.isin(["RNA", "PROT", "PROT_OLINK"])) |
                  ((c.dataset == "rat_train") & c.tissue.isin(o["tissues"]) & c.layer.isin(CORE_LAYERS))]
        focus = focus.sort_values("col_order").column_id.tolist()
        perm = c[hx & (c.layer == "RNA")].column_id.tolist()
        ra = ranked.rank_association(S, rk, c.column_id.tolist(), nperm=o["nperm"], perm_cids=perm, seed=o["seed"])
        R["rank_all"] = ra
        R["rank_focus"] = ra[ra.column_id.isin(focus)].merge(
            ranked.reference_calibration(S, ra[ra.column_id.isin(focus)], ra), on="column_id", how="left")
        T("R rank association")
        dis = ranked.disease_pathways(rinfo["full"])
        R["pw_dis"] = dis
        R["pw_store"], R["pw_store_long"] = ranked.pathway_concordance(S, dis, focus, cutoff=o["cutoff"])
        hcols = [x for x in focus if x.startswith("human_acute")]
        R["pw_pkg"], R["pw_pkg_long"] = ranked.pathway_concordance(S, dis, hcols, exercise_side="precomputed",
                                                                   cutoff=o["cutoff"])
        T("P pathway concordance")
    uni = intrinsic_universe(S) if o["universe"] == "muscle-intrinsic" else None
    R["scores"] = core.score_columns(S, core_cols, sig, cap=cap, cutoff=cutoff, universe=uni)
    R["grid"] = core.grid_long(S, core_cols, sig, cap=cap)
    R["meff"] = core.effective_n(S, sig.genes, core_cols)
    M = R["meff"]["n_genes"]
    frac = R["meff"]["meff_nyholt"] / M if M else np.nan
    R["scores"]["sign_p_eff"] = [core.deflated_sign_p(r.n_measured, r.n_opposed, frac) for r in R["scores"].itertuples()]
    T("B/C grid + cameraPR")
    R["nulls"] = nulls.run_nulls(S, core_cols, sig, B=o["nboot"], seed=o["seed"])
    if o["pool_sets"] and not metab_mode:
        grp = set(o["pool_groups"])
        genes = [g for g in sig.shown if not grp or sig.groups.get(g) in grp]
        pools = nulls.legacy_pools(S, o["pool_sets"])
        R["pool"], _ = nulls.pool_null(S, core_cols, genes, pools, np.random.default_rng(o["seed"]), B=o["nboot"])
        R["pool_genes"] = genes
    T("D calibration nulls")
    R["coverage"], R["cov_blocks"] = layers.coverage(S, sig, core_cols)
    if metab_mode:  # gene-level detection flags and RNA-vs-protein pairs do not apply to metabolites
        R["flags"], R["disc"], R["disc_genes"] = pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    else:
        R["flags"] = layers.detection_flags(S, sig.shown)
        T("A coverage + detection flags")
        pairs = layers.layer_pairs(S)
        pairs = pairs[[S.cols.loc[p, "tissue"] in o["tissues"] or S.cols.loc[p, "tissue"] in ("ADIPOSE",)
                       for p in pairs.prot]]
        R["disc"], R["disc_genes"] = layers.discordance(S, sig, pairs, cutoff=cutoff)
        T("E layer discordance")
    ev_cids = S.cols[S.cols.layer == "METAB"].column_id.tolist() if metab_mode else None
    R["ev_cols"], R["ev_rank"], R["ev_summary"], R["ev_bykind"] = everywhere.run(S, sig, cutoff=cutoff,
                                                                                min_n=o["min_n"], cids=ev_cids)
    T("F everywhere")
    R["sens"] = sensitivity(S, sig, core_cols, R, o)
    T("sensitivity")
    R["myh"] = pd.DataFrame(columns=["column_id", "MYH7", "MYH2", "MYH1", "MYH4"]) if metab_mode \
        else fibre_markers(S, core_cols)
    R["pathways"] = {}
    if metab is not None:  # 3-layer pathway panel (independent of the signature)
        for tis, sp in PATHWAY_PANELS:
            try:
                R["pathways"][(tis, sp)] = metab.pathway_panel(None, tissue=tis, species=sp)
            except Exception as e:  # metabolomics store not built
                R["pathways_error"] = f"{type(e).__name__}: {e}"
        T("3-layer pathways")
    R["timer"] = T
    return R


def sensitivity(S, sig, cids, R, o):
    """Re-derive each core-column conclusion under alternative settings; flag conclusions that flip."""
    base = R["scores"].set_index("column_id")
    cut = o["cutoff"]
    alt_cut = 0.10 if cut == 0.05 else 0.05
    uni = intrinsic_universe(S)
    musc = [] if R.get("metab_mode") else \
        [c for c in cids if S.cols.loc[c, "tissue"] in ("VL", "SKM-GN", "SKM-VL", "HEART")]
    alt_uni = core.score_columns(S, musc, sig, cutoff=cut, universe=uni).set_index("column_id") if o["universe"] == "all" \
        else core.score_columns(S, musc, sig, cutoff=cut).set_index("column_id")
    nb = nulls.run_nulls(S, cids, sig, B=200 if o["nboot"] != 200 else 1000, seed=o["seed"] + 100, kinds=("class",))
    nb = nb.set_index("column_id")
    ncl = R["nulls"][R["nulls"].null == "class"].set_index("column_id")
    rows = []
    for c in cids:
        b = base.loc[c]
        v0 = verdict(b.camera_t, b.camera_fdr, cut)
        rows.append(dict(column_id=c, setting=f"FDR cutoff {cut} vs {alt_cut}", conclusion="set-level verdict",
                         base=v0, alternative=verdict(b.camera_t, b.camera_fdr, alt_cut)))
        if c in alt_uni.index:
            a = alt_uni.loc[c]
            rows.append(dict(column_id=c, setting="universe: all genes vs muscle-intrinsic" if o["universe"] == "all"
                             else "universe: muscle-intrinsic vs all genes", conclusion="set-level verdict", base=v0,
                             alternative=verdict(a.camera_t, a.camera_fdr, cut)))
        p0, p1 = ncl.loc[c, "pct_t"], nb.loc[c, "pct_t"]
        f = lambda p: "not estimable" if not np.isfinite(p) else ("above 95th percentile" if p >= 95 else "within null")
        rows.append(dict(column_id=c, setting=f"null size {o['nboot']} vs {200 if o['nboot'] != 200 else 1000}",
                         conclusion="pathway-class null", base=f(p0), alternative=f(p1)))
    s = pd.DataFrame(rows)
    s["flipped"] = s.base != s.alternative
    return s


def fibre_markers(S, cids):
    rows = []
    for c in cids:
        f = S.frame(c).set_index("gene_symbol_human")
        rec = {"column_id": c}
        for m in ["MYH7", "MYH2", "MYH1", "MYH4"]:
            rec[m] = f.stat.get(m, np.nan)
        rows.append(rec)
    return pd.DataFrame(rows)


# ================================================================================================ report
def write(R, outdir, command="", toggles=None):
    from . import render
    o, S, sig = R["opts"], R["S"], R["sig"]
    outdir = Path(outdir)
    (outdir / "figures").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)
    cut, cap = o["cutoff"], o["cap"]
    cids = R["core_cols"]
    cw = {c: col_words(S.cols.loc[c]) for c in S.cols.column_id}
    fp = render.fmt_p
    tables = {}

    def tab(key, df):
        df.to_csv(outdir / "tables" / f"{key}.csv", index=False)
        tables[key] = df
        return f"tables/{key}.csv"

    # ------------------------------------------------------------------ tables
    sc = R["scores"].copy()
    sc.insert(1, "comparison", sc.column_id.map(cw))
    nl = R["nulls"]
    for kind, lab in [("abundance", "abund"), ("class", "class")]:
        d = nl[nl.null == kind].set_index("column_id")
        sc[f"pct_{lab}"] = sc.column_id.map(d.pct_t)
        sc[f"stratum_{lab}"] = sc.column_id.map(d.stratum)
    sc["verdict"] = [verdict(t, q, cut) for t, q in zip(sc.camera_t, sc.camera_fdr)]
    sc["early_flag"] = sc.column_id.map(S.cols.early_human_muscle).astype("boolean").fillna(False).astype(bool)
    tab("column_scores", sc)
    tab("grid_long", R["grid"].assign(comparison=R["grid"].column_id.map(cw)))
    tab("nulls", nl.drop(columns=[c for c in nl.columns if c.startswith("_")]).assign(comparison=nl.column_id.map(cw)))
    tab("mapping", sig.table)
    tab("coverage", R["coverage"].rename(columns={b: f"{'Human' if b[0] == 'human_acute' else 'Rat'} {b[1]} {b[2]}"
                                                   for b in R["cov_blocks"]}))
    tab("detection_flags", R["flags"])
    tab("layer_pairs", R["disc"])
    tab("layer_genes", R["disc_genes"])
    tab("everywhere_columns", R["ev_cols"])
    tab("everywhere_ranking", R["ev_rank"])
    tab("sensitivity", R["sens"].assign(comparison=R["sens"].column_id.map(cw)))
    tab("fibre_markers", R["myh"].assign(comparison=R["myh"].column_id.map(cw)))
    if "pool" in R:
        tab("pool_null", R["pool"].assign(comparison=R["pool"].column_id.map(cw)))

    # ------------------------------------------------------------------ figures
    def fig(key, f):
        return render.save_fig(f, outdir / "figures" / f"{key}.png")

    with render.mpl_style():
        svg = {"coverage": fig("coverage", figures.coverage(R["coverage"].head(80), R["cov_blocks"])),
               "grid": fig("grid", figures.grid(S, R["grid"], cids, sig, cap, cut)),
               "camera": fig("camera", figures.camera(S, R["scores"], cids, cut)),
               "nulls": fig("nulls", figures.nulls(S, R["scores"], nl, cids)),
               "everywhere": fig("everywhere", figures.everywhere(R["ev_rank"]))}
        d = R["disc"]
        show = []
        for pick in ([("human_acute", "EE_vs_CON_24h"), ("rat_train", "F_8w")] if len(d) else []):
            m = d[(d.dataset == pick[0]) & (d.contrast == pick[1]) & ~d.cross & d.tissue.isin(["VL", "SKM-GN"])]
            if len(m) and m.n_sig_both.iloc[0] > 0:
                show.append((m.label.iloc[0], (m.rna.iloc[0], m.prot.iloc[0])))
        if len(d):
            svg["layers"] = fig("layers", figures.discordance(d, R["disc_genes"], show))
        # time course per signature gene (rat gastrocnemius training; human vastus lateralis after one bout)
        dirs = dict(zip(sig.shown, sig.shown_dirs))
        for key, tis, sp in [("traj_rat", "SKM-GN", "rat"), ("traj_human", "VL", "human")]:
            if R.get("metab_mode"):
                break
            Tt = core.trajectory_table(S, sig.shown, tis, species=sp, cutoff=cut)
            if len(Tt):
                gs = [g for g in sig.shown if g in set(Tt.gene)]
                if len(gs) > 12:
                    mx = Tt.groupby("gene").logFC.apply(lambda v: v.abs().max())
                    gs = [g for g in gs if g in set(mx.nlargest(12).index)]
                svg[key] = fig(key, figures.trajectories(Tt, gs, dirs, tis, species=sp))
                tab(key, Tt)
        if R.get("ranked_mode"):
            ref = R["rank_all"].merge(S.cols.reset_index(drop=True)[["column_id", "dataset", "tissue", "layer", "kind"]], on="column_id")
            svg["rank_assoc"] = fig("rank_assoc", figures.rank_association(
                S, R["rank_focus"], ref[ref.kind == "reference (non-exercise)"]))
            svg["pw_counts"] = fig("pw_counts", figures.pathway_counts(S, R["pw_store"]))
            tab("rank_association_all", R["rank_all"].assign(comparison=R["rank_all"].column_id.map(cw)))
            tab("disease_pathways", R["pw_dis"])
            tab("pathway_concordance_long", R["pw_store_long"])
            tab("pathway_concordance_precomputed_long", R["pw_pkg_long"])
        if R.get("pathways"):
            from . import metab
            for (tis, sp), pw in R["pathways"].items():
                if len(pw):
                    svg[f"pathways_{tis}"] = fig(f"pathways_{tis}", metab.pathway_heatmap(pw, tis))
                    tab(f"pathways_{tis}", pw)
    R["timer"]("figures")

    ctx = context(R, svg, tables, cw, fp, command, outdir)
    if toggles:
        ctx["toggles"] = toggles
    prov = provenance(R, command)
    (outdir / "provenance.json").write_text(json.dumps(prov, indent=2, default=str), encoding="utf-8")
    ctx["provenance_json"] = json.dumps(prov, indent=2, default=str)
    render.render_report(ctx, outdir / "report.html")
    R["timer"]("render")
    return outdir / "report.html"


def provenance(R, command):
    o = dict(R["opts"])
    o.pop("quiet", None)
    sp = store.load_provenance()
    return {
        "tool": f"motrpac_probe {__version__}",
        "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
        "command": command,
        "repo_git_sha": store.git_sha(),
        "signature": {"name": R["sig"].name, "sha256_16": _hash_file(R["opts"]["_path"]),
                      "path": str(R["opts"]["_path"])},
        "store": {"contrasts_parquet_sha256_16": store.store_hash(), "built_at": sp.get("built_at"),
                  "rows": sp.get("store_rows"), "columns": sp.get("store_columns")},
        "sources": {k: sp.get(k) for k in ["MotrpacRatTraining6moData", "MotrpacHumanPreSuspensionAnalysis",
                                            "human_data_collection",
                                            "join_table_version", "gtex", "mitocarta"]},
        "legacy_scripts_sha256_16": legacy.script_hashes(),
        "python": sys.version.split()[0],
        "packages": {m.__name__: m.__version__ for m in [np, pd] + _mods()},
        "parameters": {k: (list(v) if isinstance(v, tuple) else v) for k, v in o.items() if not k.startswith("_")},
        "timings_s": {k: round(v, 2) for k, v in R["timer"].rows},
    }


def _mods():
    import jinja2
    import matplotlib
    import pyarrow
    import scipy
    return [scipy, matplotlib, jinja2, pyarrow]


# ------------------------------------------------------------------------------------------------ context
def _tint(t):
    if not np.isfinite(t):
        return ""
    a = min(abs(t), 6) / 6 * 0.55
    return f"background: rgba(31,78,121,{a:.2f})" if t > 0 else f"background: rgba(165,28,48,{a:.2f})"


def headline_html(S, sc, cut):
    """Two compact matrices: rows = tissue x layer, cols = time. Cell = opposed/measured and cameraPR t."""
    parts = []
    sc = sc.set_index("column_id")
    for ds, title in [("human_acute", "Human acute bout (vastus lateralis, exercise vs control)"),
                      ("rat_train", "Rat endurance training (trained vs sedentary)")]:
        c = S.cols.loc[[x for x in sc.index if S.cols.loc[x, "dataset"] == ds]]
        if not len(c):
            continue
        times = list(dict.fromkeys(core.short_words(r) for r in c.itertuples()))
        rows = list(dict.fromkeys(zip(c.tissue, c.layer)))
        h = [f'<table class="headline"><caption>{html.escape(title)}</caption><thead><tr><th></th>']
        for t in times:
            flag = " *" if ("15–45" in t and ds == "human_acute") else ""
            h.append(f"<th>{html.escape(t)}{flag}</th>")
        h.append("</tr></thead><tbody>")
        for tis, lay in rows:
            h.append(f"<tr><th>{html.escape(core.TISSUE_WORDS.get(tis, tis))} {core.LAYER_WORDS.get(lay, lay)}</th>")
            for t in times:
                m = c[(c.tissue == tis) & (c.layer == lay)]
                m = [x for x in m.itertuples() if core.short_words(x) == t]
                if not m:
                    h.append('<td class="na">not measured</td>')
                    continue
                r = sc.loc[m[0].column_id]
                if not r.n_measured:
                    h.append('<td class="na">no measured genes</td>')
                    continue
                star = "*" if np.isfinite(r.camera_fdr) and r.camera_fdr < cut else ""
                tt = "t NA" if not np.isfinite(r.camera_t) else f"t {r.camera_t:+.1f}{star}"
                pc = "" if not np.isfinite(r.pct_class) else f"<br><span class=\"sub\">class pct {r.pct_class:.0f}</span>"
                h.append(f'<td style="{_tint(r.camera_t if star else 0)}">{r.n_opposed}/{r.n_measured} opp.<br>'
                         f'{tt}{pc}</td>')
            h.append("</tr>")
        h.append("</tbody></table>")
        parts.append("".join(h))
    return "\n".join(parts)


def context(R, svg, tables, cw, fp, command, outdir):
    from . import render
    o, S, sig = R["opts"], R["S"], R["sig"]
    cut = o["cutoff"]
    sc = tables["column_scores"]
    mp = sig.table
    n_in, n_map = len(mp), int((mp.status == "mapped").sum())
    unm = mp[mp.status != "mapped"]
    B = o["nboot"]

    # ---------------- headline sentences (templated)
    sent = [f"{n_map} of {n_in} input rows mapped to a MoTrPAC gene; {len(sig.genes)} are counted in the statistics"
            + (f" and {len(sig.shown) - len(sig.genes)} are shown as context ({', '.join(o['context_groups'])})"
               if len(sig.shown) > len(sig.genes) else "") + "."]
    for lay in [l for l in ["RNA", "PROT", "METAB"] if (S.cols.loc[sc.column_id, "layer"] == l).any()]:
        d = sc[(S.cols.loc[sc.column_id, "layer"] == lay).to_numpy() & ~sc.early_flag.to_numpy()]
        d = d[d.camera_t.notna()]
        if not len(d):
            continue
        best = d.loc[d.camera_t.idxmax()]
        nopp = int((d.verdict == "opposed").sum())
        nsame = int((d.verdict == "same direction").sum())
        sent.append(f"{core.LAYER_WORDS[lay][0].upper() + core.LAYER_WORDS[lay][1:]}: {nopp} of {len(d)} comparisons are opposed at the set level "
                    f"(cameraPR FDR < {cut}) and {nsame} move the same way as the signature; the strongest opposition "
                    f"is {best.comparison} (t = {best.camera_t:+.2f}, {best.n_opposed}/{best.n_measured} genes opposed, "
                    f"pathway-class null percentile {best.pct_class:.0f}).")
    headline = dict(title="Headline: does the exercise response oppose the signature, by layer and time?",
                    how_to_read=(f"Each cell: signature genes moved against their disease direction / genes measured, "
                                 f"and the signed cameraPR t (+ = exercise opposes the signature; * = BH FDR < {cut} "
                                 f"across the {len(sc)} comparisons shown). Blue tint = opposed, red tint = same "
                                 f"direction as the disease; untinted = no set-level shift. 'class pct' = percentile "
                                 f"of t among {B} random gene sets with the same MitoCarta / complex / secreted "
                                 f"composition (95+ = more opposed than its pathway class). * in a header = human "
                                 f"15–45 min, read as biopsy composition (see Caveats)."),
                    table_html=headline_html(S, sc, cut), sentences=sent)
    first = [f"Coverage: {n_map}/{n_in} mapped"
             + (f"; unmapped: {', '.join(unm.iloc[:, 0].astype(str).head(12))}" + (" …" if len(unm) > 12 else "")
                if len(unm) else "") + ".",
             "Caveat: 'opposed' is a sign comparison of group-level contrasts, not a therapeutic claim; see the "
             "specificity null (a random set from the same pathway class can look just as opposed)."]

    def T(key, title, caption, frame, formats=None, max_rows=None, row_class=None):
        if frame is None or len(frame) == 0:  # never render an empty table
            frame = pd.DataFrame({"Result": ["no rows: no signature gene is measured in enough comparisons"]})
        f = outdir / "tables" / f"{key}.csv"
        if key not in tables:
            frame.to_csv(f, index=False)
        return dict(kind="table", title=title, caption=caption, csv=f"tables/{key}.csv",
                    html=render.df_to_html(frame, table_id=key, formats=formats or {}, max_rows=max_rows,
                                           row_class=row_class))

    def F(key, title, how):
        return dict(kind="figure", id=f"fig-{key}", title=title, how_to_read=how, svg=svg[key],
                    png=f"figures/{key}.png")

    num2 = "{:+.2f}".format
    sections = []
    # ---------------- A coverage
    cov = tables["coverage"]
    fl = tables["detection_flags"]
    fl_view = fl.rename(columns={"tissue": "Tissue", "gene": "Gene", "rna_columns": "RNA columns measured",
                                 "rna_baseline": "RNA baseline", "rna_baseline_tertile": "RNA baseline tertile",
                                 "prot_columns": "Protein columns measured", "prot_n_missing": "Proteomics missing values",
                                 "gtex_tpm": "GTEx v8 TPM", "gtex_tertile": "GTEx tertile",
                                 "n_collapsed_rna": "RNA features collapsed", "n_collapsed_prot": "Protein features collapsed",
                                 "baseline_unit": "Baseline unit"})
    unm_view = unm.rename(columns={"status": "Status"})
    sections.append(dict(id="coverage", title="A. Coverage and detection power", intro=[
        f"{n_map} of {n_in} rows mapped. Mapping routes: " + "; ".join(
            f"{k} {v}" for k, v in mp[mp.status == 'mapped'].map_method.value_counts().items()) + ".",
        "Detection-power flags are shown because in these data 'which layer responded' is largely which assay could "
        "see the gene (model/discordance_report.pdf): low baseline mRNA or proteomics missing values make a "
        "non-response uninformative."], items=[
        F("coverage", "Where each signature gene is measured",
          "Cell = fraction of the block's exercise comparisons in which the gene is measured (white = never). "
          "No test; this is coverage only."),
        T("coverage_by_block", "Columns measured per gene and block", "Entries are measured / available comparisons.", cov,
          max_rows=60),
        T("detection_flags", "Detection-power flags per gene and tissue",
          "RNA baseline: human log2 CPM (limma AveExpr), rat log2(sedentary mean count + 1); tertiles among all "
          "genes measured in the tissue. Proteomics missing values: rat numNAs (human not published). GTEx v8 "
          "median TPM of the matched human tissue.", fl_view, formats={"RNA baseline": "{:.1f}", "GTEx v8 TPM": "{:.1f}"},
          max_rows=80),
        T("unmapped", "Rows not used", "Every input row that was not used, with the reason.",
          unm_view if len(unm_view) else pd.DataFrame({"Status": ["all rows mapped"]}))]))

    # ---------------- B grid
    meff = R["meff"]
    sv = sc[["comparison", "n_measured", "n_opposed", "n_same", "sign_p", "sign_p_eff", "n_sig", "n_sig_opposed"]]
    sv = sv.rename(columns={"comparison": "Comparison", "n_measured": "Genes measured", "n_opposed": "Opposed",
                            "n_same": "Same direction", "sign_p": "Sign test p", "sign_p_eff": "Sign test p (effective n)",
                            "n_sig": f"Genes at BH FDR < {cut}", "n_sig_opposed": f"Opposed at BH FDR < {cut}"})
    sections.append(dict(id="grid", title="B. Agreement grid (gene by comparison)", intro=[
        f"Genes are correlated across comparisons (mean pairwise Spearman r = {meff['mean_r']:.2f} over "
        f"{meff['n_genes']} genes), so the effective number of independent genes is {meff['meff_nyholt']:.1f} (Nyholt) or "
        f"{meff['meff_liji']:.1f} (Li–Ji) of {meff['n_genes']}. Sign-test p values treat genes as independent and are "
        f"optimistic; the effective-n column rescales n and k by {meff['meff_nyholt']:.1f}/{meff['n_genes']}. cameraPR "
        "(section C) is the primary test." if meff["n_genes"] >= 2 else
        "Too few genes are measured across comparisons to estimate an effective number of genes."], items=[
        F("grid", "Signature genes against each exercise comparison",
          f"Cell = sign(exercise logFC) × disease direction × {'min(|stat|, ' + format(o['cap'], 'g') + ')' if o['cap'] else '|stat| (no cap)'}; blue = exercise moves the gene "
          f"against the disease, red = same way; dot = BH FDR < {cut} (fdr_bh, within comparison); hatched = not "
          "measured. Rows grouped by signature group."),
        dict(kind="interactive", id="grid-interactive", title="Interactive grid: filter, hover, click a gene",
             how_to_read="Same cells as the grid above (colour = agreement × |stat|, blue = opposed). Use the Tissue "
                         "and Layer menus to filter; hover a cell for logFC, stat and fdr_bh; click a gene to draw "
                         "its statistic across timepoints in both layers below (line style = group or sex).",
             spec_json=interactive.grid_spec(S, R["grid"], R["core_cols"]), fallback_svg=svg["grid"]),
        T("grid_summary", "Per comparison: genes opposed and sign test",
          "Two-sided binomial sign test of opposed vs 0.5 (descriptive; optimistic for co-regulated genes).", sv,
          formats={"Sign test p": fp, "Sign test p (effective n)": fp})]))

    # ---------------- C cameraPR
    cv = sc[["comparison", "n_measured", "camera_t", "camera_p", "camera_fdr", "n_up", "camera_t_up", "camera_fdr_up",
             "n_down", "camera_t_down", "camera_fdr_down", "verdict"]].rename(columns={
                 "comparison": "Comparison", "n_measured": "Genes", "camera_t": "t (all)", "camera_p": "p (all)",
                 "camera_fdr": "FDR (all)", "n_up": "Up genes", "camera_t_up": "t (up-in-disease half)",
                 "camera_fdr_up": "FDR (up half)", "n_down": "Down genes", "camera_t_down": "t (down-in-disease half)",
                 "camera_fdr_down": "FDR (down half)", "verdict": "Set-level verdict"})
    sections.append(dict(id="camera", title="C. Set-level opposition (cameraPR, primary evidence)", intro=[
        "cameraPR (Python port of limma, inter-gene correlation 0.01, validated against limma) tests whether the "
        "signature's statistics are shifted relative to all other genes in the same comparison. Member statistics "
        "are multiplied by −direction, so a positive t means exercise opposes the signature. The up-in-disease and "
        "down-in-disease halves are tested separately as well."], items=[
        F("camera", "Set-level opposition per comparison (cameraPR t)",
          f"Circle = whole signature (filled = BH FDR < {cut} across the comparisons shown); triangles = up-in-disease "
          "(▲) and down-in-disease (▼) halves; dashed lines = |t| = 1.96. Above zero = exercise opposes."),
        T("camera_pr", "cameraPR per comparison", "t, two-sided p and BH FDR across the comparisons shown.", cv,
          formats={"t (all)": num2, "p (all)": fp, "FDR (all)": fp, "t (up-in-disease half)": num2,
                   "FDR (up half)": fp, "t (down-in-disease half)": num2, "FDR (down half)": fp})]))

    # ---------------- D nulls
    nv = sc[["comparison", "camera_t", "pct_abund", "stratum_abund", "pct_class"]].copy()
    nb = R["nulls"]
    for kind, lab in [("abundance", "abund"), ("class", "class")]:
        d = nb[nb.null == kind].set_index("column_id")
        nv[f"median_{lab}"] = sc.column_id.map(d.null_t_median).to_numpy()
    nv = nv[["comparison", "camera_t", "median_abund", "pct_abund", "stratum_abund", "median_class", "pct_class"]]
    nv = nv.rename(columns={"comparison": "Comparison", "camera_t": "Observed t", "median_abund": "Abundance null median t",
                            "pct_abund": "Abundance null percentile", "stratum_abund": "Abundance measure",
                            "median_class": "Class null median t", "pct_class": "Class null percentile"})
    cls = S.ann.reindex(sig.genes)[nulls.CLASS_FLAGS].astype("boolean").fillna(False)
    comp = (f"Signature composition (counted genes): {int(cls.mitocarta3.sum())} MitoCarta3.0, "
            f"{int(cls.go_any_complex.sum())} GO:CC complex subunits, {int(cls.secreted_proxy.sum())} secreted/"
            f"extracellular (GO:CC proxy), of {len(sig.genes)}.")
    items = [F("nulls", "Is the opposition more than expected for genes of this abundance, or of this pathway class?",
               f"Black dot = observed signed cameraPR t; grey bar = 5–95% of {B} abundance-matched random sets; "
               f"blue bar = 5–95% of {B} sets matched on MitoCarta / GO:CC complex / secreted class; tick = null median. "
               "A dot inside the blue bar means the pathway class explains the result."),
             T("null_percentiles", "Null percentiles per comparison",
               f"Percentile = 100 × P(null t ≤ observed t), {B} draws each, drawn gene by gene within the stratum "
               "(with replacement, as in NARRATIVE R1b).", nv,
               formats={"Observed t": num2, "Abundance null median t": num2, "Class null median t": num2,
                        "Abundance null percentile": "{:.0f}", "Class null percentile": "{:.0f}"})]
    if "pool" in R:
        pv = R["pool"].assign(Comparison=R["pool"].column_id.map(cw))[["Comparison", "obs_t", "null_t_median", "pct_t",
                                                                        "pool_size"]]
        pv = pv.rename(columns={"obs_t": "Observed t (pool genes)", "null_t_median": "Pool null median t",
                                "pct_t": "Pool null percentile", "pool_size": "Pool size"})
        items.append(T("pool_null", f"Fixed-pool null for {len(R['pool_genes'])} genes ({', '.join(o['pool_sets'])})",
                       "NARRATIVE R1a: random sets of the same size from the pool genes measured in both layers of "
                       "the anchor tissue.", pv, formats={"Observed t (pool genes)": num2, "Pool null median t": num2,
                                                          "Pool null percentile": "{:.1f}"}))
    sections.append(dict(id="nulls", title="D. Calibration nulls (specificity)", intro=[comp], items=items))

    # ---------------- E layers
    if R.get("metab_mode") or not len(tables["layer_pairs"]):
        sections.append(dict(id="layers", title="E. Layer discordance", intro=[
            "Not applicable: RNA-versus-protein pairing needs gene-level signatures measured in both layers."],
            items=[]))
    else:
        d = tables["layer_pairs"]
        dsent = [layers.sentence(r) for r in d.itertuples() if not r.cross and r.tissue in ("VL", "SKM-GN", "HEART")]
        dv = d[["label", "n_genes_both", "rho_all", "n_sig_both", "rna_opposed", "prot_opposed", "both_opposed", "rho_sig",
                "class_concordant", "class_RNA_only", "class_PROT_only", "class_opposite", "class_ns"]].rename(columns={
                    "label": "Comparison", "n_genes_both": "Genes in both layers", "rho_all": "Genome-wide ρ",
                    "n_sig_both": "Signature genes in both", "rna_opposed": "RNA opposed", "prot_opposed": "Protein opposed",
                    "both_opposed": "Both opposed", "rho_sig": "Signature ρ", "class_concordant": "Concordant",
                    "class_RNA_only": "RNA only", "class_PROT_only": "Protein only", "class_opposite": "Opposite",
                    "class_ns": "Neither"})
        lg = tables["layer_genes"]
        if len(lg):
            fl2 = R["flags"].copy()
            fl2["tissue_code"] = fl2.tissue.map({"Human vastus lateralis": "VL", "Rat gastrocnemius": "SKM-GN",
                                                 "Rat vastus lateralis": "SKM-VL", "Rat heart": "HEART"})
            lg = lg[~lg.cross].merge(pd.DataFrame({"prot": lg.prot.unique()}).assign(
                tissue_code=lambda x: x.prot.map(S.cols.tissue)), on="prot")
            lg = lg.merge(fl2, left_on=["gene_symbol_human", "tissue_code"], right_on=["gene", "tissue_code"], how="left")
            lgv = lg[["pair", "gene_symbol_human", "direction", "stat_rna", "stat_prot", "layer_class", "rna_baseline_tertile",
                      "prot_n_missing", "gtex_tertile", "n_collapsed_rna", "n_collapsed_prot"]].rename(columns={
                          "pair": "Comparison", "gene_symbol_human": "Gene", "direction": "Disease direction",
                          "stat_rna": "RNA stat", "stat_prot": "Protein stat", "layer_class": "Which layer answered",
                          "rna_baseline_tertile": "RNA baseline tertile", "prot_n_missing": "Proteomics missing values",
                          "gtex_tertile": "GTEx tertile", "n_collapsed_rna": "RNA features", "n_collapsed_prot": "Protein features"})
            lgv["Which layer answered"] = lgv["Which layer answered"].map({"concordant": "both, same sign",
                                                                             "opposite": "both, opposite sign",
                                                                             "RNA_only": "RNA only", "PROT_only": "protein only",
                                                                             "ns": "neither"})
        else:
            lgv = pd.DataFrame({"Comparison": ["no signature gene measured in both layers"]})
        sections.append(dict(id="layers", title="E. Layer discordance: which layer answered, and could it have?",
                             intro=dsent[:14] + (["…"] if len(dsent) > 14 else []), items=[
            F("layers", "RNA versus protein for the signature, by timepoint",
              "Left: Spearman ρ between RNA and protein statistics (Fisher-z 95% CI; scripts/06), all genes measured in "
              "both layers (grey) vs signature genes (blue). Right: signature genes, both axes signed so that + = opposed; "
              "blue = opposed in both layers, grey = one layer, red = neither."),
            T("layer_pairs", "Paired RNA / protein summary per comparison",
              f"Layer classes use BH FDR < {cut} in each layer (scripts/08 rule).", dv,
              formats={"Genome-wide ρ": "{:.2f}", "Signature ρ": "{:.2f}"}),
            dict(kind="details", summary="Per-gene layer classes with detection-power flags", items=[
                T("layer_genes_view", "Which layer answered, per gene, with detection flags",
                  "A gene that is 'RNA only' with high proteomics missingness, or 'protein only' with low baseline mRNA, "
                  "is a detection question before it is a biological one.", lgv,
                  formats={"RNA stat": num2, "Protein stat": num2}, max_rows=300)])]))

    # ---------------- 3-layer pathways
    pitems = []
    for (tis, sp), pw in R.get("pathways", {}).items():
        if f"pathways_{tis}" not in svg:
            continue
        from .metab import PATHWAY_CAPTION
        tw = f"{'Human' if sp == 'human' else 'Rat'} {core.TISSUE_WORDS.get(tis, tis)}"
        pitems.append(F(f"pathways_{tis}", f"{tw}: metabolic pathways in RNA, protein and metabolites",
                        PATHWAY_CAPTION))
        pitems.append(dict(kind="details", summary=f"{tw}: pathway cameraPR table", items=[
            T(f"pathways_{tis}", "Pathway cameraPR t per layer and time", "t > 0 = set up with exercise.",
              tables[f"pathways_{tis}"], max_rows=200)]))
    sections.append(dict(id="pathways", title="Metabolic pathways across three layers (context, not the signature)",
                         intro=["The same pathways scored in RNA, protein and metabolites at matching timepoints. "
                                "Metabolite pools are not flux; the gene-metabolite map is many-to-many "
                                "(store/pathway_map.csv) and genes and metabolites are never matched one-to-one."]
                         if pitems else [R.get("pathways_error", "Metabolomics layer not built.")], items=pitems))

    # ---------------- F everywhere
    es, bk = R["ev_summary"], R["ev_bykind"]
    rk = R["ev_rank"]
    view = rk[rk.ranked]
    top = lambda col: view.sort_values([col, "n_measured", "sign_p"], ascending=[False, False, True]).head(20)[[
        "unit_words", "kind", "n_measured", "n_opposed", "n_same", col, "sign_p", "sign_q_bh"]].rename(columns={
            "unit_words": "Unit", "kind": "Comparison kind", "n_measured": "Cells measured", "n_opposed": "Opposed",
            "n_same": "Same direction", col: "Fraction " + ("opposed" if col == "frac_opposed" else "same"),
            "sign_p": "Sign test p", "sign_q_bh": "BH q"})
    bkv = bk.rename(columns={"kind": "Comparison kind", "units": "Units", "median_frac_opposed": "Median fraction opposed",
                             "q_opposed": "Opposed-majority units, q < 0.05", "q_same": "Same-majority units, q < 0.05"})
    sections.append(dict(id="everywhere", title="F. Everywhere: every tissue, layer and comparison", intro=[
        f"{es['n_columns']} comparison columns grouped into {es['n_units']} units (rat sexes pooled per week); "
        f"{es['n_ranked']} units with at least {es['min_n']} measured cells are ranked. Nominal sign test p < 0.05: "
        f"{es['nominal_opposed']} opposed-majority and {es['nominal_same']} same-majority units (about "
        f"{es['expected_by_chance']:.0f} expected by chance). After BH: {es['q_opposed']} opposed-majority and "
        f"{es['q_same']} same-majority units at q < 0.05.",
        "Non-exercise reference comparisons (control time course, baseline group differences, EE vs RE) are a "
        "negative control: they show how often the signature looks 'opposed' when no exercise contrast is involved."],
        items=[F("everywhere", "Opposition fraction in every ranked unit, by comparison kind",
                 "Point = one unit (size ∝ cells measured; filled = BH q < 0.05, two-sided sign test); black bar = "
                 "median. Grey = non-exercise reference comparisons. Dashed line = 0.5."),
               T("everywhere_kind", "By comparison kind", "Medians and BH-significant units per kind.", bkv,
                 formats={"Median fraction opposed": "{:.2f}"}),
               T("everywhere_opp", "Top 20 units by fraction opposed", "Sign test is descriptive (genes correlated).",
                 top("frac_opposed"), formats={"Fraction opposed": "{:.2f}", "Sign test p": fp, "BH q": fp}),
               T("everywhere_same", "Top 20 units by fraction same-direction", "", top("frac_same"),
                 formats={"Fraction same": "{:.2f}", "Sign test p": fp, "BH q": fp})]))

    # ---------------- G caveats
    early = sc[sc.early_flag]
    ecap = "; ".join(f"{r.comparison}: t = {r.camera_t:+.2f}" for r in early.itertuples() if np.isfinite(r.camera_t))
    myh = tables["fibre_markers"]
    mv = myh.assign(Comparison=myh.column_id.map(cw))[["Comparison", "MYH7", "MYH2", "MYH1", "MYH4"]].rename(columns={
        "MYH7": "MYH7 (type I)", "MYH2": "MYH2 (IIa)", "MYH1": "MYH1 (IIx)", "MYH4": "MYH4 (IIb)"})
    n_contr = int(S.ann.reindex(sig.genes).go_contractile_fiber.astype("boolean").fillna(False).sum())
    cav = fixed_caveats() + [
        f"Human 15–45 min comparisons are flagged: a shared RNA + protein drop of myofibre programs with a rise of "
        f"blood-derived programs (NARRATIVE R2) means these columns partly measure biopsy composition. This "
        f"signature there: {ecap or 'not measured'}.",
        f"{n_contr} of {len(sig.genes)} counted genes are GO:CC contractile-fibre genes. Rat training shifts fibre type "
        f"(IIb → IIa/IIx; NARRATIVE F3), so contractile genes can move with composition; the fibre-type marker "
        f"statistics are listed below.",
        f"Sign tests are optimistic: effective number of genes {meff['meff_nyholt']:.1f} of {meff['n_genes']} "
        f"(Nyholt). cameraPR corrects for inter-gene correlation and is the primary test."]
    gv = pd.DataFrame(GUARDRAILS, columns=["Safe statement", "Unsafe upgrade"])
    sections.append(dict(id="caveats", title="G. Caveats and interpretation guardrails", intro=cav, items=[
        T("guardrails", "Safe statements and the unsafe upgrades to avoid",
          "From the team's planning workbook ('Start Here', section 4), generalised to any disease signature.", gv),
        T("fibre_markers", "Fibre-type marker statistics per comparison",
          "Gene-level stat (moderated t / DESeq2 z / limma t) of the myosin heavy-chain markers; read a contractile "
          "signature against these.", mv, formats={c: num2 for c in mv.columns[1:]})]))

    # ---------------- sensitivity
    sv2 = tables["sensitivity"]
    flips = sv2[sv2.flipped]
    sens_v = flips.assign(Comparison=flips.column_id.map(cw))[["Comparison", "setting", "conclusion", "base",
                                                               "alternative"]].rename(columns={
        "setting": "Setting changed", "conclusion": "Conclusion", "base": "Default", "alternative": "Alternative"})
    ssum = sv2.groupby("setting").flipped.agg(["sum", "size"]).reset_index().rename(columns={
        "setting": "Setting changed", "sum": "Conclusions that flip", "size": "Conclusions checked"})
    sections.append(dict(id="sensitivity", title="Sensitivity of the conclusions", intro=[
        f"{int(sv2.flipped.sum())} of {len(sv2)} conclusion checks flip under an alternative setting. The stat cap "
        f"({format(o['cap'], 'g') if o['cap'] else 'none'}) only changes grid colours, never a statistic."], items=[
        T("sensitivity_summary", "Conclusions checked per setting", "", ssum),
        T("sensitivity_flips", "Conclusions that flip",
          "Listed conclusions should be reported with the setting that supports them.",
          sens_v if len(sens_v) else pd.DataFrame({"Result": ["no conclusion flips"]}), max_rows=100)]))

    if R.get("ranked_mode"):
        sections[1:1] = ranked_sections(R, tables, cw, fp, T, F, num2, cut)
    titems = [F(k, t, "Line = exercise effect per gene (log2 fold change vs sedentary rats / resting controls), band = "
                      "95% CI from the published test; blue RNA, orange protein; filled point = BH FDR < "
                      f"{cut}; dashed line = no change. Rat: female and male pooled by inverse-variance weighting. "
                      "Up to 12 genes (largest effects); every gene is in the CSV.")
              for k, t in [("traj_rat", "Each signature gene over 1–8 weeks of training (rat gastrocnemius)"),
                           ("traj_human", "Each signature gene after one bout (human vastus lateralis; EE solid, "
                                          "RE dashed)")] if k in svg]
    if titems:
        sections.insert(0, dict(id="timecourse", title="Time course of each signature gene, RNA and protein",
                                intro=["The team's rat training-course plot, generalised to any signature, both "
                                       "layers and both species."], items=titems))

    meta = [("Signature", sig.name), ("Genes counted", str(len(sig.genes))), ("Rows in file", str(n_in)),
            ("Run", dt.datetime.now().strftime("%Y-%m-%d %H:%M")), ("Store", store.store_hash()),
            ("Repo", store.git_sha()[:10])]
    glossary = GLOSSARY
    return dict(title=f"MoTrPAC probe: {sig.name}",
                subtitle="How much of this signature does the exercise response oppose, by omic layer and time?",
                meta=meta, headline=headline, first_screen=first, sections=sections, glossary=glossary,
                toggles={s: True for s in SECTIONS} | ({k: False for k in SECTIONS if k not in o["sections"]}
                                                       if o.get("sections") else {}))


TEAM_CREDIT = ("Generalised from the team's PAH blood module (`MoTrPAC Hackathon/scripts` 07.5–12 on main: GSE33463 "
               "IPAH-minus-healthy PBMC ranking vs MoTrPAC blood RNA).")


def ranked_sections(R, tables, cw, fp, T, F, num2, cut):
    """Sections for a full disease ranking: whole-ranking association and GO:BP pathway concordance."""
    S, ri = R["S"], R["rinfo"]
    ra = R["rank_focus"].copy()
    ra["Comparison"] = ra.column_id.map(cw)
    v = ra[["Comparison", "n_shared", "rho", "ci_low", "ci_high", "perm_p", "perm_fdr", "direction",
            "ref_abs_rho_median", "pct_vs_reference"]].rename(columns={
        "n_shared": "Shared genes", "rho": "Spearman ρ", "ci_low": "95% CI low", "ci_high": "95% CI high",
        "perm_p": "Permutation p", "perm_fdr": "Permutation BH q", "direction": "Direction",
        "ref_abs_rho_median": "Median |ρ|, reference contrasts", "pct_vs_reference": "|ρ| percentile vs reference"})
    ok = ra.dropna(subset=["rho"])
    lo, hi = ok.loc[ok.rho.idxmin()], ok.loc[ok.rho.idxmax()]
    n_perm = int(ra.perm_p.notna().sum())
    sents = [f"Disease ranking: {ri['n_full']:,} genes ({ri['stat_column']}), {ri['n_ranked']:,} mapped to MoTrPAC. "
             f"Derived directional signature for sections B–G: {ri['n_derived_up']} up and {ri['n_derived_down']} "
             f"down ({ri['derived_rule']}).",
             f"Most opposed comparison: {lo.Comparison} (ρ = {lo.rho:+.3f}, 95% CI {lo.ci_low:+.3f} to "
             f"{lo.ci_high:+.3f}); most similar: {hi.Comparison} (ρ = {hi.rho:+.3f}). |ρ| below about 0.1 means the "
             "exercise response shares little of the disease ranking, however small its p value."]
    sec_r = dict(id="rank", title="R. Whole-ranking association (disease ranking vs each exercise comparison)",
                 intro=[TEAM_CREDIT] + sents, items=[
        F("rank_assoc", "Does the whole disease ranking resemble or oppose each exercise comparison?",
          "Spearman ρ between the disease statistic and the exercise statistic over shared genes, with Fisher-z 95% "
          "CI (+ = exercise moves genes the same way as the disease, − = opposed). Grey points = the same ranking "
          "against non-exercise reference contrasts of that tissue and layer (the calibration)."),
        T("rank_association", "Rank association per comparison",
          f"Permutation p = two-sided gene-label shuffle ({R['opts']['nperm']:,} shuffles; computed for the {n_perm} "
          "human exercise-vs-control RNA comparisons) — it ignores gene-gene correlation, so it is exploratory, "
          "as in the team pipeline. The reference percentile compares |ρ| with the non-exercise contrasts.", v,
          formats={"Spearman ρ": "{:+.3f}", "95% CI low": "{:+.3f}", "95% CI high": "{:+.3f}", "Permutation p": fp,
                   "Permutation BH q": fp, "Median |ρ|, reference contrasts": "{:.3f}",
                   "|ρ| percentile vs reference": "{:.0f}"}, max_rows=80)])
    ps, pk, dis = R["pw_store"].copy(), R["pw_pkg"].copy(), R["pw_dis"]
    for d in (ps, pk):
        d["Comparison"] = d.column_id.map(cw)
    cols = {"shared_sets": "Shared sets", "both_bh": "BH < 0.05 on both sides", "same_direction": "Same direction",
            "opposite_direction": "Opposite direction", "rho_set_t": "Spearman ρ of set t"}
    lg = R["pw_store_long"]
    top = lg[lg.both_bh].copy() if len(lg) else lg
    if len(top):
        top["Comparison"] = top.column_id.map(cw)
        top = top.sort_values("fdr_disease")[["Comparison", "set", "t_disease", "fdr_disease", "t_exercise",
                                              "fdr_exercise", "concordance"]].rename(columns={
            "set": "GO:BP set", "t_disease": "t, disease", "fdr_disease": "BH q, disease",
            "t_exercise": "t, exercise", "fdr_exercise": "BH q, exercise", "concordance": "Direction"})
    else:
        top = pd.DataFrame({"Result": ["no GO:BP set is BH < 0.05 on both sides"]})
    sec_p = dict(id="pathways_go", title="P. Pathway concordance (GO:BP, both sides tested)", intro=[
        TEAM_CREDIT,
        f"Disease side: {len(dis):,} GO:BP sets (MotrpacHumanPreSuspensionAnalysis MOLECULAR_SIGNATURES, 10–500 genes "
        f"in the ranking) tested with cameraPR; {int((dis.fdr < 0.05).sum())} at BH < 0.05. Exercise side: the same "
        "test on every comparison (this table), and, for human comparisons, MoTrPAC's own precomputed results "
        "(second table; this is the team pipeline).",
        f"Across the comparisons shown: {int(ps.same_direction.sum())} same-direction and "
        f"{int(ps.opposite_direction.sum())} opposite-direction set × comparison rows (same test on both sides); "
        f"with MoTrPAC's precomputed results for human comparisons: {int(pk.same_direction.sum())} and "
        f"{int(pk.opposite_direction.sum())}. GO sets share genes, so rows are not independent mechanisms."],
        items=[F("pw_counts", "Pathways significant on both sides, by comparison",
                 "Bars: GO:BP sets with BH < 0.05 in both the disease ranking and the exercise comparison; red = same "
                 "direction as the disease, blue = opposite (exercise opposes). cameraPR, inter-gene correlation 0.01."),
               T("pw_store", "Per comparison (same cameraPR on both sides)", "", ps[["Comparison"] + list(cols)].rename(
                   columns=cols), formats={"Spearman ρ of set t": "{:+.3f}"}, max_rows=80),
               T("pw_pkg", "Per human comparison, MoTrPAC precomputed pathway results (team pipeline)",
                 "Exercise side = MotrpacHumanPreSuspensionAnalysis CAMERA_RESULTS (GOBP).",
                 pk[["Comparison"] + list(cols)].rename(columns=cols), formats={"Spearman ρ of set t": "{:+.3f}"}),
               dict(kind="details", summary="Jointly significant sets", items=[
                   T("pw_joint", "GO:BP sets significant on both sides", "Same cameraPR on both sides.", top,
                     formats={"t, disease": num2, "t, exercise": num2, "BH q, disease": fp, "BH q, exercise": fp},
                     max_rows=200)])])
    return [sec_r, sec_p]


GLOSSARY = [
    ("logFC", "log2 fold change of the exercise comparison (exercise minus control, or trained minus sedentary)."),
    ("stat", "gene-level test statistic of the comparison: human limma/dream moderated t; rat RNA DESeq2 Wald z; rat "
             "protein limma t. One feature per gene (max |stat| among features mapping to the gene; scripts/03_join.py)."),
    ("fdr_bh", "Benjamini–Hochberg adjusted p, recomputed within each dataset × tissue × layer × comparison over all "
               "tested features (the rat package publishes a more conservative cross-tissue BY)."),
    ("agreement", "sign(exercise logFC) × disease direction. −1 = exercise moves the gene against the disease "
                  "('opposed'); +1 = same direction."),
    ("cameraPR", "limma's pre-ranked competitive gene-set test (Wu & Smyth 2012), here a validated Python port with "
                 "inter-gene correlation 0.01; t > 0 = signature opposed."),
    ("Abundance null", "random gene sets matched gene-by-gene on abundance decile in the same comparison."),
    ("Class null", "random gene sets matched on MitoCarta3.0 / GO:CC complex / secreted composition."),
    ("Effective n", "number of independent genes implied by their correlation across comparisons (Nyholt 2004; "
                    "Li & Ji 2005)."),
    ("EE, RE, CON", "human endurance exercise, resistance exercise, and resting control groups."),
    ("EE–CON 24 h", "human delta-delta contrast: (EE at 24 h − EE pre) − (CON at 24 h − CON pre)."),
    ("F 8 wk / M 8 wk", "rat female / male, trained for 8 weeks vs sex-matched sedentary controls."),
    ("VL, gastroc.", "vastus lateralis; gastrocnemius (rat SKM-GN)."),
    ("RefMet", "Metabolomics Workbench reference nomenclature for metabolite names."),
]


def main_run(sig_path, name=None, outdir=None, command="", **kw):
    kw["_path"] = str(sig_path)
    R = compute(sig_path, name=name, **kw)
    name = name or R["sig"].name
    out = Path(outdir) if outdir else OUT / name
    p = write(R, out, command=command, toggles=None)
    return p, R
