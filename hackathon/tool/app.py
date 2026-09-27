"""motrpac_probe app: give it a disease gene/protein list, see what exercise does to those genes in MoTrPAC,
layer by layer and over time. A thin layer over the package (no analysis code here).

    pip install "marimo>=0.16"
    marimo run app.py
"""
import marimo

__generated_with = "0.16.5"
app = marimo.App(width="full", app_title="motrpac_probe", css_file="src/motrpac_probe/templates/style.css")


@app.cell
def _():
    import json
    import tempfile
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd

    from motrpac_probe import core, figures, layers, nulls, render, run, signature
    from motrpac_probe.core import LAYER_WORDS, TISSUE_WORDS, col_words
    from motrpac_probe.paths import EXAMPLES, OUT

    S = core.get_store()  # loaded once when the app starts (~4 s)
    mo.Html("<style>.markdown h1,.markdown h2,.markdown h3{font-family:var(--font-body);font-weight:600} "
            "table{font-size:13px} .scrollbox{max-height:560px;overflow:auto}</style>")
    return (EXAMPLES, LAYER_WORDS, OUT, Path, S, TISSUE_WORDS, col_words, core, figures, json, layers, mo, np,
            nulls, pd, render, run, signature, tempfile)


@app.cell
def _(EXAMPLES, json, mo):
    # ---- inputs (sidebar) --------------------------------------------------------------------------------------
    manifest = json.loads((EXAMPLES / "examples.json").read_text(encoding="utf-8"))
    examples = {e["title"]: e for e in manifest if e["file"].endswith(".csv") and e.get("kind") != "metabolite demo"}
    _args = mo.cli_args()
    _default = next((t for t, e in examples.items() if e["name"] == _args.get("example",
                     "pah_muscle_lower9_malenfant2015")), list(examples)[0])
    example = mo.ui.dropdown(options=list(examples), value=_default, label="Example signature")
    paste = mo.ui.text_area(placeholder="or paste your own, one gene per line:\nNDUFA9,-1\nLDHA,1", rows=5,
                            label="Paste genes (gene, +1 up / -1 down in disease)")
    upload = mo.ui.file(filetypes=[".csv"], label="or upload a CSV")
    species = mo.ui.radio(options={"Rat: 1-8 weeks of endurance training": "rat",
                                   "Human: one exercise bout (15 min - 24 h)": "human"},
                          value=("Human: one exercise bout (15 min - 24 h)" if _args.get("species") == "human"
                                 else "Rat: 1-8 weeks of endurance training"), label="Exercise data")
    fdr = mo.ui.slider(0.01, 0.25, step=0.01, value=0.05, label="FDR threshold (BH)", show_value=True)
    return example, examples, fdr, paste, species, upload


@app.cell
def _(S, TISSUE_WORDS, mo, species):
    _ds = "rat_train" if species.value == "rat" else "human_acute"
    _c = S.cols[S.cols.dataset == _ds]
    _both = sorted(set(_c[_c.layer == "RNA"].tissue) & set(_c[_c.layer == "PROT"].tissue)) or sorted(set(_c.tissue))
    _first = [t for t in ["SKM-GN", "VL", "HEART"] if t in _both]
    _opts = {f"{TISSUE_WORDS.get(t, t)} ({t})": t for t in _first + [t for t in _both if t not in _first]}
    tissue = mo.ui.dropdown(options=_opts, value=list(_opts)[0], label="Tissue")
    return (tissue,)


@app.cell
def _(example, fdr, mo, paste, species, tissue, upload):
    mo.sidebar([mo.md("### motrpac_probe"), example, paste, upload, mo.md("---"), species, tissue, fdr,
                mo.md("<small>Public MoTrPAC summary statistics. 'Opposed' = exercise moves the gene against its "
                      "disease direction; a hypothesis, not a treatment effect.</small>")], width="360px")
    return


@app.cell
def _(EXAMPLES, Path, example, examples, paste, signature, tempfile, upload):
    # ---- the signature: pasted text > uploaded file > example --------------------------------------------------
    if paste.value.strip():
        src = Path(tempfile.mkdtemp()) / "pasted.csv"
        _lines = [l.replace("\t", ",").replace(" ", ",").strip(",") for l in paste.value.strip().splitlines() if l.strip()]
        if not _lines[0].lower().startswith("gene"):
            _lines = ["gene_symbol,direction"] + _lines
        src.write_text("\n".join(_lines) + "\n", encoding="utf-8")
        sig_label = "pasted list"
    elif upload.value:
        src = Path(tempfile.mkdtemp()) / upload.name()
        src.write_bytes(upload.contents())
        sig_label = upload.name()
    else:
        src = EXAMPLES / examples[example.value]["file"]
        sig_label = example.value
    sig = signature.load(src)
    dirs = dict(zip(sig.shown, sig.shown_dirs))
    return dirs, sig, sig_label, src


@app.cell
def _(S, core, dirs, fdr, sig, species, tissue):
    # ---- the lookup: each gene's exercise effect over time, per layer --------------------------------------------
    T = core.trajectory_table(S, sig.shown, tissue.value, species=species.value, cutoff=fdr.value)
    lookup = core.lookup_table(T, dirs)
    return T, lookup


@app.cell
def _(T, TISSUE_WORDS, dirs, fdr, mo, np, sig, sig_label, species, tissue):
    # ---- one-line summary, always visible ----------------------------------------------------------------------
    _t = TISSUE_WORDS.get(tissue.value, tissue.value)
    _parts = []
    if len(T):
        _last = T[T.x == T.x.max()]
        _last = _last[_last.series == _last.series.iloc[0]]
        for _l, _w in [("RNA", "RNA"), ("PROT", "protein")]:
            _d = _last[_last.layer == _l]
            if len(_d):
                _opp = (np.sign(_d.logFC) * _d.gene.map(dirs) < 0)
                _parts.append(f"**{_w}**: {int(_opp.sum())} of {len(_d)} genes move against the disease "
                              f"({int((_opp & _d.significant).sum())} at FDR < {fdr.value:.2f})")
    _when = "after 8 weeks of training" if species.value == "rat" else "24 h after endurance exercise"
    mo.md(f"### {sig_label}\n{len(sig.genes)} genes mapped of {len(sig.table)} rows. "
          f"{('In ' + ('rat' if species.value == 'rat' else 'human') + ' ' + _t + ' ' + _when + ': ' + '; '.join(_parts) + '.') if _parts else 'No signature gene is measured in this tissue.'}")
    return


@app.cell
def _(T, dirs, figures, lookup, mo, render, sig, species, tissue):
    # ---- tab 1: time course (the team's rat plot, generalised) ------------------------------------------------
    if len(T):
        _genes = [g for g in sig.shown if g in set(T.gene)]
        if len(_genes) > 12:  # show the 12 genes with the largest effects; all are in the table
            _mx = T.groupby("gene").logFC.apply(lambda v: v.abs().max())
            _genes = [g for g in _genes if g in set(_mx.nlargest(12).index)]
        with render.mpl_style():
            _svg = render.fig_to_svg(figures.trajectories(T, _genes, dirs, tissue.value, species=species.value))
        tab_time = mo.vstack([
            mo.Html(_svg),
            mo.md("<small>Line = exercise effect (log2 fold change vs sedentary or control), band = 95% CI; blue RNA, "
                  "orange protein; filled point = significant at the chosen FDR; dashed line = no change. Rat: female "
                  "and male pooled. Up to 12 genes shown (largest effects); the table has every gene.</small>"),
            mo.accordion({"Table: every gene, every time point (* = significant)": mo.Html(
                f'<div class="scrollbox">{render.df_to_html(lookup, table_id="lookup")}</div>')})])
    else:
        tab_time = mo.md("None of these genes is measured in this tissue.")
    return (tab_time,)


@app.cell
def _(S, T, core, dirs, fdr, layers, mo, pd, render, species, tissue):
    # ---- tab 2: do RNA and protein agree? ---------------------------------------------------------------------
    _times = list(dict.fromkeys(T.sort_values(["series", "x"]).apply(lambda r: (r.series, r.time), axis=1))) if len(T) else []
    _rows = []
    for _s, _tm in _times:
        _d = T[(T.series == _s) & (T.time == _tm)].pivot_table(index="gene", columns="layer",
                                                                  values=["logFC", "significant"], aggfunc="first")
        if ("logFC", "RNA") not in _d or ("logFC", "PROT") not in _d:
            continue
        for _g, _r in _d.iterrows():
            _a, _b = _r[("logFC", "RNA")], _r[("logFC", "PROT")]
            if pd.isna(_a) or pd.isna(_b):
                continue
            _sa, _sb = bool(_r[("significant", "RNA")]), bool(_r[("significant", "PROT")])
            _cls = ("both, same sign" if _sa and _sb and _a * _b > 0 else "both, opposite sign" if _sa and _sb
                    else "RNA only" if _sa else "protein only" if _sb else "neither")
            _rows.append({"Time": f"{_s.replace('trained, F+M pooled', '')} {_tm}".strip(), "Gene": _g,
                          "RNA log2FC": _a, "Protein log2FC": _b, "Which layer responded": _cls,
                          "Protein vs disease": "opposed" if _b * dirs.get(_g, 0) < 0 else "same"})
    _agree = pd.DataFrame(_rows)
    _counts = (_agree.groupby(["Time", "Which layer responded"]).size().unstack(fill_value=0).reset_index()
               if len(_agree) else pd.DataFrame())
    _short = {"both, same sign": "both", "both, opposite sign": "both, opposite", "RNA only": "RNA only",
              "protein only": "protein only", "neither": "–"}
    _wide = (_agree.assign(c=_agree["Which layer responded"].map(_short))
             .pivot_table(index="Gene", columns="Time", values="c", aggfunc="first")
             .reindex(columns=list(dict.fromkeys(_agree.Time))).reset_index()) if len(_agree) else pd.DataFrame()
    if len(_wide) and tissue.value in ("VL", "SKM-GN", "SKM-VL", "HEART"):
        _f = layers.detection_flags(S, sorted(set(_agree.Gene)))
        _want = {"VL": "Human vastus lateralis", "SKM-GN": "Rat gastrocnemius", "SKM-VL": "Rat vastus lateralis",
                 "HEART": "Rat heart"}[tissue.value]
        _f = _f[_f.tissue == _want][["gene", "rna_baseline_tertile", "prot_n_missing", "gtex_tertile"]]
        _f.columns = ["Gene", "mRNA level", "Protein missing values", "GTEx abundance"]
        _wide = _wide.merge(_f, on="Gene", how="left")
    tab_layers = mo.vstack([
        mo.md(f"**Which layer responded?** Per gene and time point: significant (BH FDR < {fdr.value:.2f}) in RNA, "
              "in protein, in both, or neither. The right-hand columns say whether each assay could see the gene "
              "well: a low mRNA level or missing protein values make 'no response' uninformative."),
        mo.Html(render.df_to_html(_wide, table_id="agree_wide", formats={"Protein missing values": "{:.0f}"}))
        if len(_wide) else mo.md("No gene is measured in both RNA and protein here."),
        mo.accordion({"Counts per time point": mo.Html(render.df_to_html(_counts, table_id="layer_counts")),
                      "Every gene and time point with log2 fold changes": mo.Html(
                          f'<div class="scrollbox">{render.df_to_html(_agree, table_id="agree", formats={"RNA log2FC": "{:+.2f}", "Protein log2FC": "{:+.2f}"})}</div>')})
        if len(_agree) else mo.md("")])
    return (tab_layers,)


@app.cell
def _(S, core, fdr, mo, nulls, render, sig, species, tissue):
    # ---- tab 3: is it specific? set-level test + pathway-class null ---------------------------------------------
    _c = S.cols
    _ds = "rat_train" if species.value == "rat" else "human_acute"
    _cids = _c[(_c.dataset == _ds) & (_c.tissue == tissue.value) & _c.layer.isin(["RNA", "PROT"]) &
               _c.kind.isin(core.EXERCISE_KINDS)].sort_values("col_order").column_id.tolist()
    _sc = core.score_columns(S, _cids, sig, cutoff=fdr.value)
    _nl = nulls.run_nulls(S, _cids, sig, B=200, kinds=("class",)).set_index("column_id")
    _info = S.cols.loc[_sc.column_id]
    _sc["Row"] = [f"{'RNA' if r.layer == 'RNA' else 'protein'}, " +
                  ({"F": "female", "M": "male"}.get(r.sex, "") if r.dataset == "rat_train" else
                   {"EE": "endurance", "RE": "resistance"}.get(r.group, r.group)) for r in _info.itertuples()]
    _sc["Time"] = [(r.time[:-1] + " wk") if r.dataset == "rat_train" else
                   {"15-45m": "15–45 min", "3.5h": "3.5–4 h", "24h": "24 h"}.get(r.time, r.time) for r in _info.itertuples()]
    _pct = _sc.column_id.map(_nl.pct_t)
    _sc["cell"] = [("–" if t != t else f"{t:+.1f}{'*' if q < fdr.value else ''} ({p:.0f})")
                   for t, q, p in zip(_sc.camera_t, _sc.camera_fdr, _pct)]
    _m = _sc.pivot_table(index="Row", columns="Time", values="cell", aggfunc="first")
    _m = _m.reindex(columns=list(dict.fromkeys(_sc.Time)), index=list(dict.fromkeys(_sc.Row))).reset_index()
    _n = int(((_sc.camera_fdr < fdr.value) & (_sc.camera_t > 0)).sum())
    tab_specific = mo.vstack([
        mo.md("**Does the whole gene set shift against the disease, and is that more than any gene set of the same "
              "kind would show?** Each cell: cameraPR t (> 0 = opposed; * = BH FDR < "
              f"{fdr.value:.2f}) and, in brackets, the class-null percentile: the signature against 200 random sets "
              "with the same mitochondrial / complex / secreted make-up and abundance. Near 100 = specific to this "
              "signature; 30–70 = typical of its gene class."),
        mo.Html(render.df_to_html(_m, table_id="specific")),
        mo.md(f"{_n} of {len(_sc)} comparisons are opposed at the set level.")])
    return (tab_specific,)


@app.cell
def _(mo, render, species, tissue):
    # ---- tab 4: metabolites (pathway level, independent of the signature) -------------------------------------
    try:
        from motrpac_probe import metab
        _pw = metab.pathway_panel(None, tissue=tissue.value, species=species.value)
        with render.mpl_style():
            _svg = render.fig_to_svg(metab.pathway_heatmap(_pw, tissue.value))
        tab_metab = mo.vstack([mo.Html(_svg), mo.md(f"<small>{metab.PATHWAY_CAPTION}</small>")])
    except Exception as _e:
        tab_metab = mo.md(f"No three-layer pathway panel for this tissue ({type(_e).__name__}).")
    return (tab_metab,)


@app.cell
def _(mo, pd, render):
    from motrpac_probe.caveats import GUARDRAILS, fixed_caveats
    export = mo.ui.run_button(label="Write the full report")
    tab_about = mo.vstack([
        mo.md("\n".join(f"- {c}" for c in fixed_caveats())),
        mo.Html(render.df_to_html(pd.DataFrame(GUARDRAILS, columns=["Safe statement", "Unsafe upgrade"]),
                                  table_id="guardrails", sortable=False)),
        mo.md("The full report adds coverage, the agreement grid over all muscle and heart comparisons, abundance "
              "and class nulls, every MoTrPAC tissue, and provenance."), export])
    return export, tab_about


@app.cell
def _(mo, tab_about, tab_layers, tab_metab, tab_specific, tab_time):
    mo.ui.tabs({"Time course": tab_time, "Do RNA and protein agree?": tab_layers, "Is it specific?": tab_specific,
                "Metabolites": tab_metab, "Caveats and full report": tab_about})
    return


@app.cell
def _(OUT, export, fdr, mo, run, sig, src):
    mo.stop(not export.value)
    _p, _ = run.main_run(src, name=sig.name, outdir=OUT / f"app_{sig.name}", command="app export",
                         cutoff=fdr.value, nboot=1000, quiet=True)
    mo.md(f"Report written: [{_p}]({_p.as_uri()})")
    return


if __name__ == "__main__":
    app.run()
