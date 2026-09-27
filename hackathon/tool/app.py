"""motrpac_probe interactive app (marimo). Thin layer over the package API: no analysis code here.

Run:     marimo run app.py                      (store is preloaded at start-up)
Check:   marimo export html app.py -o /tmp/app.html [-- --example NAME | --paste FILE.csv]
         (--example: examples.json name selected by default; --paste: CSV whose first two columns are pasted)
"""
import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium", app_title="motrpac_probe", css_file="src/motrpac_probe/templates/style.css")


@app.cell
def _():
    import json
    import tempfile
    import time
    from pathlib import Path

    import marimo as mo
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    from motrpac_probe import core, figures, layers, metab, render, run
    from motrpac_probe.core import CORE_TISSUES, LAYER_WORDS, TISSUE_WORDS, col_words, short_words
    from motrpac_probe.paths import EXAMPLES, OUT

    _t = time.time()
    S = core.get_store()                 # gene store (RNA / protein / phospho)
    SM = metab.metab_store()             # + metabolomics, for the pathway panel
    load_s = time.time() - _t
    # marimo's markdown uses its own heading font; match the reports (style.css is loaded via css_file)
    mo.Html("<style>.markdown h1, .markdown h2, .markdown h3, .prose h1, .prose h2, .prose h3 "
            "{font-family: var(--font-body); font-weight: 600;}</style>")
    return (CORE_TISSUES, EXAMPLES, LAYER_WORDS, TISSUE_WORDS, OUT, Path, S, SM, col_words, core, figures, json, layers, load_s, metab, mo,
            np, pd, plt, render, run, short_words, tempfile, time)


@app.cell
def _(EXAMPLES, json, load_s, mo, pd):
    manifest = json.loads((EXAMPLES / "examples.json").read_text())
    gene_examples = {e["title"]: e for e in manifest if e.get("kind") != "metabolite demo"}
    _args = mo.cli_args()
    _by_name = {e["name"]: e["title"] for e in gene_examples.values()}
    _default = _by_name.get(_args.get("example", "pah_muscle_malenfant2015"), "(none)")
    _paste0 = ""
    if _args.get("paste"):
        _p = pd.read_csv(_args.get("paste"), dtype=str)
        _paste0 = "\n".join(f"{a},{b}" for a, b in zip(_p.iloc[:, 0], _p.iloc[:, 1]))
    sig_form = mo.md(
        "{example}\n\n{upload}\n\n{paste}"
    ).batch(
        example=mo.ui.dropdown(options=["(none)"] + list(gene_examples), value=_default, label="Example signature"),
        upload=mo.ui.file(filetypes=[".csv"], label="Upload a CSV"),
        paste=mo.ui.text_area(value=_paste0, placeholder="GENE,direction per line, e.g.\nNDUFA9,-1\nHBB,1", rows=6,
                              full_width=True, label="or paste genes"),
    ).form(submit_button_label="Use this signature", bordered=False)
    mo.vstack([
        mo.md("# motrpac_probe\nDoes exercise (MoTrPAC) move a disease signature's genes against the disease, in "
              f"RNA and in protein? Store loaded in {load_s:.1f} s.\n\n"
              "## 1. Signature\nColumns: `gene_symbol` (or `uniprot` / `ensembl` / `rat_symbol`), `direction` "
              "(+1 up in disease, −1 down), optional `group`. Upload wins over paste, paste over the example."),
        sig_form])
    return gene_examples, sig_form


@app.cell
def _(CORE_TISSUES, mo):
    opt_form = mo.md("{tissues}\n\n{cutoff} {universe} {nboot}").batch(
        tissues=mo.ui.multiselect(options=CORE_TISSUES, value=CORE_TISSUES, label="Core tissues"),
        cutoff=mo.ui.dropdown(options=["0.05", "0.10"], value="0.05", label="BH FDR cutoff"),
        universe=mo.ui.dropdown(options=["all", "muscle-intrinsic"], value="all", label="cameraPR background"),
        nboot=mo.ui.dropdown(options=["200", "1000"], value="200", label="Random sets per null"),
    ).form(submit_button_label="Apply options", bordered=False)
    mo.vstack([mo.md("## 2. Options\nDefaults are used until you press *Apply options*. 1,000 random sets per "
                     "null (the report default) takes about 10 s longer than 200."), opt_form])
    return (opt_form,)


@app.cell
def _(EXAMPLES, Path, gene_examples, mo, opt_form, run, sig_form, tempfile, time):
    # forms return None until first submitted: use the defaults shown in the widgets
    sv = sig_form.value if sig_form.value is not None else sig_form.element.value
    ov = opt_form.value if opt_form.value is not None else opt_form.element.value
    if sv["upload"]:
        src = Path(tempfile.mkdtemp()) / sv["upload"][0].name
        src.write_bytes(sv["upload"][0].contents)
        how = "uploaded file"
    elif sv["paste"].strip():
        _lines = [l.replace("\t", ",").strip() for l in sv["paste"].strip().splitlines() if l.strip()]
        if not _lines[0].lower().startswith(("gene", "uniprot", "ensembl", "rat_symbol")):
            _lines = ["gene_symbol,direction"] + _lines
        src = Path(tempfile.mkdtemp()) / "pasted.csv"
        src.write_text("\n".join(_lines) + "\n")
        how = "pasted genes"
    elif sv["example"] != "(none)":
        src = EXAMPLES / gene_examples[sv["example"]]["file"]
        how = f"example: {sv['example']}"
    else:
        src = None
    mo.stop(src is None, mo.md("Choose an example, upload a CSV or paste genes."))
    _t = time.time()
    with mo.status.spinner(title="Scoring the signature against every comparison ..."):
        R = run.compute(src, name=src.stem, cutoff=float(ov["cutoff"]), nboot=int(ov["nboot"]),
                        tissues=tuple(ov["tissues"]), universe=ov["universe"], quiet=True)
    compute_s = time.time() - _t
    mo.md(f"**Signature:** {how} (`{src.name}`) · {len(R['sig'].genes)} counted genes · "
          f"{len(R['core_cols'])} comparisons · computed in {compute_s:.1f} s")
    return (R,)


@app.cell
def _(LAYER_WORDS, R, TISSUE_WORDS, mo, render):
    _mp = R["sig"].table
    _unm = _mp[_mp.status != "mapped"]
    _cov = R["coverage"].rename(columns={b: f"{'Human' if b[0] == 'human_acute' else 'Rat'} {TISSUE_WORDS[b[1]]} "
                                            f"{LAYER_WORDS[b[2]]}" for b in R["cov_blocks"]})
    _cov = _cov.rename(columns={"gene": "Gene", "group": "Group"})
    _idcol = [c for c in ["gene_symbol", "uniprot", "ensembl", "rat_symbol"] if c in _unm.columns][:1]
    mo.vstack([
        mo.md(f"## 3. Coverage\n{int((_mp.status == 'mapped').sum())} of {len(_mp)} input rows mapped to a human "
              f"gene symbol; {len(R['sig'].genes)} are counted in the statistics (phenotype-group genes are shown "
              "but not counted). Cell = number of comparisons (timepoints × sexes or groups) in which the gene is "
              "measured, out of the total for that tissue and layer."),
        mo.accordion({f"Table: coverage of the {len(_cov)} signature genes": mo.Html(
            render.df_to_html(_cov, max_rows=80, table_id="app-cov"))}),
        mo.md(f"**Not used ({len(_unm)} rows)**") if len(_unm) else mo.md("All input rows mapped."),
        mo.Html(render.df_to_html(_unm[_idcol + ["status"]].rename(columns={"status": "Reason"}),
                                  table_id="app-unm")) if len(_unm) else mo.md(""),
    ])
    return


@app.cell
def _(R, S, figures, mo, plt, render):
    _genes = R["sig"].shown
    gene = mo.ui.dropdown(options=_genes, value=_genes[0] if _genes else None, label="Gene", searchable=True)
    with render.mpl_style():
        _f = figures.grid(S, R["grid"], R["core_cols"], R["sig"], R["opts"]["cap"], R["opts"]["cutoff"])
        _svg = render.fig_to_svg(_f)
        plt.close(_f)
    mo.vstack([mo.md("## 4. Agreement grid\nCell = agreement × |stat| capped at 4: blue = exercise moves the gene "
                     f"against the disease, red = same direction as the disease, hatched = not measured; dot = "
                     f"BH FDR < {R['opts']['cutoff']}. One row per signature gene, one column per comparison."),
               mo.Html(_svg),
               mo.md("### Gene trajectory\nPick a gene to see its exercise statistic over time in RNA and protein."),
               gene])
    return (gene,)


@app.cell
def _(R, S, col_words, gene, mo, np, plt, render, short_words):
    mo.stop(gene.value is None)
    from motrpac_probe.render import INK, LAYER_COLORS
    _L = R["grid"][R["grid"].gene == gene.value].merge(
        S.cols.reset_index(drop=True)[["column_id", "dataset", "tissue", "layer", "sex", "group", "col_order"]].rename(
            columns={"group": "arm"}),
        on="column_id").sort_values("col_order")
    _panels = list(dict.fromkeys(zip(_L.dataset, _L.tissue)))
    with render.mpl_style():
        _fig, _axes = plt.subplots(1, len(_panels), figsize=(3.3 * len(_panels) + 0.6, 4.2), squeeze=False)
        for _ax, (_ds, _tis) in zip(_axes[0], _panels):
            _d = _L[(_L.dataset == _ds) & (_L.tissue == _tis)]
            _ser = "sex" if _ds == "rat_train" else "arm"
            _xt = list(dict.fromkeys(short_words(S.cols.loc[c]).split(" ", 1)[1] for c in _d.column_id))
            for (_lay, _s), _g in _d.groupby(["layer", _ser], sort=False):
                _x = [_xt.index(short_words(S.cols.loc[c]).split(" ", 1)[1]) for c in _g.column_id]
                _ls = "-" if _s in ("F", "EE") else "--"
                _ax.plot(_x, _g.stat, _ls, color=LAYER_COLORS[_lay], lw=1.6,
                         label=f"{'protein' if _lay == 'PROT' else _lay}, {_s}")
                _sig = (_g.fdr_bh < R["opts"]["cutoff"]).to_numpy()
                _ax.scatter(np.array(_x)[_sig], _g.stat[_sig], color=LAYER_COLORS[_lay], s=30, zorder=3)
                _ax.scatter(np.array(_x)[~_sig], _g.stat[~_sig], facecolor="white", edgecolor=LAYER_COLORS[_lay],
                            s=30, zorder=3)
            _ax.axhline(0, color=INK, lw=0.6)
            _ax.set_xticks(range(len(_xt)), _xt, rotation=30, ha="right", fontsize=11)
            _ax.set_title(col_words(S.cols.loc[_d.column_id.iloc[0]]).split(", ")[0].rsplit(" ", 1)[0], fontsize=11)
            _ax.legend(fontsize=10, frameon=False, ncol=2, loc="upper left", bbox_to_anchor=(0, -0.3))
        _axes[0][0].set_ylabel("exercise stat (t or z)", fontsize=11)
        _fig.tight_layout()
        _svg = render.fig_to_svg(_fig)
        plt.close(_fig)
    _tab = _L.assign(Comparison=[col_words(S.cols.loc[c]) for c in _L.column_id])[
        ["Comparison", "logFC", "stat", "fdr_bh", "agreement"]].rename(
        columns={"stat": "stat (t or z)", "fdr_bh": "BH FDR", "agreement": "agreement (−1 = opposed)"})
    _dir = "up" if dict(zip(R["sig"].shown, R["sig"].shown_dirs))[gene.value] > 0 else "down"
    mo.vstack([mo.md(f"**{gene.value}** (listed {_dir} in disease). Filled point = BH FDR < "
                     f"{R['opts']['cutoff']}; solid = female / EE, dashed = male / RE."),
               mo.Html(_svg),
               mo.accordion({f"Table: {gene.value} in every comparison ({len(_tab)} rows)": mo.Html(
                   render.df_to_html(_tab, table_id="app-gene", formats={"BH FDR": render.fmt_p}))})])
    return


@app.cell
def _(R, S, col_words, figures, mo, plt, render):
    _sc = R["scores"].copy()
    _nl = R["nulls"]
    for _k, _lab in [("abundance", "Abundance null (percentile)"), ("class", "Class null (percentile)")]:
        _sc[_lab] = _sc.column_id.map(_nl[_nl.null == _k].set_index("column_id").pct_t)
    _sc["Comparison"] = [col_words(S.cols.loc[c]) for c in _sc.column_id]
    _v = _sc[["Comparison", "n_measured", "n_opposed", "camera_t", "camera_fdr", "Abundance null (percentile)",
              "Class null (percentile)"]].rename(columns={
        "n_measured": "Measured", "n_opposed": "Opposed", "camera_t": "cameraPR t", "camera_fdr": "FDR"})
    with render.mpl_style():
        _f1 = figures.camera(S, R["scores"], R["core_cols"], R["opts"]["cutoff"])
        _f2 = figures.nulls(S, R["scores"], _nl, R["core_cols"])
        _s1, _s2 = render.fig_to_svg(_f1), render.fig_to_svg(_f2)
        plt.close(_f1)
        plt.close(_f2)
    mo.vstack([mo.md("## 5. Set-level opposition and null percentiles\ncameraPR t > 0 = the signature as a whole "
                     "is opposed by exercise (circles; triangles = up-genes and down-genes alone). Percentile = "
                     f"where the signature's t falls among {R['opts']['nboot']} random gene sets matched on "
                     "abundance, or on MitoCarta / GO:CC complex / secreted class. A high class-null percentile "
                     "means the result is specific to these genes, not to their pathway class."),
               mo.Html(_s1), mo.Html(_s2),
               mo.accordion({f"Table: signature score per comparison ({len(_v)} rows; measured / opposed = "
                             "number of signature genes; cameraPR t > 0 = opposed)": mo.Html(render.df_to_html(
                   _v, table_id="app-cam", formats={"cameraPR t": "{:+.2f}", "FDR": render.fmt_p,
                                                    "Abundance null (percentile)": "{:.0f}",
                                                    "Class null (percentile)": "{:.0f}"}))})])
    return


@app.cell
def _(R, figures, layers, mo, plt, render):
    _d = R["disc"]
    mo.stop(not len(_d), mo.md("## 6. Layer discordance\nNo RNA–protein pairs for this signature."))
    _sents = [layers.sentence(r) for r in _d.itertuples() if not r.cross]
    with render.mpl_style():
        _f = figures.discordance(_d, R["disc_genes"], [])
        _svg = render.fig_to_svg(_f)
        plt.close(_f)
    _t = _d[["label", "n_genes_both", "rho_all", "n_sig_both", "rna_opposed", "prot_opposed", "rho_sig"]].rename(
        columns={"label": "Comparison", "n_genes_both": "Genes in both layers (all)", "rho_all": "ρ all genes",
                 "n_sig_both": "Signature genes in both layers", "rna_opposed": "Opposed in RNA",
                 "prot_opposed": "Opposed in protein", "rho_sig": "ρ signature genes"})
    mo.vstack([mo.md("## 6. Layer discordance\nSpearman ρ between the RNA and protein statistics of the same "
                     "comparison, over all genes and over the signature genes."),
               mo.md("\n".join(f"- {s}" for s in _sents[:8])),
               mo.Html(_svg),
               mo.accordion({f"Table: RNA vs protein per comparison ({len(_t)} rows)": mo.Html(render.df_to_html(
                   _t, table_id="app-disc", formats={"ρ all genes": "{:.2f}", "ρ signature genes": "{:.2f}"}))})])
    return


@app.cell
def _(mo):
    pw_tissue = mo.ui.dropdown(options={"Rat gastrocnemius (SKM-GN)": "SKM-GN", "Rat heart (HEART)": "HEART",
                                        "Human vastus lateralis (VL)": "VL"}, value="Rat gastrocnemius (SKM-GN)",
                               label="Tissue")
    return (pw_tissue,)


@app.cell
def _(SM, metab, mo, plt, pw_tissue, render):
    _tis = pw_tissue.value
    _pw = metab.pathway_panel(SM, tissue=_tis, species="human" if _tis == "VL" else "rat")
    _f = metab.pathway_heatmap(_pw, _tis)
    _w = _f.get_size_inches()[0] * 96  # show at native size (text stays >= 10 pt), scroll if wider than the page
    _svg = render.fig_to_svg(_f)
    plt.close(_f)
    mo.vstack([mo.md("## 7. Metabolic pathways in three layers\nIndependent of the signature: how each pathway's "
                     "genes (RNA, protein) and metabolites respond to exercise in one tissue."),
               pw_tissue,
               mo.Html(f'<div style="overflow-x:auto"><div style="min-width:{_w:.0f}px">{_svg}</div></div>'),
               mo.md(metab.PATHWAY_CAPTION)])
    return


@app.cell
def _(mo, pd, render):
    from motrpac_probe.caveats import GUARDRAILS, fixed_caveats
    mo.vstack([mo.md("## 8. Caveats"), mo.md("\n".join(f"- {c}" for c in fixed_caveats())),
               mo.Html(render.df_to_html(pd.DataFrame(GUARDRAILS, columns=["Safe statement", "Unsafe upgrade"]),
                                         table_id="app-guard", sortable=False))])
    return


@app.cell
def _(mo):
    export = mo.ui.run_button(label="Export report")
    mo.vstack([mo.md("## 9. Export\nWrites the full static report (all sections, tables as CSV, provenance) for "
                     "the current signature and options."), export])
    return (export,)


@app.cell
def _(OUT, R, export, mo, run):
    mo.stop(not export.value)
    with mo.status.spinner(title="Writing report ..."):
        _p = run.write(R, OUT / f"app_{R['sig'].name}", command="marimo app: Export report")
    mo.vstack([mo.md(f"Report written to `{_p}` ([open]({_p.as_uri()}); browsers may block file links from "
                     "the app, use the download instead)."),
               mo.download(data=_p.read_bytes(), filename=f"{R['sig'].name}_report.html", mimetype="text/html",
                           label="Download report.html")])
    return


if __name__ == "__main__":
    app.run()
