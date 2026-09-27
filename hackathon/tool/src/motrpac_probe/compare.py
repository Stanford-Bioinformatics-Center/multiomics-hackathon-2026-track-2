"""Several disease signatures side by side: is exercise opposition disease-specific, or a property of the
kind of genes a signature contains?

Method details: docs/METHODS.md#comparepy
"""
import datetime as dt
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import __version__, core, nulls, signature, store  # noqa: E402
from .paths import OUT  # noqa: E402

COLUMNS = [("human_acute|VL|RNA|EE_vs_CON_24h", "Human VL RNA\nEE 24 h"),
           ("human_acute|VL|RNA|RE_vs_CON_24h", "Human VL RNA\nRE 24 h"),
           ("human_acute|VL|PROT|EE_vs_CON_24h", "Human VL protein\nEE 24 h"),
           ("human_acute|VL|PROT|RE_vs_CON_24h", "Human VL protein\nRE 24 h"),
           ("rat_train|SKM-GN|PROT|F_4w", "Rat gastroc. protein\nF 4 wk"),
           ("rat_train|SKM-GN|PROT|M_4w", "Rat gastroc. protein\nM 4 wk"),
           ("rat_train|SKM-GN|PROT|F_8w", "Rat gastroc. protein\nF 8 wk"),
           ("rat_train|SKM-GN|PROT|M_8w", "Rat gastroc. protein\nM 8 wk"),
           ("rat_train|SKM-GN|RNA|F_8w", "Rat gastroc. RNA\nF 8 wk"),
           ("rat_train|SKM-GN|RNA|M_8w", "Rat gastroc. RNA\nM 8 wk"),
           ("rat_train|HEART|PROT|F_8w", "Rat heart protein\nF 8 wk"),
           ("rat_train|HEART|PROT|M_8w", "Rat heart protein\nM 8 wk")]
MITO_SHARE = 0.3   # a signature is "mitochondria-dominated" when >= 30% of counted genes are MitoCarta3.0


def compute(paths, titles=None, nboot=1000, seed=20260926, cutoff=0.05, kinds=None):
    S = core.get_store()
    cids = [c for c, _ in COLUMNS]
    rows = []
    for i, p in enumerate(paths):
        from . import ranked
        sig = ranked.load_ranked(p)[1] if ranked.is_ranked(p) else signature.load(p)
        sc = core.score_columns(S, cids, sig, cutoff=cutoff).set_index("column_id")
        nl = nulls.run_nulls(S, cids, sig, B=nboot, seed=seed)
        ann = S.ann.reindex(sig.genes)
        mito = float(ann.mitocarta3.astype("boolean").fillna(False).mean()) if len(sig.genes) else np.nan
        infl = S.gs.get("HALLMARK_INFLAMMATORY_RESPONSE", set()) | S.gs.get("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", set())
        for cid, lab in COLUMNS:
            a = nl[(nl.column_id == cid) & (nl.null == "abundance")].iloc[0]
            c = nl[(nl.column_id == cid) & (nl.null == "class")].iloc[0]
            r = sc.loc[cid]
            rows.append(dict(signature=sig.name, title=(titles or {}).get(sig.name, sig.name), n_genes=len(sig.genes),
                             kind=(kinds or {}).get(sig.name, "disease"),
                             mito_share=mito, inflam_ecm_share=float(np.mean([g in infl for g in sig.genes])) if sig.genes else np.nan,
                             column_id=cid, column=lab.replace("\n", " "), n_measured=r.n_measured,
                             n_opposed=r.n_opposed, camera_t=r.camera_t, camera_p=r.camera_p, pct_abund=a.pct_t,
                             pct_class=c.pct_t))
    df = pd.DataFrame(rows)
    q = np.full(len(df), np.nan)
    ok = df.camera_p.notna().to_numpy()
    from scipy.stats import false_discovery_control
    q[ok] = false_discovery_control(df.camera_p[ok], method="bh")
    df["camera_fdr"] = q
    return df


def figure(df):
    from .render import DIVERGING, INK
    sigs = list(dict.fromkeys(df.signature))
    titles = df.drop_duplicates("signature").set_index("signature").title
    M = df.pivot(index="signature", columns="column_id", values="camera_t").loc[sigs, [c for c, _ in COLUMNS]]
    P = df.pivot(index="signature", columns="column_id", values="pct_abund").loc[sigs, [c for c, _ in COLUMNS]]
    Q = df.pivot(index="signature", columns="column_id", values="camera_fdr").loc[sigs, [c for c, _ in COLUMNS]]
    fig, ax = plt.subplots(figsize=(1.05 * len(COLUMNS) + 3.6, 0.62 * len(sigs) + 2.2))
    cmap = DIVERGING.reversed().copy()   # blue = positive t = opposed
    cmap.set_bad("white")
    im = ax.imshow(np.ma.masked_invalid(np.clip(M.to_numpy(float), -8, 8)), cmap=cmap, vmin=-8, vmax=8, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            t, p, q = M.iat[i, j], P.iat[i, j], Q.iat[i, j]
            if not np.isfinite(t):
                ax.text(j, i, "n/a", ha="center", va="center", fontsize=9, color="#777777")
                continue
            star = "*" if np.isfinite(q) and q < 0.05 else ""
            ax.text(j, i, f"{t:+.1f}{star}\n{p:.0f}", ha="center", va="center", fontsize=10,
                    color="white" if abs(t) > 5 else INK)
    ax.set_xticks(range(len(COLUMNS)), [lab for _, lab in COLUMNS], fontsize=9.5, rotation=45, ha="right")
    ax.set_yticks(range(len(sigs)), [titles[s] for s in sigs], fontsize=10.5)
    for x in (3.5, 9.5):
        ax.axvline(x, color=INK, lw=1)
    ax.tick_params(length=0)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("cameraPR t (+ = exercise opposes)", fontsize=10)
    fig.tight_layout()
    return fig


def lesson(df, cutoff=0.05):
    """Rule-based sentences; the 'mitochondria' lesson is stated only when the numbers support it."""
    per = df.groupby("signature", sort=False).agg(title=("title", "first"), kind=("kind", "first"),
                                      mito=("mito_share", "first"),
                                      n_opp=("camera_fdr", lambda q: int(((q < cutoff) & (df.loc[q.index, "camera_t"] > 0)).sum())),
                                      n_cols=("camera_t", lambda t: int(t.notna().sum())),
                                      median_class=("pct_class", "median"))
    dis = per[per.kind == "disease"]  # controls never enter the lesson
    mito = dis[dis.mito >= MITO_SHARE]
    other = dis[dis.mito < MITO_SHARE]
    out = [f"{r.title}{'' if r.kind == 'disease' else ' (' + r.kind + ')'}: opposed at FDR < {cutoff} in {r.n_opp} of {r.n_cols} columns; MitoCarta share "
           f"{100 * r.mito:.0f}%; median pathway-class null percentile {r.median_class:.0f}."
           for r in per.itertuples()]
    if len(mito) >= 2 and len(other) >= 1:
        fm = (mito.n_opp / mito.n_cols).median()
        fo = (other.n_opp / other.n_cols).median()
        if fm >= 0.4 and fo <= fm / 2:
            out.insert(0, f"Disease signatures dominated by mitochondrial genes (MitoCarta share ≥ {100 * MITO_SHARE:.0f}%: "
                          f"{', '.join(mito.title)}) are opposed in a median {100 * fm:.0f}% of columns, the others in "
                          f"{100 * fo:.0f}%. Exercise raises mitochondrial programs; that opposition is not "
                          f"disease-specific, and the pathway-class null percentiles show it per disease.")
        else:
            out.insert(0, f"Mitochondria-dominated disease signatures are opposed in a median {100 * fm:.0f}% of columns and "
                          f"the others in {100 * fo:.0f}%; these numbers do not support a simple 'mitochondria-down "
                          "signatures are opposed, others are not' reading.")
    return out, per.reset_index()


def main_compare(paths, titles=None, outdir=None, command="", nboot=1000, kinds=None):
    from . import render
    out = Path(outdir) if outdir else OUT / "compare"
    (out / "figures").mkdir(parents=True, exist_ok=True)
    (out / "tables").mkdir(parents=True, exist_ok=True)
    df = compute(paths, titles=titles, nboot=nboot, kinds=kinds)
    df.to_csv(out / "tables" / "compare_long.csv", index=False)
    sents, per = lesson(df)
    per.to_csv(out / "tables" / "compare_per_signature.csv", index=False)
    with render.mpl_style():
        svg = render.save_fig(figure(df), out / "figures" / "compare.png")
    view = df[["title", "column", "n_measured", "n_opposed", "camera_t", "camera_fdr", "pct_abund", "pct_class",
               "mito_share"]].rename(columns={
        "title": "Signature", "column": "Comparison", "n_measured": "Genes measured", "n_opposed": "Opposed",
        "camera_t": "cameraPR t", "camera_fdr": "FDR", "pct_abund": "Abundance null percentile",
        "pct_class": "Class null percentile", "mito_share": "MitoCarta share"})
    view.to_csv(out / "tables" / "compare_table.csv", index=False)
    pv = per.rename(columns={"signature": "File", "title": "Signature", "kind": "Kind", "mito": "MitoCarta share",
                             "n_opp": "Columns opposed (FDR < 0.05)", "n_cols": "Columns measured",
                             "median_class": "Median class null percentile"})
    pv.to_csv(out / "tables" / "compare_summary.csv", index=False)
    prov = {"tool": f"motrpac_probe {__version__}", "timestamp": dt.datetime.now().isoformat(timespec="seconds"),
            "command": command, "store": store.store_hash(), "repo_git_sha": store.git_sha(),
            "signatures": [str(p) for p in paths], "nboot": nboot}
    (out / "provenance.json").write_text(json.dumps(prov, indent=2), encoding="utf-8")
    ctx = dict(
        title="MoTrPAC probe: several disease signatures side by side",
        subtitle="Is the exercise opposition specific to a disease, or to the kind of genes in its signature?",
        meta=[("Signatures", str(len(paths))), ("Store", store.store_hash())],
        headline=dict(title="Set-level opposition per disease signature, layer and time",
                      how_to_read=(f"Cell: signed cameraPR t (+ / blue = exercise opposes the disease direction; "
                                   f"* = BH FDR < 0.05 over all cells) and, below it, the percentile of t among {nboot} "
                                   "abundance-matched random gene sets. Columns: human vastus lateralis 24 h after "
                                   "one bout; rat gastrocnemius and heart after 4–8 weeks of training."),
                      table_html="", sentences=sents),
        first_screen=["The class-null column in the table asks whether a random set with the same MitoCarta / "
                      "complex / secreted composition would look just as opposed."],
        sections=[dict(id="compare", title="Figure and table", intro=[], items=[
            dict(kind="figure", id="fig-compare", title="Opposition by disease signature, layer and time",
                 how_to_read=("Rows = disease signatures and controls; blue = opposed, red = same direction as the "
                              "disease; number under t = abundance-matched null percentile; n/a = fewer than 2 "
                              "genes measured."), svg=svg, png="figures/compare.png"),
            dict(kind="table", title="Per signature", caption="", csv="tables/compare_summary.csv",
                 html=render.df_to_html(pv, table_id="compare_summary",
                                        formats={"MitoCarta share": "{:.2f}", "Median class null percentile": "{:.0f}"})),
            dict(kind="table", title="All cells", caption="Percentiles = 100 × P(null t ≤ observed t).",
                 csv="tables/compare_table.csv",
                 html=render.df_to_html(view, table_id="compare_table", formats={
                     "cameraPR t": "{:+.2f}", "FDR": render.fmt_p, "Abundance null percentile": "{:.0f}",
                     "Class null percentile": "{:.0f}", "MitoCarta share": "{:.2f}"}, max_rows=400))])],
        provenance_json=json.dumps(prov, indent=2))
    render.render_report(ctx, out / "report.html")
    return out / "report.html", df
