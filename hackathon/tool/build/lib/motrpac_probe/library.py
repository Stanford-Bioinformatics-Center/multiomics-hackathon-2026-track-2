"""MoTrPAC as a signature library: every store column exported as UP / DOWN gene sets (GMT + JSON index), and a
query that ranks all columns against a user signature.

Set definition (build):
  UP   = genes with fdr_bh < cutoff and logFC > 0 in the column, ranked by |stat| descending, first `cap` kept;
  DOWN = genes with fdr_bh < cutoff and logFC < 0, same ranking and cap.
  All 422 store columns (every kind and layer) are exported. Sets with 0 genes are not written to the GMT but are
  listed in the JSON index with n_genes = 0.

Query (per column):
  The cameraPR / sign-test scores come from everywhere.run (signed cameraPR t, + = the column moves the signature
  genes AGAINST the signature). Set overlap uses the counted signature genes split by direction (sig_up, sig_down)
  and the column's library sets (UP, DOWN, after the FDR cutoff and the cap):
    n_overlap_same    = |sig_up & UP|   + |sig_down & DOWN|
    n_overlap_opposed = |sig_up & DOWN| + |sig_down & UP|
    denominator       = |(sig_up | sig_down) | (UP | DOWN)|
    jaccard_same      = n_overlap_same / denominator
    jaccard_opposed   = n_overlap_opposed / denominator
  (0 when the denominator is 0).
"""
import datetime as dt
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from . import paths, store
from .core import col_words, get_store

GMT = "motrpac_library.gmt"
INDEX = "motrpac_library_index.json"


def _log(msg):
    print(msg, file=sys.stderr)


def set_name(c, direction):
    return "|".join(["MoTrPAC", c.species, c.tissue, c.layer, c.contrast, direction]).replace("\t", " ")


def column_sets(S, cutoff=0.05, cap=250):
    """Yields (column row, direction, genes ranked by |stat| desc and capped, n significant before cap)."""
    for c in S.cols.sort_values("col_order").itertuples():
        pos = S.col(c.column_id)
        fdr, lfc, stat, g = S.fdr[pos], S.lfc[pos], S.stat[pos], S.gcode[pos]
        for direction, m in (("UP", (fdr < cutoff) & (lfc > 0)), ("DOWN", (fdr < cutoff) & (lfc < 0))):
            idx = np.flatnonzero(m)
            idx = idx[np.argsort(-np.abs(stat[idx]), kind="stable")]
            n_sig = len(idx)
            if cap:
                idx = idx[:cap]
            yield c, direction, [str(x) for x in S.genes[g[idx]]], n_sig


def _provenance(cutoff, cap):
    p = store.load_provenance()
    return {"store_sha256_16": store.store_hash(), "store_built_at": p.get("built_at", "NA"),
            "MotrpacRatTraining6moData": p.get("MotrpacRatTraining6moData", "NA"),
            "MotrpacHumanPreSuspensionAnalysis": p.get("MotrpacHumanPreSuspensionAnalysis", "NA"),
            "join_table_version": p.get("join_table_version", "NA"),
            "built_at": dt.datetime.now().isoformat(timespec="seconds"), "cutoff": cutoff, "cap": cap}


README = """MoTrPAC signature library (motrpac_probe)
==========================================

Files
  {gmt}          gene sets, standard GMT format (one set per line:
                                 name<TAB>description<TAB>gene1<TAB>gene2...). Genes are human
                                 HGNC symbols (rat genes mapped to human orthologs as in the store).
  {index}   one JSON object per set, including empty sets (n_genes = 0, not in
                                 the GMT): set name, column_id, species, dataset, tissue, layer,
                                 contrast, label, kind, sex, time, direction, n_genes,
                                 n_significant_before_cap, cutoff, cap, provenance.

Definition
  Every MoTrPAC comparison column in the store ({ncols} columns: rat training vs sedentary, human acute
  exercise vs control and post-pre, and the non-exercise reference comparisons; RNA, protein,
  phosphosite and Olink layers) gives two sets:
    UP   = genes with BH FDR < {cutoff} and logFC > 0, ranked by |stat| descending, at most {cap} genes;
    DOWN = genes with BH FDR < {cutoff} and logFC < 0, same ranking and cap.
  Set names: MoTrPAC|species|tissue|layer|contrast|UP or DOWN (e.g. MoTrPAC|rat|SKM-GN|PROT|F_8w|UP).
  Phosphosite columns are collapsed to genes as in the store.

Loading
  GMT is the library format used by Enrichr and GSEA. Upload the .gmt file as a custom library in
  Enrichr-style tools, pass it to GSEA / fgsea / clusterProfiler::read.gmt / gseapy (gene_sets=path),
  or read it line by line (split on tabs; field 1 = name, field 2 = description, rest = genes).

Caveats
  Set membership depends on the FDR cutoff ({cutoff}) and the cap ({cap}); columns with few significant
  genes give small or empty sets, and a gene missing from a set may simply not be measured in that
  column. The full per-gene statistics are in the store (contrasts.parquet), not in the GMT.
  Built {built_at} from store {sha}.
"""


def build(cutoff=0.05, cap=250, outdir=None):
    """Write the GMT, the JSON index and README.txt; returns (gmt path, index list)."""
    t0 = time.time()
    outdir = Path(outdir or paths.LIBRARY)
    outdir.mkdir(parents=True, exist_ok=True)
    S = get_store()
    _log(f"[library] store loaded in {time.time() - t0:.1f} s")
    prov = _provenance(cutoff, cap)
    lines, index = [], []
    for c, direction, genes, n_sig in column_sets(S, cutoff, cap):
        name = set_name(c, direction)
        desc = (f"{col_words(c)} ({c.kind}) {direction.lower()}, BH FDR < {cutoff:g}, "
                f"top {cap} by |stat|").replace("\t", " ")
        if genes:
            lines.append("\t".join([name, desc] + genes))
        index.append(dict(set_name=name, description=desc, column_id=c.column_id, species=c.species,
                          dataset=c.dataset, tissue=c.tissue, layer=c.layer, contrast=c.contrast, label=c.label,
                          kind=c.kind, sex=c.sex, time=c.time, direction=direction, n_genes=len(genes),
                          n_significant_before_cap=int(n_sig), cutoff=cutoff, cap=cap, in_gmt=bool(genes),
                          provenance=prov))
    names = [r["set_name"] for r in index]
    assert len(names) == len(set(names)), "duplicate set names"
    gmt = outdir / GMT
    gmt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    (outdir / INDEX).write_text(json.dumps(index, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    (outdir / "README.txt").write_text(README.format(gmt=GMT, index=INDEX, ncols=len(S.cols), cutoff=f"{cutoff:g}",
                                                     cap=cap, built_at=prov["built_at"],
                                                     sha=prov["store_sha256_16"]), encoding="utf-8")
    _log(f"[library] build: {len(index)} sets, {len(lines)} non-empty, {time.time() - t0:.1f} s")
    return gmt, index


def summarize(index):
    n = np.array([r["n_genes"] for r in index])
    ne = n[n > 0]
    return dict(n_sets=len(n), n_nonempty=int(len(ne)), n_empty=int((n == 0).sum()),
                median_size_nonempty=float(np.median(ne)) if len(ne) else float("nan"),
                median_size_all=float(np.median(n)), n_at_cap=int((n == (index[0]["cap"] or -1)).sum()))


# ------------------------------------------------------------------------------------------------ query
HEADERS = {"words": "Comparison", "kind": "Kind", "n_measured": "Genes measured", "n_opposed": "Opposed",
           "n_same": "Same", "sign_p": "Sign test p", "camera_t": "cameraPR t", "camera_fdr": "cameraPR FDR",
           "n_overlap_opposed": "Overlap opposed", "n_overlap_same": "Overlap same",
           "jaccard_opposed": "Jaccard opposed", "jaccard_same": "Jaccard same", "n_set_up": "UP set size",
           "n_set_down": "DOWN set size", "column_id": "column_id"}


def overlaps(S, sig, cutoff=0.05, cap=250):
    su = {g for g, d in zip(sig.genes, sig.dirs) if d > 0}
    sd = {g for g, d in zip(sig.genes, sig.dirs) if d < 0}
    sets = {}
    for c, direction, genes, _ in column_sets(S, cutoff, cap):
        sets.setdefault(c.column_id, {})[direction] = set(genes)
    rows = []
    for cid, d in sets.items():
        up, dn = d["UP"], d["DOWN"]
        same = len(su & up) + len(sd & dn)
        opp = len(su & dn) + len(sd & up)
        den = len((su | sd) | (up | dn))
        rows.append(dict(column_id=cid, n_set_up=len(up), n_set_down=len(dn), n_overlap_same=same,
                         n_overlap_opposed=opp, jaccard_same=same / den if den else 0.0,
                         jaccard_opposed=opp / den if den else 0.0))
    return pd.DataFrame(rows)


def query(sig_path, top=25, cutoff=0.05, cap=250):
    """Returns (sig, opposed ranking, concordant ranking, per-column table), rankings with readable headers."""
    from . import everywhere, signature
    S = get_store()
    sig = signature.load(sig_path)
    per, _, _, _ = everywhere.run(S, sig, cutoff=cutoff)
    per = per.merge(overlaps(S, sig, cutoff, cap), on="column_id", how="left")
    keep = list(HEADERS)
    opp = per.sort_values(["camera_t", "jaccard_opposed"], ascending=[False, False], na_position="last")
    conc = per.sort_values(["camera_t", "jaccard_same"], ascending=[True, False], na_position="last")
    fmt = lambda d: d[keep].rename(columns=HEADERS).reset_index(drop=True).rename_axis("Rank").reset_index() \
        .assign(Rank=lambda x: x.Rank + 1)
    return sig, fmt(opp), fmt(conc), per


def _text_table(df, n):
    cols = ["Rank", "Comparison", "Kind", "Genes measured", "Opposed", "cameraPR t", "cameraPR FDR",
            "Jaccard opposed", "Jaccard same"]
    d = df[cols].head(n).copy()
    d["cameraPR t"] = d["cameraPR t"].map(lambda x: "NA" if pd.isna(x) else f"{x:.2f}")
    d["cameraPR FDR"] = d["cameraPR FDR"].map(lambda x: "NA" if pd.isna(x) else f"{x:.2g}")
    for c in ("Jaccard opposed", "Jaccard same"):
        d[c] = d[c].map(lambda x: f"{x:.3f}")
    return d.to_string(index=False)


METHOD = [
    "Every MoTrPAC comparison column in the store (all kinds and layers) is scored against the counted signature "
    "genes. cameraPR t is the signed pre-ranked CAMERA statistic on the column's full moderated-statistic ranking "
    "(inter-gene correlation 0.01); signature member statistics are multiplied by minus the signature direction, "
    "so a positive t means the column moves the signature genes against the signature (opposed) and a negative t "
    "means it moves them the same way (concordant). cameraPR FDR is Benjamini-Hochberg over all columns with a "
    "score. Opposed / Same count measured signature genes whose logFC sign is against / with the signature.",
    "Set overlap uses the column's library sets: UP = genes with BH FDR &lt; {cutoff} and logFC &gt; 0, DOWN = BH FDR "
    "&lt; {cutoff} and logFC &lt; 0, each ranked by |stat| and capped at {cap} genes. With sig_up / sig_down the "
    "counted signature genes by direction: Overlap same = |sig_up &cap; UP| + |sig_down &cap; DOWN|; Overlap "
    "opposed = |sig_up &cap; DOWN| + |sig_down &cap; UP|; both Jaccard indices divide by |(sig_up &cup; sig_down) "
    "&cup; (UP &cup; DOWN)|. Because the sets depend on the FDR cutoff and the cap, overlap is a coarse secondary "
    "view; the ranking uses cameraPR t (ties broken by the matching Jaccard index). Columns with fewer than 2 "
    "measured signature genes have no cameraPR t and are listed last.",
]


def query_cli(a):
    import markupsafe

    from . import render
    t0 = time.time()
    cutoff, cap = 0.05, 250
    sig, opp, conc, per = query(a.signature, top=a.top, cutoff=cutoff, cap=cap)
    t1 = time.time()
    out = Path(a.out) if a.out else paths.OUT / f"library_query_{sig.name}"
    (out / "tables").mkdir(parents=True, exist_ok=True)
    opp.to_csv(out / "tables" / "library_opposed.csv", index=False)
    conc.to_csv(out / "tables" / "library_concordant.csv", index=False)
    n_up, n_dn = int((sig.dirs > 0).sum()), int((sig.dirs < 0).sum())
    print(f"Signature {sig.name}: {len(sig.genes)} counted genes ({n_up} up, {n_dn} down); "
          f"{len(per)} MoTrPAC columns scored.")
    print(f"\nTop {a.top} opposed (cameraPR t descending):")
    print(_text_table(opp, a.top))
    print(f"\nTop {a.top} concordant (cameraPR t ascending):")
    print(_text_table(conc, a.top))

    fm = {"cameraPR t": "{:.2f}", "cameraPR FDR": render.fmt_p, "Sign test p": render.fmt_p,
          "Jaccard opposed": "{:.3f}", "Jaccard same": "{:.3f}"}
    show = [c for c in opp.columns if c != "column_id"]
    tabs = (f"<h3>Most opposed (cameraPR t descending)</h3>"
            + render.df_to_html(opp[show], table_id="library_opposed", formats=fm, max_rows=a.top)
            + f"<h3>Most concordant (cameraPR t ascending)</h3>"
            + render.df_to_html(conc[show], table_id="library_concordant", formats=fm, max_rows=a.top))
    n_opp = int(((per.camera_fdr < 0.05) & (per.camera_t > 0)).sum())
    n_same = int(((per.camera_fdr < 0.05) & (per.camera_t < 0)).sum())
    n_scored = int(per.camera_t.notna().sum())
    prov = {"command": "mprobe " + " ".join(sys.argv[1:]), "signature": str(a.signature),
            "n_counted_genes": len(sig.genes), "n_columns": len(per), **_provenance(cutoff, cap)}
    ctx = dict(
        title=f"MoTrPAC library query: {sig.name}",
        subtitle="All MoTrPAC comparison columns ranked by opposition and concordance with the signature.",
        meta=[("Signature", sig.name), ("Counted genes", f"{len(sig.genes)} ({n_up} up, {n_dn} down)"),
              ("Columns", str(len(per))), ("Store", prov["store_sha256_16"]), ("Run", prov["built_at"])],
        headline=dict(title="Ranked MoTrPAC comparisons",
                      how_to_read="Positive cameraPR t: the comparison moves the signature genes against the "
                                  "signature (opposed); negative: the same way (concordant). FDR is BH over all "
                                  "scored columns. Click a header to sort; the full rankings are in the CSV files.",
                      table_html=markupsafe.Markup(tabs),
                      sentences=[f"{n_opp} of {n_scored} scored columns are opposed at cameraPR FDR < 0.05 and "
                                 f"{n_same} are concordant."]),
        sections=[dict(id="method", title="Method", intro=[],
                       items=[{"kind": "text", "html": markupsafe.Markup(
                           "".join(f"<p>{p.format(cutoff=f'{cutoff:g}', cap=cap)}</p>" for p in METHOD))}])],
        provenance_json=json.dumps(prov, indent=2),
    )
    html = render.render_report(ctx, out / "report.html")
    (out / "provenance.json").write_text(json.dumps(prov, indent=2) + "\n", encoding="utf-8")
    probs = render.validate_report(html)
    print(f"\nwrote {out}/tables/library_opposed.csv, library_concordant.csv, report.html"
          + (f" (validation: {probs})" if probs else " (validation: ok)"))
    _log(f"[library] query: scoring {t1 - t0:.1f} s, total {time.time() - t0:.1f} s")
    return out


def build_cli(a):
    gmt, index = build(cutoff=a.cutoff, cap=a.cap)
    s = summarize(index)
    print(f"wrote {gmt}, {gmt.parent / INDEX}, {gmt.parent / 'README.txt'}")
    print(f"{s['n_sets']} sets ({s['n_nonempty']} non-empty, {s['n_empty']} empty); median size of non-empty sets "
          f"{s['median_size_nonempty']:g}; {s['n_at_cap']} sets at the cap of {a.cap}")
