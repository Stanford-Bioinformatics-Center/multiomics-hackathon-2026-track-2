"""Static gallery: copies the example run folders and the signature library into tool/site/ and writes index.html.

All links are relative, so `python -m http.server --directory tool/site` serves everything.
"""
import json
import shutil
import sys
import time
from pathlib import Path

import markupsafe
import pandas as pd

from . import paths, render
from .library import GMT, INDEX



def manifest():
    """examples/examples.json (name, file, title, kind, description, run_args); falls back to examples/*.csv."""
    m = paths.EXAMPLES / "examples.json"
    if m.exists():
        return json.loads(m.read_text())
    return [dict(name=f.stem, file=f.name, title=f.stem, kind="", description="", run_args=[])
            for f in sorted(paths.EXAMPLES.glob("*.csv"))]


def query_example():
    """The example used for the library-query demo: manifest entry with "query": true, else the first disease."""
    ex = manifest()
    for e in ex:
        if e.get("query"):
            return e["name"]
    dis = [e["name"] for e in ex if e.get("kind") == "disease"]
    return dis[0] if dis else (ex[0]["name"] if ex else "")


RUN_PARTS = ("report.html", "figures", "tables", "provenance.json")
HEADLINE = [  # (column header, comparison text in column_scores.csv)
    ("Human VL RNA EE 24 h: opposed/measured, t", "Human vastus lateralis RNA, EE–CON 24 h"),
    ("Human VL protein RE 24 h: opposed/measured, t", "Human vastus lateralis protein, RE–CON 24 h"),
    ("Rat gastrocnemius protein F 8 wk: opposed/measured, t", "Rat gastrocnemius protein, female 8 wk"),
]


def _copy(src, dst):
    if src.is_dir():
        shutil.copytree(src, dst, dirs_exist_ok=True)
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def _n_counted(run):
    m = pd.read_csv(run / "tables" / "mapping.csv", dtype=str, keep_default_na=False)
    ctx = ["phenotype"]
    pj = run / "provenance.json"
    if pj.exists():
        ctx = json.loads(pj.read_text()).get("parameters", {}).get("context_groups", ctx) or []
    grp = m["group"] if "group" in m.columns else pd.Series([""] * len(m))
    return int(((m.status == "mapped") & ~grp.isin(ctx)).sum())


def _headline(run):
    cs = pd.read_csv(run / "tables" / "column_scores.csv")
    h = {}
    for label, comp in HEADLINE:
        r = cs[cs.comparison == comp]
        if r.empty or pd.isna(r.iloc[0].n_measured):
            h[label] = None
        else:
            r = r.iloc[0]
            t = "NA" if pd.isna(r.camera_t) else f"{r.camera_t:.2f}"
            h[label] = f"{int(r.n_opposed)}/{int(r.n_measured)}, t = {t}"
    ok = cs.camera_t.notna()
    if ok.any():
        b = cs.loc[cs.camera_t[ok].abs().idxmax()]
        h["Strongest set-level comparison: t, class-null percentile"] = \
            f"{b.comparison}: t = {b.camera_t:+.2f}, pct {b.pct_class:.0f}"
    else:
        h["Strongest set-level comparison: t, class-null percentile"] = None
    scored = cs.camera_t.notna()
    h["Opposed at cameraPR FDR < 0.05 / scored"] = \
        f"{int(((cs.camera_fdr < 0.05) & (cs.camera_t > 0)).sum())}/{int(scored.sum())}"
    h["Same direction at cameraPR FDR < 0.05 / scored"] = \
        f"{int(((cs.camera_fdr < 0.05) & (cs.camera_t < 0)).sum())}/{int(scored.sum())}"
    return h


def build_site(site=None, out=None):
    t0 = time.time()
    site = Path(site or paths.SITE)
    out = Path(out or paths.OUT)
    site.mkdir(parents=True, exist_ok=True)
    rows, notes = [], []
    cmp_ = out / "compare"
    if (cmp_ / "report.html").exists():
        dst = site / "compare"
        if dst.exists():
            shutil.rmtree(dst)
        for part in RUN_PARTS:
            if (cmp_ / part).exists():
                _copy(cmp_ / part, dst / part)
        summ = pd.read_csv(cmp_ / "tables" / "compare_summary.csv")
        rows.append(dict(name="All disease signatures side by side", href="compare/report.html",
                         description="Generality view: every example signature against the same human 24 h and rat "
                                     "4–8 wk comparisons (mprobe compare).", n_genes=len(summ),
                         headline={"Signatures compared": str(len(summ))}))
    else:
        notes.append("no compare report in out/compare (run `mprobe compare`)")
    for run in sorted(out.glob("discord_*")):  # omic-discordance reports (mprobe discord)
        if not (run / "report.html").exists():
            continue
        dst = site / run.name
        if dst.exists():
            shutil.rmtree(dst)
        for part in RUN_PARTS:
            if (run / part).exists():
                _copy(run / part, dst / part)
        prm = json.loads((run / "provenance.json").read_text()).get("parameters", {})
        rows.append(dict(name=f"Layer discordance: {prm.get('species', '')} {prm.get('tissue', run.name[8:])}",
                         href=f"{run.name}/report.html",
                         description="Which omic layer responds, when, and why RNA and protein disagree "
                                     "(mprobe discord): timescale, agreement, detection power, protein-only calls, "
                                     "lag, sex, 3-layer pathways.", n_genes=None, headline={}))
    for e in manifest():
        name, desc = e["name"], (e.get("title", e["name"]) + ". " + e.get("description", "")).strip()
        run = out / name
        if not (run / "report.html").exists() or not (run / "tables" / "column_scores.csv").exists():
            notes.append(f"skipped {name}: no report.html / tables/column_scores.csv in {run}")
            continue
        dst = site / name
        if dst.exists():
            shutil.rmtree(dst)  # generated copy, rebuilt from out/
        for part in RUN_PARTS:
            if (run / part).exists():
                _copy(run / part, dst / part)
        rows.append(dict(name=name, href=f"{name}/report.html", description=desc, n_genes=_n_counted(run),
                         headline=_headline(run)))

    lib = site / "library"
    lib.mkdir(parents=True, exist_ok=True)
    lib_links = []
    for f in (GMT, INDEX, "README.txt"):
        if (paths.LIBRARY / f).exists():
            shutil.copy2(paths.LIBRARY / f, lib / f)
            lib_links.append(f)
        else:
            notes.append(f"library file missing: {paths.LIBRARY / f} (run `mprobe library build`)")
    QUERY_EXAMPLE = query_example()
    q = out / f"library_query_{QUERY_EXAMPLE}"
    qdst = lib / f"query_{QUERY_EXAMPLE}"
    if (q / "report.html").exists():
        if qdst.exists():
            shutil.rmtree(qdst)
        for part in RUN_PARTS:
            if (q / part).exists():
                _copy(q / part, qdst / part)
    else:
        notes.append(f"no library query for {QUERY_EXAMPLE} in {q} (run `mprobe library query`)")

    intro = [
        "Each example is a gene signature (direction = up or down in the condition) scored by motrpac_probe against "
        "MoTrPAC exercise and training comparisons. A positive cameraPR t means the comparison moves the signature "
        "genes against the signature (opposed); a negative t means it moves them the same way.",
        "The headline columns are read from each run's tables/column_scores.csv (the core muscle and heart "
        "comparisons); percentiles are against the pathway-class-matched random-set null. Each report holds the "
        "figures, full tables and provenance.",
    ]
    links = []
    if lib_links:
        links += [f'<a href="library/{f}">{f}</a>' for f in lib_links]
    if (qdst / "report.html").exists():
        links.append(f'<a href="library/query_{QUERY_EXAMPLE}/report.html">library query for {QUERY_EXAMPLE}</a>')
    if links:
        intro.append(markupsafe.Markup(
            "Signature library (every MoTrPAC comparison as UP and DOWN gene sets, BH FDR &lt; 0.05, top 250 by "
            "|stat|): " + ", ".join(links) + "."))
    idx = render.render_site_index(rows, site / "index.html", title="motrpac_probe: disease signatures against MoTrPAC exercise data", intro=intro)

    problems = {}
    for h in sorted(site.rglob("*.html")):
        p = render.validate_report(h)
        if p:
            problems[str(h.relative_to(site))] = p
    for n in notes:
        print(n, file=sys.stderr)
    print(f"[site] {len(rows)} examples, {len(list(site.rglob('*.html')))} HTML files, "
          f"validation problems: {sum(len(v) for v in problems.values())}; {time.time() - t0:.1f} s", file=sys.stderr)
    for k, v in problems.items():
        for x in v:
            print(f"  {k}: {x}", file=sys.stderr)
    return idx, rows, notes, problems


def build_cli(a):
    idx, rows, notes, problems = build_site()
    print(idx)
