"""mprobe command line. Every analysis is also available as a Python function (motrpac_probe.run, .library)."""
import argparse
import shlex
import sys
import time


def _store_build(a):
    from . import store
    store.build(hash_inputs=not a.no_hash)


def _csv(s):
    return tuple(x.strip() for x in s.split(",") if x.strip()) if s else ()


def _run(a):
    from . import run
    kw = dict(cutoff=a.cutoff, cap=0 if a.no_cap else a.cap, nboot=a.nboot, seed=a.seed,
              context_groups=_csv(a.context_groups) if a.context_groups is not None else None, min_n=a.min_n, universe=a.universe,
              tissues=_csv(a.tissues) or None, pool_sets=_csv(a.pool_sets), pool_groups=_csv(a.pool_groups),
              sections=_csv(a.sections) or None, quiet=a.quiet)
    kw = {k: v for k, v in kw.items() if v is not None and v != ()}
    if a.config:
        from .config import load_config
        kw = {**load_config(a.config, "run"), **kw}
    cmd = "mprobe " + " ".join(shlex.quote(x) for x in sys.argv[1:])
    p, _ = run.main_run(a.signature, name=a.name, outdir=a.out, command=cmd, **kw)
    print(p)


def _compare(a):
    import json
    from pathlib import Path
    from . import compare
    titles, kinds = {}, {}
    for p in a.signatures:  # optional titles / kinds from examples.json next to the files
        m = Path(p).resolve().parent / "examples.json"
        if m.exists():
            for e in json.loads(m.read_text()):
                titles[e["name"]] = e.get("title", e["name"])
                kinds[e["name"]] = e.get("kind", "disease")
    cmd = "mprobe " + " ".join(shlex.quote(x) for x in sys.argv[1:])
    p, _ = compare.main_compare(a.signatures, titles=titles, outdir=a.out, command=cmd, nboot=a.nboot, kinds=kinds)
    print(p)


def _discord(a):
    from . import discord
    cmd = "mprobe " + " ".join(shlex.quote(x) for x in sys.argv[1:])
    cfg = {}
    if a.config:
        from .config import load_config
        cfg = load_config(a.config, "discord")
    p, _ = discord.main_discord(a.species, a.tissue, signature=a.signature,
                                cutoff=a.cutoff if a.cutoff is not None else cfg.get("cutoff", 0.05),
                                layers_=_csv(a.layers) if a.layers else tuple(cfg.get("layers", ())) or None,
                                outdir=a.out, command=cmd, quiet=a.quiet,
                                universe=a.universe or cfg.get("universe", "all"))
    print(p)


def main(argv=None):
    p = argparse.ArgumentParser(prog="mprobe", description="Probe MoTrPAC exercise contrasts with a disease signature.")
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("store", help="build the parquet store from the hackathon CSVs")
    ss = s.add_subparsers(dest="sub", required=True)
    b = ss.add_parser("build")
    b.add_argument("--no-hash", action="store_true", help="skip sha256 of input files")
    b.set_defaults(fn=_store_build)

    r = sub.add_parser("run", help="probe MoTrPAC with one signature -> out/NAME/report.html")
    r.add_argument("--signature", required=True, help="CSV: gene_symbol|uniprot|ensembl|rat_symbol, direction, ...")
    r.add_argument("--name", help="run name (default: file stem)")
    r.add_argument("--out", help="output directory (default: tool/out/NAME)")
    r.add_argument("--config", help="YAML/JSON file with any of these options (CLI flags win)")
    r.add_argument("--cutoff", type=float, help="BH FDR cutoff for dots and verdicts (0.05)")
    r.add_argument("--cap", type=float, help="|stat| cap for grid colours (4)")
    r.add_argument("--no-cap", action="store_true", help="do not cap |stat| in the grid")
    r.add_argument("--nboot", type=int, help="random sets per null (1000)")
    r.add_argument("--seed", type=int, help="random seed (20260926)")
    r.add_argument("--context-groups", default=None,
                   help="signature groups shown but not counted (comma list; default 'phenotype')")
    r.add_argument("--tissues", default=None, help="core tissues for sections B-E")
    r.add_argument("--universe", choices=["all", "muscle-intrinsic"], default=None,
                   help="cameraPR background: all genes, or MitoCarta + contractile fibre + SR genes")
    r.add_argument("--min-n", type=int, help="min measured cells to rank a unit in section F (10)")
    r.add_argument("--pool-sets", default="", help="optional fixed-pool null: gene-set names (comma list)")
    r.add_argument("--pool-groups", default="", help="signature groups used for the fixed-pool null")
    r.add_argument("--sections", default="", help="only render these sections (comma list)")
    r.add_argument("--quiet", action="store_true")
    r.set_defaults(fn=_run)

    lib = sub.add_parser("library", help="MoTrPAC as a signature library (GMT + JSON index)")
    ls = lib.add_subparsers(dest="sub", required=True)
    lb = ls.add_parser("build", help="export every contrast as up/down gene sets")
    lb.add_argument("--cutoff", type=float, default=0.05)
    lb.add_argument("--cap", type=int, default=250, help="max genes per set, ranked by |stat| (250)")
    lb.set_defaults(fn=lambda a: __import__("motrpac_probe.library", fromlist=["x"]).build_cli(a))
    lq = ls.add_parser("query", help="rank MoTrPAC contrasts by opposition / concordance with a signature")
    lq.add_argument("--signature", required=True)
    lq.add_argument("--top", type=int, default=25)
    lq.add_argument("--out", help="output directory for CSVs + report.html (default: tool/out/library_query_NAME)")
    lq.set_defaults(fn=lambda a: __import__("motrpac_probe.library", fromlist=["x"]).query_cli(a))

    c = sub.add_parser("compare", help="several signatures side by side -> out/compare/report.html")
    c.add_argument("--signatures", nargs="+", required=True)
    c.add_argument("--out", help="output directory (default: tool/out/compare)")
    c.add_argument("--nboot", type=int, default=1000)
    c.set_defaults(fn=_compare)

    d = sub.add_parser("discord", help="omic-discordance report for one species x tissue -> out/discord_TISSUE/")
    d.add_argument("--species", choices=["rat", "human"], required=True)
    d.add_argument("--tissue", required=True, help="store tissue code, e.g. SKM-GN, HEART, VL")
    d.add_argument("--signature", help="optional signature CSV: agreement panels also restricted to it")
    d.add_argument("--layers", default=None, help="RNA,PROT,PHOSPHO,METAB")
    d.add_argument("--cutoff", type=float, help="BH FDR cutoff (0.05)")
    d.add_argument("--universe", choices=["all", "muscle-intrinsic"], default=None)
    d.add_argument("--out", help="output directory (default: tool/out/discord_TISSUE)")
    d.add_argument("--config", help="YAML/JSON file with any of these options (CLI flags win)")
    d.add_argument("--quiet", action="store_true")
    d.set_defaults(fn=_discord)

    st = sub.add_parser("site", help="static gallery of the example reports -> tool/site/")
    sts = st.add_subparsers(dest="sub", required=True)
    sts.add_parser("build").set_defaults(fn=lambda a: __import__("motrpac_probe.site", fromlist=["x"]).build_cli(a))

    a = p.parse_args(argv)
    t0 = time.time()
    a.fn(a)
    print(f"[mprobe] {a.cmd} done in {time.time() - t0:.1f} s", file=sys.stderr)


if __name__ == "__main__":
    main()
