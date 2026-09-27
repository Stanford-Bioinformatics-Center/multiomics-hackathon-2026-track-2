#!/usr/bin/env python3
"""Rebuild everything from the read-only hackathon inputs: store -> example runs -> compare -> library -> discord
-> site. Driven by examples/examples.json. Usage: python scripts/build_all.py [--skip-store] [--only examples|site]."""
import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

TOOL = Path(__file__).resolve().parents[1]
MPROBE = [sys.executable, "-m", "motrpac_probe.cli"]


def sh(args):
    t = time.time()
    print("+ mprobe " + " ".join(args), flush=True)
    subprocess.run(MPROBE + args, check=True, cwd=TOOL)
    print(f"  {time.time() - t:.1f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-store", action="store_true")
    ap.add_argument("--only", choices=["store", "examples", "site"])
    a = ap.parse_args()
    ex = json.loads((TOOL / "examples" / "examples.json").read_text(encoding="utf-8"))
    genes = [e for e in ex if e.get("kind") != "metabolite demo"]
    if not a.skip_store and a.only in (None, "store"):
        sh(["store", "build"])
        try:
            subprocess.run([sys.executable, "-c", "from motrpac_probe import store; store.build_metab()"], check=True,
                           cwd=TOOL)
        except Exception as e:  # metabolomics optional
            print(f"METAB layer not built: {e}")
    if a.only in (None, "examples"):
        for e in ex:  # gene and metabolite signatures alike (run detects metabolite files)
            sh(["run", "--quiet", "--signature", f"examples/{e['file']}", "--name", e["name"]] + e.get("run_args", []))
        sh(["compare", "--signatures"] + [f"examples/{e['file']}" for e in genes])
        sh(["library", "build"])
        q = next((e for e in genes if e.get("query")), next(e for e in genes if e.get("kind") == "disease"))
        sh(["library", "query", "--signature", f"examples/{q['file']}"])
        for sp, t in [("rat", "SKM-GN"), ("rat", "HEART"), ("human", "VL")]:
            sh(["discord", "--quiet", "--species", sp, "--tissue", t])
    if a.only in (None, "site"):
        sh(["site", "build"])


if __name__ == "__main__":
    main()
