#!/usr/bin/env python3
"""Robustness matrix for `mprobe run`: generate inputs, run them in parallel, collect results.

  python tests/matrix/run_matrix.py --workers 12 --out /scratch/mprobe_matrix     (desktop, via Slurm)

Writes tests/matrix/inputs/*.csv (deterministic, seed 7), tests/matrix/results.csv (one row per case: exit code,
runtime, validation problems, key numbers) and tests/matrix/seed_compare.csv. tests/test_matrix.py asserts on those;
tests/REPORT.md is written by write_report(). Metabolite cases live in tests/test_matrix_metab.py.
"""
import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TOOL = HERE.parents[1]
INP = HERE / "inputs"
SIZES = [5, 25, 100, 500]
POOLS = {"random": None, "mito": "mitocarta3", "secreted": "secreted_proxy", "contractile": "go_contractile_fiber",
         "ribosomal": "GOCC_RIBOSOME"}


def _universe():
    """Genes measured in human VL RNA 24 h (so 'random' means measurable), with annotations."""
    from motrpac_probe import core
    S = core.get_store()
    f = S.frame("human_acute|VL|RNA|EE_vs_CON_24h")
    g = sorted(set(f.gene_symbol_human.astype(str)))
    return S, g


def make_inputs():
    INP.mkdir(parents=True, exist_ok=True)
    S, uni = _universe()
    rng = np.random.default_rng(7)
    ann = S.ann.reindex(uni)
    cases = []

    def write(name, df, group):
        df.to_csv(INP / f"{name}.csv", index=False)
        cases.append(dict(case=name, group=group, file=f"tests/matrix/inputs/{name}.csv", args=[]))

    def dirs(n, kind):
        if kind == "up":
            return np.ones(n, int)
        if kind == "down":
            return -np.ones(n, int)
        return np.where(np.arange(n) % 2 == 0, 1, -1)

    for pool, flag in POOLS.items():
        if flag is None:
            cand = uni
        elif flag in ann.columns:
            cand = [g for g, v in zip(uni, ann[flag].fillna(False)) if v]
        else:
            cand = sorted(set(S.gs[flag]) & set(uni))
        for n in SIZES:
            k = min(n, len(cand))
            genes = list(rng.choice(cand, k, replace=False))
            for dk in (["half", "up", "down"] if pool == "random" else ["half"]):
                write(f"size{n}_{pool}_{dk}", pd.DataFrame({"gene_symbol": genes, "direction": dirs(k, dk)}),
                      "size x composition")
    # ID types (25 random genes)
    base = list(rng.choice(uni, 25, replace=False))
    d = dirs(25, "half")
    idm = pd.read_parquet(TOOL / "store" / "id_map.parquet")
    uni_map = idm[idm.id_type == "uniprot"].drop_duplicates("gene_symbol_human").set_index("gene_symbol_human").id
    ens = idm[(idm.id_type == "ensembl") & idm.id.str.startswith("ENSG")].drop_duplicates("gene_symbol_human")
    ens = ens.set_index("gene_symbol_human").id
    rat = idm[idm.id_type == "rat_symbol"].drop_duplicates("gene_symbol_human").set_index("gene_symbol_human").id
    write("id_human_symbols", pd.DataFrame({"gene_symbol": base, "direction": d}), "id types")
    write("id_uniprot_column", pd.DataFrame({"uniprot": [uni_map.get(g, "") for g in base], "direction": d}), "id types")
    write("id_uniprot_in_symbol", pd.DataFrame({"gene_symbol": [uni_map.get(g, g) for g in base], "direction": d}),
          "id types")
    write("id_ensembl", pd.DataFrame({"gene_symbol": [ens.get(g, g) for g in base], "direction": d}), "id types")
    write("id_rat_symbols", pd.DataFrame({"gene_symbol": [rat.get(g, g) for g in base], "direction": d}), "id types")
    write("id_mixed_case", pd.DataFrame({"gene_symbol": [g.lower() if i % 2 else g.title() for i, g in enumerate(base)],
                                         "direction": d}), "id types")
    write("id_duplicates", pd.DataFrame({"gene_symbol": base + base[:5] + base[5:8],
                                         "direction": list(d) + list(d[:5]) + list(-d[5:8])}), "id types")
    unk = [f"NOTAGENE{i}" for i in range(6)]
    write("id_20pct_unknown", pd.DataFrame({"gene_symbol": base[:24] + unk, "direction": list(d[:24]) + [1] * 6}),
          "id types")
    prot = set(S.frame("human_acute|VL|PROT|EE_vs_CON_24h").gene_symbol_human.astype(str))
    noprot = [g for g in uni if g not in prot]
    write("absent_from_proteomics", pd.DataFrame({"gene_symbol": list(rng.choice(noprot, 25, replace=False)),
                                                  "direction": d}), "edge cases")
    heart = set(S.frame("rat_train|HEART|PROT|F_8w").gene_symbol_human.astype(str))
    or_genes = [g for g in S.genes.astype(str) if g.startswith("OR") and g not in heart and g not in prot][:8]
    write("empty_overlap_tissue", pd.DataFrame({"gene_symbol": or_genes, "direction": dirs(len(or_genes), "half")}),
          "edge cases")
    write("weights_present", pd.DataFrame({"gene_symbol": base, "direction": d,
                                           "weight": np.round(rng.normal(0, 1, 25), 3)}), "weights")
    write("weights_absent", pd.DataFrame({"gene_symbol": base, "direction": d, "weight": [""] * 25}), "weights")
    write("direction_words", pd.DataFrame({"gene_symbol": base, "direction": np.where(d > 0, "up", "down")}),
          "id types")
    # real signatures and controls
    ex = json.loads((TOOL / "examples" / "examples.json").read_text())
    for e in ex:
        if e.get("kind") == "metabolite demo":
            continue
        cases.append(dict(case=f"example_{e['name']}", group=f"example: {e['kind']}", file=f"examples/{e['file']}",
                          args=e.get("run_args", [])))
    # seed reproducibility
    for e in ("pah_muscle_malenfant2015", "random_matched"):
        f = next(x for x in ex if x["name"] == e)
        cases.append(dict(case=f"seed2_{e}", group="seed", file=f"examples/{f['file']}",
                          args=f.get("run_args", []) + ["--seed", "12345"]))
    (HERE / "cases.json").write_text(json.dumps(cases, indent=1))
    return cases


def run_case(c, out, py):
    od = Path(out) / c["case"]
    t = time.time()
    p = subprocess.run([py, "-m", "motrpac_probe.cli", "run", "--quiet", "--signature", str(TOOL / c["file"]),
                        "--name", c["case"], "--out", str(od)] + c["args"], capture_output=True, text=True, cwd=TOOL)
    rt = time.time() - t
    rec = dict(case=c["case"], group=c["group"], exit=p.returncode, runtime_s=round(rt, 1),
               stderr_tail=p.stderr.strip().splitlines()[-1][:300] if p.stderr.strip() else "")
    if p.returncode == 0:
        from motrpac_probe import render
        rec["validation_problems"] = "; ".join(render.validate_report(od / "report.html"))
        mp = pd.read_csv(od / "tables" / "mapping.csv", dtype=str)
        sc = pd.read_csv(od / "tables" / "column_scores.csv")
        rec.update(n_rows=len(mp), n_mapped=int((mp.status == "mapped").sum()),
                   n_unmapped=int((mp.status != "mapped").sum()),
                   n_columns_no_genes=int((sc.n_measured == 0).sum()),
                   n_opposed_fdr=int(((sc.camera_fdr < 0.05) & (sc.camera_t > 0)).sum()),
                   n_same_fdr=int(((sc.camera_fdr < 0.05) & (sc.camera_t < 0)).sum()),
                   median_pct_abund=sc.pct_abund.median(), median_pct_class=sc.pct_class.median(),
                   frac_pct_class_outside_2_98=float((~sc.pct_class.dropna().between(2, 98)).mean())
                   if sc.pct_class.notna().any() else np.nan,
                   rat_gn_prot_f8_t=sc.set_index("column_id").camera_t.get("rat_train|SKM-GN|PROT|F_8w", np.nan),
                   n_flips=int(pd.read_csv(od / "tables" / "sensitivity.csv").flipped.sum()),
                   headline_has_no_genes_cell=int("no measured genes" in (od / "report.html").read_text()))
    return rec


def seed_compare(out):
    rows = []
    for e in ("pah_muscle_malenfant2015", "random_matched"):
        a = pd.read_csv(Path(out) / f"example_{e}" / "tables" / "column_scores.csv").set_index("column_id")
        b = pd.read_csv(Path(out) / f"seed2_{e}" / "tables" / "column_scores.csv").set_index("column_id")
        for col, kind in [("camera_t", "non-null"), ("n_opposed", "non-null"), ("sign_p", "non-null"),
                          ("pct_abund", "null-based"), ("pct_class", "null-based")]:
            d = (a[col] - b[col]).abs()
            rows.append(dict(signature=e, quantity=col, kind=kind, max_abs_diff=float(d.max()),
                             median_abs_diff=float(d.median())))
    df = pd.DataFrame(rows)
    df.to_csv(HERE / "seed_compare.csv", index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="/tmp/mprobe_matrix")
    a = ap.parse_args()
    cases = make_inputs()
    print(f"{len(cases)} cases, {a.workers} workers", flush=True)
    t = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        recs = list(ex.map(run_case, cases, [a.out] * len(cases), [sys.executable] * len(cases)))
    res = pd.DataFrame(recs)
    res.to_csv(HERE / "results.csv", index=False)
    seed_compare(a.out)
    print(f"done in {time.time() - t:.0f} s; failures: {int((res.exit != 0).sum())}", flush=True)
    print(res[["case", "exit", "runtime_s", "n_mapped", "n_unmapped", "n_opposed_fdr", "median_pct_class"]].to_string())


if __name__ == "__main__":
    main()
