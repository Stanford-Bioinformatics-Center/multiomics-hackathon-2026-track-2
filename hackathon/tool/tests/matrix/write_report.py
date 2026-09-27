#!/usr/bin/env python3
"""Write tests/REPORT.md from tests/matrix/results.csv + seed_compare.csv (numbers only, no free text per run)."""
from pathlib import Path

import pandas as pd

M = Path(__file__).resolve().parent
res = pd.read_csv(M / "results.csv")
seed = pd.read_csv(M / "seed_compare.csv")
FIXED = M / "FIXED.md"   # hand-maintained log of bugs the matrix found and how they were fixed


def md(df, floatfmt=".2f"):
    """Markdown table without the optional tabulate dependency."""
    f = lambda v: "" if pd.isna(v) else (format(v, floatfmt) if isinstance(v, float) else str(v).replace("|", "/"))
    out = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    out += ["| " + " | ".join(f(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(out)


ok = res[res.exit == 0]
lines = ["# Robustness matrix report", "",
         f"`tests/matrix/run_matrix.py`: {len(res)} cases, run with `mprobe run` (1,000 random sets per null) in "
         f"parallel on the desktop node. Assertions: `tests/test_matrix.py`. Metabolite cases: "
         "`tests/test_matrix_metab.py`.", "",
         "## Summary", "",
         f"- Exit 0: {len(ok)} / {len(res)}.",
         f"- Reports that validate (no broken refs, captions present, no empty tables): "
         f"{int((ok.validation_problems.fillna('') == '').sum())} / {len(ok)}.",
         f"- Runtime per run: median {res.runtime_s.median():.1f} s, max {res.runtime_s.max():.1f} s (limit 90 s).",
         ""]
g = ok.groupby("group").agg(cases=("case", "size"), median_runtime_s=("runtime_s", "median"),
                            median_class_pct=("median_pct_class", "median"),
                            cases_with_setlevel_verdicts=("n_opposed_fdr", lambda s: int((s + ok.loc[s.index, "n_same_fdr"] > 0).sum())))
lines += ["## By group", "", md(g.reset_index(), ".1f"), ""]
cols = ["case", "exit", "runtime_s", "n_rows", "n_mapped", "n_unmapped", "n_columns_no_genes", "n_opposed_fdr",
        "n_same_fdr", "median_pct_abund", "median_pct_class", "rat_gn_prot_f8_t", "n_flips"]
lines += ["## All cases", "",
          "n_opposed_fdr / n_same_fdr = core comparisons (of 52) with cameraPR BH FDR < 0.05 opposed / same "
          "direction; median_pct_* = median null percentile over the 52; rat_gn_prot_f8_t = signed cameraPR t in "
          "rat gastrocnemius protein, female 8 wk; n_flips = sensitivity-table conclusions that change under an "
          "alternative setting.", "", md(res[cols], ".2f"), ""]
lines += ["## Seed reproducibility (seed 20260926 vs 12345)", "",
          "Non-null quantities must be identical; null percentiles vary with the random draws.", "",
          md(seed, ".2f"), ""]
fail = res[res.exit != 0]
lines += ["## Failures", "", md(fail[["case", "stderr_tail"]]) if len(fail) else "None.", ""]
if FIXED.exists():
    lines += ["## Found and fixed", "", FIXED.read_text().strip(), ""]
(M.parent / "REPORT.md").write_text("\n".join(lines) + "\n")
print(M.parent / "REPORT.md")
