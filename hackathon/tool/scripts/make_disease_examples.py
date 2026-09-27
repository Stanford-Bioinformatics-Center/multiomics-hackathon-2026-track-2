#!/usr/bin/env python3
"""Generate the published non-PAH disease signatures in tool/examples/ and the
examples manifest tool/examples/examples.json.

Run from anywhere:  python tool/scripts/make_disease_examples.py
Reads only cached files in tool/examples/sources/ (see sources/RETRIEVAL.tsv):
  * Enrichr GMT libraries (downloaded 2026-09-27 from
    https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=<lib>)
  * msigdbr_cgp_subset.csv, written by tool/scripts/fetch_msigdbr_cgp.R (msigdbr 26.1.1)
No gene is typed in by hand. Every gene comes from one of those files.

Source choice (details in tool/examples/README.md):
  * Disease_Signatures_from_GEO_{up,down}_2014 was the first choice. It was rejected
    because its heart failure and septic shock sets fail a marker check (see
    sanity_report). For example, NPPA/NPPB are absent from both HF "up" sets and
    MMP8/CD177/S100A12 are absent from the septic shock "up" set, and the sets barely
    overlap MSigDB THUM_SYSTOLIC_HEART_FAILURE.
  * Disease_Perturbations_from_GEO_{up,down} and Aging_Perturbations_from_GEO_{up,down}
    (crowd-extracted GEO signatures, Wang et al. 2016 Nat Commun 7:12846, served by
    Enrichr, Kuleshov et al. 2016 NAR 44:W90) pass the same check and are used for the
    heart, blood and aging examples.
  * No explicit up/down T2D skeletal muscle signature exists in these sources. The only
    human T2D muscle entries in Disease_Perturbations are cultured myotubes (GSE12643),
    and they appear twice with up and down swapped. The T2D example is therefore
    MSigDB MOOTHA_VOXPHOS alone (direction -1).
"""
import json
import re
from pathlib import Path

import pandas as pd

TOOL = Path(__file__).resolve().parents[1]
OUT = TOOL / "examples"
SRC = OUT / "sources"

COLS = ["gene_symbol", "direction", "group", "weight", "source"]
ENRICHR = "Enrichr (Kuleshov 2016 NAR 44:W90)"
CREEDS = "crowd-extracted GEO signature (Wang 2016 Nat Commun 7:12846)"


def read_gmt(path):
    """{set name: [genes in file order]}. Enrichr lines are name<TAB><TAB>g1<TAB>g2...;
    a gene may carry ',weight', which is split off and returned separately."""
    sets, weights = {}, {}
    for line in open(path):
        parts = line.rstrip("\n").split("\t")
        genes, w = [], []
        for tok in parts[2:]:
            if not tok:
                continue
            g, _, wt = tok.partition(",")
            genes.append(g)
            w.append(float(wt) if wt else None)
        sets[parts[0]] = genes
        weights[parts[0]] = w
    return sets, weights


def gmt_pair(lib_fmt, name):
    up, _ = read_gmt(SRC / (lib_fmt % "up"))
    dn, _ = read_gmt(SRC / (lib_fmt % "down"))
    return up[name], dn[name]


ORF = re.compile(r"^(C(?:\d+|X|Y))ORF(\d+)$")


def fix_case(g):
    """Enrichr upper-cases open-reading-frame symbols (C2ORF74); HGNC writes C2orf74."""
    return ORF.sub(r"\1orf\2", g)


def frame(up, dn, group_up, group_dn, source):
    """Signature table. Genes listed in both directions are dropped (reported)."""
    up, dn = [fix_case(g) for g in up], [fix_case(g) for g in dn]
    both = set(up) & set(dn)
    up = list(dict.fromkeys(g for g in up if g not in both))
    dn = list(dict.fromkeys(g for g in dn if g not in both))
    rows = [(g, 1, group_up) for g in up] + [(g, -1, group_dn) for g in dn]
    df = pd.DataFrame(rows, columns=["gene_symbol", "direction", "group"])
    df["weight"] = ""
    df["source"] = source
    return df[COLS], sorted(both)


def write(df, name, dropped):
    df.to_csv(OUT / name, index=False)
    n = df.direction.value_counts()
    print(f"{name}: {n.get(1, 0)} up, {n.get(-1, 0)} down"
          + (f"; dropped {len(dropped)} genes listed in both directions" if dropped else ""))


# ---------------------------------------------------------------- signatures

def t2d_mootha():
    m = pd.read_csv(SRC / "msigdbr_cgp_subset.csv")
    v = m[m.gs_name == "MOOTHA_VOXPHOS"]
    pmid = v.gs_pmid.iloc[0]
    ver = v.msigdbr_version.iloc[0]
    src = (f"MSigDB C2:CGP MOOTHA_VOXPHOS (msigdbr {ver}; PMID {pmid}; Mootha 2003 Nat Genet "
           "34:267). OXPHOS gene set that GSEA found coordinately down in vastus lateralis of "
           "IGT/T2D vs NGT men (Mootha 2003 Nat Genet 34:267); all genes -1; no up set exists "
           "in the source")
    df = pd.DataFrame({"gene_symbol": v.gene_symbol.values, "direction": -1,
                       "group": "OXPHOS-down (MOOTHA_VOXPHOS)", "weight": "", "source": src})
    return df[COLS], []


def hf_hannenhalli():
    name = "Cardiomyopathy C0878544 human GSE5406 sample 88"
    up, dn = gmt_pair("Disease_Perturbations_from_GEO_%s.gmt", name)
    src = (f"{ENRICHR} Disease_Perturbations_from_GEO_up/_down '{name}'; {CREEDS}; GSE5406 "
           "(Hannenhalli 2006 Circulation 114:1269, human LV, ischemic/idiopathic "
           "cardiomyopathy vs non-failing); up/down = genes up/down in failing vs non-failing LV")
    return frame(up, dn, "HF-up", "HF-down", src)


def aging_liu():
    name = "Human biceps brachii muscles 24 years vs 70 years GDS4858 aging:33"
    up, dn = gmt_pair("Aging_Perturbations_from_GEO_%s.gmt", name)
    src = (f"{ENRICHR} Aging_Perturbations_from_GEO_up/_down '{name}'; {CREEDS}; GDS4858/GSE38718 "
           "(Liu 2013 J Gerontol A 68:1035, human biceps brachii, old 65-76 y vs young 19-28 y); "
           "up/down = genes up/down in old vs young")
    return frame(up, dn, "aging-up", "aging-down", src)


def sepsis_cvijanovich():
    name = "Septic Shock C0036983 human GSE9692 sample 307"
    up, dn = gmt_pair("Disease_Perturbations_from_GEO_%s.gmt", name)
    src = (f"{ENRICHR} Disease_Perturbations_from_GEO_up/_down '{name}'; {CREEDS}; GSE9692 "
           "(Cvijanovich 2008 Physiol Genomics 34:127, pediatric whole blood, septic shock day 1 "
           "vs normal controls); up/down = genes up/down in septic shock vs control")
    return frame(up, dn, "septic-shock-up", "septic-shock-down", src)


# ---------------------------------------------------------------- sanity check

MARKERS = {
    "heart failure": (["NPPA", "NPPB", "FRZB", "SFRP4", "ASPN"], ["SERPINA3", "FCN3", "MYH6"]),
    "septic shock": (["MMP8", "CD177", "S100A12", "ARG1", "OLFM4"], ["HLA-DRA", "PRF1", "NKG7", "KLRD1"]),
    "aging muscle": (["MT1X", "MT2A", "CDKN1A", "H19"], []),
}
CANDIDATES = [
    ("heart failure", "Disease_Signatures_from_GEO_%s_2014.gmt", "Cardiomyopathy, Dilated GSE3586"),
    ("heart failure", "Disease_Signatures_from_GEO_%s_2014.gmt", "Cardiomyopathy GSE5406"),
    ("heart failure", "Disease_Perturbations_from_GEO_%s.gmt",
     "Cardiomyopathy C0878544 human GSE5406 sample 88"),
    ("septic shock", "Disease_Signatures_from_GEO_%s_2014.gmt", "Septic Shock GSE9692"),
    ("septic shock", "Disease_Perturbations_from_GEO_%s.gmt",
     "Septic Shock C0036983 human GSE9692 sample 307"),
    ("aging muscle", "Aging_Perturbations_from_GEO_%s.gmt",
     "Human skeletal muscle (vastus lateralis) 24 years vs 71 years GDS288 aging:93"),
    ("aging muscle", "Aging_Perturbations_from_GEO_%s.gmt",
     "Human biceps brachii muscles 24 years vs 70 years GDS4858 aging:33"),
]


def sanity_report():
    """Print which textbook markers each candidate set contains, and its overlap with
    MSigDB THUM_SYSTOLIC_HEART_FAILURE_UP/DN. This is how the sources were chosen."""
    m = pd.read_csv(SRC / "msigdbr_cgp_subset.csv")
    thum_up = set(m.gene_symbol[m.gs_name == "THUM_SYSTOLIC_HEART_FAILURE_UP"])
    thum_dn = set(m.gene_symbol[m.gs_name == "THUM_SYSTOLIC_HEART_FAILURE_DN"])
    print("\nsanity check (markers expected up | down; overlap with THUM HF up/down):")
    for kind, lib, name in CANDIDATES:
        up, dn = map(set, gmt_pair(lib, name))
        mu, md = MARKERS[kind]
        extra = (f"  THUM up∩up={len(up & thum_up)} dn∩dn={len(dn & thum_dn)}"
                 if kind == "heart failure" else "")
        print(f"  [{lib.split('_%s')[0].replace('%s', '')}] {name}: n={len(up)}/{len(dn)} "
              f"up-markers {[g for g in mu if g in up]} down-markers {[g for g in md if g in dn]}"
              + extra)


# ---------------------------------------------------------------- manifest

MANIFEST = [
    dict(name="pah_muscle_malenfant2015", title="PAH, skeletal muscle proteome (Malenfant 2015)",
         kind="disease", tissue_hint="muscle",
         description="19 iTRAQ proteins plus 6 phenotype genes from vastus lateralis of IPAH patients vs controls.",
         run_args=["--pool-sets", "KEGG_OXIDATIVE_PHOSPHORYLATION,KEGG_CITRATE_CYCLE_TCA_CYCLE",
                   "--pool-groups", "OXPHOS-down,TCA/transport-down"]),
    dict(name="pah_blood_cheadle2012_eds", title="PAH, blood erythroid signature (Cheadle 2012)",
         kind="disease", tissue_hint="blood",
         description="22 erythroid development genes reported up in PAH peripheral blood (GSE33463).",
         run_args=[]),
    dict(name="hostrup2022_hiit_proteome", title="HIIT training, muscle proteome (Hostrup 2022)",
         kind="positive control", tissue_hint="muscle",
         description="Vastus lateralis proteins changed after 5 weeks of HIIT; should agree with MoTrPAC training.",
         run_args=[]),
    dict(name="random_matched", title="Random abundance-matched genes",
         kind="negative control", tissue_hint="muscle",
         description="25 random genes matched on expression decile to the PAH muscle signature; should score near the null.",
         run_args=[]),
    dict(name="random_mito9", title="Random MitoCarta genes marked down",
         kind="pathway-class control", tissue_hint="muscle",
         description="9 random mitochondrial genes all set to -1, showing what the pathway class alone scores.",
         run_args=[]),
    dict(name="tca_intermediates_demo", title="TCA intermediates (synthetic metabolite demo)",
         kind="metabolite demo", tissue_hint="muscle",
         description="Synthetic RefMet metabolite list that exercises the metabolomics path; not a published signature.",
         run_args=[]),
    dict(name="type2_diabetes_muscle_mootha2003", title="Type 2 diabetes, skeletal muscle (Mootha 2003)",
         kind="disease", tissue_hint="muscle",
         description="MSigDB MOOTHA_VOXPHOS, the OXPHOS gene set coordinately down in T2D/IGT vastus lateralis; down-only.",
         run_args=[]),
    dict(name="heart_failure_lv_hannenhalli2006", title="Heart failure, left ventricle (Hannenhalli 2006)",
         kind="disease", tissue_hint="heart",
         description="Genes up and down in failing vs non-failing human left ventricle (GSE5406), Enrichr crowd GEO signature.",
         run_args=[]),
    dict(name="aging_muscle_liu2013", title="Aging, skeletal muscle (Liu 2013)",
         kind="disease", tissue_hint="muscle",
         description="Genes up and down in old vs young human biceps brachii (GDS4858/GSE38718), Enrichr crowd GEO signature.",
         run_args=[]),
    dict(name="septic_shock_blood_cvijanovich2008", title="Septic shock, whole blood (Cvijanovich 2008)",
         kind="disease", tissue_hint="blood",
         description="Genes up and down in pediatric septic shock vs control whole blood (GSE9692), Enrichr crowd GEO signature.",
         run_args=[]),
]


def main():
    for fn, name in [(t2d_mootha, "type2_diabetes_muscle_mootha2003.csv"),
                     (hf_hannenhalli, "heart_failure_lv_hannenhalli2006.csv"),
                     (aging_liu, "aging_muscle_liu2013.csv"),
                     (sepsis_cvijanovich, "septic_shock_blood_cvijanovich2008.csv")]:
        df, dropped = fn()
        write(df, name, dropped)
    sanity_report()

    entries = []
    for e in MANIFEST:
        f = OUT / f"{e['name']}.csv"
        assert f.exists(), f
        entries.append({"name": e["name"], "file": f.name, **{k: e[k] for k in
                        ("title", "kind", "tissue_hint", "description", "run_args")}})
    listed = {e["file"] for e in entries}
    unlisted = sorted(p.name for p in OUT.glob("*.csv") if p.name not in listed)
    assert not unlisted, f"example files missing from manifest: {unlisted}"
    (OUT / "examples.json").write_text(json.dumps(entries, indent=2) + "\n")
    print(f"\nwrote examples.json ({len(entries)} entries)")


if __name__ == "__main__":
    main()
