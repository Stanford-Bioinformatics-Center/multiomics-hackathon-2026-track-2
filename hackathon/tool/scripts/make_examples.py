#!/usr/bin/env python3
"""Generate the example signature CSVs in tool/examples/.

Run from anywhere:  python tool/scripts/make_examples.py
Inputs are read (never written) from the hackathon directory (HACK = tool/..).
Output format: gene_symbol, direction (+1/-1), group, weight, source[, uniprot].
The metabolite demo uses refmet_name instead of gene_symbol.
"""
from pathlib import Path

import numpy as np
import pandas as pd

TOOL = Path(__file__).resolve().parents[1]
HACK = TOOL.parent
OUT = TOOL / "examples"
TEAMMATE = Path("/home/dani/.t3/worktrees/multiomics-hackathon-2026-track-2/"
                "t3code-0192f3ca/MoTrPAC Hackathon")

JOIN = HACK / "data" / "join_table_v2.csv"
FMAP = HACK / "data" / "raw" / "human_feature_to_gene_muscle_v1.csv"
RNA_DA = HACK / "data" / "raw" / "human_muscle_rna_DA_v1.csv"
HOSTRUP = HACK / "protonly" / "ext" / "hostrup2022_proteome_training.csv"
MITOCARTA = HACK / "model" / "ext" / "mitocarta3_human.csv"
MALENFANT_TM = TEAMMATE / "data" / "from paper" / "pah_lower_proteins_malenfant2015.csv"
BLOOD_TM = TEAMMATE / "outputs" / "PAH blood" / "pah_blood_gene_ranked.csv.gz"

COLS = ["gene_symbol", "direction", "group", "weight", "source"]

# Copied verbatim from hackathon/scripts/05_pah_grid.py (SIGNATURE); that script is
# not imported because it runs its pipeline at import time.
# (gene, PAH/control ratio, direction, group)
SIGNATURE = [
    ("NDUFA9", 0.71, -1, "OXPHOS-down"), ("UQCRC2", 0.74, -1, "OXPHOS-down"),
    ("UQCRC1", 0.74, -1, "OXPHOS-down"), ("ATP5MG", 0.70, -1, "OXPHOS-down"),
    ("ATP5F1B", 0.74, -1, "OXPHOS-down"),
    ("IDH2", 0.74, -1, "TCA/transport-down"), ("OGDH", 0.77, -1, "TCA/transport-down"),
    ("SLC25A4", 0.65, -1, "TCA/transport-down"), ("ECH1", 0.57, -1, "TCA/transport-down"),
    ("GLO1", 1.39, 1, "glycolytic-up"), ("FBP2", 1.33, 1, "glycolytic-up"),
    ("PEBP1", 1.29, 1, "contractile/other-up"), ("ATP2A1", 1.33, 1, "contractile/other-up"),
    ("MYH1", 2.42, 1, "contractile/other-up"), ("MYH2", 1.26, 1, "contractile/other-up"),
    ("MYH7", 1.26, 1, "contractile/other-up"), ("MYLPF", 1.58, 1, "contractile/other-up"),
    ("APOA1", 1.29, 1, "contractile/other-up"), ("PDLIM3", 1.58, 1, "contractile/other-up"),
    ("CS", None, -1, "phenotype"), ("TFAM", None, -1, "phenotype"),
    ("PDHA1", None, -1, "phenotype"), ("PDHB", None, -1, "phenotype"),
    ("LDHA", None, 1, "phenotype"), ("LDHB", None, 1, "phenotype"),
]
PAH_GENES = [s[0] for s in SIGNATURE]
PAH_MITO9 = ["NDUFA9", "UQCRC2", "UQCRC1", "ATP5MG", "ATP5F1B", "IDH2", "OGDH", "SLC25A4", "ECH1"]
PAPER_SYMBOL = {"ATP5L": "ATP5MG", "ATP5B": "ATP5F1B"}

EDS22 = ("ALAS2 AHSP CA1 EPB42 FECH GLRX5 GSPT1 GYPB GYPE HBA2 HBB HBD HBM HBQ1 BPGM "
         "MYL4 SELENBP1 SLC25A37 SNCA TMCC2 TSPAN5 GATA1").split()

RNA_24H_SHORT = "Endur.post_24_hr - Control.post_24_hr (delta-delta)"


def write(df, name):
    OUT.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT / name, index=False)
    return df


def join_slices():
    j = pd.read_csv(JOIN)
    h = j[j.dataset == "human_acute"]
    return j, {
        "hVL_RNA_24h": h[(h.tissue == "VL") & (h.layer == "RNA") & (h.contrast == "EE_vs_CON_24h")],
        "hVL_PROT_24h": h[(h.tissue == "VL") & (h.layer == "PROT") & (h.contrast == "EE_vs_CON_24h")],
        "rGN_PROT_F8w": j[(j.dataset == "rat_train") & (j.tissue == "SKM-GN")
                          & (j.layer == "PROT") & (j.contrast == "F_8w")],
    }


# ----------------------------------------------------------------------------- 1
def make_pah_muscle(j):
    notes = []
    fm = pd.read_csv(FMAP, usecols=["assay", "feature_id", "gene_symbol", "uniprot"])
    fm = fm[fm.assay == "prot-pr"]
    jp = j[(j.dataset == "human_acute") & (j.layer == "PROT")][["gene_symbol_human", "feature_id"]]
    cand = pd.concat([fm.rename(columns={"gene_symbol": "g"})[["g", "uniprot"]],
                      jp.rename(columns={"gene_symbol_human": "g", "feature_id": "uniprot"})])
    cand = cand.dropna().drop_duplicates()

    def pick(g):
        u = sorted(set(cand.loc[cand.g == g, "uniprot"]))
        canon = [x for x in u if "-" not in x]
        if len(canon) > 1:
            notes.append(f"{g}: several canonical accessions {canon}; left empty")
            return ""
        if canon:
            return canon[0]
        notes.append(f"{g}: no accession found")
        return ""

    rows = []
    for g, ratio, d, grp in SIGNATURE:
        pheno = grp == "phenotype"
        rows.append({
            "gene_symbol": g, "direction": d, "group": grp,
            "weight": "" if ratio is None else ratio,
            "source": ("Malenfant 2015 enzyme activity / blot (phenotype)" if pheno
                       else "Malenfant 2015 J Mol Med (iTRAQ VL, 4 IPAH vs 4 controls)"),
            "uniprot": "" if pheno else pick(g),
        })
    df = pd.DataFrame(rows, columns=COLS + ["uniprot"])

    # cross-check against teammate's table of the 9 down proteins
    if MALENFANT_TM.exists():
        tm = pd.read_csv(MALENFANT_TM)
        tm["g"] = tm.paper_symbol.replace(PAPER_SYMBOL)
        m = df.merge(tm, left_on="gene_symbol", right_on="g", how="inner")
        for _, r in m.iterrows():
            if r.uniprot != r.uniprot_accession:
                notes.append(f"UniProt disagreement {r.gene_symbol}: ours {r.uniprot} vs teammate {r.uniprot_accession}")
            if not np.isclose(float(r.weight), float(r.pah_to_control_ratio)):
                notes.append(f"ratio disagreement {r.gene_symbol}: ours {r.weight} vs teammate {r.pah_to_control_ratio}")
        notes.append(f"teammate cross-check: {len(m)}/{len(tm)} rows compared")
    else:
        notes.append("teammate Malenfant file not found; no cross-check")
    return write(df, "pah_muscle_malenfant2015.csv"), notes


# ----------------------------------------------------------------------------- 2
def make_blood():
    notes = []
    src = ("Cheadle 2012 (GSE33463) erythroid development signature, "
           "subset of the 149-gene EDS (Table S2)")
    w = {}
    if BLOOD_TM.exists():
        b = pd.read_csv(BLOOD_TM).drop_duplicates("gene_symbol").set_index("gene_symbol")
        w = {g: round(float(b.at[g, "logFC"]), 4) for g in EDS22 if g in b.index}
        src += ("; weight = IPAH-vs-healthy limma logFC in GSE33463 (30 IPAH, 41 healthy) "
                "from the team's re-analysis (pah_blood_gene_ranked)")
        notes.append(f"logFC found for {len(w)}/{len(EDS22)} genes")
    df = pd.DataFrame({"gene_symbol": EDS22, "direction": 1, "group": "erythroid (EDS subset)",
                       "weight": [w.get(g, "") for g in EDS22], "source": src}, columns=COLS)
    return write(df, "pah_blood_cheadle2012_eds.csv"), notes


# ----------------------------------------------------------------------------- 3
def make_hostrup():
    h = pd.read_csv(HOSTRUP)
    n_flag = int(h.significant_fdr05.sum())
    h = h[(h.q_value < 0.05) & h.gene_symbol.notna() & (h.gene_symbol.astype(str).str.strip() != "")]
    n_groups = len(h)
    h = h.sort_values(["q_value", "p_value"]).drop_duplicates("gene_symbol", keep="first")
    h = h[h.logFC != 0]
    d = np.sign(h.logFC).astype(int)
    df = pd.DataFrame({
        "gene_symbol": h.gene_symbol.values, "direction": d.values,
        "group": np.where(d > 0, "HIIT-up", "HIIT-down"),
        "weight": h.logFC.round(4).values,
        "source": "Hostrup 2022 eLife 11:e69802, VL proteome, n=8 men, 5 wk HIIT, q<0.05",
        "uniprot": h.leading_accession.values,
    }, columns=COLS + ["uniprot"]).sort_values(["direction", "weight"], ascending=[False, False])
    notes = [f"authors' significant_fdr05 flag: {n_flag}; q<0.05 with symbol: {n_groups} groups; "
             f"after dedup by symbol: {len(df)}"]
    return write(df, "hostrup2022_hiit_proteome.csv"), notes


# ----------------------------------------------------------------------------- 4
def make_random_matched(sl):
    rna = pd.read_csv(RNA_DA, usecols=["contrast_short", "feature_id", "AveExpr"])
    rna = rna[rna.contrast_short == RNA_24H_SHORT].drop_duplicates("feature_id")
    col = sl["hVL_RNA_24h"][["gene_symbol_human", "feature_id"]]
    ab = col.merge(rna[["feature_id", "AveExpr"]], on="feature_id", how="inner")
    ab = ab.dropna(subset=["AveExpr"]).drop_duplicates("gene_symbol_human")
    ab["decile"] = pd.qcut(ab.AveExpr, 10, labels=False)
    dec = ab.set_index("gene_symbol_human").decile
    pah_dec = [int(dec[g]) for g in PAH_GENES if g in dec.index]
    fallback = int(np.round(np.median(pah_dec)))
    missing = [g for g in PAH_GENES if g not in dec.index]

    rng = np.random.default_rng(20260927)
    pool = ab[~ab.gene_symbol_human.isin(PAH_GENES)]
    used, picks = set(), []
    for g, _, d, _ in SIGNATURE:
        k = int(dec[g]) if g in dec.index else fallback
        c = sorted(set(pool.loc[pool.decile == k, "gene_symbol_human"]) - used)
        r = c[rng.integers(len(c))]
        used.add(r)
        picks.append((r, d, g, k))
    df = pd.DataFrame({
        "gene_symbol": [p[0] for p in picks], "direction": [p[1] for p in picks],
        "group": "random (abundance-matched)", "weight": "",
        "source": ("random draw, seed 20260927, matched on human VL RNA AveExpr decile "
                   "to pah_muscle_malenfant2015"),
    }, columns=COLS)
    notes = [f"{len(ab)} genes with AveExpr in hVL RNA EE_vs_CON_24h; PAH genes not measured "
             f"there: {missing} -> decile {fallback} (median PAH decile)",
             "pairs (random <- PAH, decile): " + ", ".join(f"{p[0]}<-{p[2]}({p[3]})" for p in picks)]
    return write(df, "random_matched.csv"), notes


# ----------------------------------------------------------------------------- 5
def make_random_mito9(sl):
    mc = set(pd.read_csv(MITOCARTA).Symbol.dropna())
    keep = mc.copy()
    for k in ["hVL_RNA_24h", "hVL_PROT_24h", "rGN_PROT_F8w"]:
        keep &= set(sl[k].gene_symbol_human)
    keep -= set(PAH_MITO9) | {"CS", "TFAM", "PDHA1", "PDHB"}
    pool = sorted(keep)
    rng = np.random.default_rng(20260928)
    pick = list(rng.choice(pool, 9, replace=False))
    df = pd.DataFrame({
        "gene_symbol": pick, "direction": -1, "group": "random MitoCarta", "weight": "",
        "source": ("random draw of 9 MitoCarta3.0 genes, seed 20260928; all marked -1 "
                   "(down in 'disease') on purpose"),
    }, columns=COLS)
    return write(df, "random_mito9.csv"), [f"eligible MitoCarta pool: {len(pool)} genes"]


# ----------------------------------------------------------------------------- 6
def make_tca():
    # RefMet names verified against MotrpacRatTraining6moData::METAB_FEATURE_ID_MAP
    # (metabolite_refmet): alpha-ketoglutarate = "Oxoglutaric acid", acetylcarnitine = "CAR(2:0)".
    mets = [("Citric acid", -1), ("Isocitric acid", -1), ("cis-Aconitic acid", -1),
            ("Oxoglutaric acid", -1), ("Succinic acid", -1), ("Fumaric acid", -1),
            ("Malic acid", -1), ("Pyruvic acid", -1), ("Lactic acid", 1), ("CAR(2:0)", -1)]
    df = pd.DataFrame({
        "refmet_name": [m[0] for m in mets], "direction": [m[1] for m in mets],
        "group": "TCA intermediates (synthetic demo)", "weight": "",
        "source": "SYNTHETIC demo, not a published signature",
    }, columns=["refmet_name", "direction", "group", "weight", "source"])
    return write(df, "tca_intermediates_demo.csv"), []


def main():
    j, sl = join_slices()
    results = {
        "pah_muscle_malenfant2015.csv": make_pah_muscle(j),
        "pah_blood_cheadle2012_eds.csv": make_blood(),
        "hostrup2022_hiit_proteome.csv": make_hostrup(),
        "random_matched.csv": make_random_matched(sl),
        "random_mito9.csv": make_random_mito9(sl),
        "tca_intermediates_demo.csv": make_tca(),
    }
    for name, (df, notes) in results.items():
        up, dn = int((df.direction == 1).sum()), int((df.direction == -1).sum())
        line = f"{name}: n={len(df)} (+1: {up}, -1: {dn})"
        if "gene_symbol" in df:
            g = set(df.gene_symbol)
            line += "; in join: " + ", ".join(
                f"{k}={len(g & set(v.gene_symbol_human))}" for k, v in sl.items())
        print(line)
        for n in notes:
            print("    " + n)


if __name__ == "__main__":
    main()
