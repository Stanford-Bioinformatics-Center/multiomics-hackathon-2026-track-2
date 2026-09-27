"""Build and load the MoTrPAC contrast store (tool/store/*.parquet). Schema: tool/store/SCHEMA.md.

Rows = one gene per comparison column (dataset x tissue x layer x contrast), exactly as scripts/03_join.py and
scripts/10_everywhere.py define them:
  * columns present in data/join_table_v2.csv are copied from it verbatim (source = "join_table_v2");
  * every other column (19 rat tissues, human muscle/adipose/blood, all contrasts incl. non-exercise references,
    PHOSPHO) is collapsed from the raw DA tables with the `bh` and `collapse` functions lifted from
    10_everywhere.py (source = "alltissue_v1"): BH over all tested features before collapsing, rat RNA on/off
    filter, max-|stat| feature per human gene.
"""
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

from . import legacy
from .paths import (ANNOT, COLUMNS, CONTRASTS, EVERY_RAW, FIGTAB, HACK, JOIN, MODEL_EXT, PROVENANCE, RAW, STORE)

GENESETS = STORE / "genesets.parquet"
IDMAP = STORE / "id_map.parquet"

H_TIME_CODE = {"during_20_min": "during20min", "during_40_min": "during40min", "post_10_min": "10min",
               "post_15_30_45_min": "15to45min", "post_3.5_4_hr": "3.5to4h", "post_24_hr": "24h",
               "pre_exercise": "pre"}
H_TISSUE = {"muscle": "VL", "adipose": "ADIPOSE", "blood": "BLOOD"}
JOIN_TISSUE_LAYERS = {("rat_train", "SKM-GN", "RNA"), ("rat_train", "SKM-GN", "PROT"), ("rat_train", "SKM-VL", "RNA"),
                      ("rat_train", "HEART", "RNA"), ("rat_train", "HEART", "PROT")}
# matched human tissue in GTEx v8 (CFDE) for each store tissue; None = no reasonable GTEx match
GTEX_TISSUE = {"VL": "Muscle - Skeletal", "SKM-GN": "Muscle - Skeletal", "SKM-VL": "Muscle - Skeletal",
               "HEART": "Heart - Left Ventricle", "WAT-SC": "Adipose - Subcutaneous",
               "ADIPOSE": "Adipose - Subcutaneous", "LIVER": "Liver", "KIDNEY": "Kidney - Cortex", "LUNG": "Lung",
               "CORTEX": "Brain - Cortex", "ADRNL": "Adrenal Gland", "BLOOD": "Whole Blood",
               "COLON": "Colon - Transverse", "HIPPOC": "Brain - Hippocampus", "HYPOTH": "Brain - Hypothalamus",
               "OVARY": "Ovary", "SMLINT": "Small Intestine - Terminal Ileum", "SPLEEN": "Spleen", "TESTES": "Testis",
               "BAT": None, "VENACV": None}
CAT_COLS = ["column_id", "dataset", "species", "tissue", "layer", "contrast", "gene_symbol_human", "feature_id",
            "stat_type", "source"]
ROW_COLS = ["column_id", "dataset", "species", "tissue", "layer", "contrast", "gene_symbol_human", "feature_id",
            "logFC", "stat", "stat_type", "p", "fdr", "fdr_bh", "n_collapsed", "training_fdr", "baseline_expr",
            "prot_n_missing", "source"]


def human_contrast_code(short):
    """'Endur.post_24_hr - Control.post_24_hr (delta-delta)' -> 'EE_vs_CON_24h' (03_join.py names for its six)."""
    lhs, rhs = [s.strip() for s in short.replace(" (delta-delta)", "").split(" - ")]
    g1, t1 = lhs.split(".", 1)
    g2, t2 = rhs.split(".", 1)
    G = legacy.GRP
    g1, g2, c1, c2 = G[g1], G[g2], H_TIME_CODE[t1], H_TIME_CODE[t2]
    if c1 == c2 == "pre":
        return f"{g1}_vs_{g2}_pre"
    if g1 == g2:
        return f"{g1}_post_vs_pre_{c1}"
    return f"{g1}_vs_{g2}_{c1}"


def column_id(dataset, tissue, layer, contrast):
    return f"{dataset}|{tissue}|{layer}|{contrast}"


def sha256(path, n=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while b := f.read(n):
            h.update(b)
    return h.hexdigest()


def _finish(d, source):
    d = d.copy()
    d["species"] = np.where(d.dataset == "human_acute", "human", "rat")
    d["column_id"] = d.dataset + "|" + d.tissue + "|" + d.layer + "|" + d.contrast
    d["source"] = source
    for c in ROW_COLS:
        if c not in d:
            d[c] = np.nan
    return d[ROW_COLS]


# ------------------------------------------------------------------------------------------------- rat
def _rat_fmap():
    rf2g = pd.read_csv(RAW / "rat_alltissue_feature_to_gene_v1.csv", usecols=["feature_ID", "gene_symbol"],
                       low_memory=False).dropna().drop_duplicates()
    orth = pd.read_csv(RAW / "rat_to_human_gene_v1.csv", usecols=["RAT_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL"],
                       low_memory=False).dropna().drop_duplicates()
    fmap = rf2g.merge(orth, left_on="gene_symbol", right_on="RAT_SYMBOL")[["feature_ID", "HUMAN_ORTHOLOG_SYMBOL"]]
    return fmap.drop_duplicates()


def build_rat(log):
    """All rat tissues x RNA/PROT/PHOSPHO from the alltissue dumps (10_everywhere.py logic, all genes).
    Returns (rows for columns NOT in the join table, rows for the join-table tissue x layers for cross-checking),
    and the per-feature covariates needed to annotate the join rows."""
    fmap = _rat_fmap()
    keep, overlap, cov = [], [], []
    for assay, layer, statcol in [("trnscrpt", "RNA", "zscore"), ("prot", "PROT", "tscore"),
                                  ("phospho", "PHOSPHO", "tscore")]:
        extra = (["comparison_average_intensity", "reference_average_intensity"] if layer == "RNA" else ["numNAs"])
        d = pd.read_csv(RAW / f"rat_alltissue_{assay}_DA_v1.csv", low_memory=False,
                        usecols=["tissue", "feature_ID", "sex", "comparison_group", "logFC", statcol, "p_value",
                                 "adj_p_value", "selection_fdr"] + extra)
        d["contrast"] = d.sex.str[0].str.upper() + "_" + d.comparison_group
        d["fdr_bh"] = d.groupby(["tissue", "contrast"]).p_value.transform(legacy.bh)
        if layer == "RNA":
            low = d[["comparison_average_intensity", "reference_average_intensity"]].min(axis=1)
            d = d[~((low < 1) | (d.logFC.abs() > 5))]
            d["baseline_expr"] = np.log2(d.reference_average_intensity + 1)
        else:
            d["prot_n_missing"] = d.numNAs
        cov.append(d[["tissue", "feature_ID", "contrast"] + (["baseline_expr"] if layer == "RNA" else ["prot_n_missing"])]
                   .assign(layer=layer))
        d = d.merge(fmap, on="feature_ID")
        d = d.rename(columns={"feature_ID": "feature_id", statcol: "stat", "p_value": "p", "adj_p_value": "fdr",
                              "selection_fdr": "training_fdr", "HUMAN_ORTHOLOG_SYMBOL": "gene_symbol_human"})
        d = d.assign(dataset="rat_train", layer=layer, stat_type="zscore" if statcol == "zscore" else "t")
        d = legacy.collapse(d.dropna(subset=["stat"]).rename(columns={"gene_symbol_human": "gene"}),
                            ["dataset", "tissue", "layer", "contrast", "gene"]).rename(columns={"gene": "gene_symbol_human"})
        is_join = np.array([(ds, t, l) in JOIN_TISSUE_LAYERS for ds, t, l in zip(d.dataset, d.tissue, d.layer)])
        keep.append(_finish(d[~is_join], "alltissue_v1"))
        overlap.append(d[is_join])
        log(f"rat {layer}: {d.tissue.nunique()} tissues, {len(d)} gene rows ({(~is_join).sum()} kept, "
            f"{is_join.sum()} overlap join)")
        del d
    return pd.concat(keep, ignore_index=True), pd.concat(overlap, ignore_index=True), pd.concat(cov, ignore_index=True)


# ------------------------------------------------------------------------------------------------- human
HSRC = [("muscle", "RNA", "transcript-rna-seq", RAW / "human_muscle_rna_DA_v1.csv"),
        ("muscle", "PROT", "prot-pr", RAW / "human_muscle_prot_DA_v1.csv"),
        ("muscle", "PHOSPHO", "prot-ph", RAW / "human_muscle_phospho_DA_v1.csv"),
        ("adipose", "RNA", "transcript-rna-seq", EVERY_RAW / "human_adipose_rna_DA_v1.csv.gz"),
        ("adipose", "PROT", "prot-pr", EVERY_RAW / "human_adipose_prot_DA_v1.csv.gz"),
        ("adipose", "PHOSPHO", "prot-ph", EVERY_RAW / "human_adipose_phospho_DA_v1.csv.gz"),
        ("blood", "RNA", "transcript-rna-seq", EVERY_RAW / "human_blood_rna_DA_v1.csv.gz"),
        ("blood", "PROT_OLINK", "prot-ol", EVERY_RAW / "human_blood_olink_DA_v1.csv.gz")]


def build_human(log):
    hf2g = pd.read_csv(EVERY_RAW / "human_feature_to_gene_all_v1.csv.gz", low_memory=False)
    hf2g = hf2g.dropna(subset=["gene_symbol"]).drop_duplicates(["assay", "feature_id"])
    join_short = set(legacy.H_CONTRASTS)
    keep, overlap, cov, cmeta = [], [], [], []
    for tissue, layer, assay, f in HSRC:
        head = pd.read_csv(f, nrows=0).columns
        extra = [c for c in ["AveExpr"] if c in head]
        d = pd.read_csv(f, usecols=["contrast_short", "feature_id", "logFC", "t", "p_value", "adj_p_value"] + extra,
                        low_memory=False)
        d["fdr_bh"] = d.groupby("contrast_short").p_value.transform(legacy.bh)
        T = H_TISSUE[tissue]
        cmeta.append(pd.DataFrame({"contrast_raw": d.contrast_short.unique()}).assign(tissue=T, layer=layer))
        d["contrast"] = d.contrast_short.map({c: human_contrast_code(c) for c in d.contrast_short.unique()})
        if "AveExpr" in d:
            d["baseline_expr"] = d.AveExpr
            cov.append(d[["contrast", "feature_id", "baseline_expr"]].assign(tissue=T, layer=layer))
        sym = hf2g[hf2g.assay == assay].set_index("feature_id").gene_symbol
        d["gene"] = d.feature_id.map(sym)
        d = d.dropna(subset=["gene"]).rename(columns={"t": "stat", "p_value": "p", "adj_p_value": "fdr"})
        d = d.assign(dataset="human_acute", tissue=T, layer=layer, stat_type="t_moderated")
        d = legacy.collapse(d.dropna(subset=["stat"]), ["dataset", "tissue", "layer", "contrast", "gene"])
        d = d.rename(columns={"gene": "gene_symbol_human"})
        is_join = (d.tissue == "VL") & d.layer.isin(["RNA", "PROT"]) & d.contrast_short.isin(join_short)
        keep.append(_finish(d[~is_join], "alltissue_v1"))
        overlap.append(d[is_join])
        log(f"human {tissue} {layer}: {d.contrast.nunique()} contrasts, {len(d)} gene rows ({(~is_join).sum()} kept)")
    return (pd.concat(keep, ignore_index=True), pd.concat(overlap, ignore_index=True),
            pd.concat(cov, ignore_index=True), pd.concat(cmeta, ignore_index=True))


# ------------------------------------------------------------------------------------------------- join
def build_join(rat_cov, hum_cov, log):
    j = pd.read_csv(JOIN, low_memory=False)
    n_csv = len(j)
    j["contrast"] = j.contrast.astype(str)
    # covariates for the join rows: human AveExpr by (contrast, feature); rat baseline count / numNAs by
    # (tissue, layer, feature, contrast)
    h = j.dataset == "human_acute"
    hc = hum_cov[(hum_cov.tissue == "VL")].rename(columns={"baseline_expr": "b"})
    jh = j[h].merge(hc[["layer", "contrast", "feature_id", "b"]].drop_duplicates(["layer", "contrast", "feature_id"]),
                    on=["layer", "contrast", "feature_id"], how="left")
    jh["baseline_expr"] = jh.b
    rc = rat_cov.rename(columns={"feature_ID": "feature_id"})
    rc = rc.drop_duplicates(["tissue", "layer", "feature_id", "contrast"])
    jr = j[~h].merge(rc, on=["tissue", "layer", "feature_id", "contrast"], how="left")
    out = pd.concat([jh.drop(columns="b"), jr], ignore_index=True)
    assert len(out) == n_csv, (len(out), n_csv)
    log(f"join_table_v2: {n_csv} rows copied")
    return _finish(out, "join_table_v2"), n_csv


def crosscheck(join_rows, overlap_rows):
    """Compare the join table with the same columns rebuilt from the alltissue dumps."""
    key = ["dataset", "tissue", "layer", "contrast", "gene_symbol_human"]
    o = overlap_rows.copy()
    o["tissue"] = o.tissue.replace({"muscle": "VL"})
    m = join_rows[key + ["stat", "fdr_bh"]].merge(o[key + ["stat", "fdr_bh"]], on=key, how="outer",
                                                   suffixes=("_join", "_rebuilt"), indicator=True)
    both = m[m._merge == "both"]
    return dict(n_join=int((m._merge != "right_only").sum()), n_rebuilt=int((m._merge != "left_only").sum()),
                n_both=len(both), n_stat_identical=int(np.isclose(both.stat_join, both.stat_rebuilt).sum()),
                n_fdr_bh_identical=int(np.isclose(both.fdr_bh_join, both.fdr_bh_rebuilt).sum()),
                n_join_only=int((m._merge == "left_only").sum()), n_rebuilt_only=int((m._merge == "right_only").sum()))


# ------------------------------------------------------------------------------------------------- columns
def column_index(rows, hmeta):
    g = rows.groupby("column_id", observed=True).agg(
        dataset=("dataset", "first"), species=("species", "first"), tissue=("tissue", "first"),
        layer=("layer", "first"), contrast=("contrast", "first"), source=("source", "first"),
        n_genes=("gene_symbol_human", "size"), n_fdrbh05=("fdr_bh", lambda s: int((s < 0.05).sum())),
        stat_type=("stat_type", "first")).reset_index()
    for c in ["dataset", "species", "tissue", "layer", "contrast", "source", "stat_type"]:
        g[c] = g[c].astype(str)
    raw = {(t, l, human_contrast_code(c)): c for c, t, l in zip(hmeta.contrast_raw, hmeta.tissue, hmeta.layer)}
    raw.update({("VL", l, v): k for k, v in legacy.H_CONTRASTS.items() for l in ("RNA", "PROT")})
    recs = []
    for r in g.itertuples():
        if r.dataset == "human_acute":
            short = raw[(r.tissue, r.layer, r.contrast)]
            lab, kind, key = legacy.parse_human(short)
            grp = r.contrast.split("_")[0]
            time = lab.rsplit(" ", 1)[-1]
            sex, week, tord = "", "", key[3]
        else:
            short, lab, kind = "", r.contrast.replace("_", " "), "training vs sedentary"
            grp, sex, week = "trained", r.contrast[0], r.contrast[2:]
            time, tord = week, int(week[:-1])
        recs.append(dict(column_id=r.column_id, contrast_raw=short, label=lab, kind=kind, group=grp, sex=sex,
                         time=time, time_order=tord,
                         early_human_muscle=bool(r.dataset == "human_acute" and r.tissue == "VL"
                                                 and "15to45min" in r.contrast),
                         gtex_tissue=GTEX_TISSUE.get(r.tissue) or ""))
    g = g.merge(pd.DataFrame(recs), on="column_id")
    order = {t: i for i, t in enumerate(["VL", "ADIPOSE", "BLOOD"] + legacy.RAT_ORDER)}
    lorder = {l: i for i, l in enumerate(["RNA", "PROT", "PHOSPHO", "PROT_OLINK"])}
    kord = {"exercise vs control": 0, "exercise post−pre": 1, "reference (non-exercise)": 2, "training vs sedentary": 0}
    g["_o"] = list(zip(g.dataset.map({"human_acute": 0, "rat_train": 1}), g.tissue.map(order), g.layer.map(lorder),
                       g.kind.map(kord), g.group, g.sex, g.time_order, g.contrast))
    g = g.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    g["col_order"] = np.arange(len(g))
    return g


# ------------------------------------------------------------------------------------------------- annotations
def build_annotations(genes):
    ann = pd.DataFrame({"gene_symbol_human": sorted(genes)})
    # GTEx v8 (gencode v26) uses some retired symbols (ATP5B, ATP5L, ...): match on Ensembl gene id first
    # (store symbol -> ENSG via the MoTrPAC human feature map and the RGD ortholog table), then on symbol.
    gtex = pd.read_csv(MODEL_EXT / "gtex_v8_gene_median_tpm.gct.gz", sep="\t", skiprows=2)
    gtex["ens"] = gtex.Name.str.replace(r"\.\d+$", "", regex=True)
    tcols = sorted({v for v in GTEX_TISSUE.values() if v})
    by_ens = gtex.groupby("ens")[tcols].max()
    by_sym = gtex.groupby("Description")[tcols].max()
    idm = build_idmap()
    ens = idm[(idm.id_type == "ensembl") & idm.id.str.startswith("ENSG")]
    ens = ens[ens.id.isin(by_ens.index)].drop_duplicates("gene_symbol_human").set_index("gene_symbol_human").id
    for t in tcols:
        v = ann.gene_symbol_human.map(ens).map(by_ens[t])
        ann[f"gtex_tpm|{t}"] = v.fillna(ann.gene_symbol_human.map(by_sym[t])).astype("float32")
    # GO:CC flags exactly as scripts/08_discordance_model.py (msigdbr drops terms > 2000 genes, hence the proxies)
    go = pd.read_csv(MODEL_EXT / "gocc_msigdbr.csv")
    gsets = go.groupby("gs_name").gene_symbol.apply(set)

    def members(*names):
        return set().union(*[gsets[n] for n in names])
    extracell = members("GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE", "GOCC_BASEMENT_MEMBRANE", "GOCC_COLLAGEN_TRIMER",
                        "GOCC_BLOOD_MICROPARTICLE", "GOCC_PLATELET_ALPHA_GRANULE_LUMEN",
                        *[n for n in gsets.index if n.endswith("LIPOPROTEIN_PARTICLE")])
    complex_terms = [n for n in gsets.index if "COMPLEX" in n]
    flags = {"mitocarta3": set(pd.read_csv(MODEL_EXT / "mitocarta3_human.csv").Symbol.dropna()),
             "go_mitochondrion": gsets["GOCC_MITOCHONDRION"],
             "stable_complex": members("GOCC_RESPIRATORY_CHAIN_COMPLEX", "GOCC_RIBOSOME", "GOCC_PROTEASOME_COMPLEX",
                                       "GOCC_SPLICEOSOMAL_COMPLEX"),
             "go_any_complex": members(*complex_terms),
             "secreted_proxy": extracell,
             "go_contractile_fiber": gsets["GOCC_CONTRACTILE_MUSCLE_FIBER"]}
    for k, s in flags.items():
        ann[k] = ann.gene_symbol_human.isin(s)
    return ann


def build_genesets():
    parts = [pd.read_csv(RAW / "genesets_v1.csv").assign(source="data/raw/genesets_v1.csv (msigdbr)"),
             pd.read_csv(FIGTAB / "genesets_extra.csv").assign(source="fig/tables/genesets_extra.csv (msigdbr)"),
             pd.read_csv(FIGTAB / "genesets_pkg.csv").assign(source="MotrpacHumanPreSuspensionAnalysis 2.0.8"),
             pd.read_csv(MODEL_EXT / "gocc_msigdbr.csv").assign(source="model/ext/gocc_msigdbr.csv (msigdbr GO:CC)")]
    g = pd.concat([p[["gs_name", "gene_symbol", "source"]] for p in parts], ignore_index=True)
    g = pd.concat([g, pd.DataFrame({"gs_name": "MYH_SET", "gene_symbol": ["MYH1", "MYH2", "MYH4", "MYH7"],
                                    "source": "scripts/05_pah_grid.py"})], ignore_index=True)
    return g.drop_duplicates(["gs_name", "gene_symbol"]).reset_index(drop=True)


def build_idmap():
    hf = pd.read_csv(EVERY_RAW / "human_feature_to_gene_all_v1.csv.gz", low_memory=False,
                     usecols=["assay", "feature_id", "gene_symbol", "uniprot"]).dropna(subset=["gene_symbol"])
    hm = pd.read_csv(RAW / "human_feature_to_gene_muscle_v1.csv", low_memory=False,
                     usecols=["assay", "feature_id", "gene_symbol", "uniprot"]).dropna(subset=["gene_symbol"])
    h = pd.concat([hf, hm])
    h = h[h.assay.isin(["prot-pr", "prot-ol", "prot-ph"])]
    uni = pd.concat([pd.DataFrame({"id": h.uniprot, "gene_symbol_human": h.gene_symbol}),
                     pd.DataFrame({"id": h.feature_id[h.assay == "prot-pr"], "gene_symbol_human": h.gene_symbol[h.assay == "prot-pr"]})])
    uni = uni.dropna()
    uni["id"] = uni.id.astype(str)
    uni = pd.concat([uni, uni.assign(id=uni.id.str.replace(r"-\d+$", "", regex=True))])
    uni = uni.drop_duplicates().assign(id_type="uniprot")
    orth = pd.read_csv(RAW / "rat_to_human_gene_v1.csv", low_memory=False,
                       usecols=["RAT_SYMBOL", "RAT_ENSEMBL_ID", "HUMAN_ORTHOLOG_SYMBOL", "HUMAN_ORTHOLOG_ENSEMBL_ID"])
    orth = orth.dropna(subset=["RAT_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL"]).drop_duplicates()
    rat = pd.DataFrame({"id": orth.RAT_SYMBOL, "gene_symbol_human": orth.HUMAN_ORTHOLOG_SYMBOL, "id_type": "rat_symbol"})
    rens = orth.dropna(subset=["RAT_ENSEMBL_ID"])
    rens = pd.DataFrame({"id": rens.RAT_ENSEMBL_ID, "gene_symbol_human": rens.HUMAN_ORTHOLOG_SYMBOL,
                         "id_type": "ensembl"})
    hr = pd.read_csv(EVERY_RAW / "human_feature_to_gene_all_v1.csv.gz", low_memory=False,
                     usecols=["assay", "feature_id", "gene_symbol", "ensembl_gene"]).dropna(subset=["gene_symbol"])
    hr = hr[hr.assay == "transcript-rna-seq"]
    hens = pd.concat([pd.DataFrame({"id": hr.feature_id, "gene_symbol_human": hr.gene_symbol}),
                      pd.DataFrame({"id": hr.ensembl_gene, "gene_symbol_human": hr.gene_symbol}),
                      pd.DataFrame({"id": orth.HUMAN_ORTHOLOG_ENSEMBL_ID, "gene_symbol_human": orth.HUMAN_ORTHOLOG_SYMBOL})])
    hens = hens.dropna()
    hens["id"] = hens.id.astype(str).str.replace(r"\.\d+$", "", regex=True)
    hens = hens.drop_duplicates().assign(id_type="ensembl")
    return pd.concat([uni, rat, rens, hens], ignore_index=True)[["id_type", "id", "gene_symbol_human"]].drop_duplicates()


# ------------------------------------------------------------------------------------------------- provenance
def git_sha(path=HACK):
    try:
        return subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "NA"


def r_pkg_version(pkg):
    for lib in [os.path.expanduser("~/miniconda3/envs/motrpac/lib/R/library")]:
        f = os.path.join(lib, pkg, "DESCRIPTION")
        if os.path.exists(f):
            for line in open(f):
                if line.startswith("Version:"):
                    return line.split(":", 1)[1].strip()
    return "NA"


def build(log=print, hash_inputs=True):
    STORE.mkdir(parents=True, exist_ok=True)
    t0 = dt.datetime.now()
    rat_rows, rat_overlap, rat_cov, = build_rat(log)
    hum_rows, hum_overlap, hum_cov, hmeta = build_human(log)
    join_rows, n_csv = build_join(rat_cov, hum_cov, log)
    xc = crosscheck(join_rows, pd.concat([rat_overlap, hum_overlap], ignore_index=True))
    log(f"cross-check join vs rebuilt: {xc}")
    rows = pd.concat([join_rows, rat_rows, hum_rows], ignore_index=True)
    del rat_rows, hum_rows, rat_overlap, hum_overlap
    for c in CAT_COLS:
        rows[c] = rows[c].astype(str).astype("category")
    for c in ["baseline_expr", "prot_n_missing", "training_fdr"]:
        rows[c] = rows[c].astype("float32")
    rows["n_collapsed"] = rows.n_collapsed.astype("int16")
    cols = column_index(rows, hmeta)
    rows = rows.sort_values(["column_id", "gene_symbol_human"]).reset_index(drop=True)
    rows.to_parquet(CONTRASTS, compression="zstd", index=False)
    cols.to_parquet(COLUMNS, index=False)
    ann = build_annotations(set(rows.gene_symbol_human.astype(str)))
    ann.to_parquet(ANNOT, index=False)
    build_genesets().to_parquet(GENESETS, index=False)
    build_idmap().to_parquet(IDMAP, index=False)

    inputs = [JOIN, RAW / "rat_alltissue_trnscrpt_DA_v1.csv", RAW / "rat_alltissue_prot_DA_v1.csv",
              RAW / "rat_alltissue_phospho_DA_v1.csv", RAW / "rat_alltissue_feature_to_gene_v1.csv",
              RAW / "rat_to_human_gene_v1.csv", RAW / "human_feature_to_gene_muscle_v1.csv",
              RAW / "genesets_v1.csv", FIGTAB / "genesets_extra.csv", FIGTAB / "genesets_pkg.csv",
              EVERY_RAW / "human_feature_to_gene_all_v1.csv.gz", MODEL_EXT / "gtex_v8_gene_median_tpm.gct.gz",
              MODEL_EXT / "mitocarta3_human.csv", MODEL_EXT / "gocc_msigdbr.csv"] + [s[3] for s in HSRC]
    prov = [("built_at", t0.isoformat(timespec="seconds")),
            ("build_seconds", f"{(dt.datetime.now() - t0).total_seconds():.0f}"),
            ("repo_git_sha", git_sha()),
            ("python", sys.version.split()[0]), ("pandas", pd.__version__), ("numpy", np.__version__),
            ("MotrpacRatTraining6moData", r_pkg_version("MotrpacRatTraining6moData")),
            ("MotrpacHumanPreSuspensionAnalysis", r_pkg_version("MotrpacHumanPreSuspensionAnalysis")),
            ("MotrpacRatTraining6mo", r_pkg_version("MotrpacRatTraining6mo")),
            ("join_table_version", "join_table_v2 (scripts/03_join.py v2: rat RNA on/off rows dropped)"),
            ("join_table_rows_csv", str(n_csv)),
            ("join_table_rows_store", str(int((rows.source == "join_table_v2").sum()))),
            ("store_rows", str(len(rows))), ("store_columns", str(len(cols))),
            ("crosscheck_join_vs_rebuilt", json.dumps(xc)),
            ("gtex", "GTEx v8 gene median TPM (model/ext/gtex_v8_gene_median_tpm.gct.gz)"),
            ("mitocarta", "MitoCarta3.0 human (model/ext/mitocarta3_human.csv)"),
            ("legacy_script_sha256_16", json.dumps(legacy.script_hashes()))]
    for p in inputs:
        prov.append((f"input:{os.path.relpath(p, HACK)}",
                     json.dumps({"bytes": os.path.getsize(p), "sha256": sha256(p) if hash_inputs else "skipped"})))
    pd.DataFrame(prov, columns=["key", "value"]).to_parquet(PROVENANCE, index=False)
    log(f"wrote {CONTRASTS} ({len(rows)} rows, {len(cols)} columns) in {(dt.datetime.now() - t0).total_seconds():.0f} s")
    return rows, cols


# ------------------------------------------------------------------------------------------------- load
_CACHE = {}


def load(columns=None):
    """Store rows as a DataFrame (categoricals kept). `columns`: subset of row fields to read."""
    key = tuple(columns) if columns else None
    if key not in _CACHE:
        _CACHE[key] = pd.read_parquet(CONTRASTS, columns=list(columns) if columns else None)
    return _CACHE[key]


def load_columns():
    if "cols" not in _CACHE:
        _CACHE["cols"] = pd.read_parquet(COLUMNS)
    return _CACHE["cols"]


def load_annotations():
    if "ann" not in _CACHE:
        _CACHE["ann"] = pd.read_parquet(ANNOT).set_index("gene_symbol_human")
    return _CACHE["ann"]


def load_genesets():
    if "gs" not in _CACHE:
        g = pd.read_parquet(GENESETS)
        _CACHE["gs"] = {k: set(v) for k, v in g.groupby("gs_name").gene_symbol}
    return _CACHE["gs"]


def load_idmap():
    if "id" not in _CACHE:
        _CACHE["id"] = pd.read_parquet(IDMAP)
    return _CACHE["id"]


def load_provenance():
    if "prov" not in _CACHE:
        _CACHE["prov"] = dict(pd.read_parquet(PROVENANCE).itertuples(index=False))
    return _CACHE["prov"]


def store_hash():
    return sha256(CONTRASTS)[:16]


# ------------------------------------------------------------------------------------------------- METAB layer
# Metabolomics as a third layer (schema: store/SCHEMA.md "## METAB layer"). Same row schema as contrasts.parquet;
# for layer == "METAB" the gene_symbol_human column holds the RefMet NAME of the metabolite (not a gene), and
# feature_id holds "<platform>:<feature_ID>" (rat) or the human RefMet feature_id with its platform prefix.
RAW_METAB = STORE / "raw_metab"
METAB = STORE / "metab.parquet"
METAB_FEATURES = STORE / "metab_features.parquet"
METAB_COLUMNS = STORE / "metab_columns.parquet"
PATHWAY_MAP = STORE / "pathway_map.csv"
METAB_SOURCE = "metab_v1"
ISTD_TAG = "[iSTD]"


def _metab_features(rat_fm, hum_fm):
    """One row per RefMet name used in the METAB rows (union rat/human); first non-missing annotation wins."""
    r = rat_fm.rename(columns={"metabolite_refmet": "refmet_name"}).assign(species="rat")
    h = hum_fm.dropna(subset=["refmet_name"]).copy()
    h["kegg_id"] = h.kegg_id.fillna(h.f2g_kegg_id)  # RefMet download first, MoTrPAC HUMAN_FEATURE_TO_GENE second
    h["species"] = "human"
    keep = ["refmet_name", "species", "refmet_id", "super_class", "main_class", "sub_class", "kegg_id", "hmdb_id",
            "pubchem_cid", "inchi_key", "class_source"]
    u = pd.concat([r[keep], h[keep]], ignore_index=True)
    u = u[~u.refmet_name.str.contains(ISTD_TAG, regex=False)]
    f = u.groupby("refmet_name", sort=True).agg(
        species=("species", lambda s: "/".join(sorted(set(s)))),
        **{c: (c, "first") for c in keep if c not in ("refmet_name", "species")}).reset_index()
    f["pubchem_cid"] = f.pubchem_cid.map(lambda v: "" if pd.isna(v) else str(int(float(v))))
    f = f.rename(columns={"kegg_id": "kegg", "hmdb_id": "hmdb", "pubchem_cid": "pubchem"})
    # notation-insensitive key: 'CAR(16:0(OH))' (rat, RefMet 2021 notation) == 'CAR 16:0;OH' (human, current RefMet)
    f["name_key"] = f.refmet_name.str.lower().str.replace(r"[\s();]", "", regex=True)
    f["class_source"] = f.class_source.fillna("none")
    cols = ["refmet_name", "super_class", "main_class", "sub_class", "kegg", "hmdb", "pubchem", "class_source",
            "refmet_id", "inchi_key", "species", "name_key"]
    return f[cols]


def build_metab(log=print):
    """METAB rows (store/metab.parquet), feature annotation (metab_features.parquet) and a column index
    (metab_columns.parquet, same fields as columns.parquet). Rules:
      * fdr_bh = BH within dataset x tissue x contrast over all tested features before collapsing;
      * rat: one row per RefMet name per column = the platform feature with max |t| (legacy.collapse),
        n_collapsed = number of platform features for that name; internal standards ('[iSTD]') are tested in the
        BH but not kept as rows; rat PLASMA stays 'PLASMA';
      * human: feature_id is already a RefMet name (one lowest-CV platform per metabolite, so n_collapsed = 1);
        the two non-RefMet conventional-assay features (KET, NEFA) are dropped; contrasts coded with
        human_contrast_code()."""
    rat_fm = pd.read_csv(RAW_METAB / "rat_metab_feature_map_v1.csv.gz", low_memory=False)
    hum_fm = pd.read_csv(RAW_METAB / "human_metab_feature_map_v1.csv.gz", low_memory=False)

    # ---- rat
    d = pd.read_csv(RAW_METAB / "rat_metab_DA_v1.csv.gz", low_memory=False,
                    usecols=["tissue", "dataset", "feature_ID", "metabolite_refmet", "sex", "comparison_group", "logFC",
                             "tscore", "p_value", "adj_p_value", "selection_fdr", "reference_average_intensity"])
    n_raw = len(d)
    d["contrast"] = d.sex.str[0].str.upper() + "_" + d.comparison_group
    d["fdr_bh"] = d.groupby(["tissue", "contrast"]).p_value.transform(legacy.bh)
    istd = d.metabolite_refmet.str.contains(ISTD_TAG, regex=False)
    d = d[~istd]
    d = d.rename(columns={"tscore": "stat", "p_value": "p", "adj_p_value": "fdr", "selection_fdr": "training_fdr",
                          "reference_average_intensity": "baseline_expr", "metabolite_refmet": "gene"})
    d["feature_id"] = d.dataset + ":" + d.feature_ID
    d = d.assign(dataset="rat_train", layer="METAB", stat_type="t")
    d = legacy.collapse(d.dropna(subset=["stat"]), ["dataset", "tissue", "layer", "contrast", "gene"])
    rat = _finish(d.rename(columns={"gene": "gene_symbol_human"}), METAB_SOURCE)
    log(f"rat METAB: {n_raw} DA rows ({int(istd.sum())} internal-standard rows dropped after BH) -> {len(rat)} "
        f"RefMet rows, {rat.tissue.nunique()} tissues")

    # ---- human
    h = pd.read_csv(RAW_METAB / "human_metab_DA_v1.csv.gz", low_memory=False,
                    usecols=["tissue", "platform", "contrast_short", "feature_id", "logFC", "t", "p_value",
                             "adj_p_value", "AveExpr"])
    h["fdr_bh"] = h.groupby(["tissue", "contrast_short"]).p_value.transform(legacy.bh)
    refmet = set(hum_fm.refmet_name.dropna())
    nonref = ~h.feature_id.isin(refmet)
    log(f"human METAB: dropping {int(nonref.sum())} rows of non-RefMet features {sorted(h.feature_id[nonref].unique())}")
    h = h[~nonref]
    h["T"] = h.tissue.map(H_TISSUE)
    hmeta = h[["contrast_short", "T"]].drop_duplicates().rename(columns={"contrast_short": "contrast_raw",
                                                                         "T": "tissue"}).assign(layer="METAB")
    h["contrast"] = h.contrast_short.map({c: human_contrast_code(c) for c in h.contrast_short.unique()})
    h = h.rename(columns={"t": "stat", "p_value": "p", "adj_p_value": "fdr", "AveExpr": "baseline_expr"})
    h["gene"] = h.feature_id
    h["feature_id"] = h.platform + ":" + h.feature_id
    h = h.drop(columns="tissue").rename(columns={"T": "tissue"})
    h = h.assign(dataset="human_acute", layer="METAB", stat_type="t")
    h = legacy.collapse(h.dropna(subset=["stat"]), ["dataset", "tissue", "layer", "contrast", "gene"])
    hum = _finish(h.rename(columns={"gene": "gene_symbol_human"}), METAB_SOURCE)
    log(f"human METAB: {len(hum)} RefMet rows, {hum.column_id.nunique()} columns")

    rows = pd.concat([rat, hum], ignore_index=True)
    for c in CAT_COLS:
        rows[c] = rows[c].astype(str).astype("category")
    for c in ["baseline_expr", "prot_n_missing", "training_fdr"]:
        rows[c] = rows[c].astype("float32")
    rows["n_collapsed"] = rows.n_collapsed.astype("int16")
    rows = rows.sort_values(["column_id", "gene_symbol_human"]).reset_index(drop=True)

    cols = column_index(rows, hmeta)
    # column_index() has no order for METAB / PLASMA: order here, after every column of columns.parquet
    tord = {t: i for i, t in enumerate(["VL", "ADIPOSE", "BLOOD"] + legacy.RAT_ORDER + ["PLASMA"])}
    kord = {"exercise vs control": 0, "exercise post−pre": 1, "reference (non-exercise)": 2, "training vs sedentary": 0}
    cols["_o"] = list(zip(cols.dataset.map({"human_acute": 0, "rat_train": 1}), cols.tissue.map(tord),
                          cols.kind.map(kord), cols.group, cols.sex, cols.time_order, cols.contrast))
    cols = cols.sort_values("_o").drop(columns="_o").reset_index(drop=True)
    base = int(pd.read_parquet(COLUMNS, columns=["col_order"]).col_order.max()) + 1 if COLUMNS.exists() else 0
    cols["col_order"] = base + np.arange(len(cols))

    feats = _metab_features(rat_fm, hum_fm)
    missing = set(rows.gene_symbol_human.astype(str)) - set(feats.refmet_name)
    assert not missing, f"METAB names without a feature row: {sorted(missing)[:10]}"
    feats = feats[feats.refmet_name.isin(set(rows.gene_symbol_human.astype(str)))].reset_index(drop=True)
    rows.to_parquet(METAB, compression="zstd", index=False)
    feats.to_parquet(METAB_FEATURES, index=False)
    cols.to_parquet(METAB_COLUMNS, index=False)
    log(f"wrote {METAB} ({len(rows)} rows, {len(cols)} columns), {METAB_FEATURES} ({len(feats)} names)")
    return rows, cols, feats


def load_metab():
    """(METAB rows, METAB column index). Rows have the contrasts.parquet schema; gene_symbol_human = RefMet name."""
    if "metab" not in _CACHE:
        _CACHE["metab"] = (pd.read_parquet(METAB), pd.read_parquet(METAB_COLUMNS))
    return _CACHE["metab"]


def load_metab_features():
    if "metab_feat" not in _CACHE:
        _CACHE["metab_feat"] = pd.read_parquet(METAB_FEATURES)
    return _CACHE["metab_feat"]


def load_pathway_map():
    """Curated many-to-many pathway map (store/pathway_map.csv): pathway, side (gene|metabolite), member, rule,
    source, measured_in. Genes and metabolites are never joined one-to-one."""
    if "pwmap" not in _CACHE:
        _CACHE["pwmap"] = pd.read_csv(PATHWAY_MAP, keep_default_na=False)
    return _CACHE["pwmap"]
