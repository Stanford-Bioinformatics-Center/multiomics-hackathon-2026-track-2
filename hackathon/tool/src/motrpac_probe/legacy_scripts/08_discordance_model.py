#!/usr/bin/env python3
"""08_discordance_model.py — when do RNA and protein training responses agree, and can gene properties
predict it?  Rat SKM-GN and HEART, trained vs sedentary; train in one tissue, test in the other.

Run from hackathon/:  python3 scripts/08_discordance_model.py   (~6 min on danilogin CPU; deterministic)
Then:                 python3 model/make_report_pdf.py           (-> model/discordance_report.pdf)
Inputs (read-only):
  data/join_table_v2.csv                      gene-level rat_train effects (03_join.py)
  data/raw/rat_{skmgn,heart}_{trnscrpt,prot}_DA_v1.csv   baseline counts / proteomics NAs per feature
  data/raw/rat_skmgn_phospho_DA_v1.csv, rat_feature_to_gene_v1.csv, rat_to_human_gene_v1.csv
  model/ext/gocc_msigdbr.csv                  msigdbr 26.1.1 C5 GO:CC (exported with the motrpac conda env)
  model/ext/mitocarta3_human.csv              MitoCarta3.0 human (Broad, downloaded 2026-09-26)
  model/ext/gtex_v8_gene_median_tpm.gct.gz    GTEx v8 median TPM per tissue (downloaded 2026-09-26)
Outputs: model/tables/*.csv, model/catalog_v1.csv, model/fig/*.png

Order of analysis:
  1. Stratified Spearman(RNA stat, PROT stat) per tissue x contrast within feature strata (no model).
  2. PRIMARY target: continuous residual of PROT stat regressed on RNA stat (OLS per tissue x contrast, all
     genes measured in both layers, no threshold). Primary model Ridge (L2 linear), secondary HGB.
  3. SECONDARY targets: class labels (gene x tissue x contrast, fdr_bh<0.05 in >=1 layer):
       concordant both sig same sign; opposite both sig opposite sign; RNA_only / PROT_only one layer sig.
     Concordant < 50 in every tissue x contrast -> binary collapse:
       Task A  concordant vs discordant (RNA_only + PROT_only + opposite)
       Task B  RNA_only vs PROT_only (which layer carries the discordant response)
       Task C  3-class, run only if concordant >= 50 in both tissues (it is not: SKM-GN 8w pooled = 46)
     Primary model L2 logistic regression, secondary HGB; baselines majority-class and mitochondrial-only.
Sexes are pooled (F_Xw + M_Xw) per time point; sex is NOT a feature; 8w is repeated within males only.
"""
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.impute import SimpleImputer
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import balanced_accuracy_score, r2_score, roc_auc_score
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=UserWarning)
RNG = 0
RAW, EXT, OUT = "data/raw", "model/ext", "model"
TAB, FIG = f"{OUT}/tables", f"{OUT}/fig"
os.makedirs(TAB, exist_ok=True)
os.makedirs(FIG, exist_ok=True)
TISSUES = ["SKM-GN", "HEART"]
WEEKS = ["1w", "2w", "4w", "8w"]
CONTRASTS = [f"{s}_{wk}" for s in "FM" for wk in WEEKS]
MIN_CLASS = 30         # binary tasks: each class needs >= 30 genes
MIN_CONCORDANT_3CLASS = 50
XT = [("SKM-GN", "HEART"), ("HEART", "SKM-GN")]

# ------------------------------------------------------------------ labels
j = pd.read_csv("data/join_table_v2.csv")
j = j[(j.dataset == "rat_train") & j.tissue.isin(TISSUES)]
w = j.pivot_table(index=["tissue", "contrast", "gene_symbol_human"], columns="layer",
                  values=["stat", "logFC", "fdr_bh", "n_collapsed"], aggfunc="first")
w.columns = [f"{a}_{b}" for a, b in w.columns]
fid = j.pivot_table(index=["tissue", "contrast", "gene_symbol_human"], columns="layer",
                    values="feature_id", aggfunc="first")
w = w.join(fid.rename(columns={"RNA": "feature_id_RNA", "PROT": "feature_id_PROT"}))
w = w.dropna(subset=["stat_RNA", "stat_PROT"]).reset_index()
w["sex"] = w.contrast.str[0].map({"F": "female", "M": "male"})
w["week"] = w.contrast.str[2:]
rs, ps = w.fdr_bh_RNA < .05, w.fdr_bh_PROT < .05
same = np.sign(w.stat_RNA) == np.sign(w.stat_PROT)
w["label"] = np.select([rs & ps & same, rs & ps & ~same, rs & ~ps, ~rs & ps],
                       ["concordant", "opposite", "RNA_only", "PROT_only"], "ns")
w["label_binary"] = np.where(w.label == "ns", "ns",
                             np.where(w.label == "concordant", "concordant", "discordant"))

# continuous target: residual of PROT stat on RNA stat, per tissue x contrast
w["resid_prot_on_rna"] = np.nan
for _, idx in w.groupby(["tissue", "contrast"]).groups.items():
    x, y = w.loc[idx, "stat_RNA"].values, w.loc[idx, "stat_PROT"].values
    b1, b0 = np.polyfit(x, y, 1)
    w.loc[idx, "resid_prot_on_rna"] = y - (b0 + b1 * x)


# ------------------------------------------------------------------ features
# baseline RNA abundance: sedentary-control mean normalized count of the RNA feature kept in the join
# (reference_average_intensity in the DESeq2 DA table); proteomics missingness: numNAs of the PROT feature.
# Proteomics "average intensity" columns are TMT log2 ratios to a pooled reference, not abundance -> unused.
base = []
for tissue, tc in [("SKM-GN", "skmgn"), ("HEART", "heart")]:
    r = pd.read_csv(f"{RAW}/rat_{tc}_trnscrpt_DA_v1.csv",
                    usecols=["feature_ID", "sex", "comparison_group", "reference_average_intensity"])
    r = r.rename(columns={"feature_ID": "feature_id_RNA", "comparison_group": "week",
                          "reference_average_intensity": "rna_baseline_count"})
    p = pd.read_csv(f"{RAW}/rat_{tc}_prot_DA_v1.csv", usecols=["feature_ID", "sex", "comparison_group", "numNAs"])
    p = p.rename(columns={"feature_ID": "feature_id_PROT", "comparison_group": "week", "numNAs": "prot_n_missing"})
    base.append((tissue, r.drop_duplicates(["feature_id_RNA", "sex", "week"]),
                 p.drop_duplicates(["feature_id_PROT", "sex", "week"])))
parts = []
for tissue, r, p in base:
    t = w[w.tissue == tissue].merge(r, on=["feature_id_RNA", "sex", "week"], how="left")
    parts.append(t.merge(p, on=["feature_id_PROT", "sex", "week"], how="left"))
w = pd.concat(parts, ignore_index=True)
w["rna_log2_baseline"] = np.log2(w.rna_baseline_count + 1)

# GO:CC (msigdbr drops terms >2000 genes, so 'extracellular region' and 'protein-containing complex' are
# absent -> proxies below, stated in RESULTS.md)
go = pd.read_csv(f"{EXT}/gocc_msigdbr.csv")
gsets = go.groupby("gs_name").gene_symbol.apply(set)
def members(*names):
    return set().union(*[gsets[n] for n in names])
extracell = members("GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE", "GOCC_BASEMENT_MEMBRANE", "GOCC_COLLAGEN_TRIMER",
                    "GOCC_BLOOD_MICROPARTICLE", "GOCC_PLATELET_ALPHA_GRANULE_LUMEN",
                    *[n for n in gsets.index if n.endswith("LIPOPROTEIN_PARTICLE")])
complex_terms = [n for n in gsets.index if "COMPLEX" in n]
FLAGS = {
    "go_mitochondrion": gsets["GOCC_MITOCHONDRION"],
    "go_resp_chain": gsets["GOCC_RESPIRATORY_CHAIN_COMPLEX"],
    "go_ribosome": gsets["GOCC_RIBOSOME"],
    "go_proteasome": gsets["GOCC_PROTEASOME_COMPLEX"],
    "go_spliceosome": gsets["GOCC_SPLICEOSOMAL_COMPLEX"],
    "go_any_complex": members(*complex_terms),
    "go_ecm": gsets["GOCC_EXTERNAL_ENCAPSULATING_STRUCTURE"],
    "go_extracellular": extracell,
    "go_contractile_fiber": gsets["GOCC_CONTRACTILE_MUSCLE_FIBER"],
}
mc = pd.read_csv(f"{EXT}/mitocarta3_human.csv")
FLAGS["mitocarta3"] = set(mc.Symbol.dropna())
for k, s in FLAGS.items():
    w[k] = w.gene_symbol_human.isin(s).astype(int)
n_terms = go.groupby("gene_symbol").gs_name.nunique()
w["n_gocc_terms"] = w.gene_symbol_human.map(n_terms).fillna(0)

# CV group = smallest GO:CC complex term (5-300 genes) the gene belongs to, else the gene itself
csize = go[go.gs_name.isin(complex_terms)].groupby("gs_name").gene_symbol.nunique()
cterm = go[go.gs_name.isin(csize[(csize >= 5) & (csize <= 300)].index)].copy()
cterm["size"] = cterm.gs_name.map(csize)
fam = cterm.sort_values(["size", "gs_name"]).drop_duplicates("gene_symbol").set_index("gene_symbol").gs_name
w["cv_group"] = w.gene_symbol_human.map(fam).fillna("gene:" + w.gene_symbol_human)

# phosphosites detected in rat SKM-GN phosphoproteome (only phospho table available; used for both tissues)
ph = pd.read_csv(f"{RAW}/rat_skmgn_phospho_DA_v1.csv", usecols=["feature_ID"]).drop_duplicates()
f2g = pd.read_csv(f"{RAW}/rat_feature_to_gene_v1.csv", usecols=["feature_ID", "gene_symbol"]).dropna().drop_duplicates()
orth = pd.read_csv(f"{RAW}/rat_to_human_gene_v1.csv", usecols=["RAT_SYMBOL", "HUMAN_ORTHOLOG_SYMBOL"]).dropna().drop_duplicates()
phg = ph.merge(f2g, on="feature_ID").merge(orth, left_on="gene_symbol", right_on="RAT_SYMBOL")
nsite = phg.groupby("HUMAN_ORTHOLOG_SYMBOL").feature_ID.nunique()
skm_prot_genes = set(j[(j.tissue == "SKM-GN") & (j.layer == "PROT")].gene_symbol_human)
w["n_phosphosites_skmgn"] = w.gene_symbol_human.map(nsite)
w.loc[w.n_phosphosites_skmgn.isna() & w.gene_symbol_human.isin(skm_prot_genes), "n_phosphosites_skmgn"] = 0
# genes not in the SKM-GN proteome stay NA (not measurable there)


# GTEx v8 median TPM in the matching human tissue (orthologous gene symbol; max over duplicated symbols)
gtex = pd.read_csv(f"{EXT}/gtex_v8_gene_median_tpm.gct.gz", sep="\t", skiprows=2,
                   usecols=["Description", "Muscle - Skeletal", "Heart - Left Ventricle"])
gtex = gtex.groupby("Description").max()
GTEX_COL = {"SKM-GN": "Muscle - Skeletal", "HEART": "Heart - Left Ventricle"}
w["gtex_log2_tpm"] = np.nan
for t, col in GTEX_COL.items():
    m = w.tissue == t
    w.loc[m, "gtex_log2_tpm"] = np.log2(w.loc[m, "gene_symbol_human"].map(gtex[col]) + 1)
# stable-complex subunit = member of respiratory chain / ribosome / proteasome / spliceosome (GO:CC)
w["stable_complex"] = w[["go_resp_chain", "go_ribosome", "go_proteasome", "go_spliceosome"]].max(axis=1)

FEATURES = ["rna_log2_baseline", "gtex_log2_tpm", "prot_n_missing", "n_collapsed_RNA", "n_collapsed_PROT",
            "go_mitochondrion", "mitocarta3", "go_resp_chain", "go_ribosome", "go_proteasome",
            "go_spliceosome", "go_any_complex", "go_ecm", "go_extracellular", "go_contractile_fiber",
            "n_phosphosites_skmgn", "n_gocc_terms"]
feat_cov = w.groupby("tissue")[FEATURES].apply(lambda d: d.notna().mean()).T.round(3)
feat_cov.to_csv(f"{TAB}/feature_coverage.csv")

counts = w.groupby(["tissue", "contrast"]).label.value_counts().unstack().fillna(0).astype(int)
counts = counts[["concordant", "RNA_only", "PROT_only", "opposite", "ns"]]
counts["n_both_layers"] = counts.sum(axis=1)
counts.to_csv(f"{TAB}/class_counts.csv")
pooled = w.groupby(["tissue", "week"]).label.value_counts().unstack().fillna(0).astype(int)
pooled = pooled[["concordant", "RNA_only", "PROT_only", "opposite", "ns"]]
pooled.to_csv(f"{TAB}/class_counts_sexpooled.csv")
print(counts.to_string())

# ================================================================== 1. stratified Spearman (no model)
# rho(RNA stat, PROT stat) over ALL genes measured in both layers, within each stratum; 95% percentile
# bootstrap over genes (B=300). delta = level minus reference level, from paired bootstrap iterations.
B = 300
w["abund_tertile"] = np.nan
for _, idx in w.groupby(["tissue", "contrast"]).groups.items():
    v = w.loc[idx, "gtex_log2_tpm"]
    ok = v.notna()
    w.loc[v[ok].index, "abund_tertile"] = pd.qcut(v[ok].rank(method="first"), 3, labels=[1, 2, 3]).astype(float)
w["phospho_bin"] = pd.cut(w.n_phosphosites_skmgn, [-0.5, 0.5, 3.5, np.inf], labels=["0", "1-3", ">3"]).astype(object)
STRATA = {  # name -> (column, ordered levels, labels, reference level)
    "mitochondrial (MitoCarta3.0)": ("mitocarta3", [0, 1], ["no", "yes"], 0),
    "stable-complex subunit (GO:CC)": ("stable_complex", [0, 1], ["no", "yes"], 0),
    "secreted/extracellular (GO:CC proxy)": ("go_extracellular", [0, 1], ["no", "yes"], 0),
    "GTEx abundance tertile (matched tissue)": ("abund_tertile", [1.0, 2.0, 3.0], ["low", "mid", "high"], 1.0),
    "phosphosites in rat SKM-GN": ("phospho_bin", ["0", "1-3", ">3"], ["0", "1-3", ">3"], "0"),
}

def rank_rho(x, y):
    return np.corrcoef(pd.Series(x).rank().values, pd.Series(y).rank().values)[0, 1]

def boot_rho(x, y, rng):
    x, y = np.asarray(x), np.asarray(y)
    out = np.empty(B)
    for b in range(B):
        i = rng.integers(0, len(x), len(x))
        out[b] = rank_rho(x[i], y[i])
    return out

srows = []
rng = np.random.default_rng(RNG)
for (t, c), g in w.groupby(["tissue", "contrast"]):
    bs_all = boot_rho(g.stat_RNA, g.stat_PROT, rng)
    srows.append(dict(tissue=t, contrast=c, stratum="all genes", level="all", n=len(g),
                      rho=rank_rho(g.stat_RNA, g.stat_PROT), lo=np.percentile(bs_all, 2.5), hi=np.percentile(bs_all, 97.5)))
    for sname, (col, levels, labs, ref) in STRATA.items():
        boots = {}
        for lev, lab in zip(levels, labs):
            s = g[g[col] == lev]
            if len(s) < 10:
                srows.append(dict(tissue=t, contrast=c, stratum=sname, level=lab, n=len(s)))
                continue
            boots[lev] = boot_rho(s.stat_RNA, s.stat_PROT, rng)
            srows.append(dict(tissue=t, contrast=c, stratum=sname, level=lab, n=len(s),
                              rho=rank_rho(s.stat_RNA, s.stat_PROT),
                              lo=np.percentile(boots[lev], 2.5), hi=np.percentile(boots[lev], 97.5)))
        for lev, lab in zip(levels, labs):
            if lev != ref and lev in boots and ref in boots:
                d = boots[lev] - boots[ref]
                srows[-len(levels) + levels.index(lev)].update(
                    delta_vs_ref=srows[-len(levels) + levels.index(lev)]["rho"] - srows[-len(levels) + levels.index(ref)]["rho"],
                    delta_lo=np.percentile(d, 2.5), delta_hi=np.percentile(d, 97.5))
strat = pd.DataFrame(srows)
strat.to_csv(f"{TAB}/stratified_rho.csv", index=False)
print(strat[strat.contrast.str.endswith("8w")].round(3).to_string())

# ================================================================== 2. PRIMARY: continuous residual
def ridge():
    return make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler(), Ridge(1.0))
def hgb_reg():
    return HistGradientBoostingRegressor(learning_rate=0.05, max_iter=200, max_leaf_nodes=15,
                                         min_samples_leaf=20, l2_regularization=1.0, random_state=RNG)
REGS = {"Ridge": ridge, "HGB": hgb_reg}

def boot_reg(y, p, n=300):
    rng = np.random.default_rng(RNG); y, p = np.asarray(y), np.asarray(p); out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y)); out.append(r2_score(y[i], p[i]))
    return np.percentile(out, [2.5, 97.5])

rrows, reg_fitted, resid_pred = [], {}, {}
for week in WEEKS + ["8w_male_only"]:
    if week == "8w_male_only":
        dd = {t: w[(w.tissue == t) & (w.contrast == "M_8w")] for t in TISSUES}
    else:
        dd = {t: w[(w.tissue == t) & (w.week == week)] for t in TISSUES}
    for mname, make in REGS.items():
        for t in TISSUES:
            d = dd[t]; pred = np.zeros(len(d))
            for tr, te in GroupKFold(n_splits=5).split(d, groups=d.cv_group):
                pred[te] = make().fit(d.iloc[tr][FEATURES], d.resid_prot_on_rna.iloc[tr]).predict(d.iloc[te][FEATURES])
            lo, hi = boot_reg(d.resid_prot_on_rna, pred)
            rrows.append(dict(week=week, model=mname, eval=f"CV {t}", n=len(d), R2=r2_score(d.resid_prot_on_rna, pred),
                              R2_lo=lo, R2_hi=hi, spearman=spearmanr(d.resid_prot_on_rna, pred)[0]))
        for tr, te in XT:
            m = make().fit(dd[tr][FEATURES], dd[tr].resid_prot_on_rna)
            pred = m.predict(dd[te][FEATURES])
            lo, hi = boot_reg(dd[te].resid_prot_on_rna, pred)
            rrows.append(dict(week=week, model=mname, eval=f"{tr} -> {te}", n=len(dd[te]),
                              R2=r2_score(dd[te].resid_prot_on_rna, pred), R2_lo=lo, R2_hi=hi,
                              spearman=spearmanr(dd[te].resid_prot_on_rna, pred)[0]))
            reg_fitted[(week, mname, tr)] = m
            if week in WEEKS:
                resid_pred[(week, mname, te)] = pd.Series(pred, index=dd[te].index)
rmetrics = pd.DataFrame(rrows)
rmetrics.to_csv(f"{TAB}/metrics_residual.csv", index=False)
print(rmetrics[rmetrics.week.str.startswith("8w")].round(3).to_string())

# residual: permutation importance (drop in cross-tissue R2) + Ridge standardized coefficients, 8w
rimp = []
for mname in REGS:
    for tr, te in XT:
        d = w[(w.tissue == te) & (w.week == "8w")]
        pi = permutation_importance(reg_fitted[("8w", mname, tr)], d[FEATURES], d.resid_prot_on_rna,
                                    scoring="r2", n_repeats=10, random_state=RNG)
        for f, mu, sd in zip(FEATURES, pi.importances_mean, pi.importances_std):
            rimp.append(dict(model=mname, train=tr, test=te, feature=f, dR2_mean=mu, dR2_sd=sd))
rimp = pd.DataFrame(rimp)
rimp.to_csv(f"{TAB}/residual_permutation_importance_8w.csv", index=False)
coef = []
for tr in TISSUES:
    pipe = reg_fitted[("8w", "Ridge", tr)]
    names = list(FEATURES) + [f"{FEATURES[i]}_missing" for i in pipe[0].indicator_.features_]
    for n_, b in zip(names, pipe[-1].coef_):
        coef.append(dict(train=tr, feature=n_, std_coef=b))
pd.DataFrame(coef).to_csv(f"{TAB}/residual_ridge_coefficients_8w.csv", index=False)

# ================================================================== 3. SECONDARY: class labels
def hgb():
    return HistGradientBoostingClassifier(learning_rate=0.05, max_iter=200, max_leaf_nodes=15,
                                          min_samples_leaf=20, l2_regularization=1.0,
                                          class_weight="balanced", random_state=RNG)
def logreg():
    return make_pipeline(SimpleImputer(strategy="median", add_indicator=True), StandardScaler(),
                         LogisticRegression(C=1.0, class_weight="balanced", max_iter=5000))
MODELS = {"LogReg_L2": (logreg, FEATURES), "HGB": (hgb, FEATURES), "mito_only": (logreg, ["go_mitochondrion"])}
TASKS = {"A": ("concordant", ["concordant", "RNA_only", "PROT_only", "opposite"]),
         "B": ("RNA_only", ["RNA_only", "PROT_only"]),
         "C": (None, ["concordant", "RNA_only", "PROT_only"])}

def task_data(df, task):
    pos, keep = TASKS[task]
    d = df[df.label.isin(keep)].copy()
    d["y"] = (d.label == pos).astype(int) if pos else d.label.map({c: i for i, c in enumerate(keep)})
    return d

def score(y, p):
    if p.ndim == 1:
        return roc_auc_score(y, p), balanced_accuracy_score(y, (p >= 0.5).astype(int))
    return roc_auc_score(y, p, multi_class="ovr", average="macro"), balanced_accuracy_score(y, p.argmax(1))

def boot_ci(y, p, n=300):
    rng = np.random.default_rng(RNG)
    y = np.asarray(y); out = []
    for _ in range(n):
        i = rng.integers(0, len(y), len(y))
        if len(np.unique(y[i])) < len(np.unique(y)):
            continue
        out.append(score(y[i], p[i])[0])
    return np.percentile(out, [2.5, 97.5])

def proba(m, X):
    p = m.predict_proba(X)
    return p[:, 1] if p.shape[1] == 2 else p

def run_cv(d, mname):
    make, cols = MODELS[mname]
    oof = np.zeros((len(d), d.y.nunique())) if d.y.nunique() > 2 else np.zeros(len(d))
    cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=RNG)
    for tr, te in cv.split(d, d.y, d.cv_group):
        m = make().fit(d.iloc[tr][cols], d.y.iloc[tr])
        oof[te] = proba(m, d.iloc[te][cols])
    return oof

def runnable(d, task):
    if task == "C":
        if (d.label == "concordant").sum() < MIN_CONCORDANT_3CLASS:
            return False, f"skipped: concordant <{MIN_CONCORDANT_3CLASS} (binary collapse applies)"
        nclass = 3
    else:
        nclass = 2
    if d.y.nunique() < nclass or d.y.value_counts().min() < MIN_CLASS:
        return False, f"skipped: a class has <{MIN_CLASS} genes"
    return True, ""

rows, oof_store, fitted = [], {}, {}
for week in WEEKS + ["8w_male_only"]:
    for task in ["A", "B", "C"]:
        if task == "C" and week != "8w":
            continue
        sel = (lambda t: (w.tissue == t) & (w.contrast == "M_8w")) if week == "8w_male_only" else \
              (lambda t: (w.tissue == t) & (w.week == week))
        dd = {t: task_data(w[sel(t)], task) for t in TISSUES}
        ok = {t: runnable(dd[t], task) for t in TISSUES}
        for mname in MODELS:
            base_row = dict(week=week, task=task, model=mname)
            for t in TISSUES:
                mn = int(dd[t].y.value_counts().min()) if len(dd[t]) else 0
                if not ok[t][0]:
                    rows.append({**base_row, "eval": f"CV {t}", "n": len(dd[t]), "min_class_n": mn, "note": ok[t][1]})
                    continue
                p = run_cv(dd[t], mname)
                auc, bac = score(dd[t].y.values, p); lo, hi = boot_ci(dd[t].y.values, p)
                rows.append({**base_row, "eval": f"CV {t}", "n": len(dd[t]), "min_class_n": mn,
                             "AUC": auc, "AUC_lo": lo, "AUC_hi": hi, "bal_acc": bac})
                oof_store[(week, task, mname, t)] = (dd[t].index, p)
            for tr, te in XT:
                ev = f"{tr} -> {te}"
                if not (ok[tr][0] and ok[te][0]):
                    rows.append({**base_row, "eval": ev, "n": len(dd[te]),
                                 "note": (ok[tr][1] or ok[te][1]) + " (train or test)"})
                    continue
                make, cols = MODELS[mname]
                m = make().fit(dd[tr][cols], dd[tr].y)
                p = proba(m, dd[te][cols])
                auc, bac = score(dd[te].y.values, p); lo, hi = boot_ci(dd[te].y.values, p)
                rows.append({**base_row, "eval": ev, "n": len(dd[te]), "min_class_n": int(dd[te].y.value_counts().min()),
                             "AUC": auc, "AUC_lo": lo, "AUC_hi": hi, "bal_acc": bac})
                fitted[(week, task, mname, tr)] = m
        for t in TISSUES:  # majority-class baseline: constant score -> AUC 0.5, balanced accuracy 0.5
            if len(dd[t]):
                rows.append(dict(week=week, task=task, model="majority", eval=f"any, test {t}", n=len(dd[t]),
                                 AUC=0.5, bal_acc=0.5, note="constant prediction; 0.5 by construction"))
metrics = pd.DataFrame(rows)
metrics.to_csv(f"{TAB}/metrics_classification.csv", index=False)
print(metrics[metrics.week.str.startswith("8w")].round(3).to_string())

# permutation importance (drop in cross-tissue AUC), LogReg (primary) and HGB, 8w
imp_rows = []
for task in ["A", "B"]:
    for mname in ["LogReg_L2", "HGB"]:
        for tr, te in XT:
            key = ("8w", task, mname, tr)
            if key not in fitted:
                continue
            d = task_data(w[(w.tissue == te) & (w.week == "8w")], task)
            pi = permutation_importance(fitted[key], d[FEATURES], d.y, scoring="roc_auc", n_repeats=20, random_state=RNG)
            for f, mu, sd in zip(FEATURES, pi.importances_mean, pi.importances_std):
                imp_rows.append(dict(task=task, model=mname, train=tr, test=te, feature=f, dAUC_mean=mu, dAUC_sd=sd))
imp = pd.DataFrame(imp_rows)
imp.to_csv(f"{TAB}/permutation_importance_8w.csv", index=False)
top5 = {t: imp[(imp.task == t) & (imp.model == "LogReg_L2")].groupby("feature").dAUC_mean.mean()
             .sort_values(ascending=False).head(5).index.tolist() for t in ["A", "B"]}
lcoef = []
for task in ["A", "B"]:
    for tr in TISSUES:
        pipe = fitted[("8w", task, "LogReg_L2", tr)]
        names = list(FEATURES) + [f"{FEATURES[i]}_missing" for i in pipe[0].indicator_.features_]
        for n_, b in zip(names, pipe[-1].coef_[0]):
            lcoef.append(dict(task=task, train=tr, feature=n_, std_coef=b))
pd.DataFrame(lcoef).to_csv(f"{TAB}/logreg_coefficients_8w.csv", index=False)

# partial dependence of the top-5 (LogReg importance) features, both models, both training tissues
pdp = {}
for task, feats in top5.items():
    for mname in ["LogReg_L2", "HGB"]:
        for tr in TISSUES:
            m = fitted[("8w", task, mname, tr)]
            d = task_data(w[(w.tissue == tr) & (w.week == "8w")], task)
            for f in feats:
                vals = d[f].dropna()
                grid = np.unique(vals) if vals.nunique() <= 10 else np.unique(np.quantile(vals, np.linspace(.05, .95, 20)))
                X = d[FEATURES].copy(); ys = []
                for g in grid:
                    X[f] = g
                    ys.append(m.predict_proba(X)[:, 1].mean())
                pdp[(task, mname, tr, f)] = (grid, np.array(ys))
pd.DataFrame([dict(task=k[0], model=k[1], train=k[2], feature=k[3], x=x, mean_pred=y)
              for k, (xs, ys) in pdp.items() for x, y in zip(xs, ys)]).to_csv(f"{TAB}/partial_dependence_8w.csv", index=False)

# raw RNA-only fraction by gene property, tissue x sex, 8w (no model)
d8 = w[(w.week == "8w") & w.label.isin(["RNA_only", "PROT_only"])].copy()
d8["y"] = (d8.label == "RNA_only").astype(int)
d8["prot_miss"] = (d8.prot_n_missing > 0).astype(int)
PROPS = [("RNA baseline: top vs bottom quartile", "rna_log2_baseline"), ("GTEx TPM: top vs bottom quartile", "gtex_log2_tpm"),
         ("proteomics missing values >0", "prot_miss"), ("MitoCarta3.0", "mitocarta3"),
         ("stable-complex subunit", "stable_complex"), ("GO:CC extracellular proxy", "go_extracellular")]
f6 = []
for t in TISSUES:
    for sx in ["female", "male"]:
        g = d8[(d8.tissue == t) & (d8.sex == sx)]
        for lab, col in PROPS:
            if col in ("rna_log2_baseline", "gtex_log2_tpm"):
                q1, q3 = g[col].quantile([.25, .75])
                a, b = g[g[col] >= q3], g[g[col] <= q1]
            else:
                a, b = g[g[col] == 1], g[g[col] == 0]
            f6.append(dict(tissue=t, sex=sx, property=lab, rate_with=a.y.mean(), n_with=len(a),
                           rate_without=b.y.mean(), n_without=len(b), n_rna_only=int(g.y.sum()), n_prot_only=int((1 - g.y).sum())))
f6 = pd.DataFrame(f6)
f6.to_csv(f"{TAB}/rna_only_rate_by_property_8w.csv", index=False)

# ================================================================== catalog
cat = w.copy()
for task, col in [("A", "p_concordant"), ("B", "p_RNA_only")]:
    for mname, sfx in [("LogReg_L2", "lr"), ("HGB", "hgb")]:
        cat[f"{col}_xtissue_{sfx}"] = np.nan
        cat[f"{col}_cv_oof_{sfx}"] = np.nan
        for week in WEEKS:
            for tr, te in XT:
                key = (week, task, mname, tr)
                if key in fitted:
                    idx = cat[(cat.tissue == te) & (cat.week == week)].index
                    cat.loc[idx, f"{col}_xtissue_{sfx}"] = fitted[key].predict_proba(cat.loc[idx, FEATURES])[:, 1]
            for t in TISSUES:
                k = (week, task, mname, t)
                if k in oof_store:
                    ix, p = oof_store[k]
                    cat.loc[ix, f"{col}_cv_oof_{sfx}"] = p
for mname, sfx in [("Ridge", "ridge"), ("HGB", "hgb")]:
    cat[f"resid_pred_xtissue_{sfx}"] = np.nan
    for (week, mn, te), s in resid_pred.items():
        if mn == mname:
            cat.loc[s.index, f"resid_pred_xtissue_{sfx}"] = s
pcols = [c for c in cat.columns if c.startswith(("p_concordant_", "p_RNA_only_", "resid_pred_"))]
catcols = ["tissue", "contrast", "sex", "week", "gene_symbol_human", "feature_id_RNA", "feature_id_PROT",
           "logFC_RNA", "stat_RNA", "fdr_bh_RNA", "logFC_PROT", "stat_PROT", "fdr_bh_PROT", "label", "label_binary",
           "resid_prot_on_rna", "cv_group", "rna_baseline_count", "stable_complex", "abund_tertile", "phospho_bin"] + \
          FEATURES + pcols
cat[catcols].sort_values(["tissue", "contrast", "gene_symbol_human"]).to_csv(f"{OUT}/catalog_v1.csv", index=False)
print(f"catalog rows: {len(cat)}")

# ================================================================== figures
SURF, INK, INK2, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
S1, S2, S3, S4 = "#2a78d6", "#eb6834", "#1baf7a", "#4a3aa7"
C = {"concordant": S1, "RNA_only": S2, "PROT_only": S3, "opposite": INK2}
MC = {"LogReg_L2": S1, "Ridge": S1, "HGB": S2, "mito_only": S3}
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "axes.edgecolor": INK2,
                     "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "text.color": INK,
                     "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True})
EVALS = ["CV SKM-GN", "CV HEART", "SKM-GN -> HEART", "HEART -> SKM-GN"]
ECOL = dict(zip(EVALS, [S1, S2, S3, S4]))
TNAME = {"A": "A: concordant vs discordant", "B": "B: RNA-only vs PROT-only"}

# fig1 class counts
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, t in zip(axes, TISSUES):
    c = counts.loc[t].reindex(CONTRASTS); x = np.arange(len(c)); bottom = np.zeros(len(c))
    for k in ["concordant", "RNA_only", "PROT_only", "opposite"]:
        ax.bar(x, c[k], bottom=bottom, color=C[k], width=0.7, label=k, edgecolor=SURF, linewidth=1)
        bottom += c[k].values
    for xi, b in zip(x, bottom):
        ax.text(xi, b + 8, int(b), ha="center", va="bottom", fontsize=7, color=INK2)
    ax.set_xticks(x); ax.set_xticklabels(c.index)
    ax.set_title(f"rat {t}: genes significant in >=1 layer (BH FDR<0.05)", fontsize=9, loc="left")
axes[0].set_ylabel("genes"); axes[1].legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(f"{FIG}/fig1_class_counts.png", dpi=150); plt.close(fig)

# stratified_rho: rows = strata, cols = tissue; x = contrast; one series per stratum level
LEVCOL = [S1, S2, S3]
fig, axes = plt.subplots(len(STRATA), 2, figsize=(11, 13.5), sharex=True, sharey=True)
for r, (sname, (col, levels, labs, ref)) in enumerate(STRATA.items()):
    for cix, t in enumerate(TISSUES):
        ax = axes[r, cix]
        a = strat[(strat.tissue == t) & (strat.stratum == "all genes")].set_index("contrast").reindex(CONTRASTS)
        ax.plot(range(len(CONTRASTS)), a.rho, color="#9a9893", lw=1, ls="--", label="all genes", zorder=1)
        for k, lab in enumerate(labs):
            q = strat[(strat.tissue == t) & (strat.stratum == sname) & (strat.level == lab)].set_index("contrast").reindex(CONTRASTS)
            x = np.arange(len(CONTRASTS)) + (k - (len(labs) - 1) / 2) * 0.18
            ax.errorbar(x, q.rho, yerr=[q.rho - q.lo, q.hi - q.rho], fmt="o", ms=5, color=LEVCOL[k], elinewidth=1.2,
                        capsize=0, label=f"{lab} (median n={int(q.n.median())})", zorder=3)
        ax.axhline(0, color=INK2, lw=0.8)
        ax.axvline(3.5, color=GRID, lw=1.5)
        ax.set_title(f"{t} | {sname}", fontsize=8.5, loc="left")
        ax.legend(frameon=False, fontsize=7, loc="upper left", ncol=2)
        if cix == 0:
            ax.set_ylabel("Spearman rho\n(RNA stat, PROT stat)", fontsize=8)
for ax in axes[-1]:
    ax.set_xticks(range(len(CONTRASTS))); ax.set_xticklabels(CONTRASTS)
fig.tight_layout(); fig.savefig(f"{FIG}/stratified_rho.png", dpi=150); plt.close(fig)

# fig3 residual R2 by week (Ridge primary, HGB secondary)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, mname in zip(axes, ["Ridge", "HGB"]):
    s = rmetrics[(rmetrics.model == mname) & rmetrics.week.isin(WEEKS)]
    for k, e in enumerate(EVALS):
        q = s[s["eval"] == e].set_index("week").reindex(WEEKS)
        x = np.arange(len(WEEKS)) + (k - 1.5) * 0.08
        ax.errorbar(x, q.R2, yerr=[q.R2 - q.R2_lo, q.R2_hi - q.R2], fmt="-o", ms=5, lw=1.5, color=ECOL[e],
                    elinewidth=1, capsize=0, label=e)
    ax.axhline(0, color=INK2, lw=1, ls="--")
    ax.set_xticks(range(len(WEEKS))); ax.set_xticklabels(WEEKS); ax.set_xlabel("weeks of training (sexes pooled)")
    ax.set_title(f"residual PROT-on-RNA stat, {mname}" + (" (primary)" if mname == "Ridge" else " (secondary)"),
                 fontsize=9, loc="left")
axes[0].set_ylabel("out-of-sample R² (95% bootstrap CI)"); axes[0].legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig3_residual_r2_by_week.png", dpi=150); plt.close(fig)

# fig4 residual importance, 8w
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
for ax, mname in zip(axes, ["Ridge", "HGB"]):
    s = rimp[rimp.model == mname]
    order = s.groupby("feature").dR2_mean.mean().sort_values().index
    y = np.arange(len(order))
    for i, (tr, te) in enumerate(XT):
        q = s[s.train == tr].set_index("feature").reindex(order)
        ax.errorbar(q.dR2_mean, y + (i - .5) * 0.3, xerr=q.dR2_sd, fmt="o", ms=5, color=[S1, S2][i],
                    elinewidth=1.2, capsize=0, label=f"train {tr}, test {te}")
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("drop in test R² when feature is permuted")
    ax.set_title(f"residual target, {mname}, 8w, cross-tissue", fontsize=9, loc="left")
axes[1].legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(f"{FIG}/fig4_residual_importance_8w.png", dpi=150); plt.close(fig)

# fig5 classification AUC at 8w (pooled) and 8w males only
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
for r, week in enumerate(["8w", "8w_male_only"]):
    for cix, task in enumerate(["A", "B"]):
        ax = axes[r, cix]
        m8 = metrics[(metrics.week == week) & (metrics.task == task)]
        for i, mname in enumerate(["LogReg_L2", "HGB", "mito_only"]):
            s = m8[m8.model == mname].set_index("eval").reindex(EVALS)
            x = np.arange(len(EVALS)) + (i - 1) * 0.22
            ax.errorbar(x, s.AUC, yerr=[s.AUC - s.AUC_lo, s.AUC_hi - s.AUC], fmt="o", ms=6, color=MC[mname],
                        elinewidth=1.5, capsize=0, label=mname)
        ax.axhline(0.5, color=INK2, lw=1, ls="--", label="majority / chance")
        ax.set_xticks(range(len(EVALS))); ax.set_xticklabels([e.replace(" -> ", "\n-> ") for e in EVALS], fontsize=8)
        ax.set_title(f"{TNAME[task]} | {'8w, sexes pooled' if week == '8w' else '8w, males only'}", fontsize=9, loc="left")
        ax.set_ylim(0.3, 1.0)
    axes[r, 0].set_ylabel("ROC AUC (95% bootstrap CI)")
axes[0, 0].legend(frameon=False, fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig5_auc_8w.png", dpi=150); plt.close(fig)

# fig6 permutation importance for classification: mean over the two cross-tissue directions, bar = range
fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
for ax, task in zip(axes, ["A", "B"]):
    s = imp[imp.task == task]
    order = s[s.model == "LogReg_L2"].groupby("feature").dAUC_mean.mean().sort_values().index
    y = np.arange(len(order))
    for i, mname in enumerate(["LogReg_L2", "HGB"]):
        g = s[s.model == mname].groupby("feature").dAUC_mean.agg(["mean", "min", "max"]).reindex(order)
        ax.errorbar(g["mean"], y + (i - .5) * 0.3, xerr=[g["mean"] - g["min"], g["max"] - g["mean"]], fmt="o", ms=5,
                    color=MC[mname], elinewidth=1.2, capsize=0, label=mname)
    ax.axvline(0, color=INK2, lw=1)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=8)
    ax.set_xlabel("drop in cross-tissue test AUC when permuted\n(dot = mean of both directions, bar = range)")
    ax.set_title(f"{TNAME[task]}, 8w", fontsize=9, loc="left")
axes[1].legend(frameon=False, fontsize=7, loc="lower right")
fig.tight_layout(); fig.savefig(f"{FIG}/fig6_permutation_importance_8w.png", dpi=150); plt.close(fig)

# fig7a/b partial dependence
for task in ["A", "B"]:
    feats = top5[task]
    fig, axes = plt.subplots(1, 5, figsize=(13, 3.1), sharey=True)
    for ax, f in zip(axes, feats):
        for mname in ["LogReg_L2", "HGB"]:
            for tr, ls in zip(TISSUES, ["-", "--"]):
                xs, ys = pdp[(task, mname, tr, f)]
                ax.plot(xs, ys, ls, marker="o" if len(xs) <= 10 else None, ms=5, lw=2, color=MC[mname],
                        label=f"{mname}, trained on {tr}")
        ax.set_xlabel(f, fontsize=8)
    axes[0].set_ylabel("mean P(concordant)" if task == "A" else "mean P(RNA-only)")
    axes[0].legend(frameon=False, fontsize=6.5)
    fig.suptitle(f"{TNAME[task]}: partial dependence at 8w on the top-5 features (LogReg permutation importance)",
                 fontsize=9, x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(f"{FIG}/fig7{'ab'[task == 'B']}_partial_dependence_task{task}_8w.png", dpi=150); plt.close(fig)

# fig8 classification AUC over training duration (LogReg primary)
fig, axes = plt.subplots(1, 2, figsize=(10, 3.6), sharey=True)
for ax, task in zip(axes, ["A", "B"]):
    s = metrics[(metrics.task == task) & (metrics.model == "LogReg_L2") & metrics.week.isin(WEEKS)]
    for k, e in enumerate(EVALS):
        q = s[s["eval"] == e].set_index("week").reindex(WEEKS)
        x = np.arange(len(WEEKS)) + (k - 1.5) * 0.08
        ax.errorbar(x, q.AUC, yerr=[q.AUC - q.AUC_lo, q.AUC_hi - q.AUC], fmt="-o", ms=5, lw=1.5, color=ECOL[e],
                    elinewidth=1, capsize=0, label=e)
    ax.axhline(0.5, color=INK2, lw=1, ls="--")
    ax.set_xticks(range(len(WEEKS))); ax.set_xticklabels(WEEKS); ax.set_xlabel("weeks of training (sexes pooled)")
    ax.set_title(TNAME[task] + " (LogReg)", fontsize=9, loc="left"); ax.set_ylim(0.3, 1.0)
axes[0].set_ylabel("ROC AUC (95% bootstrap CI)"); axes[1].legend(frameon=False, fontsize=7, loc="lower left")
fig.tight_layout(); fig.savefig(f"{FIG}/fig8_auc_by_week.png", dpi=150); plt.close(fig)

# fig9 raw RNA-only fraction by property, tissue x sex
fig, axes = plt.subplots(1, 4, figsize=(14, 3.8), sharey=True)
for ax, ((t, sx), g) in zip(axes, f6.groupby(["tissue", "sex"], sort=False)):
    for yi, r in enumerate(g.itertuples()):
        ax.plot([r.rate_with, r.rate_without], [yi, yi], color=GRID, lw=2, zorder=2)
        ax.scatter(r.rate_with, yi, s=40, color=S2, zorder=3, label="has property / top quartile" if yi == 0 else None)
        ax.scatter(r.rate_without, yi, s=40, color="#9a9893", zorder=3, label="lacks property / bottom quartile" if yi == 0 else None)
        ax.text(1.04, yi, f"n={r.n_with}", va="center", ha="left", fontsize=7, color=INK2)
    ax.set_yticks(range(len(g))); ax.set_yticklabels(g.property, fontsize=8); ax.set_xlim(0, 1)
    ax.set_title(f"{t} {sx}, 8w\n{g.n_rna_only.iloc[0]} RNA-only / {g.n_prot_only.iloc[0]} PROT-only genes", fontsize=8, loc="left")
    ax.set_xlabel("fraction RNA-only", fontsize=8)
axes[0].invert_yaxis()
axes[0].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(0, -0.22), ncol=2)
fig.tight_layout(); fig.savefig(f"{FIG}/fig9_rna_only_rate_by_property_8w.png", dpi=150); plt.close(fig)
print("done")
