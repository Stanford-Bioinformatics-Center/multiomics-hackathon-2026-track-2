"""Scoring shared by every command: does exercise move each signature gene against the disease (agreement,
sign test), and is the whole set shifted (signed cameraPR, + = opposed)?

Method details: docs/METHODS.md#corepy
"""
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import binomtest, false_discovery_control, spearmanr

from . import legacy, store

CAP = 4.0
INTER_GENE_COR = 0.01
CORE_TISSUES = ["VL", "SKM-GN", "SKM-VL", "HEART"]
CORE_LAYERS = ["RNA", "PROT"]
EXERCISE_KINDS = ["exercise vs control", "training vs sedentary"]
TISSUE_WORDS = {"VL": "vastus lateralis", "SKM-GN": "gastrocnemius", "SKM-VL": "vastus lateralis",
                "HEART": "heart", "ADIPOSE": "subcutaneous adipose", "BLOOD": "blood", "WAT-SC": "subcutaneous WAT",
                "BAT": "brown adipose", "LIVER": "liver", "KIDNEY": "kidney", "LUNG": "lung", "CORTEX": "cortex",
                "ADRNL": "adrenal", "COLON": "colon", "HIPPOC": "hippocampus", "HYPOTH": "hypothalamus",
                "OVARY": "ovary", "SMLINT": "small intestine", "SPLEEN": "spleen", "TESTES": "testes",
                "VENACV": "vena cava", "PLASMA": "plasma"}
LAYER_WORDS = {"RNA": "RNA", "PROT": "protein", "PHOSPHO": "phosphosite", "PROT_OLINK": "plasma protein (Olink)",
               "METAB": "metabolite"}
H_TIME_WORDS = {"15to45min": "15–45 min", "3.5to4h": "3.5–4 h", "24h": "24 h", "10min": "10 min",
                "during20min": "20 min during", "during40min": "40 min during", "pre": "baseline"}


def bh(p):
    return legacy.bh(p)


def col_words(c):
    """Human-readable label for a column index row (e.g. 'Rat gastrocnemius protein, female 8 wk')."""
    sp = "Human" if c.dataset == "human_acute" else "Rat"
    tis = TISSUE_WORDS.get(c.tissue, c.tissue)
    lay = LAYER_WORDS.get(c.layer, c.layer)
    if c.dataset == "rat_train":
        when = f"{'female' if c.sex == 'F' else 'male'} {c.time[:-1]} wk"
    else:
        when = c.label.replace("−", "–")
        for k, v in {"15-45m": "15–45 min", "3.5h": "3.5–4 h", "24h": "24 h", "10m": "10 min",
                     "20m-during": "20 min during", "40m-during": "40 min during"}.items():
            if when.endswith(k):
                when = when[: -len(k)] + v
    return f"{sp} {tis} {lay}, {when}"


def short_words(c):
    """Compact column label for figure axes: 'EE 24 h' / 'F 8 wk'."""
    if c.dataset == "rat_train":
        return f"{c.sex} {c.time[:-1]} wk"
    return col_words(c).split(", ", 1)[1].replace("EE–CON ", "EE ").replace("RE–CON ", "RE ")


class Store:
    """Array view of the parquet store: one contiguous slice of rows per comparison column."""

    def __init__(self, extra_rows=None):
        r = store.load()
        if extra_rows is not None:
            r = pd.concat([r, extra_rows], ignore_index=True)
        self.rows = r
        self.cols = store.load_columns().set_index("column_id", drop=False)
        g = r.gene_symbol_human
        self.genes = g.cat.categories
        self.gcode = g.cat.codes.to_numpy()
        self.stat = r.stat.to_numpy(float)
        self.lfc = r.logFC.to_numpy(float)
        self.fdr = r.fdr_bh.to_numpy(float)
        self.base = r.baseline_expr.to_numpy(float)
        self.miss = r.prot_n_missing.to_numpy(float)
        self.ncoll = r.n_collapsed.to_numpy(float)
        self.pos = r.groupby("column_id", observed=True).indices
        self.ann = store.load_annotations()
        self.gs = store.load_genesets()

    def codes(self, genes):
        return self.genes.get_indexer(pd.Index(genes))

    def col(self, cid):
        return self.pos[cid]

    def frame(self, cid):
        return self.rows.iloc[self.pos[cid]]

    def select(self, tissues=None, layers=None, kinds=None, datasets=None):
        c = self.cols
        m = np.ones(len(c), bool)
        for field_, vals in [("tissue", tissues), ("layer", layers), ("kind", kinds), ("dataset", datasets)]:
            if vals is not None:
                m &= c[field_].isin(vals).to_numpy()
        return c[m].sort_values("col_order")


_STORE = {}


def get_store():
    if "s" not in _STORE:
        _STORE["s"] = Store()
    return _STORE["s"]


# --------------------------------------------------------------------------------------------- per column
def camera_signed(stat, member, signs):
    """cameraPR t of `member` after multiplying member stats by `signs` (so + = opposed). (t, p) or NaN."""
    m = int(member.sum())
    if m < 2 or m >= len(stat) - 2:
        return np.nan, np.nan
    s = stat.copy()
    s[member] = s[member] * signs
    t, _, p = legacy.camera_pr(s, member, INTER_GENE_COR)
    return float(t), float(p)


def camera_vec(stat, sets, signs=None, inter_gene_cor=INTER_GENE_COR):
    """Vectorised cameraPR t for B random sets (rows of `sets`, integer positions) with member stats multiplied by
    `signs` (length m). Algebraically identical to legacy.camera_pr applied to each flipped vector."""
    G = len(stat)
    B, m = sets.shape
    if m < 2:
        return np.full(B, np.nan)
    x = stat[sets]
    xs = x * signs if signs is not None else x
    tot = stat.sum() - x.sum(1) + xs.sum(1)
    mean_ = tot / G
    var_ = ((stat ** 2).sum() - G * mean_ ** 2) / (G - 1)
    delta = G / (G - m) * (xs.mean(1) - mean_)
    vp = ((G - 1) * var_ - delta ** 2 * m * (G - m) / G) / (G - 2)
    vif = 1 + (m - 1) * inter_gene_cor
    return delta / np.sqrt(vp * (vif / m + 1 / (G - m)))


@dataclass
class Sig:
    """A mapped signature: counted genes (used for statistics) and shown genes (incl. context groups)."""
    name: str
    table: pd.DataFrame            # full input table with mapping status
    genes: list                    # counted human symbols (unique)
    dirs: np.ndarray               # +1 / -1 aligned with genes
    shown: list = field(default_factory=list)       # counted + context genes, display order
    shown_dirs: np.ndarray = None
    groups: dict = field(default_factory=dict)      # gene -> group label


def score_column(S, cid, sig, cap=CAP, cutoff=0.05, universe=None):
    """Counts, sign test and signed cameraPR for one column. `universe`: optional set of genes restricting the
    ranking (e.g. muscle-intrinsic)."""
    pos = S.col(cid)
    g = S.gcode[pos]
    stat, lfc, fdr = S.stat[pos], S.lfc[pos], S.fdr[pos]
    if universe is not None:
        keep = np.isin(g, S.codes(sorted(universe)))
        g, stat, lfc, fdr = g[keep], stat[keep], lfc[keep], fdr[keep]
    dmap = np.zeros(len(S.genes) + 1)
    codes = S.codes(sig.genes)
    ok = codes >= 0
    dmap[codes[ok]] = sig.dirs[ok]
    d = dmap[g]
    member = d != 0
    agree = pd.Series(np.where(member, np.sign(lfc) * d, np.nan))
    n, k, p = legacy.sign_test(agree)
    a = agree.dropna()
    a = a[a != 0]
    out = dict(column_id=cid, n_universe=len(g), n_measured=n, n_opposed=k, n_same=int((a > 0).sum()),
               frac_opposed=k / n if n else np.nan, sign_p=p,
               n_sig=int((member & (fdr < cutoff)).sum()),
               n_sig_opposed=int((member & (fdr < cutoff) & (np.sign(lfc) * d < 0)).sum()),
               n_sig_same=int((member & (fdr < cutoff) & (np.sign(lfc) * d > 0)).sum()))
    signs = -d[member]
    out["camera_t"], out["camera_p"] = camera_signed(stat, member, signs)
    up, dn = member & (d > 0), member & (d < 0)
    out["n_up"], out["n_down"] = int(up.sum()), int(dn.sum())
    out["camera_t_up"], out["camera_p_up"] = camera_signed(stat, up, -np.ones(up.sum()))
    out["camera_t_down"], out["camera_p_down"] = camera_signed(stat, dn, np.ones(dn.sum()))
    return out


SCORE_FIELDS = ["column_id", "n_universe", "n_measured", "n_opposed", "n_same", "frac_opposed", "sign_p", "n_sig",
                "n_sig_opposed", "n_sig_same", "camera_t", "camera_p", "n_up", "n_down", "camera_t_up", "camera_p_up",
                "camera_t_down", "camera_p_down"]


def score_columns(S, cids, sig, **kw):
    df = pd.DataFrame([score_column(S, c, sig, **kw) for c in cids], columns=SCORE_FIELDS)
    for col in ["camera_p", "camera_p_up", "camera_p_down"]:
        q = np.full(len(df), np.nan)
        ok = df[col].notna().to_numpy()
        if ok.any():
            q[ok] = false_discovery_control(df.loc[ok, col], method="bh")
        df[col.replace("_p", "_fdr")] = q
    return df


def grid_long(S, cids, sig, cap=CAP):
    """Gene x column long table over the shown genes (NA rows kept)."""
    rows = []
    codes = dict(zip(sig.shown, S.codes(sig.shown)))
    dirs = dict(zip(sig.shown, sig.shown_dirs))
    for cid in cids:
        pos = S.col(cid)
        g = S.gcode[pos]
        where = {c: i for i, c in enumerate(g)} if len(sig.shown) > 40 else None
        for gene in sig.shown:
            c = codes[gene]
            i = (where.get(c) if where is not None else (np.flatnonzero(g == c)[0] if (g == c).any() else None)) \
                if c >= 0 else None
            if i is None:
                rows.append(dict(column_id=cid, gene=gene, direction=dirs[gene], logFC=np.nan, stat=np.nan,
                                 fdr_bh=np.nan, n_collapsed=np.nan))
            else:
                p = pos[i]
                rows.append(dict(column_id=cid, gene=gene, direction=dirs[gene], logFC=S.lfc[p], stat=S.stat[p],
                                 fdr_bh=S.fdr[p], n_collapsed=S.ncoll[p]))
    L = pd.DataFrame(rows)
    L["agreement"] = np.sign(L.logFC) * L.direction
    L["cell"] = L.agreement * np.minimum(L.stat.abs(), cap if cap else np.inf)
    L["group"] = L.gene.map(sig.groups)
    return L


# --------------------------------------------------------------------------------------------- effective n
def effective_n(S, genes, cids, min_periods=10):
    """NARRATIVE R3: Spearman correlation of the genes' stats across columns -> Nyholt and Li-Ji effective
    number of independent genes. Genes measured in < min_periods columns are dropped; NaN correlations -> 0."""
    mat = {}
    for cid in cids:
        f = S.frame(cid)
        mat[cid] = f.set_index("gene_symbol_human").stat.reindex(genes).to_numpy()
    mat = pd.DataFrame(mat, index=genes)
    mat = mat[mat.notna().sum(axis=1) >= min_periods]
    M = len(mat)
    if M < 2:
        return dict(n_genes=M, meff_nyholt=float(M), meff_liji=float(M), mean_r=np.nan, top_eigen_frac=np.nan)
    corr = mat.T.corr(method="spearman", min_periods=min_periods).fillna(0.0)
    lam = np.clip(np.linalg.eigvalsh(corr.to_numpy()), 0, None)
    meff_nyholt = 1 + (M - 1) * (1 - np.var(lam, ddof=1) / M)
    meff_liji = float(np.sum((lam >= 1) + (lam - np.floor(lam))))
    off = corr.to_numpy()[np.triu_indices(M, 1)]
    return dict(n_genes=M, meff_nyholt=float(meff_nyholt), meff_liji=meff_liji, mean_r=float(off.mean()),
                min_r=float(off.min()), max_r=float(off.max()), top_eigen_frac=float(lam.max() / lam.sum()))


def deflated_sign_p(n, k, frac):
    """Sign test with n and k scaled by the effective-number fraction (Meff / M)."""
    if not n or not np.isfinite(frac):
        return np.nan
    ne = max(1, int(round(n * frac)))
    ke = int(round(k * frac))
    return binomtest(min(ke, ne), ne, 0.5).pvalue


def spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 5:
        return np.nan, np.nan, int(ok.sum())
    r, p = spearmanr(x[ok], y[ok])
    return float(r), float(p), int(ok.sum())


# --------------------------------------------------------------------------------------------- time course
RAT_WEEKS = [("1w", 1), ("2w", 2), ("4w", 4), ("8w", 8)]
HUMAN_TIMES = [("15to45min", "15–45 min"), ("3.5to4h", "3.5–4 h"), ("24h", "24 h")]


def trajectory_table(S, genes, tissue="SKM-GN", layers=("RNA", "PROT"), species=None, cutoff=0.05):
    """Each gene's exercise effect over time, per layer, with a 95% CI (SE = |logFC / stat| from the published test).
    Rat: trained vs sedentary at 1/2/4/8 weeks, female and male pooled by inverse-variance weighting (sedentary = 0).
    Human: exercise vs control at 15-45 min, 3.5-4 h and 24 h, endurance (EE) and resistance (RE).
    `significant` = BH FDR < cutoff in any contributing comparison."""
    rows = []
    rat = (species or ("human" if tissue == "VL" else "rat")) == "rat"
    for layer in layers:
        if rat:
            points = [(w, x, [f"rat_train|{tissue}|{layer}|{s}_{w}" for s in "FM"], "trained, F+M pooled")
                      for w, x in RAT_WEEKS]
        else:
            points = [(lab, i, [f"human_acute|{tissue}|{layer}|{g}_vs_CON_{t}"], g)
                      for g in ("EE", "RE") for i, (t, lab) in enumerate(HUMAN_TIMES)]
        for xlab, x, cids, series in points:
            parts = [S.frame(c) for c in cids if c in S.pos]
            if not parts:
                continue
            d = pd.concat(parts)
            d = d[d.gene_symbol_human.isin(genes)][["gene_symbol_human", "logFC", "stat", "fdr_bh"]].copy()
            d["se"] = (d.logFC / d.stat).abs()
            for g, dd in d.dropna(subset=["se"]).groupby("gene_symbol_human", observed=True):
                w = 1 / dd.se.clip(lower=1e-6) ** 2
                m, se = float((w * dd.logFC).sum() / w.sum()), float(1 / np.sqrt(w.sum()))
                rows.append(dict(gene=str(g), layer=layer, x=x, time=xlab, series=series, logFC=m, se=se,
                                 lo=m - 1.96 * se, hi=m + 1.96 * se, fdr_min=float(dd.fdr_bh.min()),
                                 significant=bool((dd.fdr_bh < cutoff).any())))
    return pd.DataFrame(rows)


def lookup_table(T, dirs):
    """Readable gene x time table from trajectory_table: '+0.46*' (* = significant at the chosen FDR), plus the
    direction relative to the disease at the last time point (opposed / same as disease)."""
    if not len(T):
        return pd.DataFrame()
    t = T.copy()
    t["cell"] = [f"{v:+.2f}{'*' if s else ''}" for v, s in zip(t.logFC, t.significant)]
    t["col"] = (t.series.replace({"trained, F+M pooled": ""}) + " " + t.time).str.strip()
    order = list(dict.fromkeys(t.sort_values(["series", "x"]).col))
    w = t.pivot_table(index=["gene", "layer"], columns="col", values="cell", aggfunc="first")[order].reset_index()
    last = t.sort_values(["series", "x"]).groupby(["gene", "layer"]).tail(1).set_index(["gene", "layer"])
    lab = []
    for g, l in zip(w.gene, w.layer):
        v = last.loc[(g, l)]
        opp = np.sign(v.logFC) * dirs.get(g, 0) < 0
        lab.append(("opposed" if opp else "same as disease") + ("" if v.significant else " (n.s.)"))
    w[f"Direction at {order[-1]}"] = lab
    w["layer"] = w.layer.map(LAYER_WORDS).fillna(w.layer)
    return w.rename(columns={"gene": "Gene", "layer": "Layer"})
