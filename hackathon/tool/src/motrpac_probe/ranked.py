"""Full disease rankings (gene + t), generalised from the team's PAH blood pipeline: whole-ranking Spearman
association and GO:BP pathway concordance against every MoTrPAC comparison.

Method details: docs/METHODS.md#rankedpy
"""
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import false_discovery_control, norm, rankdata, t as tdist

from . import signature, store
from .core import INTER_GENE_COR, Sig

STAT_COLS = ["t", "stat", "z", "zscore", "score", "logfc"]
MIN_SET, MAX_SET = 10, 500
GOBP = store.STORE / "gobp_pkg.csv.gz"
CAMERA_PKG = store.STORE / "camera_pkg_gobp.csv.gz"
H_CODE = {"muscle": "VL", "blood": "BLOOD", "adipose": "ADIPOSE"}
ASSAY_LAYER = {"transcript-rna-seq": "RNA", "prot-pr": "PROT", "prot-ol": "PROT_OLINK"}


def is_ranked(path):
    cols = {c.strip().lower() for c in pd.read_csv(path, nrows=0).columns}
    return "direction" not in cols and bool(cols & set(STAT_COLS))


def load_ranked(path, name=None, universe=None, cutoff=0.05, cap=250, exact_symbols=False):
    """-> (ranking mapped to store symbols, derived directional Sig, info). info["full"] = the disease's full
    ranking in its own symbols (disease-side pathway universe)."""
    t = signature.read_table(path)
    sc = next(c for c in STAT_COLS if c in t.columns)
    t["stat_value"] = pd.to_numeric(t[sc], errors="coerce")
    t["direction"] = np.sign(t.stat_value).fillna(0).astype(int).astype(str)
    from .core import get_store
    universe = get_store().genes if universe is None else universe
    # duplicates in a ranking are resolved by max |stat| (the team's collapse rule), not by the list rules
    t = t.reindex(t.stat_value.abs().sort_values(ascending=False).index).reset_index(drop=True)
    for c in ("gene_symbol", "uniprot", "ensembl", "rat_symbol"):
        if c in t.columns:
            t = t[t[c].isna() | ~t[c].duplicated()]
    # the disease's own full ranking (its own gene universe) is kept for the disease-side pathway test
    full = pd.DataFrame({"gene": t.get("gene_symbol", pd.Series(dtype=str)), "stat": t.stat_value})
    full = full.dropna().drop_duplicates("gene").reset_index(drop=True)
    m = signature.map_genes(t, universe)
    if exact_symbols and "gene_symbol" in m.columns:  # the team pipeline: exact symbol match only
        m.loc[m.map_method != "human symbol", ["gene", "status"]] = [np.nan, "unmapped: exact-symbol mode"]
    ok = m[m.status == "mapped"].copy()
    fdr_col = next((c for c in ["adj_p_value", "fdr", "padj", "q_value", "adj.p.val"] if c in m.columns), None)
    rk = pd.DataFrame({"gene": ok.gene, "stat": ok.stat_value,
                       "logFC": pd.to_numeric(ok.get("logfc", pd.Series(np.nan, index=ok.index)), errors="coerce"),
                       "fdr": pd.to_numeric(ok[fdr_col], errors="coerce") if fdr_col else np.nan})
    rk = rk.dropna(subset=["stat"]).reset_index(drop=True)
    # derived directional signature for the gene-set sections
    pick = rk[rk.fdr < cutoff] if fdr_col else rk
    up = pick[pick.stat > 0].nlargest(cap, "stat")
    dn = pick[pick.stat < 0].nsmallest(cap, "stat")
    d = pd.concat([up, dn])
    genes, dirs = d.gene.tolist(), np.sign(d.stat).astype(int).to_numpy()
    sig = Sig(name=name or str(path).rsplit("/", 1)[-1].split(".")[0], table=m, genes=genes, dirs=dirs,
              shown=list(genes), shown_dirs=dirs.copy(),
              groups={g: ("up in disease (top by t)" if s > 0 else "down in disease (top by t)")
                      for g, s in zip(genes, dirs)})
    info = dict(stat_column=sc, fdr_column=fdr_col or "", n_rows=len(t), n_ranked=len(rk),
                n_derived_up=len(up), n_derived_down=len(dn), full=full, n_full=len(full),
                derived_rule=(f"{fdr_col} < {cutoff}" if fdr_col else "no FDR column: top genes by |stat|")
                + f", at most {cap} per direction by |{sc}|")
    return rk, sig, info


# ------------------------------------------------------------------------------------------------ 1. ranks
def _col_stats(S, cid):
    pos = S.col(cid)
    return pd.Series(S.stat[pos], index=S.genes[S.gcode[pos]].astype(str))


def rank_association(S, rk, cids, nperm=10000, perm_cids=None, seed=20260926, chunk=500):
    """Spearman rho per column; gene-label permutation p for columns in perm_cids (default all)."""
    rng = np.random.default_rng(seed)
    dvec = np.full(len(S.genes) + 1, np.nan)  # disease stat looked up by store gene code
    codes = S.codes(rk.gene)
    dvec[codes[codes >= 0]] = rk.stat.to_numpy(float)[codes >= 0]
    perm_cids = set(cids if perm_cids is None else perm_cids)
    rows = []
    for cid in cids:
        pos = S.col(cid)
        dv, ev = dvec[S.gcode[pos]], S.stat[pos]
        ok = np.isfinite(dv) & np.isfinite(ev)
        n = int(ok.sum())
        rec = dict(column_id=cid, n_shared=n)
        if n < 10:
            rows.append(rec)
            continue
        x = rankdata(dv[ok])
        y = rankdata(ev[ok])
        x = (x - x.mean()) / x.std()
        y = (y - y.mean()) / y.std()
        rho = float(np.mean(x * y))
        z, se = np.arctanh(np.clip(rho, -0.999999, 0.999999)), 1 / np.sqrt(n - 3)
        rec.update(rho=rho, ci_low=float(np.tanh(z - 1.96 * se)), ci_high=float(np.tanh(z + 1.96 * se)))
        if cid in perm_cids and nperm:
            hits, done = 0, 0
            while done < nperm:
                k = min(chunk, nperm - done)
                null = rng.permuted(np.broadcast_to(x, (k, n)), axis=1) @ y / n  # shuffled gene labels
                hits += int((np.abs(null) >= abs(rho)).sum())
                done += k
            rec.update(perm_p=(1 + hits) / (nperm + 1), nperm=nperm)
        rows.append(rec)
    df = pd.DataFrame(rows)
    for c in ["perm_p"]:
        if c in df and df[c].notna().any():
            ok = df[c].notna()
            df.loc[ok, "perm_fdr"] = false_discovery_control(df.loc[ok, c], method="bh")
    df["direction"] = np.where(df.get("rho", np.nan) > 0, "same as disease", "opposed")
    return df


def reference_calibration(S, ra, all_ra):
    """For each exercise column: percentile of |rho| among non-exercise reference columns of the same dataset x
    tissue x layer (NaN when fewer than 3 reference columns)."""
    ref = all_ra.merge(S.cols.reset_index(drop=True)[["column_id", "dataset", "tissue", "layer", "kind"]], on="column_id")
    ref = ref[ref.kind == "reference (non-exercise)"]
    out = []
    for r in ra.merge(S.cols.reset_index(drop=True)[["column_id", "dataset", "tissue", "layer"]], on="column_id").itertuples():
        pool = ref[(ref.dataset == r.dataset) & (ref.tissue == r.tissue) & (ref.layer == r.layer)].rho.dropna()
        out.append(dict(column_id=r.column_id, n_reference=len(pool),
                        ref_abs_rho_median=float(pool.abs().median()) if len(pool) else np.nan,
                        pct_vs_reference=100 * float((pool.abs() <= abs(r.rho)).mean()) if len(pool) >= 3 else np.nan))
    return pd.DataFrame(out)


# ------------------------------------------------------------------------------------------------ 2. pathways
def gobp_sets():
    g = pd.read_csv(GOBP)
    return g.groupby("gs_name").gene_symbol.apply(lambda s: sorted(set(s))).to_dict()


_MEMB = {}


def _membership(sets):
    """Global sparse set x gene membership matrix (built once per set collection)."""
    key = id(sets)
    if key not in _MEMB:
        vocab = sorted({g for mem in sets.values() for g in mem})
        vi = {g: i for i, g in enumerate(vocab)}
        names = list(sets)
        rows = [r for r, nm in enumerate(names) for _ in sets[nm]]
        cols = [vi[g] for nm in names for g in sets[nm]]
        M = sparse.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(names), len(vocab)))
        M.data[:] = 1
        _MEMB[key] = (np.array(names), vi, M.tocsc())
    return _MEMB[key]


def camera_sets(stat, genes, sets, min_size=MIN_SET, max_size=MAX_SET):
    """cameraPR (limma port, inter.gene.cor 0.01) for many sets at once. stat/genes: the ranking universe."""
    G = len(stat)
    names_all, vi, Mg = _membership(sets)
    pos = np.array([i for i, g in enumerate(genes) if g in vi], int)
    vcols = np.array([vi[genes[i]] for i in pos], int)
    M = sparse.csr_matrix((Mg.shape[0], G))
    if len(pos):
        sub = Mg[:, vcols].tocoo()
        M = sparse.csr_matrix((sub.data, (sub.row, pos[sub.col])), shape=(Mg.shape[0], G))
    m_all = np.asarray(M.sum(1)).ravel()
    keep = (m_all >= min_size) & (m_all <= max_size)
    if not keep.any():
        return pd.DataFrame(columns=["set", "n", "t", "p", "fdr", "direction"])
    M, names, m = M[keep], names_all[keep], m_all[keep]
    mean_stat, var_stat = stat.mean(), stat.var(ddof=1)
    delta = G / (G - m) * (M @ stat / m - mean_stat)
    var_pooled = ((G - 1) * var_stat - delta ** 2 * m * (G - m) / G) / (G - 2)
    vif = 1 + (m - 1) * INTER_GENE_COR
    t = delta / np.sqrt(var_pooled * (vif / m + 1 / (G - m)))
    p = 2 * tdist.sf(np.abs(t), G - 2)
    return pd.DataFrame({"set": names, "n": m.astype(int), "t": t, "p": p,
                         "fdr": false_discovery_control(p, method="bh"),
                         "direction": np.where(t > 0, "Up", "Down")})


def disease_pathways(rk, sets=None):
    sets = gobp_sets() if sets is None else sets
    return camera_sets(rk.stat.to_numpy(float), rk.gene.tolist(), sets)


def exercise_pathways_store(S, cid, sets):
    e = _col_stats(S, cid)
    e = e[~e.index.duplicated()]
    return camera_sets(e.to_numpy(float), e.index.tolist(), sets)


def exercise_pathways_precomputed(S, cids):
    """Package CAMERA_RESULTS (GOBP) for human columns, keyed by column_id."""
    from .store import ensure_store, human_contrast_code
    ensure_store()
    c = pd.read_csv(CAMERA_PKG)
    c["column_id"] = ("human_acute|" + c.tissue.map(H_CODE) + "|" + c.assay.map(ASSAY_LAYER) + "|"
                      + c.contrast_short.map(lambda s: human_contrast_code(s)))
    c = c[c.column_id.isin(cids)]
    return {cid: pd.DataFrame({"set": d.set, "n": d.set_size, "t": d.t, "p": d.p_value, "fdr": d.adj_p_value,
                               "direction": d.direction}) for cid, d in c.groupby("column_id")}


def pathway_concordance(S, dis, cids, sets=None, exercise_side="store", cutoff=0.05):
    """Join disease-side and exercise-side set results per column. Returns (per-column summary, joined long table)."""
    sets = gobp_sets() if sets is None else sets
    ex = exercise_pathways_precomputed(S, cids) if exercise_side == "precomputed" else \
        {cid: exercise_pathways_store(S, cid, sets) for cid in cids}
    summ, long = [], []
    for cid in cids:
        e = ex.get(cid)
        if e is None or not len(e):
            continue
        j = dis.merge(e, on="set", suffixes=("_disease", "_exercise"))
        both = (j.fdr_disease < cutoff) & (j.fdr_exercise < cutoff)
        same = both & (np.sign(j.t_disease) == np.sign(j.t_exercise))
        rho = float(pd.Series(j.t_disease.to_numpy()).corr(pd.Series(j.t_exercise.to_numpy()), method="spearman")) \
            if len(j) > 10 else np.nan
        summ.append(dict(column_id=cid, exercise_side=exercise_side, shared_sets=len(j), both_bh=int(both.sum()),
                         same_direction=int(same.sum()), opposite_direction=int((both & ~same).sum()),
                         rho_set_t=rho))
        long.append(j.assign(column_id=cid, both_bh=both, concordance=np.where(
            ~both, "", np.where(same, "same direction", "opposite direction"))))
    return pd.DataFrame(summ), (pd.concat(long, ignore_index=True) if long else pd.DataFrame())


def signed_z(p, direction):
    """The team's pah_camera_z: sign(direction) x qnorm(p / 2, lower = FALSE)."""
    return np.where(np.asarray(direction) == "Up", 1, -1) * norm.isf(np.maximum(p, np.finfo(float).tiny) / 2)
