"""Calibration nulls for the opposition score (section D).

(i)  abundance-matched: each counted signature gene measured in the column is replaced by a random non-signature
     gene from the same abundance decile of that column (NARRATIVE R1b: decile of the column's own baseline, i.e.
     human AveExpr / rat sedentary mean count; columns without a baseline use GTEx v8 median TPM of the matched
     human tissue; genes without either form their own stratum). Drawn with replacement, gene by gene, in the
     same order as scripts/05_pah_grid.py R1b, and given the signature's direction vector.
(ii) pathway-class-matched: same, but the stratum is the gene's class = (MitoCarta3.0 member, GO:CC complex
     subunit, secreted/extracellular GO:CC proxy; the three flags of scripts/08_discordance_model.py) crossed with
     the abundance tertile of (i), so it is also abundance-matched. Answers "is it your signature, or its pathway
     class?".
(iii) pool null (optional, NARRATIVE R1a): random sets of the same size drawn without replacement from a fixed
     gene pool (e.g. KEGG OXPHOS + TCA genes measured in both layers), one draw list per dataset.
Scores: signed cameraPR t (primary) and n opposed. Percentile = 100 x P(null <= observed) (as in R1a); empirical
p = (1 + #null >= observed) / (B + 1) (as in R1b).
"""
import numpy as np
import pandas as pd

from .core import camera_vec

CLASS_FLAGS = ["mitocarta3", "go_any_complex", "secreted_proxy"]


def _strata_abundance(S, cid, pos):
    b = S.base[pos]
    kind = "baseline (column)"
    if np.isfinite(b).mean() < 0.5:
        t = S.cols.loc[cid, "gtex_tissue"]
        if t:
            genes = S.genes[S.gcode[pos]]
            b = S.ann[f"gtex_tpm|{t}"].reindex(genes).to_numpy(float)
            b = np.log2(b + 1)
            kind = f"GTEx v8 TPM ({t})"
        else:
            return np.zeros(len(pos), int), "none (unmatched)"
    lab = np.full(len(pos), -1)
    ok = np.isfinite(b)
    if ok.sum() >= 10:
        lab[ok] = pd.qcut(b[ok], 10, labels=False, duplicates="drop")
    return lab, kind


def _strata_class(S, pos):
    genes = S.genes[S.gcode[pos]]
    a = S.ann.reindex(genes)[CLASS_FLAGS].fillna(False).to_numpy(bool)
    return a[:, 0] * 4 + a[:, 1] * 2 + a[:, 2] * 1


def matched_draws(strata, member_idx, rng, B):
    """For each member (in order) draw B random non-member positions from the member's stratum (with
    replacement), vectorised with the same RNG stream as the per-(boot, gene) loop of R1b."""
    nonmem = np.ones(len(strata), bool)
    nonmem[member_idx] = False
    buckets, sizes = [], []
    for i in member_idx:
        b = np.flatnonzero(nonmem & (strata == strata[i]))
        if len(b) == 0:  # empty stratum: fall back to all non-members
            b = np.flatnonzero(nonmem)
        buckets.append(b)
        sizes.append(len(b))
    r = rng.integers(0, np.array(sizes), size=(B, len(member_idx)))
    return np.column_stack([buckets[k][r[:, k]] for k in range(len(member_idx))]) if len(member_idx) else r


def column_null(S, cid, sig, kind, rng, B=1000):
    """Null distribution for one column. Returns dict with observed and null summaries."""
    pos = S.col(cid)
    g = S.gcode[pos]
    stat, lfc = S.stat[pos], S.lfc[pos]
    dmap = np.zeros(len(S.genes) + 1)
    codes = S.codes(sig.genes)
    ok = codes >= 0
    dmap[codes[ok]] = sig.dirs[ok]
    d = dmap[g]
    member_idx = np.flatnonzero(d != 0)
    m = len(member_idx)
    out = dict(column_id=cid, null=kind, n_members=m)
    if m < 2:
        return {**out, "stratum": "", "obs_t": np.nan, "null_t_median": np.nan, "null_t_q05": np.nan,
                "null_t_q95": np.nan, "pct_t": np.nan, "p_emp_t": np.nan, "obs_opposed": np.nan,
                "null_opposed_median": np.nan, "pct_opposed": np.nan, "p_emp_opposed": np.nan}
    if kind == "abundance":
        strata, skind = _strata_abundance(S, cid, pos)
    else:  # class x abundance tertile, so the class null is also abundance-matched
        ab, akind = _strata_abundance(S, cid, pos)
        tert = np.where(ab >= 0, np.minimum(ab // 4, 2), 3)  # deciles 0-3 / 4-7 / 8-9 -> 0 / 1 / 2; NA -> 3
        strata = _strata_class(S, pos) * 4 + tert
        skind = f"MitoCarta x GO:CC complex x secreted x abundance tertile ({akind})"
    draws = matched_draws(strata, member_idx, rng, B)
    signs = -d[member_idx]
    null_t = camera_vec(stat, draws, signs)
    obs_t = camera_vec(stat, member_idx[None, :], signs)[0]
    null_k = (np.sign(lfc[draws]) * d[member_idx] < 0).sum(1)
    obs_k = int((np.sign(lfc[member_idx]) * d[member_idx] < 0).sum())
    return {**out, "stratum": skind, "obs_t": obs_t, "null_t_median": np.median(null_t),
            "null_t_q05": np.percentile(null_t, 5), "null_t_q95": np.percentile(null_t, 95),
            "pct_t": 100 * np.mean(null_t <= obs_t), "p_emp_t": (1 + (null_t >= obs_t).sum()) / (B + 1),
            "obs_opposed": obs_k, "null_opposed_median": float(np.median(null_k)),
            "null_opposed_q95": float(np.percentile(null_k, 95)),
            "pct_opposed": 100 * np.mean(null_k <= obs_k), "p_emp_opposed": (1 + (null_k >= obs_k).sum()) / (B + 1),
            "_null_t": null_t}


def run_nulls(S, cids, sig, B=1000, seed=20260926, kinds=("abundance", "class")):
    """Both nulls for every column; each null type has its own RNG stream (seed, seed+1)."""
    out = []
    for j, kind in enumerate(kinds):
        rng = np.random.default_rng(seed + j)
        for cid in cids:
            out.append(column_null(S, cid, sig, kind, rng, B))
    return pd.DataFrame(out)


def pool_null(S, cids, genes, pool_by_dataset, rng, B=1000, draws=None):
    """NARRATIVE R1a: random len(genes)-gene sets from a fixed pool (per dataset), all genes 'down' (-1).
    `draws` may be passed in (dict dataset -> list of arrays) to share draws across calls."""
    k = len(genes)
    if draws is None:
        draws = {ds: [rng.choice(pool, k, replace=False) for _ in range(B)] for ds, pool in pool_by_dataset.items()}
    rows = []
    for cid in cids:
        ds = S.cols.loc[cid, "dataset"]
        pos = S.col(cid)
        g = S.gcode[pos]
        stat = S.stat[pos]
        idx_of = {c: i for i, c in enumerate(g)}

        def t_of(gs):
            ii = [idx_of[c] for c in S.codes(list(gs)) if c in idx_of]
            if len(ii) < 2:
                return np.nan
            return camera_vec(stat, np.array(ii)[None, :])[0]
        obs = t_of(genes)
        null = np.array([t_of(s) for s in draws[ds]])
        rows.append(dict(column_id=cid, obs_t=obs, null_t_median=np.nanmedian(null),
                         pct_t=100 * np.nanmean(null <= obs), pool_size=len(pool_by_dataset[ds])))
    return pd.DataFrame(rows), draws


def legacy_pools(S, sets, anchors=(("human_acute", "VL"), ("rat_train", "SKM-GN"))):
    """Pool = union of gene sets intersected with genes measured in BOTH RNA and PROT of the dataset's anchor
    tissue (any contrast), sorted, exactly as scripts/05_pah_grid.py R1a."""
    union = set().union(*[S.gs[s] for s in sets])
    r = S.rows
    out = {}
    for ds, t in anchors:
        def genes(layer):
            m = (r.dataset == ds) & (r.tissue == t) & (r.layer == layer)
            return set(r.gene_symbol_human[m].astype(str))
        out[ds] = sorted(union & genes("RNA") & genes("PROT"))
    return out
