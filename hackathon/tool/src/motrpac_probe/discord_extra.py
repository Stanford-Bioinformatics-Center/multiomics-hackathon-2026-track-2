"""Discordance panels e and f: protein-without-RNA calls from summary statistics, and whether early RNA
predicts later protein (held-out delta R2 with a permutation null).

Method details: docs/METHODS.md#discord_extrapy
"""
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from scipy.stats import norm  # noqa: E402

from .render import ACCENT, INK, LAYER_COLORS, MUTED, mpl_style  # noqa: E402

__all__ = ["protein_without_rna", "lag_panel"]

SEXES = ["F", "M"]
WEEKS = ["1w", "2w", "4w", "8w"]
CELLS8 = [(s, w) for s in SEXES for w in WEEKS]          # column order of every (genes x 8) matrix
LATE_IX = [2, 3, 6, 7]                                    # F_4w, F_8w, M_4w, M_8w inside CELLS8
Z95 = float(norm.ppf(0.95))                               # one-sided 5% / 90% two-sided bound
CV_REPEATS = 5                                            # lag: held-out dR2 also averaged over this many splits
SEX_WORDS = {"F": "female", "M": "male", "pooled": "both sexes"}

# ---- reference numbers, copied from protonly/RESULTS.md section 1 (sample-level analysis, lenient = lcb_all8)
PROTONLY_REF = {
    "SKM-GN": dict(n_repl=601, n_repl_measurable=561, n_supported=374, n_lenient=116, n_strict=62,
                   n_every=12, n_opposite=16, n_undetermined=55, null_shuffle_mean=0.1,
                   mirror_repl=242, mirror_lenient=24, mirror_strict=16, conc_expected_lenient=136,
                   conc_excess_strict=34, boot_median_bounded=0.23, boot_median_supported=0.69),
    "HEART": dict(n_repl=595, n_repl_measurable=571, n_supported=353, n_lenient=91, n_strict=42,
                  n_every=4, n_opposite=33, n_undetermined=94, null_shuffle_mean=0.06,
                  mirror_repl=312, mirror_lenient=68, mirror_strict=46, conc_expected_lenient=108,
                  conc_excess_strict=22, boot_median_bounded=0.31, boot_median_supported=0.75),
}
# ---- copied from lag/RESULTS.md (join_table_v2, raw logFC, pooled sexes). held-out = mean of 20 x 10-fold CV.
LAG_REF = {
    "SKM-GN": {("fixed", "1w"): 0.0039, ("fixed", "2w"): 0.0034, ("fixed", "4w"): 0.011, ("fixed", "8w"): 0.022,
               ("conditional", "1w"): -2.2e-4, ("conditional", "2w"): -1.9e-4, ("conditional", "4w"): 0.0021},
    "HEART": {("fixed", "1w"): 3.2e-4, ("fixed", "2w"): -1.7e-4, ("fixed", "4w"): 0.0021, ("fixed", "8w"): 0.0034,
              ("conditional", "1w"): -2.4e-4, ("conditional", "2w"): -2.5e-4, ("conditional", "4w"): 3.3e-4},
}
# in-sample grid (M0 = protein at w), pooled, raw
LAG_REF_GRID = {"SKM-GN": {"1w": 0.0065, "2w": 0.0037, "4w": 0.0055},
                "HEART": {"1w": 2.8e-4, "2w": 1.9e-4, "4w": 7.0e-4}}


# ============================================================================================ shared helpers
def _layer(S, tissue, layer):
    """(n_all_genes x 8) logFC, SE and baseline arrays over CELLS8, or None if a column is missing."""
    c = S.cols
    c = c[(c.dataset == "rat_train") & (c.tissue == tissue) & (c.layer == layer)]
    by = {str(k): cid for k, cid in zip(c.contrast, c.column_id)}
    G = len(S.genes)
    lfc, se, base = (np.full((G, 8), np.nan) for _ in range(3))
    for j, (s, w) in enumerate(CELLS8):
        cid = by.get(f"{s}_{w}")
        if cid is None:
            return None
        pos = S.col(cid)
        g = S.gcode[pos]
        ok = g >= 0
        l_, t_ = S.lfc[pos][ok], S.stat[pos][ok]
        with np.errstate(all="ignore"):
            e = np.abs(l_ / t_)
        e[~np.isfinite(e) | (e <= 0)] = np.nan
        lfc[g[ok], j] = l_
        se[g[ok], j] = e
        base[g[ok], j] = S.base[pos][ok]
    return dict(lfc=lfc, se=se, base=base, cols=[by[f"{s}_{w}"] for s, w in CELLS8])


def _bh(p):
    """Benjamini-Hochberg along the last axis (no NaN allowed)."""
    p = np.asarray(p, float)
    n = p.shape[-1]
    o = np.argsort(p, axis=-1)
    ps = np.take_along_axis(p, o, -1) * n / np.arange(1, n + 1)
    ps = np.minimum.accumulate(ps[..., ::-1], axis=-1)[..., ::-1]
    q = np.empty_like(ps)
    np.put_along_axis(q, o, np.minimum(ps, 1.0), -1)
    return q


def _late(L4, S4):
    """Per-sex late effects (mean of 4w and 8w) and their SE; L4/S4 (..., 4) in order F4, F8, M4, M8."""
    LF = (L4[..., 0] + L4[..., 1]) / 2
    LM = (L4[..., 2] + L4[..., 3]) / 2
    SF = np.sqrt(S4[..., 0] ** 2 + S4[..., 1] ** 2) / 2
    SM = np.sqrt(S4[..., 2] ** 2 + S4[..., 3] ** 2) / 2
    return LF, SF, LM, SM


def _responders(L4, S4, cutoff):
    """Replicated responder rule on the four late cells (see protein_without_rna docstring)."""
    sg = np.sign(L4)
    same = np.all(sg == sg[..., :1], -1) & (sg[..., 0] != 0)
    zs = (L4 / S4).sum(-1) / 2.0                           # Stouffer of the four signed z
    q = _bh(2 * norm.sf(np.abs(zs)))
    LF, SF, LM, SM = _late(L4, S4)
    L = (LF + LM) / 2
    SE = np.sqrt(SF ** 2 + SM ** 2) / 2
    sigma = np.where(same, sg[..., 0], np.sign(L))
    repl = same & (q < cutoff) & (sigma * LF / SF > Z95) & (sigma * LM / SM > Z95)
    lcb = np.maximum(np.abs(L) - Z95 * SE, 0.0)
    return dict(sigma=sigma, repl=repl, L=L, SE=SE, lcb=lcb, q=q)


def _evidence(sigma, lcb, L4, S4):
    """Classify the other layer's pooled late effect in direction sigma against margin lcb (lenient) / lcb/2."""
    LF, SF, LM, SM = _late(L4, S4)
    x = sigma * (LF + LM) / 2
    s = np.sqrt(SF ** 2 + SM ** 2) / 2
    sup = x / s > Z95
    opp = (-x / s > Z95) & ~sup
    ub = x + Z95 * s
    free = ~sup & ~opp
    return dict(sup=sup, opp=opp, lenient=free & (ub < lcb), strict=free & (ub < 0.5 * lcb), x=x, s=s, ub=ub)


def _p_bounded(theta, s, m):
    """P(not significant either way and upper bound < m) when x ~ N(theta, s)."""
    upper = np.minimum(m - Z95 * s, Z95 * s)
    lower = -Z95 * s
    p = norm.cdf((upper - theta) / s) - norm.cdf((lower - theta) / s)
    return np.where(upper > lower, np.clip(p, 0, 1), 0.0)


def _empty(reason):
    return dict(tables={}, figures={}, numbers={}, captions={}, notes=[f"not available: {reason}"])


def _fmt(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "NA"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    ax = abs(x)
    if ax == 0:
        return "0"
    if ax < 1e-3:
        return f"{x:.1e}"
    if ax < 0.1:
        return f"{x:.4f}"
    return f"{x:.2f}"


# ============================================================================================ panel 1
def protein_without_rna(S, tissue, cutoff=0.05, nboot=500, nnull=200, seed=1, min_count=10):
    """Protein training responses with no detectable RNA response, from summary statistics (rat, one tissue).

    Rules (all one-sided tests at 5%, bounds = one-sided 95% = two-sided 90%):
      late effect per sex  = mean(4w, 8w logFC); SE = sqrt(se4^2 + se8^2) / 2 (4w and 8w are separate animal
                             groups; the shared sedentary control covariance is ignored)
      pooled late effect L = mean over sexes, SE = sqrt(SE_F^2 + SE_M^2) / 2
      replicated responder = four sex x {4w, 8w} logFC share a sign, BH q < cutoff on the Stouffer z of the four
                             signed z (logFC/SE, two-sided p, BH over the tissue's universe), and each sex's late
                             effect one-sided p < 0.05 in that sign
      RNA measurable       = median over the 8 RNA columns of (2^baseline_expr - 1) >= min_count
      classes (RNA pooled late effect x in the protein direction, SE s): supported x/s > z95; opposite -x/s > z95;
                             protein-only lenient: neither and x + z95 s < lcb (lcb = max(|L_P| - z95 SE_P, 0));
                             strict: x + z95 s < lcb / 2; undetermined otherwise
      bounded at every sex x week: lenient, and each of the 8 RNA cells' upper bound < lcb
    Returns dict(tables, figures, numbers, captions, notes); tables: classes, nulls, mirror, stability,
    illustrations; figures: calls_vs_null, trajectories.
    """
    t0 = time.time()
    P, R = _layer(S, tissue, "PROT"), _layer(S, tissue, "RNA")
    if P is None or R is None:
        return _empty(f"{tissue} needs all 8 rat_train RNA and PROT columns (F/M x 1w/2w/4w/8w)")
    rng = np.random.default_rng(seed)

    # universe: genes with finite late logFC + SE in both layers
    uni = np.all(np.isfinite(P["lfc"][:, LATE_IX]) & np.isfinite(P["se"][:, LATE_IX]) &
                 np.isfinite(R["lfc"][:, LATE_IX]) & np.isfinite(R["se"][:, LATE_IX]), 1)
    codes = np.flatnonzero(uni)
    if len(codes) < 50:
        return _empty(f"{tissue}: only {len(codes)} genes measured in both layers")
    PL, PS = P["lfc"][codes][:, LATE_IX], P["se"][codes][:, LATE_IX]
    RL, RS = R["lfc"][codes][:, LATE_IX], R["se"][codes][:, LATE_IX]
    RL8, RS8 = R["lfc"][codes], R["se"][codes]
    PL8, PS8 = P["lfc"][codes], P["se"][codes]
    with np.errstate(all="ignore"):
        count = np.nanmedian(2.0 ** R["base"][codes] - 1, axis=1)
    measurable = count >= min_count
    genes = np.asarray(S.genes)[codes]
    G = len(codes)

    # ---------------------------------------------------------------- observed calls
    pr = _responders(PL, PS, cutoff)
    ev = _evidence(pr["sigma"], pr["lcb"], RL, RS)
    resp = pr["repl"]
    den = resp & measurable
    with np.errstate(all="ignore"):
        cell_ub = pr["sigma"][:, None] * RL8 + Z95 * RS8
    every = den & ev["lenient"] & np.all(np.nan_to_num(cell_ub, nan=np.inf) < pr["lcb"][:, None], 1)
    n = dict(universe=G, repl=int(resp.sum()), repl_measurable=int(den.sum()),
             too_low=int((resp & ~measurable).sum()),
             supported=int((den & ev["sup"]).sum()), opposite=int((den & ev["opp"]).sum()),
             lenient=int((den & ev["lenient"]).sum()), strict=int((den & ev["strict"]).sum()),
             every=int(every.sum()))
    n["undetermined"] = n["repl_measurable"] - n["supported"] - n["opposite"] - n["lenient"]
    frac = (lambda k: k / n["repl_measurable"] if n["repl_measurable"] else np.nan)

    # ---------------------------------------------------------------- mirror: RNA responders, protein bounded
    rr = _responders(RL, RS, cutoff)
    rden = rr["repl"] & measurable
    evp = _evidence(rr["sigma"], rr["lcb"], PL, PS)
    m = dict(repl=int(rden.sum()), supported=int((rden & evp["sup"]).sum()),
             opposite=int((rden & evp["opp"]).sum()), lenient=int((rden & evp["lenient"]).sum()),
             strict=int((rden & evp["strict"]).sum()))
    m["undetermined"] = m["repl"] - m["supported"] - m["opposite"] - m["lenient"]

    # ---------------------------------------------------------------- null (a): label-shuffle analogue
    null_repl, null_len, null_str = [], [], []
    chunk = 25
    for b0 in range(0, nnull, chunk):
        B = min(chunk, nnull - b0)
        sim = rng.standard_normal((B, G, 4)) * PS[None]
        r0 = _responders(sim, PS[None], cutoff)
        e0 = _evidence(r0["sigma"], r0["lcb"], RL[None], RS[None])
        d0 = r0["repl"] & measurable[None]
        null_repl += list(d0.sum(1))
        null_len += list((d0 & e0["lenient"]).sum(1))
        null_str += list((d0 & e0["strict"]).sum(1))
    null_repl, null_len, null_str = map(np.asarray, (null_repl, null_len, null_str))

    # ---------------------------------------------------------------- null (b): concordance null
    xP, xR = pr["L"][den], (RL[den][:, [0, 1]].mean(1) + RL[den][:, [2, 3]].mean(1)) / 2
    sP = pr["SE"][den]
    conc_rows, beta = [], {}
    if den.sum() >= 5:
        vP = np.var(xP, ddof=1)
        rel = (vP - np.mean(sP ** 2)) / vP
        ols = np.cov(xP, xR, ddof=1)[0, 1] / vP
        beta["all"] = ols / rel
        sup_d = ev["sup"][den]
        if sup_d.sum() >= 5:
            vPs = np.var(xP[sup_d], ddof=1)
            rels = (vPs - np.mean(sP[sup_d] ** 2)) / vPs
            beta["supported"] = (np.cov(xP[sup_d], xR[sup_d], ddof=1)[0, 1] / vPs) / rels
        else:
            beta["supported"] = np.nan
        th_s = ev["s"][den]
        for bname, b_ in beta.items():
            theta = b_ * np.abs(xP)                        # true RNA effect in the protein direction
            for cname, mg, obs in [("lenient", pr["lcb"][den], n["lenient"]),
                                   ("strict", 0.5 * pr["lcb"][den], n["strict"])]:
                pb = np.nan_to_num(_p_bounded(theta, th_s, mg)) if np.isfinite(b_) else np.full(len(theta), np.nan)
                ex = float(pb.sum())
                sd = float(np.sqrt(np.sum(pb * (1 - pb))))
                conc_rows.append(dict(null=f"concordance, beta {bname}", call=cname, beta=b_, observed=obs,
                                      expected=ex, null_sd=sd, excess=obs - ex,
                                      z=(obs - ex) / sd if sd > 0 else np.nan))
        beta["ols"], beta["reliability"] = ols, rel
    conc = pd.DataFrame(conc_rows)

    # ---------------------------------------------------------------- null (c): parametric bootstrap
    f_repl, f_sup, f_len, f_str = (np.zeros(G) for _ in range(4))
    for b0 in range(0, nboot, chunk):
        B = min(chunk, nboot - b0)
        bP = PL[None] + rng.standard_normal((B, G, 4)) * PS[None]
        bR = RL[None] + rng.standard_normal((B, G, 4)) * RS[None]
        rb = _responders(bP, PS[None], cutoff)
        eb = _evidence(rb["sigma"], rb["lcb"], bR, RS[None])
        # a call counts only if the resampled responder keeps the observed direction
        db = rb["repl"] & measurable[None] & (rb["sigma"] == pr["sigma"][None])
        f_repl += db.sum(0)
        f_sup += (db & eb["sup"]).sum(0)
        f_len += (db & eb["lenient"]).sum(0)
        f_str += (db & eb["strict"]).sum(0)
    f_repl, f_sup, f_len, f_str = (f / max(nboot, 1) for f in (f_repl, f_sup, f_len, f_str))

    obs_sup, obs_len, obs_str = den & ev["sup"], den & ev["lenient"], den & ev["strict"]
    med = (lambda f, msk: float(np.median(f[msk])) if msk.any() else np.nan)
    stab = pd.DataFrame([
        dict(observed_class="RNA-supported", n_genes=int(obs_sup.sum()), call_frequency_of="RNA-supported",
             median_frequency=med(f_sup, obs_sup), n_frequency_ge_0_5=int((f_sup[obs_sup] >= 0.5).sum())),
        dict(observed_class="protein-only lenient", n_genes=int(obs_len.sum()),
             call_frequency_of="protein-only lenient", median_frequency=med(f_len, obs_len),
             n_frequency_ge_0_5=int((f_len[obs_len] >= 0.5).sum())),
        dict(observed_class="protein-only strict", n_genes=int(obs_str.sum()),
             call_frequency_of="protein-only strict", median_frequency=med(f_str, obs_str),
             n_frequency_ge_0_5=int((f_str[obs_str] >= 0.5).sum())),
        dict(observed_class="replicated responder (RNA measurable)", n_genes=int(den.sum()),
             call_frequency_of="replicated responder", median_frequency=med(f_repl, den),
             n_frequency_ge_0_5=int((f_repl[den] >= 0.5).sum())),
    ])

    # ---------------------------------------------------------------- illustrations
    idx = np.flatnonzero(obs_str)
    idx = idx[np.lexsort((-np.abs(pr["L"][idx]), -f_str[idx]))][:12]
    ill = pd.DataFrame(dict(gene=genes[idx], protein_direction=np.where(pr["sigma"][idx] > 0, "up", "down"),
                            boot_freq_strict=f_str[idx], boot_freq_lenient=f_len[idx],
                            protein_late=pr["L"][idx], protein_late_lcb=pr["lcb"][idx],
                            rna_late=(RL[idx][:, [0, 1]].mean(1) + RL[idx][:, [2, 3]].mean(1)) / 2,
                            rna_late_ub_in_protein_direction=ev["ub"][idx]))
    for j, (s, w) in enumerate(CELLS8):
        ill[f"rna_{s}_{w}"] = RL8[idx, j]
    for j, (s, w) in enumerate(CELLS8):
        ill[f"prot_{s}_{w}"] = PL8[idx, j]

    # ---------------------------------------------------------------- tables
    ref = PROTONLY_REF.get(tissue)
    cls_rows = [("replicated protein responders", n["repl"], None, "n_repl"),
                ("  RNA measurable (denominator)", n["repl_measurable"], None, "n_repl_measurable"),
                ("  RNA too low to call", n["too_low"], None, None),
                ("RNA-supported", n["supported"], frac(n["supported"]), "n_supported"),
                ("RNA significantly opposite", n["opposite"], frac(n["opposite"]), "n_opposite"),
                ("protein-only, lenient", n["lenient"], frac(n["lenient"]), "n_lenient"),
                ("protein-only, strict (subset of lenient)", n["strict"], frac(n["strict"]), "n_strict"),
                ("undetermined", n["undetermined"], frac(n["undetermined"]), "n_undetermined"),
                ("bounded at every sex x week", n["every"], frac(n["every"]), "n_every")]
    classes = pd.DataFrame([dict(call=a, n_genes=b, fraction_of_denominator=c,
                                 original_sample_level=(ref.get(k) if (ref and k) else np.nan))
                            for a, b, c, k in cls_rows])
    nulls = pd.DataFrame(
        [dict(null="label-shuffle analogue (protein logFC ~ N(0, SE))", call=c, beta=np.nan, observed=o,
              expected=float(v.mean()), null_sd=float(v.std(ddof=1)) if len(v) > 1 else np.nan,
              excess=o - float(v.mean()), z=np.nan)
         for c, o, v in [("replicated responders", n["repl_measurable"], null_repl),
                         ("lenient", n["lenient"], null_len), ("strict", n["strict"], null_str)]]
        + ([] if conc.empty else conc.to_dict("records")))
    mirror = pd.DataFrame([dict(call=a, n_genes=b, fraction=(b / m["repl"] if m["repl"] else np.nan))
                           for a, b in [("replicated RNA responders (RNA measurable)", m["repl"]),
                                        ("protein-supported", m["supported"]),
                                        ("protein significantly opposite", m["opposite"]),
                                        ("RNA-only (protein bounded), lenient", m["lenient"]),
                                        ("RNA-only (protein bounded), strict", m["strict"]),
                                        ("undetermined", m["undetermined"])]])

    # ---------------------------------------------------------------- numbers
    ce = {(r["null"], r["call"]): r for r in conc_rows}
    ba, bs = "concordance, beta all", "concordance, beta supported"
    numbers = dict(
        n_universe=G, n_repl=n["repl"], n_repl_measurable=n["repl_measurable"], n_too_low=n["too_low"],
        n_supported=n["supported"], n_opposite=n["opposite"], n_lenient=n["lenient"], n_strict=n["strict"],
        n_undetermined=n["undetermined"], n_every=n["every"],
        frac_lenient=frac(n["lenient"]), frac_strict=frac(n["strict"]),
        null_shuffle_mean_repl=float(null_repl.mean()), null_shuffle_mean_lenient=float(null_len.mean()),
        null_shuffle_mean_strict=float(null_str.mean()),
        beta_eiv_all=beta.get("all", np.nan), beta_eiv_supported=beta.get("supported", np.nan),
        beta_ols_all=beta.get("ols", np.nan), protein_reliability=beta.get("reliability", np.nan),
        conc_expected_lenient=ce.get((ba, "lenient"), {}).get("expected", np.nan),
        conc_expected_strict=ce.get((ba, "strict"), {}).get("expected", np.nan),
        conc_excess_lenient=ce.get((ba, "lenient"), {}).get("excess", np.nan),
        conc_excess_strict=ce.get((ba, "strict"), {}).get("excess", np.nan),
        conc_sup_excess_lenient=ce.get((bs, "lenient"), {}).get("excess", np.nan),
        conc_sup_excess_strict=ce.get((bs, "strict"), {}).get("excess", np.nan),
        mirror_repl=m["repl"], mirror_lenient=m["lenient"], mirror_strict=m["strict"],
        boot_median_lenient=med(f_len, obs_len), boot_median_strict=med(f_str, obs_str),
        boot_median_supported=med(f_sup, obs_sup),
        n_lenient_boot_ge_0_5=int((f_len[obs_len] >= 0.5).sum()), nboot=nboot, nnull=nnull)

    # ---------------------------------------------------------------- figures
    figures = {"calls_vs_null": _fig_calls(numbers, tissue)}
    if len(ill):
        figures["trajectories"] = _fig_traj(ill, tissue)
    numbers["runtime_s"] = round(time.time() - t0, 2)

    # ---------------------------------------------------------------- notes
    notes = [
        "Summary-statistics re-implementation of protonly/PLAN.md: every effect is a trained-vs-sedentary group "
        "contrast; SE = |logFC / stat| (DESeq2 Wald z for RNA, limma t for protein); normal quantiles replace "
        "PLAN's Welch/Satterthwaite t intervals.",
        "Late effect per sex = mean of the 4w and 8w logFC with SE = sqrt(se4^2 + se8^2)/2. The 4w and 8w groups "
        "are different animals but share the sedentary control, whose covariance is ignored (PLAN computed it "
        "exactly from sample-level data), so late-effect SEs are approximate (likely too small).",
        "Responder test: Stouffer z of the four signed sex x week z (sum/2), two-sided, BH over the "
        f"{G} genes with late logFC in both layers; PLAN used a t test of the pooled late contrast.",
        f"RNA measurable = sedentary mean normalized count >= {min_count} (count = 2^baseline_expr - 1, median "
        "over the 8 RNA columns). PLAN used mean log2 CPM >= 2, which is not in the summary tables.",
        "RNA evidence uses the RNA late effect (4w/8w, both sexes) only; PLAN's primary used all 8 RNA sex x week "
        "cells and also counted any single significant RNA cell as support (closest to PLAN's 'lcb_late' "
        "sensitivity variant, which pooled RNA over the same 4w/8w window).",
        "Genes come from the store's per-gene rows (join_table_v2 collapsing to human symbols); PLAN chose "
        "protein features effect-blind (fewest missing values) and 1:1 orthologs, so gene sets differ.",
        "Label-shuffle analogue: each gene's four late protein logFC redrawn from N(0, SE) "
        f"({nnull} times) with the observed RNA; the full rule is re-applied. It calibrates the responder step "
        "only.",
        "Concordance null: RNA late effect ~ N(beta * |protein late effect|, SE_RNA) in the protein direction; "
        "beta = errors-in-variables slope over the replicated responders (OLS slope of RNA on protein late effect "
        "divided by the protein reliability (var(L_P) - mean(SE_P^2)) / var(L_P)); expected counts are the exact "
        "expectation of that simulation (sum of per-gene probabilities), not a Monte Carlo mean. 'beta "
        "supported' uses the RNA-supported responders only. PLAN used a through-origin slope and the all-8 RNA "
        "window.",
        f"Stability: parametric bootstrap, all late logFC redrawn from N(observed, SE) {nboot} times and the calls "
        "recomputed (PLAN resampled animals). A call counts only if the gene stays a responder in the same "
        "direction.",
    ]
    if ref:
        notes.append(
            f"Comparison with protonly/RESULTS.md ({tissue}, sample-level): replicated responders "
            f"{n['repl']} here vs {ref['n_repl']} there ({n['repl'] - ref['n_repl']:+d}); RNA-measurable "
            f"{n['repl_measurable']} vs {ref['n_repl_measurable']}; supported {n['supported']} vs "
            f"{ref['n_supported']}; lenient/strict {n['lenient']}/{n['strict']} vs "
            f"{ref['n_lenient']}/{ref['n_strict']} ({n['lenient'] - ref['n_lenient']:+d}/"
            f"{n['strict'] - ref['n_strict']:+d}); every sex x week {n['every']} vs {ref['n_every']}; opposite "
            f"{n['opposite']} vs {ref['n_opposite']}; mirror lenient/strict {m['lenient']}/{m['strict']} of "
            f"{m['repl']} vs {ref['mirror_lenient']}/{ref['mirror_strict']} of {ref['mirror_repl']}; median "
            f"bootstrap frequency lenient {_fmt(numbers['boot_median_lenient'])} vs {ref['boot_median_bounded']}, "
            f"supported {_fmt(numbers['boot_median_supported'])} vs {ref['boot_median_supported']}. The "
            "rules above differ from the sample-level ones, so the counts are expected to differ; nothing was "
            "tuned to match.")
    notes.append("Report counts, not gene lists: individual protein-only calls sit near the margin (see the "
                 "stability table); the illustration genes are the most stable strict calls, shown as examples.")

    captions = {
        "classes": (f"{tissue}: replicated protein responders by RNA evidence",
                    "Counts of genes; fractions are of responders whose RNA is measurable. Lenient = the RNA upper "
                    "90% bound lies below the protein lower 90% bound; strict = below half of it. "
                    "'original_sample_level' = protonly/RESULTS.md (different, sample-level method)."),
        "nulls": (f"{tissue}: how many calls would chance or a uniformly scaled-down RNA response produce?",
                  "Label-shuffle rows: mean calls when protein has no training effect (should be ~0). "
                  "Concordance rows: expected protein-only calls if every responder's RNA moved beta times its "
                  "protein effect; 'excess' = observed - expected, the part a uniform attenuation cannot explain."),
        "mirror": (f"{tissue}: the mirror - RNA responders whose protein is bounded",
                   "Same rules with the layers swapped. If RNA-only calls are as common as protein-only calls, "
                   "the asymmetry is measurement, not biology."),
        "stability": (f"{tissue}: how often each call recurs when the logFCs are redrawn within their SE",
                      "Median over the genes observed in a class of the fraction of bootstrap draws giving the "
                      "same call. Low values mean the count is meaningful but the individual names are not."),
        "illustrations": (f"{tissue}: the most stable strict protein-only calls (illustrations)",
                          "Top genes by bootstrap frequency of the strict call, with their per sex x week RNA and "
                          "protein logFC. Examples, not an enrichment test."),
        "calls_vs_null": (f"{tissue}: protein-only calls vs what the nulls expect",
                          "Blue bars: observed protein-only calls. Dark grey: expected under the concordance null "
                          "(RNA = beta x protein, observed with its own noise). Light: label-shuffle null. "
                          "The gap between dark and grey is the excess."),
        "trajectories": (f"{tissue}: RNA and protein trajectories of the most stable strict protein-only genes",
                         "Each panel is one gene: log2 fold change vs sedentary at 1, 2, 4 and 8 weeks. Blue = RNA, "
                         "orange = protein; solid = female, dashed = male. Protein moves while RNA stays near 0."),
    }
    return dict(tables=dict(classes=classes, nulls=nulls, mirror=mirror, stability=stab, illustrations=ill),
                figures=figures, numbers=numbers, captions=captions, notes=notes)


def _fig_calls(nb, tissue):
    with mpl_style():
        fig, ax = plt.subplots(figsize=(7.5, 3.2))
        groups = ["lenient", "strict"]
        obs = [nb["n_lenient"], nb["n_strict"]]
        conc = [nb["conc_expected_lenient"], nb["conc_expected_strict"]]
        shuf = [nb["null_shuffle_mean_lenient"], nb["null_shuffle_mean_strict"]]
        y = np.arange(len(groups))
        h = 0.26
        for off, vals, col, lab in [(-h, obs, ACCENT, "observed"),
                                    (0, conc, MUTED, "concordance null (expected)"),
                                    (h, shuf, "#bdbdbd", "label-shuffle null (mean)")]:
            vv = np.nan_to_num(np.asarray(vals, float))
            ax.barh(y + off, vv, height=h * 0.9, color=col, label=lab)
            for yi, v in zip(y + off, vals):
                v = float(v)
                lab_v = "NA" if not np.isfinite(v) else (f"{v:.0f}" if abs(v) >= 10 else f"{v:.2f}")
                ax.text((0 if not np.isfinite(v) else v) + max(max(obs), 1) * 0.01, yi, lab_v,
                        va="center", fontsize=10, color=INK)
        ax.set_yticks(y, [f"protein-only\n{g}" for g in groups])
        ax.invert_yaxis()
        ax.set_xlabel(f"genes (of {nb['n_repl_measurable']} replicated responders with measurable RNA)")
        ax.legend(loc="lower right", fontsize=10)
        ax.set_xlim(0, max([max(obs)] + [v for v in conc if np.isfinite(v)] + [1]) * 1.18)
        fig.tight_layout()
    return fig


def _fig_traj(ill, tissue):
    k = len(ill)
    ncol = 4 if k > 6 else max(k, 1)
    nrow = int(np.ceil(k / ncol))
    x = np.arange(4)
    with mpl_style():
        fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.5 * nrow + 0.6), sharex=True, sharey=True,
                                 squeeze=False)
        for i, ax in enumerate(axes.flat):
            if i >= k:
                ax.set_visible(False)
                continue
            r = ill.iloc[i]
            ax.axhline(0, color=MUTED, lw=0.8)
            for lay, pre in [("RNA", "rna"), ("PROT", "prot")]:
                for s, ls in [("F", "-"), ("M", "--")]:
                    ax.plot(x, [r[f"{pre}_{s}_{w}"] for w in WEEKS], ls=ls, marker="o", ms=4,
                            color=LAYER_COLORS[lay], lw=1.6)
            ax.set_title(f"{r.gene} (protein {r.protein_direction})", fontsize=11)
            ax.set_xticks(x, ["1", "2", "4", "8"])
            if i % ncol == 0:
                ax.set_ylabel("log2 FC", fontsize=11)
            if i >= k - ncol:
                ax.set_xlabel("weeks of training", fontsize=11)
            ax.tick_params(labelsize=10)
        handles = [Line2D([], [], color=LAYER_COLORS["RNA"], lw=2, label="RNA"),
                   Line2D([], [], color=LAYER_COLORS["PROT"], lw=2, label="protein"),
                   Line2D([], [], color=INK, ls="-", label="female"),
                   Line2D([], [], color=INK, ls="--", label="male")]
        fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=11, frameon=False,
                   bbox_to_anchor=(0.5, 1.0))
        fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


# ============================================================================================ panel 2
def _lag_data(S, tissue):
    P, R = _layer(S, tissue, "PROT"), _layer(S, tissue, "RNA")
    if P is None or R is None:
        return None
    out = {}
    for si, s in enumerate(SEXES):
        j = slice(4 * si, 4 * si + 4)
        ok = np.all(np.isfinite(P["lfc"][:, j]) & np.isfinite(R["lfc"][:, j]), 1)
        out[s] = dict(codes=np.flatnonzero(ok), prot=P["lfc"][ok, j], rna=R["lfc"][ok, j])
    F, M = out["F"], out["M"]
    out["pooled"] = dict(codes=np.r_[F["codes"], M["codes"]], prot=np.vstack([F["prot"], M["prot"]]),
                         rna=np.vstack([F["rna"], M["rna"]]),
                         sexM=np.r_[np.zeros(len(F["codes"])), np.ones(len(M["codes"]))])
    return out


def _perm_rows(d, sex, rng):
    """Row permutation that shuffles genes. Pooled (as 07_lag.perm_index): one gene permutation applied to both
    sexes for genes present in both, so a gene's F and M rows stay paired; single-sex genes shuffle within sex."""
    n = len(d["codes"])
    if sex != "pooled":
        return rng.permutation(n)
    out = np.arange(n)
    nF = int((d["sexM"] == 0).sum())
    cF, cM = d["codes"][:nF], d["codes"][nF:]
    both = np.intersect1d(cF, cM)
    iF = np.searchsorted(cF, both)                           # codes are sorted within each sex
    iM = nF + np.searchsorted(cM, both)
    pg = rng.permutation(len(both))
    out[iF], out[iM] = iF[pg], iM[pg]
    for ix in (np.setdiff1d(np.arange(nF), iF), np.setdiff1d(np.arange(nF, n), iM)):
        out[ix] = ix[rng.permutation(len(ix))]
    return out


def _col(d, name):
    lay, w = name.split("_")
    return d[lay][:, WEEKS.index(w)]


def _X(d, cols, sex):
    parts = [np.ones(len(d["prot"]))] + [_col(d, c) for c in cols]
    if sex == "pooled":
        parts.append(d["sexM"])
    return np.column_stack(parts)


def _oof_sse(X, y, fold, k):
    sse = 0.0
    for f in range(k):
        te = fold == f
        b, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        sse += float(((y[te] - X[te] @ b) ** 2).sum())
    return sse


class _FoldCV:
    """Out-of-fold SSE of y ~ X0 + x for many x with fixed X0, y and folds (train normal equations per fold,
    obtained as full-data sums minus test-fold sums). Identical to per-fold least squares for full-rank X."""

    def __init__(self, X0, y, fold, k):
        self.order = np.argsort(fold, kind="stable")
        X0, y = X0[self.order], y[self.order]
        self.X0, self.y, self.k = X0, y, k
        self.bounds = np.r_[0, np.cumsum(np.bincount(fold, minlength=k))]
        self.sl = [slice(self.bounds[f], self.bounds[f + 1]) for f in range(k)]
        A, c = X0.T @ X0, X0.T @ y
        self.A0 = [A - X0[s].T @ X0[s] for s in self.sl]
        self.c0 = [c - X0[s].T @ y[s] for s in self.sl]

    def sse(self, x):
        # einsum instead of BLAS calls: these products are tiny, and multithreaded BLAS on a busy node is ~10x
        # slower than a plain loop for them
        x = x[self.order]
        X0, y = self.X0, self.y
        p = X0.shape[1]
        Xx_all, xx_all, xy_all = np.einsum("np,n->p", X0, x), np.einsum("n,n->", x, x), np.einsum("n,n->", x, y)
        A = np.empty((p + 1, p + 1))
        out = 0.0
        for f, s in enumerate(self.sl):
            X0t, xt, yt = X0[s], x[s], y[s]
            Xx = Xx_all - np.einsum("np,n->p", X0t, xt)
            A[:p, :p], A[:p, p], A[p, :p], A[p, p] = self.A0[f], Xx, Xx, xx_all - np.einsum("n,n->", xt, xt)
            b = np.linalg.solve(A, np.r_[self.c0[f], xy_all - np.einsum("n,n->", xt, yt)])
            r = yt - np.einsum("np,p->n", X0t, b[:p]) - xt * b[p]
            out += float(np.einsum("n,n->", r, r))
        return out


class _Resid:
    """Residualise vectors on a fixed design (orthonormal basis computed once; einsum, no BLAS)."""

    def __init__(self, X):
        self.Q = np.linalg.qr(X)[0]

    def __call__(self, v):
        return v - np.einsum("np,p->n", self.Q, np.einsum("np,n->p", self.Q, v))


def lag_panel(S, tissue, nperm=200, seed=1, folds=5):
    """Does RNA at week w add held-out information about the 8-week protein response?

    Target y = 8-week protein logFC (07_lag.py modelled logFC; stat was used there only for calibration).
    Rows = genes with RNA and protein logFC at all four weeks in that sex; 'pooled' stacks both sexes with a sex
    covariate and keeps a gene's two rows in the same fold / permutation.
    Model families (M1 = M0 + RNA(w)):
      grid        : M0 = protein(w)                  w = 1w, 2w, 4w
      fixed       : M0 = protein(2w)                 w = 1w, 2w, 4w, 8w (8w = same-time reference)
      conditional : M0 = protein(2w) + RNA(8w)       w = 1w, 2w, 4w   (early RNA beyond 8-week RNA)
    Held-out dR2 = out-of-fold R2(M1) - R2(M0), each relative to the intercept(+sex) model on the same folds
    (k folds over genes). Permutation null: RNA(w) permuted across genes nperm times (same folds).
    Returns dict(tables, figures, numbers, captions, notes); tables: models, n_genes; figures: delta_r2.
    """
    t0 = time.time()
    D = _lag_data(S, tissue)
    if D is None:
        return _empty(f"{tissue} needs all 8 rat_train RNA and PROT columns (F/M x 1w/2w/4w/8w)")
    if min(len(D[s]["codes"]) for s in SEXES) < 10 * folds:
        return _empty(f"{tissue}: too few genes measured in both layers at every week")
    rng = np.random.default_rng(seed)
    specs = [("grid", w, [f"prot_{w}"]) for w in ["1w", "2w", "4w"]] + \
            [("fixed", w, ["prot_2w"]) for w in WEEKS] + \
            [("conditional", w, ["prot_2w", "rna_8w"]) for w in ["1w", "2w", "4w"]]
    rows, cache = [], {}
    for sex in SEXES + ["pooled"]:
        d = D[sex]
        ug, inv = np.unique(d["codes"], return_inverse=True)
        ng = len(ug)
        fold = (rng.permutation(ng) % folds)[inv]              # folds over genes: F and M rows stay together
        extra_folds = [(rng.permutation(ng) % folds)[inv] for _ in range(CV_REPEATS - 1)]
        sse_null_x = [_oof_sse(_X(d, [], sex), _col(d, "prot_8w"), f_, folds) for f_ in extra_folds]
        y = _col(d, "prot_8w")
        Xn = _X(d, [], sex)
        sse_null = _oof_sse(Xn, y, fold, folds)
        sst = float(((y - y.mean()) ** 2).sum())
        perms = [_perm_rows(d, sex, rng) for _ in range(nperm)]
        for fam, w, base in specs:
            key = (sex, tuple(base), w)
            if key in cache:                               # grid 2w == fixed 2w
                rows.append({**cache[key], "family": fam})
                continue
            X0 = _X(d, base, sex)
            x = _col(d, f"rna_{w}")
            X1 = np.column_stack([X0, x])
            b1, *_ = np.linalg.lstsq(X1, y, rcond=None)
            r0 = 1 - float(((y - X0 @ np.linalg.lstsq(X0, y, rcond=None)[0]) ** 2).sum()) / sst
            r1 = 1 - float(((y - X1 @ b1) ** 2).sum()) / sst
            sse0 = _oof_sse(X0, y, fold, folds)
            cv = _FoldCV(X0, y, fold, folds)
            oof = (sse0 - cv.sse(x)) / sse_null
            oof_rep = [oof] + [(_oof_sse(X0, y, f_, folds) - _oof_sse(X1, y, f_, folds)) / sn_
                               for f_, sn_ in zip(extra_folds, sse_null_x)]
            res0 = _Resid(X0)
            ry = res0(y)

            def dr2_in(xv):
                rx = res0(xv)
                return float(np.einsum("n,n->", ry, rx) ** 2 / np.einsum("n,n->", rx, rx) / sst)

            null_in, null_oof = np.empty(nperm), np.empty(nperm)
            for i, ix in enumerate(perms):
                xp = x[ix]
                null_in[i] = dr2_in(xp)
                null_oof[i] = (sse0 - cv.sse(xp)) / sse_null
            r = dict(family=fam, rna_week=w, sex=sex, model=f"prot_8w ~ {' + '.join(base)} (+ rna_{w})",
                     n_genes=ng, n_rows=len(y), R2_M0=r0, R2_M1=r1, dR2_in_sample=r1 - r0, beta_rna=float(b1[-1]),
                     dR2_heldout=oof, dR2_heldout_mean_repeats=float(np.mean(oof_rep)),
                     null_heldout_q05=float(np.quantile(null_oof, 0.05)),
                     null_heldout_q50=float(np.quantile(null_oof, 0.5)),
                     null_heldout_q95=float(np.quantile(null_oof, 0.95)),
                     perm_p_heldout=(1 + int((null_oof >= oof).sum())) / (1 + nperm),
                     percentile_in_null=float((null_oof < oof).mean() * 100),
                     null_in_sample_q95=float(np.quantile(null_in, 0.95)),
                     perm_p_in_sample=(1 + int((null_in >= r1 - r0).sum())) / (1 + nperm))
            cache[key] = r
            rows.append(r)
    fam_order = {"grid": 0, "fixed": 1, "conditional": 2}
    tab = pd.DataFrame(rows)
    tab["_o"] = tab.family.map(fam_order)
    tab["_s"] = tab.sex.map({"F": 0, "M": 1, "pooled": 2})
    tab = tab.sort_values(["_o", "_s", "rna_week"]).drop(columns=["_o", "_s"]).reset_index(drop=True)
    ref = LAG_REF.get(tissue, {})
    refg = LAG_REF_GRID.get(tissue, {})
    tab["original_heldout_pooled"] = [ref.get((f, w), np.nan) if s == "pooled" else np.nan
                                      for f, w, s in zip(tab.family, tab.rna_week, tab.sex)]
    tab["original_in_sample_pooled"] = [refg.get(w, np.nan) if (s == "pooled" and f == "grid") else np.nan
                                        for f, w, s in zip(tab.family, tab.rna_week, tab.sex)]
    ngt = pd.DataFrame([dict(sex=s, n_genes=len(np.unique(D[s]["codes"])), n_rows=len(D[s]["prot"]))
                        for s in SEXES + ["pooled"]])

    pick = (lambda f, w, s="pooled": tab[(tab.family == f) & (tab.rna_week == w) & (tab.sex == s)].iloc[0])
    numbers = {}
    for f, ws in [("grid", ["1w", "2w", "4w"]), ("fixed", WEEKS), ("conditional", ["1w", "2w", "4w"])]:
        for w in ws:
            r = pick(f, w)
            numbers[f"{f}_{w}_heldout_dR2"] = r.dR2_heldout
            numbers[f"{f}_{w}_perm_p"] = r.perm_p_heldout
    numbers["n_genes_pooled"] = int(len(np.unique(D["pooled"]["codes"])))
    numbers["nperm"], numbers["folds"] = nperm, folds

    figures = {"delta_r2": _fig_lag(tab, tissue)}
    numbers["runtime_s"] = round(time.time() - t0, 2)

    notes = [
        "Target = 8-week protein logFC, predictors = logFC (as in scripts/07_lag.py; raw scale, no rank-INT). "
        "Every logFC is a trained-vs-sedentary group contrast from different animals per week, so this is a "
        "between-gene association, not a within-animal trajectory.",
        f"dR2_heldout, the permutation null and p use one {folds}-fold split over genes; "
        f"dR2_heldout_mean_repeats averages {CV_REPEATS} splits (07_lag: 20 x 10-fold) and is the estimate to quote. "
        f"The null permutes RNA(w) across genes {nperm} times on the same folds (07_lag: 1000). "
        "p = (1 + #null >= observed) / (1 + nperm), so the smallest reportable p is "
        f"{1 / (1 + nperm):.4f}.",
        "Pooled = all female rows + all male rows with a sex covariate (as 07_lag); genes in both sexes are "
        "permuted jointly, single-sex genes within sex. Cluster-robust CIs for beta are not computed here "
        "(no statsmodels); use the permutation p.",
        "All four weeks share the sedentary control, so part of any cross-week association is shared-control "
        "noise; M0 contains early protein, which absorbs that term for protein, not for RNA.",
        "The 8w row of the 'fixed' family is the same-time reference: if RNA led protein, early RNA should add "
        "more than 8-week RNA; the conditional family asks whether early RNA adds anything once 8-week RNA is in.",
    ]
    if ref:
        parts = []
        for (f, w), v in ref.items():
            parts.append(f"{f} {w}: {_fmt(pick(f, w).dR2_heldout_mean_repeats)} vs {_fmt(v)}")
        g = "; ".join(f"grid {w} in-sample: {_fmt(pick('grid', w).dR2_in_sample)} vs {_fmt(v)}"
                      for w, v in refg.items())
        notes.append(f"Comparison with lag/RESULTS.md ({tissue}, pooled, raw; held-out dR2 here (mean of "
                     f"{CV_REPEATS} x {folds}-fold) vs there (mean of 20 x 10-fold)): "
                     + "; ".join(parts) + f". {g}. Store rows = join_table_v2 as in the original; in-sample values match; held-out "
                     "differences come from the random fold assignment (5- vs 10-fold, 5 vs 20 repeats).")
    captions = {
        "models": (f"{tissue}: does RNA add information about 8-week protein?",
                   "One row per model x sex. dR2_heldout = extra fraction of between-gene variance in 8-week "
                   "protein logFC explained on held-out genes when RNA(w) is added; null_heldout_q05/q95 = the "
                   "same quantity with RNA permuted across genes; perm_p_heldout = permutation p. "
                   "'original_*' columns quote lag/RESULTS.md."),
        "n_genes": (f"{tissue}: genes used", "Genes with RNA and protein logFC at all four weeks (per sex)."),
        "delta_r2": (f"{tissue}: extra variance in 8-week protein explained by RNA at each week, on held-out genes",
                     "Dots = held-out dR2 (in % of variance, one 5-fold split) when RNA at that week is added; grey bars = 5-95% "
                     "of the same number with RNA shuffled across genes (the no-information null). A dot inside "
                     "its bar means no detectable information. Colours = what the model already contains."),
    }
    return dict(tables=dict(models=tab, n_genes=ngt), figures=figures, numbers=numbers, captions=captions,
                notes=notes)


def _fig_lag(tab, tissue):
    fams = [("grid", "over protein at the same week", "#7a5195", -0.22),
            ("fixed", "over protein at 2 wk", ACCENT, 0.0),
            ("conditional", "over protein at 2 wk + RNA at 8 wk", "#6b8e23", 0.22)]
    xpos = {w: i for i, w in enumerate(WEEKS)}
    with mpl_style():
        fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True)
        for ax, sex in zip(axes, ["F", "M", "pooled"]):
            ax.axhline(0, color=MUTED, lw=0.8)
            for fam, lab, col, off in fams:
                t = tab[(tab.family == fam) & (tab.sex == sex)]
                xs = np.array([xpos[w] for w in t.rna_week]) + off
                lo, hi = t.null_heldout_q05.to_numpy() * 100, t.null_heldout_q95.to_numpy() * 100
                ax.bar(xs, hi - lo, bottom=lo, width=0.16, color="#c4c4c4", edgecolor="none", zorder=1)
                ax.plot(xs, t.dR2_heldout * 100, "o", color=col, ms=7, zorder=3, label=f"+ RNA(w) {lab}")
            ax.set_xticks(range(4), ["1", "2", "4", "8"])
            ax.set_xlabel("week of the added RNA", fontsize=11)
            ax.set_title(SEX_WORDS[sex].capitalize() if sex != "pooled" else "Both sexes (pooled)", fontsize=12)
            ax.tick_params(labelsize=10)
        axes[0].set_ylabel("held-out ΔR² (% of variance)", fontsize=11)
        handles = [Line2D([], [], ls="", marker="o", ms=7, color=c, label=f"+ RNA(w) {lab}")
                   for _, lab, c, _ in fams] + [Patch(color="#c4c4c4", label="permutation null, 5–95%")]
        fig.legend(handles=handles, loc="lower center", ncol=2, fontsize=10, frameon=False,
                   bbox_to_anchor=(0.5, -0.02))
        fig.tight_layout(rect=(0, 0.13, 1, 1))
    return fig
