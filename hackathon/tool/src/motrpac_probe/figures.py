"""Matplotlib figures for the run report. Each function returns a Figure; captions live in report.py."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from .core import short_words  # noqa: E402

try:
    from .render import ACCENT, DIVERGING, INK, LAYER_COLORS, MUTED
except Exception:  # render layer not importable: fall back to the same palette
    from matplotlib.colors import LinearSegmentedColormap
    ACCENT, INK, MUTED = "#2f5f8a", "#111111", "#8a8a8a"
    DIVERGING = LinearSegmentedColormap.from_list("div", ["#1f4e79", "#f7f7f7", "#a51c30"])
    LAYER_COLORS = {"RNA": "#2f5f8a", "PROT": "#c26a2e", "PHOSPHO": "#6b8e23", "METAB": "#7a5195"}

BLUE_OPP = DIVERGING  # negative (blue) = exercise opposes the signature in the grid
FS = 11


def _blocks(S, cids):
    """Contiguous (start, end, label) blocks of columns sharing dataset x tissue x layer."""
    c = S.cols.loc[list(cids)]
    key = list(zip(c.dataset, c.tissue, c.layer))
    out, s = [], 0
    for i in range(1, len(key) + 1):
        if i == len(key) or key[i] != key[s]:
            ds, t, l = key[s]
            sp = "Human" if ds == "human_acute" else "Rat"
            tw = {"VL": "VL", "SKM-GN": "gastroc.", "SKM-VL": "VL", "HEART": "heart"}.get(t, t)
            out.append((s, i, f"{sp} {tw}\n{'protein' if l == 'PROT' else l}"))
            s = i
    return out


def grid(S, L, cids, sig, cap, cutoff, max_genes=60):
    genes = sig.shown[:max_genes]
    cidx = {c: i for i, c in enumerate(cids)}
    gidx = {g: i for i, g in enumerate(genes)}
    M = np.full((len(genes), len(cids)), np.nan)
    D = np.zeros_like(M, bool)
    for r in L[L.gene.isin(gidx)].itertuples():
        if np.isfinite(r.stat):
            M[gidx[r.gene], cidx[r.column_id]] = r.cell
            D[gidx[r.gene], cidx[r.column_id]] = r.fdr_bh < cutoff
    cw = 0.2
    fig, ax = plt.subplots(figsize=(max(7, len(cids) * cw + 3.2), max(3, len(genes) * cw + 2.6)))
    vmax = cap if cap else np.nanmax(np.abs(M)) if np.isfinite(M).any() else 4
    cmap = BLUE_OPP.copy()
    cmap.set_bad("white")
    im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, vmin=-vmax, vmax=vmax, aspect="equal", interpolation="none")
    yy, xx = np.where(np.isnan(M))
    for y, x in zip(yy, xx):
        ax.add_patch(Rectangle((x - .5, y - .5), 1, 1, fill=False, hatch="////", edgecolor="#d0d0d0", lw=0))
    yy, xx = np.where(D)
    ax.scatter(xx, yy, s=10, color=INK, zorder=3)
    ax.set_xticks(np.arange(-.5, len(cids)), minor=True)
    ax.set_yticks(np.arange(-.5, len(genes)), minor=True)
    ax.grid(which="minor", color="white", lw=0.8)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    for a, b, lab in _blocks(S, cids):
        if a:
            ax.axvline(a - .5, color=INK, lw=1.0)
        ax.text((a + b - 1) / 2, -1.0, lab, ha="center", va="bottom", fontsize=9)
    grp = [sig.groups.get(g, "") for g in genes]
    for i in range(1, len(genes)):
        if grp[i] != grp[i - 1]:
            ax.axhline(i - .5, color=INK, lw=0.8)
    c = S.cols.loc[list(cids)]
    ax.set_xticks(range(len(cids)), [short_words(r) for r in c.itertuples()], rotation=90, fontsize=9)
    dirs = dict(zip(sig.shown, sig.shown_dirs))
    ax.set_yticks(range(len(genes)), [f"{g} ({'up' if dirs[g] > 0 else 'down'})" for g in genes], fontsize=9)
    cb = fig.colorbar(im, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("agreement × |stat| (capped)\nblue = opposes, red = same", fontsize=9)
    cb.ax.tick_params(labelsize=9)
    fig.tight_layout()
    return fig


def camera(S, sc, cids, cutoff=0.05):
    c = S.cols.loc[list(cids)]
    x = np.arange(len(cids))
    s = sc.set_index("column_id").loc[list(cids)]
    fig, ax = plt.subplots(figsize=(max(8, len(cids) * 0.2 + 2.5), 4.2))
    ax.axhline(0, color=INK, lw=0.8)
    for v in (1.96, -1.96):
        ax.axhline(v, color=MUTED, lw=0.6, ls="--")
    for layer, col in LAYER_COLORS.items():
        m = (c.layer == layer).to_numpy()
        if not m.any():
            continue
        sigm = m & (s.camera_fdr < cutoff).to_numpy()
        ax.scatter(x[m & ~sigm], s.camera_t[m & ~sigm], s=36, facecolor="white", edgecolor=col, lw=1.4, zorder=3,
                   label=f"{'protein' if layer == 'PROT' else layer}, FDR ≥ {cutoff}")
        ax.scatter(x[sigm], s.camera_t[sigm], s=36, color=col, zorder=3,
                   label=f"{'protein' if layer == 'PROT' else layer}, FDR < {cutoff}")
        ax.scatter(x[m], s.camera_t_up[m], marker="^", s=14, color=col, alpha=0.55, lw=0, zorder=2)
        ax.scatter(x[m], s.camera_t_down[m], marker="v", s=14, color=col, alpha=0.55, lw=0, zorder=2)
    for a, b, lab in _blocks(S, cids):
        if a:
            ax.axvline(a - .5, color=INK, lw=0.6)
        ax.text((a + b - 1) / 2, 1.01, lab.replace("\n", " "), transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=8.5)
    ax.set_xticks(x, [short_words(r) for r in c.itertuples()], rotation=90, fontsize=9)
    ax.set_ylabel("cameraPR t (+ = exercise opposes)", fontsize=FS)
    ax.set_xlim(-.7, len(cids) - .3)
    ax.legend(fontsize=8.5, ncol=4, loc="upper left", bbox_to_anchor=(0, -0.42), frameon=False)
    fig.tight_layout()
    return fig


def nulls(S, sc, nl, cids):
    c = S.cols.loc[list(cids)]
    x = np.arange(len(cids))
    fig, ax = plt.subplots(figsize=(max(8, len(cids) * 0.2 + 2.5), 4.0))
    ax.axhline(0, color=INK, lw=0.8)
    for off, (kind, col, lab) in zip([-0.18, 0.18], [("abundance", "#9a9a9a", "abundance-matched null"),
                                                     ("class", ACCENT, "pathway-class-matched null")]):
        d = nl[nl.null == kind].set_index("column_id").reindex(list(cids))
        ax.vlines(x + off, d.null_t_q05, d.null_t_q95, color=col, lw=3, alpha=0.6, label=f"{lab}, 5–95%")
        ax.scatter(x + off, d.null_t_median, marker="_", s=40, color=col, zorder=3)
    s = sc.set_index("column_id").loc[list(cids)]
    ax.scatter(x, s.camera_t, s=22, color=INK, zorder=4, label="signature (observed)")
    for a, b, lab in _blocks(S, cids):
        if a:
            ax.axvline(a - .5, color=INK, lw=0.6)
        ax.text((a + b - 1) / 2, 1.01, lab.replace("\n", " "), transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=8.5)
    ax.set_xticks(x, [short_words(r) for r in c.itertuples()], rotation=90, fontsize=9)
    ax.set_ylabel("cameraPR t (+ = opposes)", fontsize=FS)
    ax.set_xlim(-.7, len(cids) - .3)
    ax.legend(fontsize=8.5, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.42), frameon=False)
    fig.tight_layout()
    return fig


def discordance(ds, G, show_pairs):
    fig, axes = plt.subplots(1, 1 + len(show_pairs), figsize=(4.4 * (1 + len(show_pairs)), 4.2))
    axes = np.atleast_1d(axes)
    ax = axes[0]
    d = ds[~ds.cross].reset_index(drop=True)
    x = np.arange(len(d))
    err = lambda lo, hi: [(lo).fillna(0).clip(lower=0), (hi).fillna(0).clip(lower=0)]  # rho = +-1 -> CI bound inside
    ax.errorbar(x - 0.15, d.rho_all, yerr=err(d.rho_all - d.rho_all_lo, d.rho_all_hi - d.rho_all), fmt="o", ms=4,
                color=MUTED, lw=1, label="all genes")
    ax.errorbar(x + 0.15, d.rho_sig, yerr=err(d.rho_sig - d.rho_sig_lo, d.rho_sig_hi - d.rho_sig),
                fmt="o", ms=4, color=ACCENT, lw=1, label="signature genes")
    ax.axhline(0, color=INK, lw=0.8)
    ax.set_xticks(x, [f"{'H' if r.dataset == 'human_acute' else 'R'} {r.tissue} {r.contrast}".replace("_vs_CON_", " ")
                      for r in d.itertuples()], rotation=90, fontsize=7.5)
    ax.set_ylabel("Spearman ρ, RNA vs protein stat", fontsize=FS)
    ax.set_title("Agreement per timepoint", loc="left", fontsize=FS)
    ax.legend(fontsize=8.5, frameon=False)
    for ax, (label, pr) in zip(axes[1:], show_pairs):
        g = G[(G.rna == pr[0]) & (G.prot == pr[1])]
        col = np.where(g.rna_opposed & g.prot_opposed, ACCENT, np.where(g.rna_opposed | g.prot_opposed, "#9a9a9a", "#a51c30"))
        ax.scatter(g.stat_rna * -g.direction, g.stat_prot * -g.direction, c=col, s=26, zorder=3)
        lab = g.assign(_m=(g.stat_rna.abs() + g.stat_prot.abs())).nlargest(25, "_m")  # label the 25 strongest
        for r in lab.itertuples():
            ax.annotate(r.gene_symbol_human, (r.stat_rna * -r.direction, r.stat_prot * -r.direction), fontsize=7.5,
                        xytext=(2, 2), textcoords="offset points")
        ax.axhline(0, color=INK, lw=0.6)
        ax.axvline(0, color=INK, lw=0.6)
        ax.set_xlabel("RNA stat × −direction (+ = opposed)", fontsize=10)
        ax.set_ylabel("protein stat × −direction", fontsize=10)
        ax.set_title(label, loc="left", fontsize=10)
    fig.tight_layout()
    return fig


def everywhere(rank):
    v = rank[rank.ranked].copy()
    kinds = ["exercise vs control", "exercise post−pre", "training vs sedentary", "reference (non-exercise)"]
    kinds = [k for k in kinds if k in set(v.kind)]
    fig, ax = plt.subplots(figsize=(9, 1.0 + 0.9 * len(kinds)))
    rng = np.random.default_rng(0)
    for i, k in enumerate(kinds):
        d = v[v.kind == k]
        y = i + rng.uniform(-0.25, 0.25, len(d))
        sig_ = (d.sign_q_bh < 0.05).to_numpy()
        col = MUTED if k.startswith("reference") else ACCENT
        ax.scatter(d.frac_opposed[~sig_], y[~sig_], s=14 + d.n_measured[~sig_] / 2, facecolor="white", edgecolor=col,
                   lw=1)
        ax.scatter(d.frac_opposed[sig_], y[sig_], s=14 + d.n_measured[sig_] / 2, color=col)
        ax.plot([d.frac_opposed.median()] * 2, [i - 0.35, i + 0.35], color=INK, lw=2)
    ax.axvline(0.5, color=INK, lw=0.6, ls="--")
    ax.set_yticks(range(len(kinds)), [f"{k} ({(v.kind == k).sum()})" for k in kinds], fontsize=10)
    ax.set_xlabel("fraction of measured signature genes opposed (per unit)", fontsize=FS)
    ax.set_xlim(-0.02, 1.02)
    ax.invert_yaxis()
    fig.tight_layout()
    return fig


def coverage(cov, blocks):
    M = np.array([[int(v.split("/")[0]) / max(1, int(v.split("/")[1])) for v in cov[b]] for b in blocks]).T
    fig, ax = plt.subplots(figsize=(1.6 + 0.62 * len(blocks), 1.2 + 0.2 * len(cov)))
    ax.imshow(M, cmap=matplotlib.colors.LinearSegmentedColormap.from_list("cov", ["#ffffff", ACCENT]), vmin=0, vmax=1,
              aspect="auto", interpolation="none")
    ax.set_xticks(range(len(blocks)), [f"{'Human' if d == 'human_acute' else 'Rat'} {t} {'protein' if l == 'PROT' else l}"
                                       for d, t, l in blocks], rotation=60, ha="right", fontsize=9)
    ax.set_yticks(range(len(cov)), cov.gene, fontsize=8.5 if len(cov) > 30 else 9.5)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    fig.tight_layout()
    return fig


def rank_association(S, ra, ref=None):
    """Spearman rho (disease ranking vs exercise statistic) per column with Fisher 95% CI; grey = reference
    (non-exercise) columns of the same tissue x layer."""
    d = ra.dropna(subset=["rho"]).merge(S.cols.reset_index(drop=True)[["column_id", "dataset", "tissue", "layer", "col_order"]],
                                        on="column_id").sort_values("col_order").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(7, 0.24 * len(d) + 2.5), 4.2))
    x = np.arange(len(d))
    for layer, col in LAYER_COLORS.items():
        m = (d.layer.str.startswith(layer)).to_numpy()
        if m.any():
            ax.errorbar(x[m], d.rho[m], yerr=[(d.rho - d.ci_low)[m].clip(lower=0), (d.ci_high - d.rho)[m].clip(lower=0)],
                        fmt="o", ms=4.5, color=col, lw=1, label="protein" if layer == "PROT" else layer)
    if ref is not None and len(ref):
        for i, r in enumerate(d.itertuples()):
            rr = ref[(ref.dataset == r.dataset) & (ref.tissue == r.tissue) & (ref.layer == r.layer)].rho.dropna()
            ax.scatter(np.full(len(rr), i), rr, s=8, color=MUTED, alpha=0.5, lw=0, zorder=1)
    ax.axhline(0, color=INK, lw=0.8)
    for a, b, lab in _blocks(S, d.column_id):
        if a:
            ax.axvline(a - .5, color=INK, lw=0.6)
        ax.text((a + b - 1) / 2, 1.01, lab.replace("\n", " "), transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=8.5)
    ax.set_xticks(x, [short_words(r) for r in S.cols.loc[d.column_id].itertuples()], rotation=90, fontsize=9)
    ax.set_ylabel("Spearman ρ (+ = same as disease, − = opposed)", fontsize=10)
    ax.set_xlim(-.7, len(d) - .3)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.42), ncol=4)
    fig.tight_layout()
    return fig


def pathway_counts(S, sm):
    """Jointly significant GO:BP sets per column, split into same and opposite direction."""
    d = sm.merge(S.cols.reset_index(drop=True)[["column_id", "col_order"]], on="column_id").sort_values("col_order").reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(max(7, 0.24 * len(d) + 2.5), 3.8))
    x = np.arange(len(d))
    ax.bar(x, d.same_direction, color="#a51c30", label="same direction as disease")
    ax.bar(x, -d.opposite_direction, color="#1f4e79", label="opposite direction (opposed)")
    ax.axhline(0, color=INK, lw=0.8)
    for a, b, lab in _blocks(S, d.column_id):
        if a:
            ax.axvline(a - .5, color=INK, lw=0.6)
        ax.text((a + b - 1) / 2, 1.01, lab.replace("\n", " "), transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=8.5)
    ax.set_xticks(x, [short_words(r) for r in S.cols.loc[d.column_id].itertuples()], rotation=90, fontsize=9)
    ax.set_ylabel("GO:BP sets, BH < 0.05 on both sides", fontsize=10)
    ax.set_xlim(-.7, len(d) - .3)
    ax.legend(fontsize=8.5, frameon=False, loc="upper left", bbox_to_anchor=(0, -0.42), ncol=2)
    fig.tight_layout()
    return fig
