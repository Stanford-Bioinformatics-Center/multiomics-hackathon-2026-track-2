"""Plot paper-selected pathways over official EE/RE control-adjusted CAMERA data."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, linkage

MODULE = Path(__file__).resolve().parents[1]
OUT = MODULE / "outputs"
EE = ["during_20_min", "during_40_min", "post_10_min",
      "post_15_30_45_min", "post_3.5_4_hr", "post_24_hr"]
COLS = [("EE", t) for t in EE] + [("RE", t) for t in EE[2:]]
TIME_LABELS = {
    "during_20_min": "during 20 min", "during_40_min": "during 40 min",
    "post_10_min": "post 10 min", "post_15_30_45_min": "post 15/30/45 min",
    "post_3.5_4_hr": "post 3.5/4 hr", "post_24_hr": "post 24 hr",
}
TIME_COLORS = dict(zip(EE, ["#d4c883", "#b3bd94", "#c6c4d0",
                            "#bcb0c2", "#aa9aac", "#92879e"]))
MODALITY_COLORS = {"EE": "#c9ac9a", "RE": "#8caea8"}


def size_for_q(q):
    # Circle area represents the published BH-adjusted CAMERA P value.
    return 20 + 250 * np.minimum(-np.log10(np.maximum(q, 1e-6)), 6) / 6


def draw(mapped, stem):
    data = mapped
    sets = data["set"].drop_duplicates().tolist()
    wide = data.pivot(index="set", columns=["modality", "timepoint"],
                      values="z.std").reindex(index=sets, columns=COLS)
    if wide.isna().any().any():
        raise ValueError("Plot would have missing set/contrast cells")
    # Cluster only rows. Chronological columns stay fixed, separated by modality.
    tree = linkage(wide.to_numpy(), method="average", metric="euclidean",
                   optimal_ordering=True)
    ordering = dendrogram(tree, no_plot=True, orientation="left")["leaves"]
    rows = [sets[i] for i in ordering]
    wide = wide.loc[rows]
    qwide = data.pivot(index="set", columns=["modality", "timepoint"],
                       values="adj_p_value").reindex(index=rows, columns=COLS)

    fig = plt.figure(figsize=(18.2, 6.7), dpi=170,
                     facecolor="white")
    gs = GridSpec(3, 4, figure=fig, width_ratios=[1.2, 10, 3.1, 3.7],
                  height_ratios=[0.34, 0.34, 5], wspace=0.05, hspace=0.035)
    dend_ax = fig.add_subplot(gs[2, 0])
    strip1 = fig.add_subplot(gs[0, 1])
    strip2 = fig.add_subplot(gs[1, 1])
    ax = fig.add_subplot(gs[2, 1])
    legend_ax = fig.add_subplot(gs[:, 3])
    legend_ax.axis("off")

    dendrogram(tree, orientation="left", ax=dend_ax, no_labels=True,
               color_threshold=0, above_threshold_color="#73808a")
    dend_ax.invert_yaxis()
    dend_ax.axis("off")

    for strip, kind in ((strip1, "modality"), (strip2, "timepoint")):
        strip.set_xlim(-0.5, len(COLS) - 0.5)
        strip.set_ylim(0, 1)
        for j, (modality, time) in enumerate(COLS):
            color = MODALITY_COLORS[modality] if kind == "modality" else TIME_COLORS[time]
            strip.add_patch(Rectangle((j - .48, .06), .96, .88,
                                      facecolor=color, edgecolor="#dddddd", lw=.6))
        strip.set_xticks([])
        strip.set_yticks([])
        for spine in strip.spines.values():
            spine.set_visible(False)

    norm = TwoSlopeNorm(vmin=-2, vcenter=0, vmax=5.2)
    cmap = plt.get_cmap("coolwarm")
    for i, name in enumerate(rows):
        for j, col in enumerate(COLS):
            z = wide.loc[name, col]
            q = qwide.loc[name, col]
            ax.scatter(j, i, s=size_for_q(q), c=[cmap(norm(z))],
                       edgecolor="#776f73", linewidth=.45, zorder=3)
    ax.set_xlim(-.5, len(COLS) - .5)
    ax.set_ylim(len(rows) - .5, -.5)
    ax.set_aspect("equal")
    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels([TIME_LABELS[t] for _, t in COLS], rotation=55,
                       ha="right", fontsize=8)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows, fontsize=10)
    ax.yaxis.tick_right()
    ax.tick_params(axis="both", length=0, pad=6)
    ax.set_xticks(np.arange(-.5, len(COLS), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="#dedede", linewidth=.8)
    ax.tick_params(which="minor", length=0)
    ax.axvline(5.5, color="#6d777c", lw=2.2)
    for spine in ax.spines.values():
        spine.set_color("#777777")

    legend_ax.text(.05, .94, "Modality", fontsize=11, fontweight="bold",
                   transform=legend_ax.transAxes)
    for y, (name, label) in zip([.88, .83], [("EE", "Endurance"),
                                            ("RE", "Resistance")]):
        legend_ax.add_patch(Rectangle((.05, y), .075, .04,
                                      color=MODALITY_COLORS[name],
                                      transform=legend_ax.transAxes))
        legend_ax.text(.15, y+.02, f"{name}  {label}", va="center", fontsize=9,
                       transform=legend_ax.transAxes)
    legend_ax.text(.05, .73, "Time point", fontsize=11, fontweight="bold",
                   transform=legend_ax.transAxes)
    for k, t in enumerate(EE):
        y = .68 - k*.052
        legend_ax.add_patch(Rectangle((.05, y), .075, .035,
                                      color=TIME_COLORS[t],
                                      transform=legend_ax.transAxes))
        legend_ax.text(.15, y+.018, TIME_LABELS[t], va="center", fontsize=8.5,
                       transform=legend_ax.transAxes)
    cax = legend_ax.inset_axes([.05, .19, .08, .16])
    colorbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=cax)
    colorbar.set_ticks([-2, 0, 2, 5])
    colorbar.ax.tick_params(labelsize=8)
    legend_ax.text(.05, .37, "CAMERA Z-score", fontsize=10, fontweight="bold",
                   transform=legend_ax.transAxes)
    legend_ax.text(.39, .33, "BH-adjusted P", fontsize=10, fontweight="bold",
                   transform=legend_ax.transAxes)
    for k, q in enumerate([1, .05, .001, .00001]):
        y = .275 - k*.065
        legend_ax.scatter(.45, y, s=size_for_q(q), c="#ad9d9d",
                          edgecolor="#776f73", linewidth=.45,
                          transform=legend_ax.transAxes)
        legend_ax.text(.59, y, f"{q:g}", va="center", fontsize=8.5,
                       transform=legend_ax.transAxes)

    title = "Table 2 BioCarta pathways in MoTrPAC blood RNA"
    subtitle = "All six BioCarta pathways reported in Table 2"
    fig.suptitle(title + "\n" + subtitle, x=.36, y=.985, fontsize=14)
    fig.text(.36, .015, "MoTrPAC: exercise minus matched control; paper: disease minus healthy. "
             "Scores come from different tests and cohorts.", ha="center", fontsize=8,
             color="#505050")
    fig.subplots_adjust(left=.045, right=.98, top=.84, bottom=.25)
    fig.savefig(OUT / f"{stem}.png", dpi=200, facecolor="white")
    fig.savefig(OUT / f"{stem}.svg", facecolor="white")
    plt.close(fig)
    print("Wrote", stem, "with row order:", ", ".join(rows))


def main():
    mapped = pd.read_csv(OUT / "paper_motrpac_biocarta_mapped.csv")
    draw(mapped, "biocarta_all_six_pathways")


if __name__ == "__main__":
    main()
