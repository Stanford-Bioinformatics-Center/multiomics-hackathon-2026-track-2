"""Plot resting levels of the top PAH-versus-healthy metabolites in every group.

Each row is one of the 12 strongest resting PAH-versus-healthy results from
script 03. Each dot is a group's mean log2 level minus the healthy mean, so 0
is the healthy group. If PAH drives a difference, only the PAH dot moves away
from 0. If the sampling setting drives it, the three Cath groups move together.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from common import OUTPUT

# Group: (label, colour, marker). Validated categorical slots (blue, orange,
# aqua, violet); marker shape and fill carry the same identity for CVD/print.
STYLE = {
    "LowRisk": ("Low-risk SSc (non-invasive)", "#2a78d6", "o", "none"),
    "Normal Pressures": ("Normal-pressure SSc (Cath)", "#eb6834", "s", None),
    "Borderline Pressures": ("Borderline SSc (Cath)", "#1baf7a", "D", None),
    "PAH": ("SSc-PAH (Cath)", "#4a3aa7", "^", None),
}
INK, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"

res = pd.read_csv(OUTPUT / "03_resting_feature_tests.csv")
top = (res[res.q_welch_PAH_vs_Healthy < 0.05]
       .sort_values("q_welch_PAH_vs_Healthy")
       .drop_duplicates("refmet_name")
       .head(12)
       .iloc[::-1])  # strongest result at the top of the plot

fig, ax = plt.subplots(figsize=(8, 6.2), facecolor=SURFACE)
ax.set_facecolor(SURFACE)
offsets = [0.24, 0.08, -0.08, -0.24]
for (group, (label, colour, marker, face)), dy in zip(STYLE.items(), offsets):
    x = top[f"mean_log2_{group}"] - top.mean_log2_Healthy
    ax.scatter(x, [i + dy for i in range(len(top))], s=46, marker=marker, color=colour,
               facecolors=face if face else colour, linewidths=1.6, label=label, zorder=3)
ax.axvline(0, color=MUTED, linewidth=1.2, zorder=2)
ax.text(0.03, len(top) - 0.35, "Healthy\n(non-invasive)", color=MUTED, fontsize=8.5, va="top")
ax.set_yticks(range(len(top)))
ax.set_yticklabels(top.refmet_name, color=INK, fontsize=9.5)
ax.set_xlabel("Mean resting level minus healthy mean (log2)", color=INK, fontsize=10)
ax.grid(axis="x", color=GRID, linewidth=0.8)
ax.set_axisbelow(True)
for side in ["top", "right", "left"]:
    ax.spines[side].set_visible(False)
ax.spines["bottom"].set_color(MUTED)
ax.tick_params(colors=MUTED, length=0)
ax.set_title("Catheter-lab groups share the shift, with or without PAH",
             loc="left", color=INK, fontsize=12, fontweight="bold", pad=26)
ax.text(0, 1.015, "Top 12 resting PAH-vs-healthy results, ST000763 (script 03)",
        transform=ax.transAxes, color=MUTED, fontsize=9.5)
ax.legend(loc="lower center", bbox_to_anchor=(0.45, -0.24), ncol=2, frameon=False,
          fontsize=9, labelcolor=INK)
fig.tight_layout()
for ext in ("png", "svg"):
    fig.savefig(OUTPUT / f"07_resting_setting_effect.{ext}", dpi=200, facecolor=SURFACE)
top.iloc[::-1][["refmet_name", "q_welch_PAH_vs_Healthy"] +
               [f"mean_log2_{g}" for g in ["Healthy"] + list(STYLE)]].to_csv(
    OUTPUT / "07_resting_setting_effect_data.csv", index=False)
print("Saved output/07_resting_setting_effect.png and .svg")
