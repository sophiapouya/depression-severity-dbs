"""
Standalone patient legend for SEEG and DBS figures (PCA and Lasso scatter plots + histograms).
Exports a tight SVG containing only the legend: colored dashed line + patient label.
Colors match plt.cm.tab10(np.linspace(0, 1, 7)) over the ALL_SUBJS order.

Run from project root: python paper_figures/patient_legend.py
"""
import os
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
})

ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
LABELS = [s.replace("DBSTRD", "TRD") for s in ALL_SUBJS]
COLORS = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))

handles = [
    mlines.Line2D([], [], color=c, linestyle="--", linewidth=1.2)
    for c in COLORS
]

fig, ax = plt.subplots(figsize=(0.1, 0.1))
ax.set_visible(False)

leg = fig.legend(
    handles, LABELS,
    loc="center",
    fontsize=5,
    frameon=False,
    handlelength=1.6,
    handletextpad=0.4,
    labelspacing=0.35,
)

fig.canvas.draw()
bb = leg.get_window_extent(renderer=fig.canvas.get_renderer())
fig.set_size_inches(bb.width / fig.dpi, bb.height / fig.dpi)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "patient_legend")
os.makedirs(OUT_DIR, exist_ok=True)

svg_path = os.path.join(OUT_DIR, "patient_legend.svg")
fig.savefig(svg_path, format="svg", bbox_inches="tight", transparent=True)
print(f"Saved: {svg_path}")

png_path = os.path.join(OUT_DIR, "patient_legend.png")
fig.savefig(png_path, dpi=300, bbox_inches="tight", transparent=True)
print(f"Saved: {png_path}")
plt.show()
