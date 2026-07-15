"""
Boxplot of per-patient true NRMSE across DBS contact conditions.
No title. Star = significant (p < 0.05), circle = not.
Run from project root: python paper_figures/pca_boxplot_dbs.py
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ── style ──────────────────────────────────────────────────────────────────────
FS = 5
matplotlib.rcParams.update({
    "font.size": FS,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": FS,
    "ytick.labelsize": FS,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "lines.linewidth": 0.75,
})

# ── constants ─────────────────────────────────────────────────────────────────
ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
CONDITIONS = ["ALL", "LEFT", "RIGHT", "LSCC", "RSCC", "LVCVS", "RVCVS", "SCC", "VCVS"]
SIG_THRESHOLD = 0.05

DATA_DIR      = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/pca_permutation_test"
SEEG_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "permutation_results", "pca_seeg")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_boxplot_dbs")
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
records = []
for condition in CONDITIONS:
    folder = os.path.join(DATA_DIR, condition)
    for fpath in sorted(glob.glob(os.path.join(folder, "*_permutation_summary.csv"))):
        df = pd.read_csv(fpath)
        for _, row in df.iterrows():
            records.append({
                "condition": condition,
                "patient_id": row["patient_id"],
                "true_nrmse": row["true_nrmse"],
                "significant": row["perm_p_value"] < SIG_THRESHOLD,
            })

# SEEG ALL condition
for fpath in sorted(glob.glob(os.path.join(SEEG_DATA_DIR, "*_OLS_LOO_seeg_ALL_permutation_summary.csv"))):
    df = pd.read_csv(fpath)
    for _, row in df.iterrows():
        records.append({
            "condition":   "SEEG ALL",
            "patient_id":  row["patient_id"],
            "true_nrmse":  row["true_nrmse"],
            "significant": row["perm_p_value"] < SIG_THRESHOLD,
        })

data = pd.DataFrame(records)
patients = sorted(data["patient_id"].unique())
patient_colors = {
    p: plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))[ALL_SUBJS.index(p)]
    for p in patients if p in ALL_SUBJS
}

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
ALL_CONDITIONS = CONDITIONS + ["SEEG ALL"]
n_cond = len(ALL_CONDITIONS)

x_positions = [i * 0.6 for i in range(n_cond)]

fig, ax = plt.subplots(figsize=(105 * MM, 45 * MM))

box_data = [data[data["condition"] == c]["true_nrmse"].dropna().values for c in ALL_CONDITIONS]
ax.boxplot(
    box_data,
    positions=x_positions,
    widths=0.35,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(color="black", linewidth=0.8),
    boxprops=dict(facecolor="lightgray", alpha=0.6, linewidth=0.5),
    whiskerprops=dict(color="gray", linewidth=0.5),
    capprops=dict(color="gray", linewidth=0.5),
)

rng = np.random.default_rng(42)
for i, condition in enumerate(ALL_CONDITIONS):
    subset = data[data["condition"] == condition]
    n = len(subset)
    jitter = rng.uniform(-0.08, 0.08, size=n)
    for j, (_, row) in enumerate(subset.iterrows()):
        x = x_positions[i] + jitter[j]
        y = row["true_nrmse"]
        color = patient_colors.get(row["patient_id"], "gray")
        marker = "*" if row["significant"] else "o"
        size = 6 if row["significant"] else 2
        ax.plot(x, y, marker=marker, color=color, markersize=size,
                markeredgecolor="none", zorder=5)

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.6, alpha=0.6)

ax.set_xticks(x_positions)
CONDITION_LABELS = ["ALL", "LEFT", "RIGHT", "L-SCC", "R-SCC", "L-VCVS", "R-VCVS", "SCC", "VCVS", "SEEG\nALL"]
ax.set_xticklabels(CONDITION_LABELS, fontsize=5, rotation=0, ha="center")
ax.set_ylim(top=1.45)
ax.set_ylabel("True NRMSE", fontsize=5)
ax.set_xlim(x_positions[0] - 0.4, x_positions[-1] + 0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout(pad=0.3)

out_path = os.path.join(OUT_DIR, "pca_boxplot_dbs.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
svg_path = os.path.join(OUT_DIR, "pca_boxplot_dbs.svg")
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path}")
plt.show()
