"""
Pairwise paired t-test p-value grid for DBS PCA NRMSE across conditions.
Upper triangle shows p-values (color-coded by significance).
Diagonal and lower triangle are blank.
Run from project root: python paper_figures/pca_ttest_grid_dbs.py
"""
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

# ── style ──────────────────────────────────────────────────────────────────────
FS = 6
matplotlib.rcParams.update({
    "font.size": FS,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
})

# ── constants ──────────────────────────────────────────────────────────────────
CONDITIONS = ["ALL", "LEFT", "RIGHT", "LSCC", "RSCC", "LVCVS", "RVCVS", "SCC", "VCVS", "SEEG ALL"]
CONDITION_LABELS = ["ALL", "LEFT", "RIGHT", "L-SCC", "R-SCC", "L-VCVS", "R-VCVS", "SCC", "VCVS", "SEEG ALL"]
SIG_THRESHOLD = 0.05

DATA_DIR      = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/pca_permutation_test"
SEEG_DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                              "permutation_results", "pca_seeg")
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_ttest_grid_dbs")
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ──────────────────────────────────────────────────────────────────
records = []
for condition in CONDITIONS[:-1]:  # DBS conditions
    folder = os.path.join(DATA_DIR, condition)
    for fpath in sorted(glob.glob(os.path.join(folder, "*_permutation_summary.csv"))):
        df = pd.read_csv(fpath)
        for _, row in df.iterrows():
            records.append({
                "condition": condition,
                "patient_id": row["patient_id"],
                "true_nrmse": row["true_nrmse"],
            })

# SEEG ALL condition
for fpath in sorted(glob.glob(os.path.join(SEEG_DATA_DIR, "*_OLS_LOO_seeg_ALL_permutation_summary.csv"))):
    df = pd.read_csv(fpath)
    for _, row in df.iterrows():
        records.append({
            "condition":  "SEEG ALL",
            "patient_id": row["patient_id"],
            "true_nrmse": row["true_nrmse"],
        })

data = pd.DataFrame(records)
wide = data.pivot_table(index="patient_id", columns="condition", values="true_nrmse")
wide = wide.reindex(columns=CONDITIONS)

# ── compute upper-triangle p-values ───────────────────────────────────────────
n = len(CONDITIONS)
p_matrix = np.full((n, n), np.nan)

for i in range(n):
    for j in range(i + 1, n):
        paired = wide[[CONDITIONS[i], CONDITIONS[j]]].dropna()
        if len(paired) >= 2:
            _, p_val = ttest_rel(paired[CONDITIONS[i]].values, paired[CONDITIONS[j]].values)
            p_matrix[i, j] = p_val

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(100 * MM, 100 * MM))
ax.set_aspect("equal")
ax.axis("off")

cell = 1.0  # each cell is 1 unit wide/tall

for i in range(n):
    for j in range(n):
        x = j * cell
        # invert y so row 0 is at the top
        y = (n - 1 - i) * cell

        if i == j:
            # diagonal: condition label
            ax.text(
                x + cell / 2, y + cell / 2,
                CONDITION_LABELS[i],
                ha="center", va="center",
                fontsize=5.5, fontweight="bold",
            )

        elif j > i:
            # upper triangle: colored cell with p-value
            p = p_matrix[i, j]
            if np.isnan(p):
                color = "#dddddd"
            elif p < SIG_THRESHOLD:
                # scale from light to dark salmon based on how significant
                intensity = max(0.0, 1.0 - p / SIG_THRESHOLD)
                color = (1.0, 1.0 - 0.45 * intensity, 1.0 - 0.45 * intensity)
            else:
                color = "#f0f0f0"

            rect = plt.Rectangle(
                (x, y), cell, cell,
                facecolor=color, edgecolor="white", linewidth=0.5,
            )
            ax.add_patch(rect)

            if not np.isnan(p):
                if p < 0.001:
                    p_str = "p<0.001"
                elif p < 0.01:
                    p_str = f"p={p:.3f}"
                else:
                    p_str = f"p={p:.2f}"
                star = "*" if p < SIG_THRESHOLD else ""
                ax.text(
                    x + cell / 2, y + cell / 2,
                    f"{p_str}{star}",
                    ha="center", va="center",
                    fontsize=4.5,
                    color="black",
                )

        # lower triangle: leave blank (no patch, no text)

# ── condition labels on axes ───────────────────────────────────────────────────
for i, label in enumerate(CONDITION_LABELS):
    # row labels on the left
    ax.text(
        -0.05, (n - 1 - i) * cell + cell / 2,
        label, ha="right", va="center", fontsize=5,
    )
    # column labels on the top
    ax.text(
        i * cell + cell / 2, n * cell + 0.05,
        label, ha="center", va="bottom", fontsize=5,
    )

ax.set_xlim(-0.5, n * cell)
ax.set_ylim(-0.2, n * cell + 0.3)

fig.tight_layout(pad=0.3)

# ── save ───────────────────────────────────────────────────────────────────────
out_png = os.path.join(OUT_DIR, "pca_ttest_grid_dbs.png")
out_svg = os.path.join(OUT_DIR, "pca_ttest_grid_dbs.svg")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"Saved: {out_png}")
fig.savefig(out_svg, format="svg", bbox_inches="tight")
print(f"Saved: {out_svg}")
plt.show()
