"""
Permutation histogram for PCA-OLS DBS results.
Pools all 1000-iteration permutation NRMSE values from
permutation_results/pca_dbs/ into one histogram, then marks
each patient's true NRMSE as a vertical line.

Zero NRMSE values are filtered out (they arise when every
LOO fold happens to predict exactly, which causes the numerator
of get_nrmse to be 0 while the denominator is non-zero —
a numerical artifact, not a meaningful result).

Run from project root: python paper_figures/perm_hist_pca_dbs.py
"""
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

PERM_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "permutation_results", "pca_dbs",
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "perm_hist_pca_dbs")
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
pooled_nrmse = []
patient_stats = {}  # subj -> {"nrmse": float, "p": float}

for subj in ALL_SUBJS:
    nrmse_csv = os.path.join(PERM_DIR, f"{subj}_OLS_LOO_ALL_permutation_nrmse.csv")
    summary_csv = os.path.join(PERM_DIR, f"{subj}_OLS_LOO_ALL_permutation_summary.csv")

    df = pd.read_csv(nrmse_csv)
    summary = pd.read_csv(summary_csv)

    perm_vals = df.loc[df["type"] == "permuted", "nrmse"].values
    true_vals = df.loc[df["type"] == "true", "nrmse"].values

    # filter zeros (numerical artifact, not meaningful)
    perm_vals = perm_vals[perm_vals != 0]
    true_vals = true_vals[true_vals != 0]

    pooled_nrmse.extend(perm_vals)
    if len(true_vals) > 0:
        p_val = float(summary["perm_p_value"].iloc[0])
        patient_stats[subj] = {"nrmse": float(true_vals[0]), "p": p_val}
    else:
        print(f"Warning: no valid true NRMSE for {subj}")

pooled_nrmse = np.array(pooled_nrmse)
print(f"Pooled permutation NRMSEs: {len(pooled_nrmse)} values across {len(ALL_SUBJS)} patients")

# ── plot ───────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4.5))
colors = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))

ax.hist(pooled_nrmse, bins=50, color="steelblue", alpha=0.6, edgecolor="none")

for color, subj in zip(colors, ALL_SUBJS):
    if subj not in patient_stats:
        continue
    true_val = patient_stats[subj]["nrmse"]
    p_val = patient_stats[subj]["p"]
    p_str = f"p<0.001" if p_val < 0.001 else f"p={p_val:.3f}"
    ax.axvline(true_val, color=color, linewidth=1.5, linestyle="--",
               label=f"{subj.replace('DBSTRD', 'TRD')}  {p_str}")

ax.set_xlabel("Permutation NRMSE", fontsize=6)
ax.set_ylabel("Count", fontsize=6)
leg = ax.legend(loc="upper right", fontsize=5, frameon=False)
sig_flags = [
    patient_stats.get(subj, {}).get("p", 1.0) < 0.05
    for subj in ALL_SUBJS
    if subj in patient_stats
]
for text, is_sig in zip(leg.get_texts(), sig_flags):
    if is_sig:
        text.set_fontweight("bold")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out_path = os.path.join(OUT_DIR, "perm_hist_pca_dbs.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
plt.show()
