"""
Permutation histogram for PCA-OLS SEEG results (greymatter filtered).
Reads chunk files from permutation_results/pca_seeg_gm/, merges on the fly,
pools all 1000-iteration permutation NRMSE values into one histogram, then marks
each patient's true NRMSE as a vertical line with p-value.
Bold labels indicate p < 0.05 (per-patient permutation test).

Run from project root: python paper_figures/perm_hist_pca_seeg_gm.py
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
N_CHUNKS = 20

PERM_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "permutation_results", "pca_seeg_gm",
)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "perm_hist_pca_seeg_gm")
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
pooled_nrmse = []
patient_stats = {}

for subj in ALL_SUBJS:
    true_csv = os.path.join(PERM_DIR, f"{subj}_true_decoder.csv")
    true_df = pd.read_csv(true_csv)
    true_nrmse = float(true_df["true_nrmse"].iloc[0])

    chunk_dfs = []
    missing = []
    for i in range(N_CHUNKS):
        chunk_csv = os.path.join(PERM_DIR, f"{subj}_chunk{i:03d}.csv")
        if os.path.exists(chunk_csv):
            chunk_dfs.append(pd.read_csv(chunk_csv))
        else:
            missing.append(i)
    if missing:
        print(f"WARNING {subj}: missing chunks {missing}")

    perm_vals = pd.concat(chunk_dfs, ignore_index=True)["nrmse"].values
    p_value = np.mean(perm_vals <= true_nrmse)

    print(f"{subj}: {len(perm_vals)} perms, true NRMSE={true_nrmse:.4f}, p={p_value:.4f}")
    pooled_nrmse.extend(perm_vals)
    patient_stats[subj] = {"nrmse": true_nrmse, "p": p_value}

pooled_nrmse = np.array(pooled_nrmse)
print(f"Pooled permutation NRMSEs: {len(pooled_nrmse)} values across {len(ALL_SUBJS)} patients")

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(45 * MM, 30 * MM))
colors = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))

ax.hist(pooled_nrmse, bins=50, color="lightpink", alpha=0.8, edgecolor="none")
ax.set_xlim(0, 4.5)
ax.set_xticks(np.arange(0, 5.0, 0.5))

for color, subj in zip(colors, ALL_SUBJS):
    true_val = patient_stats[subj]["nrmse"]
    p_val = patient_stats[subj]["p"]
    p_str = "p<0.001" if p_val < 0.001 else f"p={p_val:.3f}"
    ax.axvline(true_val, color=color, linewidth=1.5, linestyle="--",
               label=p_str)

ax.set_xlabel("Permutation NRMSE", fontsize=6)
ax.set_ylabel("Count", fontsize=6)
leg = ax.legend(loc="upper right", fontsize=5, frameon=False)
sig_flags = [patient_stats[subj]["p"] < 0.05 for subj in ALL_SUBJS]
for text, is_sig in zip(leg.get_texts(), sig_flags):
    if is_sig:
        text.set_fontweight("bold")

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out_path = os.path.join(OUT_DIR, "perm_hist_pca_seeg_gm.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
svg_path = os.path.join(OUT_DIR, "perm_hist_pca_seeg_gm.svg")
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path}")
plt.show()
