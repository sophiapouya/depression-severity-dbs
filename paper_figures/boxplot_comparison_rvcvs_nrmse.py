"""
Comparison boxplot of per-patient NRMSE across four model/feature conditions:
  1. RVCVS PCA              (permutation_results/pca_dbs_rvcvs_only/RVCVS/)
  2. RVCVS sharpwave Ridge  (permutation_results/ridge_sharpwave_ldlpfc_rvcvs/)  [Top 3]
  3. RVCVS+L-DLPFC PCA      (permutation_results/dbs_seeg_pca_permutation_test/dbs_RVCVS_seeg_dlpfc_left/)
  4. RVCVS+L-DLPFC sharpwave Ridge (permutation_results/ridge_sharpwave_ldlpfc_rvcvs/)  [Top 3]

Conditions 1 and 4 load from *_permutation_summary.csv files.
Condition 3 loads from chunk CSVs + true_decoder.csv (PCA permutation format).

Star = significant patient (perm p < 0.05), circle = not significant.

Run from project root: python paper_figures/boxplot_comparison_rvcvs_nrmse.py
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
# Full list kept for color consistency with all other paper_figures scripts
ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
PLOT_SUBJS = [s for s in ALL_SUBJS if s != "DBSTRD011"]
SIG_THRESHOLD = 0.05

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Conditions with "pattern" load from *_permutation_summary.csv files.
# Conditions with "chunk_dir" load from PCA chunk + true_decoder format.
CONDITIONS = [
    {
        "label":   "R-VCVS\nAll",
        "pattern": os.path.join(BASE, "permutation_results", "pca_dbs_rvcvs_only",
                                "RVCVS", "*_permutation_summary.csv"),
    },
    {
        "label":   "R-VCVS\nTop 3",
        "pattern": os.path.join(BASE, "permutation_results", "ridge_sharpwave_rvcvs",
                                "*_permutation_summary.csv"),
    },
    {
        "label":   "R-VCVS\n+ L-dlPFC\nAll",
        "chunk_dir": os.path.join(BASE, "permutation_results", "dbs_seeg_pca_permutation_test",
                                  "dbs_RVCVS_seeg_dlpfc_left", "chunks"),
    },
    {
        "label":   "R-VCVS\n+ L-dlPFC\nTop 3",
        "pattern": os.path.join(BASE, "permutation_results", "ridge_sharpwave_ldlpfc_rvcvs",
                                "*_permutation_summary.csv"),
    },
]

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "boxplot_comparison_rvcvs_nrmse")
os.makedirs(OUT_DIR, exist_ok=True)


def load_chunk_condition(chunk_dir, label):
    """Load true NRMSE and compute permutation p-value from PCA chunk format."""
    rows = []
    for subj in PLOT_SUBJS:
        true_csv = os.path.join(chunk_dir, f"{subj}_true_decoder.csv")
        if not os.path.exists(true_csv):
            print(f"WARNING: missing true_decoder for {subj} in {chunk_dir}")
            continue
        true_df    = pd.read_csv(true_csv)
        true_nrmse = float(true_df["true_nrmse"].iloc[0])

        chunk_files = sorted(glob.glob(os.path.join(chunk_dir, f"{subj}_chunk*.csv")))
        if not chunk_files:
            print(f"WARNING: no chunks for {subj} in {chunk_dir}")
            continue
        perm_nrmse = pd.concat([pd.read_csv(f) for f in chunk_files])["nrmse"].values
        n_perm     = len(perm_nrmse)
        perm_p     = (np.sum(perm_nrmse <= true_nrmse) + 1) / (n_perm + 1)
        rows.append({
            "condition":   label,
            "patient_id":  subj,
            "nrmse":       true_nrmse,
            "significant": perm_p < SIG_THRESHOLD,
        })
    return rows


# ── load data ─────────────────────────────────────────────────────────────────
records = []
for cond in CONDITIONS:
    if "chunk_dir" in cond:
        records.extend(load_chunk_condition(cond["chunk_dir"], cond["label"]))
    else:
        files = sorted(glob.glob(cond["pattern"]))
        if not files:
            print(f"WARNING: no files found for {cond['label']} at {cond['pattern']}")
        for fpath in files:
            df = pd.read_csv(fpath)
            for _, row in df.iterrows():
                records.append({
                    "condition":   cond["label"],
                    "patient_id":  row["patient_id"],
                    "nrmse":       row["true_nrmse"],
                    "significant": row["perm_p_value"] < SIG_THRESHOLD,
                })

data = pd.DataFrame(records)
data = data[data["patient_id"] != "DBSTRD011"]
print(data.groupby("condition")[["nrmse", "significant"]].agg(["count", "mean"]))

# Colors anchored to the full 7-patient list so they match all other paper figures
_tab10 = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))
patient_colors = {p: _tab10[ALL_SUBJS.index(p)] for p in ALL_SUBJS}

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
cond_labels = [c["label"] for c in CONDITIONS]
x_positions = [i * 0.4 for i in range(len(CONDITIONS))]

fig, ax = plt.subplots(figsize=(44 * MM, 45 * MM))

box_data = [
    data[data["condition"] == label]["nrmse"].dropna().values
    for label in cond_labels
]
ax.boxplot(
    box_data,
    positions=x_positions,
    widths=0.25,
    patch_artist=True,
    showfliers=False,
    medianprops=dict(color="black", linewidth=0.8),
    boxprops=dict(facecolor="lightgray", alpha=0.6, linewidth=0.5),
    whiskerprops=dict(color="gray", linewidth=0.5),
    capprops=dict(color="gray", linewidth=0.5),
)

# evenly-spaced tiny jitter so all 7 patients sit close to centre
jitter_offsets = np.linspace(-0.04, 0.04, len(PLOT_SUBJS))
for i, label in enumerate(cond_labels):
    subset = data[data["condition"] == label].reset_index(drop=True)
    for j, (_, row) in enumerate(subset.iterrows()):
        color  = patient_colors.get(row["patient_id"], "gray")
        marker = "*" if row["significant"] else "o"
        size   = 5 if row["significant"] else 3
        ax.plot(x_positions[i] + jitter_offsets[j % 7], row["nrmse"],
                marker=marker, color=color, markersize=size,
                markeredgecolor="none", zorder=5)

ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.6, alpha=0.6)

ax.set_xticks(x_positions)
ax.set_xticklabels(cond_labels, fontsize=4, rotation=0, ha="center")
ax.set_ylabel("NRMSE", fontsize=5)
ax.set_xlim(x_positions[0] - 0.25, x_positions[-1] + 0.25)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

fig.tight_layout(pad=0.3)

for ext in ("png", "svg"):
    out_path = os.path.join(OUT_DIR, f"boxplot_comparison_rvcvs_nrmse.{ext}")
    fig.savefig(out_path, dpi=300, bbox_inches="tight",
                format=ext if ext == "svg" else None)
    print(f"Saved: {out_path}")

plt.show()
