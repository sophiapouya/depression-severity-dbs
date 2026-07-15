"""
Permutation histograms for DBS + patient-specific best 2-SEEG-region combo,
for two DBS conditions: RVCVS and SCC.

Best combo = the L+R SEEG region pair with lowest true NRMSE for each patient.
Reads pre-computed chunk files from permutation_results/dbs_seeg_pca_permutation_test/.
Pools all permutation NRMSE values into one histogram, marks each patient's
true NRMSE as a colored vertical dashed line. Bold label = p < 0.05.

Run from project root: python paper_figures/perm_hist_seeg_dbs_best.py
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
ALL_SUBJS  = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
              "DBSTRD010", "DBSTRD011", "DBSTRD014"]
PLOT_SUBJS = [s for s in ALL_SUBJS if s != "DBSTRD011"]
L_REGIONS  = ["Lacc", "Lamy", "Ldlpfc", "Lofc", "Lvmpfc"]
R_REGIONS  = ["Racc", "Ramy", "Rdlpfc", "Rofc", "Rvmpfc"]

BASE     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEG_DBS = os.path.join(BASE, "permutation_results", "dbs_seeg_pca_permutation_test")

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "perm_hist_seeg_dbs_best")
os.makedirs(OUT_DIR, exist_ok=True)

_tab10 = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))
patient_colors = {p: _tab10[ALL_SUBJS.index(p)] for p in ALL_SUBJS}


def get_best_combo_dir(prefix, patient_id):
    """Return chunks dir for the combo with lowest true NRMSE for this patient."""
    best_nrmse = float("inf")
    best_dir   = None
    for d in os.listdir(SEEG_DBS):
        if not d.startswith(f"dbs_{prefix}_seeg_"):
            continue
        suffix = d.split(f"dbs_{prefix}_seeg_")[1]
        parts  = suffix.split("_")
        if len(parts) != 2 or parts[0] not in L_REGIONS or parts[1] not in R_REGIONS:
            continue
        csv = os.path.join(SEEG_DBS, d, "chunks", f"{patient_id}_true_decoder.csv")
        if not os.path.exists(csv):
            continue
        v = float(pd.read_csv(csv)["true_nrmse"].iloc[0])
        if v < best_nrmse:
            best_nrmse = v
            best_dir   = os.path.join(SEEG_DBS, d, "chunks")
            best_combo = suffix
    return best_dir, best_combo if best_dir else None


def run_condition(dbs_prefix, fname):
    pooled_nrmse  = []
    patient_stats = {}

    for pid in PLOT_SUBJS:
        chunks_dir, combo = get_best_combo_dir(dbs_prefix, pid)
        if chunks_dir is None:
            print(f"WARNING: no best combo for {pid} / {dbs_prefix}")
            continue

        true_csv   = os.path.join(chunks_dir, f"{pid}_true_decoder.csv")
        true_nrmse = float(pd.read_csv(true_csv)["true_nrmse"].iloc[0])

        chunk_files = sorted(glob.glob(os.path.join(chunks_dir, f"{pid}_chunk*.csv")))
        if not chunk_files:
            print(f"WARNING: no chunk files for {pid} in {chunks_dir}")
            continue
        perm_vals = pd.concat([pd.read_csv(f) for f in chunk_files])["nrmse"].values
        p_value   = (np.sum(perm_vals <= true_nrmse) + 1) / (len(perm_vals) + 1)

        print(f"{pid} ({dbs_prefix}+{combo}): {len(perm_vals)} perms, "
              f"true NRMSE={true_nrmse:.4f}, p={p_value:.4f}")
        pooled_nrmse.extend(perm_vals)
        patient_stats[pid] = {"nrmse": true_nrmse, "p": p_value, "combo": combo}

    pooled_nrmse = np.array(pooled_nrmse)

    # ── plot ──────────────────────────────────────────────────────────────────
    MM = 1 / 25.4
    fig, ax = plt.subplots(figsize=(45 * MM, 30 * MM))

    ax.hist(pooled_nrmse, bins=50, color="lightpink", alpha=0.8, edgecolor="none")
    ax.set_xlim(0, 4.5)
    ax.set_xticks(np.arange(0, 5.0, 0.5))

    for pid in PLOT_SUBJS:
        if pid not in patient_stats:
            continue
        color    = patient_colors[pid]
        true_val = patient_stats[pid]["nrmse"]
        p_val    = patient_stats[pid]["p"]
        p_str    = "p<0.001" if p_val < 0.001 else f"p={p_val:.3f}"
        ax.axvline(true_val, color=color, linewidth=1.5, linestyle="--", label=p_str)

    ax.set_xlabel("Permutation NRMSE", fontsize=6)
    ax.set_ylabel("Count", fontsize=6)
    leg = ax.legend(loc="upper right", fontsize=5, frameon=False)
    for text, pid in zip(leg.get_texts(), PLOT_SUBJS):
        if pid in patient_stats and patient_stats[pid]["p"] < 0.05:
            text.set_fontweight("bold")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout(pad=0.3)
    for ext in ("png", "svg"):
        out = os.path.join(OUT_DIR, f"{fname}.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight",
                    format=ext if ext == "svg" else None)
        print(f"Saved: {out}")
    plt.close(fig)


run_condition("RVCVS", "perm_hist_rvcvs_best")
run_condition("SCC",   "perm_hist_scc_best")
plt.show()
