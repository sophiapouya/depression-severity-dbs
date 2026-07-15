"""
Heatmaps for DBS × SEEG combined PCA permutation results.
Loads directly from chunks in permutation_results/dbs_seeg_pca_permutation_test/.
Produces 4 figures matching dbs_seeg_perm_heatmaps.py style:
  1. Per-subject heatmap: 7-panel (one per patient), each DBS × SEEG p-value
  2. Mean NRMSE heatmap:        DBS × SEEG
  3. Mean perm p-value heatmap: DBS × SEEG
  4. Fraction significant:      DBS × SEEG
X-axis = all 35 SEEG conditions (10 unilateral + 25 bilateral).
Y-axis = 7 DBS conditions.
Run from project root: python paper_figures/pca_perm_heatmaps_dbs_seeg.py
"""
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import pandas as pd
import seaborn as sns

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
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

# ── colormaps ──────────────────────────────────────────────────────────────────
CMAP_G2P = LinearSegmentedColormap.from_list("green_pink", [
    "#1a7a4a", "#95d5b2", "#fce4ec", "#c2185b",
])
CMAP_P2G = CMAP_G2P.reversed()

# ── constants ──────────────────────────────────────────────────────────────────
ALL_SUBJS      = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
                  "DBSTRD010", "DBSTRD011", "DBSTRD014"]
PATIENT_LABELS = [p.replace("DBSTRD", "TRD") for p in ALL_SUBJS]

DBS_REGIONS    = ["ALL", "LSCC", "RSCC", "LVCVS", "RVCVS", "SCC", "VCVS"]
DBS_LABELS     = {"ALL": "ALL", "LSCC": "L-SCC", "RSCC": "R-SCC",
                  "LVCVS": "L-VCVS", "RVCVS": "R-VCVS", "SCC": "SCC", "VCVS": "VCVS"}
DBS_ROW_LABELS = [DBS_LABELS[d] for d in DBS_REGIONS]

SEEG_BASE = ["acc", "amy", "dlpfc", "ofc", "vmpfc"]
SEEG_UNI  = [f"{r}_{s}" for r in SEEG_BASE for s in ("left", "right")]   # 10
SEEG_BI   = [f"L{r1}_R{r2}" for r1 in SEEG_BASE for r2 in SEEG_BASE]    # 25
SEEG_ALL  = SEEG_UNI + SEEG_BI                                            # 35

SEEG_UNI_LABELS = [("L-" if "_left" in s else "R-") + s.split("_")[0] for s in SEEG_UNI]
SEEG_BI_LABELS  = [f"L-{r1}/R-{r2}" for r1 in SEEG_BASE for r2 in SEEG_BASE]
SEEG_ALL_LABELS = SEEG_UNI_LABELS + SEEG_BI_LABELS

N_UNI = len(SEEG_UNI)   # 10  (divider line position)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR     = os.path.join(PROJECT_ROOT, "permutation_results", "dbs_seeg_pca_permutation_test")
OUT_DIR      = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_perm_heatmaps_dbs_seeg")
CACHE_CSV    = os.path.join(OUT_DIR, "_all_results_cache.csv")
os.makedirs(OUT_DIR, exist_ok=True)

MM = 1 / 25.4


def save(fig, name):
    fig.savefig(os.path.join(OUT_DIR, f"{name}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(OUT_DIR, f"{name}.svg"), format="svg", bbox_inches="tight")
    print(f"Saved: {name}")


# ── load / build data ──────────────────────────────────────────────────────────
if os.path.exists(CACHE_CSV):
    print(f"Loading cached results from {CACHE_CSV}")
    df = pd.read_csv(CACHE_CSV)
else:
    print("Building results from chunks (this may take a few minutes)...")
    records = []
    cond_dirs = sorted(glob.glob(os.path.join(DATA_DIR, "*/")))
    for k, cond_dir in enumerate(cond_dirs):
        cond_label = os.path.basename(cond_dir.rstrip("/"))
        m = re.match(r"dbs_(.+)_seeg_(.+)", cond_label)
        if not m:
            continue
        dbs_region = m.group(1)
        seeg_combo = m.group(2)
        chunks_dir = os.path.join(cond_dir, "chunks")

        true_files = glob.glob(os.path.join(chunks_dir, "*_true_decoder.csv"))
        if not true_files:
            continue
        true_data = {}
        for tf in true_files:
            pid = os.path.basename(tf).replace("_true_decoder.csv", "")
            row = pd.read_csv(tf).iloc[0]
            true_data[pid] = {
                "true_nrmse": float(row["true_nrmse"]),
                "r_val":      float(row["r_val"]),
                "r_squared":  float(row["r_squared"]),
            }

        chunk_files = glob.glob(os.path.join(chunks_dir, "*_chunk*.csv"))
        if chunk_files:
            all_chunks = pd.concat(
                [pd.read_csv(f, usecols=["patient_id", "nrmse"]) for f in chunk_files],
                ignore_index=True,
            )
            perm_groups = {pid: grp["nrmse"].values
                           for pid, grp in all_chunks.groupby("patient_id")}
        else:
            perm_groups = {}

        for pid, td in true_data.items():
            perm  = perm_groups.get(pid)
            p_val = float(np.mean(perm <= td["true_nrmse"])) if perm is not None else np.nan
            records.append({
                "patient_id":   pid,
                "dbs_region":   dbs_region,
                "seeg_combo":   seeg_combo,
                "true_nrmse":   td["true_nrmse"],
                "r_val":        td["r_val"],
                "r_squared":    td["r_squared"],
                "perm_p_value": p_val,
                "sig_0.05":     float(p_val < 0.05) if not np.isnan(p_val) else np.nan,
            })

        if (k + 1) % 50 == 0:
            print(f"  {k+1}/{len(cond_dirs)} conditions processed...")

    df = pd.DataFrame(records)
    df.to_csv(CACHE_CSV, index=False)
    print(f"Cached to {CACHE_CSV}")

def add_divider(ax):
    """Vertical line between unilateral and bilateral SEEG sections."""
    ax.axvline(N_UNI, color="black", linewidth=1.0)


# ── figure 1: per-subject (7 panels, one per patient) ─────────────────────────
fig1, axes1 = plt.subplots(
    1, len(ALL_SUBJS),
    figsize=(160 * MM, 38 * MM),
    sharey=True,
)

for ax, pid, plabel in zip(axes1, ALL_SUBJS, PATIENT_LABELS):
    sub = df[df["patient_id"] == pid]
    pivot = sub.pivot_table(
        index="dbs_region", columns="seeg_combo",
        values="perm_p_value", aggfunc="first",
    ).reindex(index=DBS_REGIONS, columns=SEEG_ALL)

    pivot_sig = sub.pivot_table(
        index="dbs_region", columns="seeg_combo",
        values="sig_0.05", aggfunc="first",
    ).reindex(index=DBS_REGIONS, columns=SEEG_ALL)
    ann = pivot_sig.map(lambda v: "*" if v == 1.0 else "")

    sns.heatmap(
        pivot, ax=ax,
        cmap=CMAP_G2P, vmin=0, vmax=0.2,
        annot=ann, annot_kws={"fontsize": 3, "color": "black", "fontweight": "bold"},
        fmt="",
        xticklabels=False, yticklabels=(DBS_ROW_LABELS if ax is axes1[0] else False),
        linewidths=0.15, linecolor="white",
        cbar=False,
    )
    add_divider(ax)
    ax.set_title(plabel, fontsize=5, pad=2)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(length=0)

axes1[0].set_yticklabels(axes1[0].get_yticklabels(), fontsize=5, rotation=0)

cbar1 = fig1.colorbar(
    axes1[-1].collections[0], ax=axes1,
    orientation="vertical", shrink=0.8, pad=0.02,
)
cbar1.set_label("Permutation p-value", fontsize=5, labelpad=4)
cbar1.ax.tick_params(labelsize=4, length=2)
cbar1.outline.set_visible(False)

fig1.tight_layout(pad=0.3, w_pad=0.2)
save(fig1, "dbs_seeg_pca_subject_heatmap")


# ── figures 2–4: summary heatmaps (DBS × SEEG) ────────────────────────────────
def summary_heatmap(pivot_data, cmap, cbar_label, fmt, vmin=None, vmax=None):
    fig, ax = plt.subplots(figsize=(120 * MM, 42 * MM))
    sns.heatmap(
        pivot_data, ax=ax,
        cmap=cmap, vmin=vmin, vmax=vmax,
        annot=True, fmt=fmt,
        annot_kws={"fontsize": 3},
        linewidths=0.2, linecolor="white",
        xticklabels=SEEG_ALL_LABELS, yticklabels=DBS_ROW_LABELS,
        cbar=False,
    )
    add_divider(ax)
    ax.set_xticklabels(SEEG_ALL_LABELS, fontsize=3.5, rotation=45, ha="right")
    ax.set_yticklabels(DBS_ROW_LABELS, fontsize=5, rotation=0)
    ax.set_xlabel("SEEG Region", fontsize=5)
    ax.set_ylabel("DBS Condition", fontsize=5)
    ax.tick_params(length=0)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("bottom", size="6%", pad=0.5)
    cbar = fig.colorbar(ax.collections[0], cax=cax, orientation="horizontal")
    cax.xaxis.set_ticks_position("bottom")
    cax.xaxis.set_label_position("bottom")
    cbar.set_label(cbar_label, fontsize=5, labelpad=3)
    cax.tick_params(labelsize=4, length=2)
    cbar.outline.set_visible(False)

    fig.tight_layout(pad=0.3)
    return fig


def make_pivot(col, agg="mean"):
    return df.pivot_table(
        index="dbs_region", columns="seeg_combo",
        values=col, aggfunc=agg,
    ).reindex(index=DBS_REGIONS, columns=SEEG_ALL)


fig2 = summary_heatmap(make_pivot("true_nrmse"),   CMAP_G2P, "Mean NRMSE",            ".2f")
fig3 = summary_heatmap(make_pivot("perm_p_value"), CMAP_G2P, "Mean perm p-value",     ".2f", vmin=0, vmax=0.25)
fig4 = summary_heatmap(make_pivot("sig_0.05"),     CMAP_P2G, "Fraction sig (p<0.05)", ".2f", vmin=0, vmax=1)

save(fig2, "dbs_seeg_pca_mean_nrmse")
save(fig3, "dbs_seeg_pca_mean_pvalue")
save(fig4, "dbs_seeg_pca_frac_sig")

# ── figure 5: per-patient perm p-value heatmaps (consistent power-norm scale) ─
# PowerNorm gamma<1 stretches the green (low p) end so small differences are
# visible, while keeping an identical scale across all patients.
from matplotlib.colors import PowerNorm

P5_VMAX = 0.25
P5_NORM = PowerNorm(gamma=0.4, vmin=0, vmax=P5_VMAX)

n_patients = len(ALL_SUBJS)
fig5, axes5 = plt.subplots(
    1, n_patients,
    figsize=(120 * MM * n_patients, 52 * MM),
    sharey=True,
)

for ax, pid, plabel in zip(axes5, ALL_SUBJS, PATIENT_LABELS):
    sub = df[df["patient_id"] == pid]
    pivot = sub.pivot_table(
        index="dbs_region", columns="seeg_combo",
        values="perm_p_value", aggfunc="first",
    ).reindex(index=DBS_REGIONS, columns=SEEG_ALL)

    sns.heatmap(
        pivot, ax=ax,
        cmap=CMAP_G2P, norm=P5_NORM,
        annot=True, fmt=".2f",
        annot_kws={"fontsize": 2.5},
        linewidths=0.2, linecolor="white",
        xticklabels=SEEG_ALL_LABELS,
        yticklabels=(DBS_ROW_LABELS if ax is axes5[0] else False),
        cbar=False,
        mask=pivot.isna(),
    )
    add_divider(ax)
    ax.set_title(plabel, fontsize=5, pad=2)
    ax.set_xticklabels(SEEG_ALL_LABELS, fontsize=2.5, rotation=45, ha="right")
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(length=0)

axes5[0].set_yticklabels(axes5[0].get_yticklabels(), fontsize=5, rotation=0)

# single shared colorbar with ticks at meaningful p-value thresholds
cbar5 = fig5.colorbar(
    axes5[-1].collections[0], ax=axes5,
    orientation="vertical", shrink=0.8, pad=0.02,
    norm=P5_NORM,
)
cbar5.set_label("Permutation p-value", fontsize=5, labelpad=4)
cbar5.set_ticks([0, 0.01, 0.05, 0.10, 0.25])
cbar5.set_ticklabels(["0", "0.01", "0.05", "0.10", "0.25"])
cbar5.ax.tick_params(labelsize=4, length=2)
cbar5.outline.set_visible(False)

fig5.tight_layout(pad=0.3, w_pad=0.2)
save(fig5, "dbs_seeg_pca_per_patient_pvalue")

plt.show()
