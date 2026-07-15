"""
Horizontal bar plot: number of PCs explaining 90% of variance per patient.
Produces two separate figures — one for DBS, one for SEEG.
Run from project root: python paper_figures/pca_components_90pct.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from config import FEATURES_DBS, FEATURES_SEEG_GM

# ── style ─────────────────────────────────────────────────────────────────────
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
})

# ── constants ─────────────────────────────────────────────────────────────────
VARIANCE_THRESHOLD = 0.90
NON_FEATURE_COLS = {"patient_id", "session_name", "catdi_score", "time"}
MODES = {"dbs": str(FEATURES_DBS), "seeg": str(FEATURES_SEEG_GM)}

ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
PATIENT_COLORS = {
    subj: plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))[i]
    for i, subj in enumerate(ALL_SUBJS)
}

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_components")
os.makedirs(OUT_DIR, exist_ok=True)


def compute_pcs(csv_path):
    df = pd.read_csv(csv_path)
    results = []
    for patient in df["patient_id"].unique():
        patient_df = df[df["patient_id"] == patient].copy()
        feature_df = patient_df.drop(columns=[c for c in NON_FEATURE_COLS if c in patient_df.columns])

        cutoff = int(np.ceil(0.70 * feature_df.shape[0]))
        feature_df = feature_df[feature_df.columns[feature_df.notna().sum() >= cutoff]]
        feature_df = feature_df.fillna(feature_df.mean())

        X = StandardScaler().fit_transform(feature_df)
        pca = PCA()
        pca.fit(X)

        cumvar = np.cumsum(pca.explained_variance_ratio_)
        n_pcs = int(np.argmax(cumvar >= VARIANCE_THRESHOLD) + 1)
        results.append({"patient": patient, "n_pcs": n_pcs})

    return pd.DataFrame(results).sort_values("patient", ascending=False)


def make_plot(results_df, mode):
    fig, ax = plt.subplots(figsize=(6, len(results_df) * 0.55 + 0.8))

    display_labels = [p.replace("DBSTRD", "TRD") for p in results_df["patient"]]
    bar_colors = [PATIENT_COLORS.get(p, "steelblue") for p in results_df["patient"]]
    bars = ax.barh(
        display_labels,
        results_df["n_pcs"],
        color=bar_colors,
        height=0.6,
    )

    for bar, val in zip(bars, results_df["n_pcs"]):
        ax.text(
            bar.get_width() + 0.3, bar.get_y() + bar.get_height() / 2,
            str(int(val)),
            va="center", ha="left", fontsize=5,
        )

    ax.set_yticklabels(ax.get_yticklabels(), fontsize=6)
    ax.set_xlabel(f"PCs explaining {int(VARIANCE_THRESHOLD * 100)}% variance", fontsize=6)
    ax.set_xlim(0, results_df["n_pcs"].max() + 5)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.18)

    out_path = os.path.join(OUT_DIR, f"pca_components_90pct_{mode}.png")
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    print(f"Saved: {out_path}")
    svg_path = os.path.join(OUT_DIR, f"pca_components_90pct_{mode}.svg")
    fig.savefig(svg_path, format="svg", bbox_inches="tight")
    print(f"Saved: {svg_path}")
    plt.show()


for mode, csv_path in MODES.items():
    make_plot(compute_pcs(csv_path), mode)
