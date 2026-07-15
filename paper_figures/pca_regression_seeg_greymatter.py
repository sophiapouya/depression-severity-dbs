"""
SEEG PCA regression LOO CV (greymatter filtered) — pooled scatter plot, patients in different colors.
PCA retains components explaining 90% of variance; OLS regression on PC scores.
Run from project root: python paper_figures/pca_regression_seeg_greymatter.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from config import FEATURES_SEEG

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
NON_FEATURE_COLS = ["patient_id", "catdi_score", "session_name", "time"]
VARIANCE_THRESHOLD = 0.90

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_regression_seeg")
os.makedirs(OUT_DIR, exist_ok=True)


# ── LOO CV per patient ────────────────────────────────────────────────────────
def run_patient(patient_df, patient_id):
    catdi_scores = patient_df["catdi_score"].reset_index(drop=True)
    feature_df = patient_df.drop(columns=[c for c in NON_FEATURE_COLS if c in patient_df.columns])
    feature_df = feature_df.reset_index(drop=True)

    # 70% presence rule
    cutoff = int(np.ceil(0.70 * feature_df.shape[0]))
    feature_df = feature_df[feature_df.columns[feature_df.notna().sum() >= cutoff]]

    measured_all, predicted_all = [], []

    for tr_idx, te_idx in LeaveOneOut().split(feature_df):
        data_tr = feature_df.iloc[tr_idx]
        data_te = feature_df.iloc[te_idx]
        y_tr = catdi_scores.iloc[tr_idx]
        y_te = catdi_scores.iloc[te_idx]

        # fill NaNs with training means
        train_mean = data_tr.mean()
        data_tr = data_tr.fillna(train_mean)
        data_te = data_te.fillna(train_mean)

        # standardize using training stats
        scaler = StandardScaler()
        scaler.fit(data_tr)
        X_tr = scaler.transform(data_tr)
        X_te = scaler.transform(data_te)

        # PCA: retain components explaining 90% variance, fit on train only
        pca = PCA(VARIANCE_THRESHOLD)
        pca.fit(X_tr)
        pc_tr = pca.transform(X_tr)
        pc_te = pca.transform(X_te)

        # OLS regression
        model = LinearRegression()
        model.fit(pc_tr, y_tr)
        predicted_all.extend(model.predict(pc_te))
        measured_all.extend(y_te)

    print(f"{patient_id}: {len(measured_all)} sessions")
    return np.array(measured_all, dtype=float), np.array(predicted_all, dtype=float)


# ── run all patients ──────────────────────────────────────────────────────────
features_path = str(FEATURES_SEEG).replace(".csv", "_greymatter.csv")
all_df = pd.read_csv(features_path)
all_df = all_df.drop(columns=["time"], errors="ignore")

results = []
for patient in all_df["patient_id"].unique():
    patient_df = all_df[all_df["patient_id"] == patient].copy()
    measured, predicted = run_patient(patient_df, patient)
    results.append({"patient_id": patient, "measured": measured, "predicted": predicted})

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(33 * MM, 33 * MM))
colors = plt.cm.tab10(np.linspace(0, 1, len(results)))

all_vals = np.concatenate([np.concatenate([r["measured"], r["predicted"]]) for r in results])
lim = (float(np.floor(all_vals.min() / 5) * 5), 100)

for color, r in zip(colors, results):
    ax.scatter(r["measured"], r["predicted"],
               color=color, s=5, alpha=0.9, linewidths=0, clip_on=False)

ax.plot(lim, lim, linestyle="--", color="black", linewidth=0.75)
ax.set_xlim(*lim)
ax.set_ylim(*lim)
ticks = np.arange(0, 101, 20)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Measured Depression Severity", fontsize=5)
ax.set_ylabel("Predicted Depression Severity", fontsize=5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out_path = os.path.join(OUT_DIR, "pca_regression_seeg_all_patients_greymatter.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
svg_path = os.path.join(OUT_DIR, "pca_regression_seeg_all_patients_greymatter.svg")
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path}")
plt.show()
