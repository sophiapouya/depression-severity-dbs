"""
SEEG Lasso LOO CV (greymatter filtered) — pooled scatter plot, patients in different colors.
Region selected most often in the inner loop is shown per patient in the legend.
Run from project root: python paper_figures/lasso_seeg_greymatter.py
"""
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from sklearn.linear_model import Lasso
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from config import BASE_DIR_SEEG, CATDI_SCORES
from src.postprocessing_functions import get_nrmse

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
BASE_DIR = str(BASE_DIR_SEEG)
CATDI_FILE = str(CATDI_SCORES)
ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
REGIONS = ["acc", "amy", "ofc", "vmpfc", "dlpfc"]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
L1_REGS = np.around(np.arange(0.1, 1.1, 0.1), 1)
OUTLIER_SD = 4

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "lasso_seeg")
os.makedirs(OUT_DIR, exist_ok=True)


# ── LOO CV per patient ────────────────────────────────────────────────────────
def run_patient(subj):
    power_df = pd.read_csv(os.path.join(
        BASE_DIR, subj, "bipolar_channels_greymatter", "power_bipolar_greymatter",
        f"{subj}_bipolar_power_greymatter.csv",
    ))

    catdi = pd.read_excel(CATDI_FILE, sheet_name=subj)
    if subj in ["DBSTRD011", "DBSTRD014"]:
        catdi["Name"] = catdi["Name"].astype(str).str.split("task-").str[-1]
    scores_map = dict(zip(catdi["Name"].astype(str), catdi["Result"]))
    power_df["scores"] = power_df["session"].astype(str).map(scores_map)
    power_df = power_df.dropna(subset=["scores"])

    meta = pd.read_csv(os.path.join(BASE_DIR, subj, f"{subj}_metadata_greymatter.csv"))
    meta["ch_name"] = meta["channel_name"]
    power_df = pd.merge(power_df, meta[["ch_name", "channel_region"]], how="left", on="ch_name")

    unique_sessions = power_df["session"].unique()
    measured_all, predicted_all, region_selections = [], [], []

    for tr_idx, te_idx in LeaveOneOut().split(unique_sessions):
        tr_sess = unique_sessions[tr_idx]
        te_sess = unique_sessions[te_idx]
        data_tr = power_df[power_df["session"].isin(tr_sess)]
        data_te = power_df[power_df["session"].isin(te_sess)]

        # inner LOO: select best region and alpha
        region_scores, region_alphas = {}, {}
        for region in REGIONS:
            reg_df = data_tr[data_tr["channel_region"] == region]
            if reg_df.empty:
                continue
            pivot = reg_df.pivot_table(index="session", columns="ch_name", values=BANDS)
            reg_scores = reg_df.groupby("session")["scores"].first().loc[pivot.index]

            alpha_nrmse = []
            for l1 in L1_REGS:
                inner_pred, inner_true = [], []
                for ii_tr, ii_te in LeaveOneOut().split(pivot):
                    X_tr = pivot.iloc[ii_tr].values
                    X_te = pivot.iloc[ii_te].values
                    y_tr = reg_scores.iloc[ii_tr].values
                    y_te = reg_scores.iloc[ii_te].values

                    mu = np.mean(X_tr, axis=0)
                    sd = np.std(X_tr, axis=0)
                    sd[sd == 0] = 1
                    X_tr = np.where(np.abs(X_tr - mu) > OUTLIER_SD * sd, mu, X_tr)
                    X_te = np.where(np.abs(X_te - mu) > OUTLIER_SD * sd, mu, X_te)

                    m = make_pipeline(StandardScaler(), Lasso(alpha=l1, random_state=0, max_iter=2000))
                    m.fit(X_tr, y_tr)
                    inner_pred.extend(m.predict(X_te))
                    inner_true.extend(y_te)

                alpha_nrmse.append(get_nrmse(inner_true, inner_pred))

            best_idx = int(np.argmin(alpha_nrmse))
            region_scores[region] = alpha_nrmse[best_idx]
            region_alphas[region] = L1_REGS[best_idx]

        if not region_scores:
            continue
        best_region = min(region_scores, key=region_scores.get)
        best_alpha = region_alphas[best_region]
        region_selections.append(best_region)

        # outer fold
        tr_reg = data_tr[data_tr["channel_region"] == best_region]
        te_reg = data_te[data_te["channel_region"] == best_region]
        pivot_tr = tr_reg.pivot_table(index="session", columns="ch_name", values=BANDS)
        pivot_te = te_reg.pivot_table(
            index="session", columns="ch_name", values=BANDS
        ).reindex(columns=pivot_tr.columns)
        y_tr = tr_reg.groupby("session")["scores"].first().loc[pivot_tr.index]
        y_te = te_reg.groupby("session")["scores"].first().loc[pivot_te.index]

        X_tr = pivot_tr.values
        X_te = pivot_te.values
        mu = np.mean(X_tr, axis=0)
        sd = np.std(X_tr, axis=0)
        sd[sd == 0] = 1
        X_tr = np.where(np.abs(X_tr - mu) > OUTLIER_SD * sd, mu, X_tr)
        X_te = np.where(np.abs(X_te - mu) > OUTLIER_SD * sd, mu, X_te)

        model = make_pipeline(StandardScaler(), Lasso(alpha=best_alpha, random_state=42, max_iter=10000))
        model.fit(X_tr, y_tr)
        predicted_all.extend(model.predict(X_te))
        measured_all.extend(y_te)

    modal_region = Counter(region_selections).most_common(1)[0][0] if region_selections else "n/a"
    print(f"{subj}: {len(measured_all)} sessions, modal region = {modal_region}")
    return np.array(measured_all, dtype=float), np.array(predicted_all, dtype=float), modal_region


# ── run all patients ──────────────────────────────────────────────────────────
print("Running Lasso LOO CV (this will take a while)...")
results = []
for subj in ALL_SUBJS:
    measured, predicted, modal_region = run_patient(subj)
    results.append({"patient_id": subj, "measured": measured,
                    "predicted": predicted, "region": modal_region})

# ── plot ───────────────────────────────────────────────────────────────────────
MM = 1 / 25.4
fig, ax = plt.subplots(figsize=(33 * MM, 33 * MM))
colors = plt.cm.tab10(np.linspace(0, 1, len(results)))

all_vals = np.concatenate([np.concatenate([r["measured"], r["predicted"]]) for r in results])
lim = (float(np.floor(all_vals.min() / 5) * 5), float(np.ceil(all_vals.max() / 5) * 5))

for color, r in zip(colors, results):
    ax.scatter(r["measured"], r["predicted"],
               color=color, s=5, alpha=0.9, linewidths=0)

ax.plot(lim, lim, linestyle="--", color="black", linewidth=0.75)
ax.set_xlim(*lim)
ax.set_ylim(*lim)
ticks = np.arange(int(lim[0]), int(lim[1]) + 1, 20)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_xlabel("Measured Depression Severity", fontsize=5)
ax.set_ylabel("Predicted Depression Severity", fontsize=5)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

out_path = os.path.join(OUT_DIR, "lasso_seeg_all_patients_greymatter.png")
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved: {out_path}")
svg_path = os.path.join(OUT_DIR, "lasso_seeg_all_patients_greymatter.svg")
fig.savefig(svg_path, format="svg", bbox_inches="tight")
print(f"Saved: {svg_path}")
plt.show()
