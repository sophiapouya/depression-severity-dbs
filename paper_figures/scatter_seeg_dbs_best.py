"""
Scatter plots (measured vs predicted CATDI) for DBS + patient-specific best
2-SEEG-region combo, for two DBS conditions: RVCVS and SCC.

Best combo = the L+R SEEG region pair with lowest true NRMSE for each patient.
Re-runs the true PCA-OLS LOO-CV decoder to get per-session predictions.
Patients colored with tab10 anchored to the 7-patient list (DBSTRD011 excluded).

Run from project root: python paper_figures/scatter_seeg_dbs_best.py
"""
import glob
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

from config import FEATURES_DBS, FEATURES_SEEG_GM

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
ID_COLS    = ["patient_id", "session_name", "catdi_score"]
VARIANCE_THRESHOLD = 0.90

BASE     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEEG_DBS = os.path.join(BASE, "permutation_results", "dbs_seeg_pca_permutation_test")

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "scatter_seeg_dbs_best")
os.makedirs(OUT_DIR, exist_ok=True)

_tab10 = plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))
patient_colors = {p: _tab10[ALL_SUBJS.index(p)] for p in ALL_SUBJS}

# ── helpers ───────────────────────────────────────────────────────────────────
def subset_dbs_columns(df, region):
    df = df.drop(columns=["time"], errors="ignore").copy()
    if region == "ALL":
        return df
    if region == "SCC":
        return df.drop(columns=[c for c in df.columns if "VCVS" in c])
    if region == "VCVS":
        return df.drop(columns=[c for c in df.columns if "SCC" in c])
    keep = [c for c in df.columns if region in c] + ID_COLS
    return df[[c for c in keep if c in df.columns]]


def get_best_combo(prefix, patient_id):
    """Return (left_region, right_region) with lowest true NRMSE for this patient."""
    best_nrmse = float("inf")
    best_left = best_right = None
    for d in os.listdir(SEEG_DBS):
        if not d.startswith(f"dbs_{prefix}_seeg_"):
            continue
        suffix = d.split(f"dbs_{prefix}_seeg_")[1]
        parts = suffix.split("_")
        if len(parts) != 2 or parts[0] not in L_REGIONS or parts[1] not in R_REGIONS:
            continue
        csv = os.path.join(SEEG_DBS, d, "chunks", f"{patient_id}_true_decoder.csv")
        if not os.path.exists(csv):
            continue
        v = float(pd.read_csv(csv)["true_nrmse"].iloc[0])
        if v < best_nrmse:
            best_nrmse = v
            best_left  = parts[0][1:]  # strip L
            best_right = parts[1][1:]  # strip R
    return best_left, best_right


def build_allowed_channels(meta_df, left_region, right_region):
    meta_df = meta_df.copy()
    meta_df["channel_name"]   = meta_df["channel_name"].str.strip()
    meta_df["channel_region"] = meta_df["channel_region"].str.strip().str.lower()
    return set(
        row["channel_name"] for _, row in meta_df.iterrows()
        if row["channel_region"] == left_region.lower()
        and row["channel_name"][0].upper() == "L"
    ) | set(
        row["channel_name"] for _, row in meta_df.iterrows()
        if row["channel_region"] == right_region.lower()
        and row["channel_name"][0].upper() == "R"
    )


def run_true_decoder(patient_df_raw):
    catdi_scores = patient_df_raw["catdi_score"].copy()
    feature_df   = patient_df_raw.drop(columns=ID_COLS)

    cutoff     = int(np.ceil(0.70 * feature_df.shape[0]))
    feature_df = feature_df[feature_df.columns[feature_df.notna().sum() >= cutoff]]

    measured_all, predicted_all = [], []
    for tr_idx, te_idx in LeaveOneOut().split(feature_df):
        data_tr = feature_df.iloc[tr_idx].copy()
        data_te = feature_df.iloc[te_idx].copy()
        y_tr    = catdi_scores.iloc[tr_idx].values

        train_mean = data_tr.mean()
        data_tr    = data_tr.fillna(train_mean)
        data_te    = data_te.fillna(train_mean)

        scaler = StandardScaler()
        X_tr   = scaler.fit_transform(data_tr)
        X_te   = scaler.transform(data_te)

        pca    = PCA(VARIANCE_THRESHOLD)
        pc_tr  = pca.fit_transform(X_tr)
        pc_te  = pca.transform(X_te)

        model  = LinearRegression()
        model.fit(pc_tr, y_tr)
        predicted_all.extend(model.predict(pc_te))
        measured_all.extend(catdi_scores.iloc[te_idx].values)

    return np.array(measured_all, dtype=float), np.array(predicted_all, dtype=float)


# ── load raw feature files once ────────────────────────────────────────────────
dbs_df  = pd.read_csv(str(FEATURES_DBS))
seeg_df = pd.read_csv(str(FEATURES_SEEG_GM))


def run_condition(dbs_prefix, title, fname):
    results = []
    for pid in PLOT_SUBJS:
        left_r, right_r = get_best_combo(dbs_prefix, pid)
        if left_r is None:
            print(f"WARNING: no best combo found for {pid} / {dbs_prefix}")
            continue

        meta_path = os.path.join(str(BASE), "..", "neuralData", "seegData", pid,
                                 f"{pid}_metadata_greymatter.csv")
        # try multiple possible paths
        for candidate in [
            os.path.join(os.path.dirname(str(FEATURES_SEEG_GM)), pid,
                         f"{pid}_metadata_greymatter.csv"),
        ]:
            if os.path.exists(candidate):
                meta_path = candidate
                break

        if not os.path.exists(meta_path):
            print(f"WARNING: metadata not found for {pid}: {meta_path}")
            continue

        meta_df  = pd.read_csv(meta_path)
        channels = build_allowed_channels(meta_df, left_r, right_r)

        dbs_pat  = subset_dbs_columns(dbs_df, dbs_prefix)
        dbs_pat  = dbs_pat[dbs_pat["patient_id"] == pid].copy()

        seeg_pat = seeg_df[seeg_df["patient_id"] == pid].copy()
        seeg_cols = [c for c in seeg_pat.columns
                     if c not in ID_COLS and any(c.startswith(ch + "_") for ch in channels)]
        seeg_pat = seeg_pat[ID_COLS + seeg_cols]

        merged = pd.merge(dbs_pat, seeg_pat, how="inner", on=ID_COLS)
        if merged.empty:
            print(f"WARNING: empty merge for {pid} / {dbs_prefix}")
            continue

        measured, predicted = run_true_decoder(merged)
        r, _ = pearsonr(measured, predicted)
        print(f"{pid} ({dbs_prefix}+L{left_r}+R{right_r}): "
              f"{len(measured)} sessions, r={r:.3f}")
        results.append({"patient_id": pid, "measured": measured, "predicted": predicted})

    # ── plot ──────────────────────────────────────────────────────────────────
    MM = 1 / 25.4
    fig, ax = plt.subplots(figsize=(33 * MM, 33 * MM))

    all_vals = np.concatenate([np.concatenate([r["measured"], r["predicted"]])
                                for r in results])
    lim = (float(np.floor(all_vals.min() / 5) * 5) - 2, 100)

    for r in results:
        color = patient_colors.get(r["patient_id"], "gray")
        ax.scatter(r["measured"], r["predicted"],
                   color=color, s=4, alpha=0.9, linewidths=0, clip_on=False)

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

    fig.tight_layout(pad=0.3)
    for ext in ("png", "svg"):
        out = os.path.join(OUT_DIR, f"{fname}.{ext}")
        fig.savefig(out, dpi=300, bbox_inches="tight",
                    format=ext if ext == "svg" else None)
        print(f"Saved: {out}")
    plt.close(fig)


run_condition("RVCVS", "R-VCVS + best SEEG combo", "scatter_rvcvs_best")
run_condition("SCC",   "SCC + best SEEG combo",    "scatter_scc_best")
plt.show()
