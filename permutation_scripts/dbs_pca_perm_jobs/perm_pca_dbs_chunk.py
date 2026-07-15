"""
DBS PCA permutation test — chunk-based for SLURM array jobs.

1000 permutations split into N_CHUNKS. Each chunk runs its slice in parallel
via joblib. Chunk 0 also runs and saves the true decoder result.

Probe types: ALL, LEFT, RIGHT, LSCC, RSCC, RVCVS, LVCVS, SCC, VCVS

Output goes to:
  permutation_results/pca_dbs/{probe_type}/
    {patient_id}_true_decoder.csv   (chunk 0 only)
    {patient_id}_chunk{chunk_id:03d}.csv

Run from project root:
  python dbs_pca_perm_jobs/perm_pca_dbs_chunk.py \
      --patient_id DBSTRD001 --probe_type ALL --chunk_id 0 --n_chunks 20
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from config import FEATURES_DBS
from src.postprocessing_functions import get_nrmse

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id",  type=str, required=True)
parser.add_argument("--probe_type",  type=str, required=True,
                    choices=["ALL", "LEFT", "RIGHT", "LSCC", "RSCC", "RVCVS", "LVCVS", "SCC", "VCVS"])
parser.add_argument("--chunk_id",    type=int, required=True)
parser.add_argument("--n_chunks",    type=int, default=20)
parser.add_argument("--n_permutations", type=int, default=1000)
parser.add_argument("--n_jobs",      type=int, default=8)
args = parser.parse_args()

patient_id   = args.patient_id
probe_type   = args.probe_type
chunk_id     = args.chunk_id
n_chunks     = args.n_chunks
n_permutations = args.n_permutations
n_jobs       = args.n_jobs

chunk_size  = int(np.ceil(n_permutations / n_chunks))
chunk_start = chunk_id * chunk_size
chunk_end   = min(chunk_start + chunk_size, n_permutations)
perm_indices = list(range(chunk_start, chunk_end))

print(f"Patient={patient_id}, probe={probe_type}, chunk={chunk_id}/{n_chunks}, "
      f"perms {chunk_start}-{chunk_end-1}, n_jobs={n_jobs}")

VARIANCE_THRESHOLD = 0.90
ID_COLUMNS  = ["patient_id", "session_name", "catdi_score"]
RANDOM_SEED = 42

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out_dir = os.path.join(PROJECT_ROOT, "permutation_results", "pca_dbs", probe_type)
os.makedirs(out_dir, exist_ok=True)

# ── load and filter data ───────────────────────────────────────────────────────
all_df = pd.read_csv(str(FEATURES_DBS))
all_df = all_df.drop(columns=["time"], errors="ignore")

patient_df_raw = all_df[all_df["patient_id"] == patient_id].copy()
if patient_df_raw.empty:
    raise ValueError(f"No data found for patient {patient_id}")

# Apply probe filtering (same logic as pca_regression.py)
if probe_type == "ALL":
    pass
elif probe_type == "LEFT":
    drop = [c for c in patient_df_raw.columns if ("RS" in c) or ("RV" in c)]
    patient_df_raw = patient_df_raw.drop(columns=drop)
elif probe_type == "RIGHT":
    drop = [c for c in patient_df_raw.columns if ("LS" in c) or ("LV" in c)]
    patient_df_raw = patient_df_raw.drop(columns=drop)
elif probe_type == "SCC":
    drop = [c for c in patient_df_raw.columns if "VCVS" in c]
    patient_df_raw = patient_df_raw.drop(columns=drop)
elif probe_type == "VCVS":
    drop = [c for c in patient_df_raw.columns if "SCC" in c]
    patient_df_raw = patient_df_raw.drop(columns=drop)
else:
    # LSCC, RSCC, RVCVS, LVCVS
    keep = [c for c in patient_df_raw.columns if probe_type in c] + ID_COLUMNS
    patient_df_raw = patient_df_raw[[c for c in keep if c in patient_df_raw.columns]]

print(f"  Feature columns after probe filter: "
      f"{len([c for c in patient_df_raw.columns if c not in ID_COLUMNS])}")

# ── core decoder ───────────────────────────────────────────────────────────────
def run_pca_decoder_once(df, catdi_override=None):
    df = df.copy()
    if catdi_override is not None:
        df["catdi_score"] = catdi_override

    catdi_scores = df["catdi_score"].copy()
    feature_df   = df.drop(columns=[c for c in ID_COLUMNS if c in df.columns])

    cutoff = int(np.ceil(0.70 * feature_df.shape[0]))
    keep   = feature_df.columns[feature_df.notna().sum() >= cutoff]
    feature_df = feature_df[keep]

    if feature_df.shape[1] == 0:
        return np.array([]), np.array([]), np.inf

    measured_all, predicted_all = [], []

    for tr_idx, te_idx in LeaveOneOut().split(feature_df):
        data_tr = feature_df.iloc[tr_idx].copy()
        data_te = feature_df.iloc[te_idx].copy()
        y_tr    = catdi_scores.iloc[tr_idx].values
        y_te    = catdi_scores.iloc[te_idx].values

        train_mean = data_tr.mean()
        data_tr = data_tr.fillna(train_mean)
        data_te = data_te.fillna(train_mean)

        scaler = StandardScaler()
        scaler.fit(data_tr)
        X_tr = scaler.transform(data_tr)
        X_te = scaler.transform(data_te)

        pca = PCA(VARIANCE_THRESHOLD)
        pca.fit(X_tr)
        pc_tr = pca.transform(X_tr)
        pc_te = pca.transform(X_te)

        model = LinearRegression()
        model.fit(pc_tr, y_tr)
        predicted_all.extend(model.predict(pc_te))
        measured_all.extend(y_te)

    if not measured_all:
        return np.array([]), np.array([]), np.inf

    measured_all  = np.array(measured_all)
    predicted_all = np.array(predicted_all)
    nrmse = get_nrmse(measured_all, predicted_all)
    return measured_all, predicted_all, nrmse

# ── true decoder (chunk 0 only) ────────────────────────────────────────────────
if chunk_id == 0:
    print("Running true decoder...")
    true_measured, true_predicted, true_nrmse = run_pca_decoder_once(patient_df_raw)

    if len(true_measured) == 0:
        raise ValueError("True decoder produced no results")

    true_r, true_pearson_p = pearsonr(true_measured, true_predicted)
    true_mse = mean_squared_error(true_measured, true_predicted)
    true_r2  = r2_score(true_measured, true_predicted)
    print(f"True NRMSE={true_nrmse:.4f}, R={true_r:.4f}")

    pd.DataFrame([{
        "patient_id":  patient_id,
        "probe_type":  probe_type,
        "true_nrmse":  true_nrmse,
        "r_val":       true_r,
        "pearson_p":   true_pearson_p,
        "mse":         true_mse,
        "r_squared":   true_r2,
        "n_sessions":  len(patient_df_raw),
    }]).to_csv(os.path.join(out_dir, f"{patient_id}_true_decoder.csv"), index=False)

# ── permutations ───────────────────────────────────────────────────────────────
def run_one_permutation(perm_idx):
    rng = np.random.default_rng(RANDOM_SEED * 10000 + perm_idx)
    shuffled = rng.permutation(patient_df_raw["catdi_score"].to_numpy())
    _, _, perm_nrmse = run_pca_decoder_once(patient_df_raw, catdi_override=shuffled)
    if (perm_idx + 1) % 10 == 0:
        print(f"  Done permutation {perm_idx}")
    return perm_idx, perm_nrmse

print(f"Running {len(perm_indices)} permutations with {n_jobs} workers...")
results = Parallel(n_jobs=n_jobs, backend="loky")(
    delayed(run_one_permutation)(i) for i in perm_indices
)

pd.DataFrame([
    {"patient_id": patient_id, "probe_type": probe_type, "iteration": idx, "nrmse": nrmse}
    for idx, nrmse in results
]).to_csv(os.path.join(out_dir, f"{patient_id}_chunk{chunk_id:03d}.csv"), index=False)

print(f"Saved chunk {chunk_id} to {out_dir}/")
