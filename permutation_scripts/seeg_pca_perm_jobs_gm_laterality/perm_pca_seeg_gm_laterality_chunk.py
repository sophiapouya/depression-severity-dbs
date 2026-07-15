"""
SEEG PCA permutation test (greymatter filtered) with hemisphere support.

Identical to seeg_pca_perm_jobs_gm but adds --hemisphere (all | left | right).
Hemisphere filtering uses each patient's metadata_greymatter.csv:
  first letter of channel_name → L = left, R = right.
  hemisphere=all keeps every channel (same behaviour as seeg_pca_perm_jobs_gm).

Output:
  {BASE_DIR_SEEG}/seeg_pca_permutation_test_gm_laterality/{hemisphere}/chunks/
    {patient_id}_true_decoder.csv   (chunk 0 only)
    {patient_id}_chunk{id:03d}.csv

Run from project root:
  python seeg_pca_perm_jobs_gm_laterality/perm_pca_seeg_gm_laterality_chunk.py \\
      --patient_id DBSTRD001 --hemisphere left --chunk_id 0 --n_chunks 20
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

from config import BASE_DIR_SEEG, FEATURES_SEEG_GM
from src.postprocessing_functions import get_nrmse

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id",     type=str, required=True)
parser.add_argument("--hemisphere",     type=str, required=True,
                    choices=["all", "left", "right"])
parser.add_argument("--chunk_id",       type=int, required=True)
parser.add_argument("--n_chunks",       type=int, default=20)
parser.add_argument("--n_permutations", type=int, default=1000)
parser.add_argument("--n_jobs",         type=int, default=8)
args = parser.parse_args()

patient_id = args.patient_id
hemisphere = args.hemisphere
chunk_id   = args.chunk_id
n_chunks   = args.n_chunks
n_jobs     = args.n_jobs

chunk_size   = int(np.ceil(args.n_permutations / n_chunks))
chunk_start  = chunk_id * chunk_size
chunk_end    = min(chunk_start + chunk_size, args.n_permutations)
perm_indices = list(range(chunk_start, chunk_end))

print(f"Patient={patient_id}, hemisphere={hemisphere}, chunk={chunk_id}/{n_chunks}, "
      f"perms {chunk_start}-{chunk_end-1}, n_jobs={n_jobs}")

VARIANCE_THRESHOLD = 0.90
ID_COLUMNS = ["patient_id", "session_name", "catdi_score"]
RANDOM_SEED = 42

base_dir   = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "seeg_pca_permutation_test_gm_laterality",
                          hemisphere, "chunks")
os.makedirs(chunks_dir, exist_ok=True)

# ── load feature data ──────────────────────────────────────────────────────────
seeg_df        = pd.read_csv(str(FEATURES_SEEG_GM))
patient_df_raw = seeg_df[seeg_df["patient_id"] == patient_id].copy()
if patient_df_raw.empty:
    raise ValueError(f"No data found for patient {patient_id}")

# ── hemisphere filtering via metadata ─────────────────────────────────────────
if hemisphere != "all":
    meta_path = os.path.join(base_dir, patient_id,
                             f"{patient_id}_metadata_greymatter.csv")
    if not os.path.exists(meta_path):
        raise FileNotFoundError(f"Metadata not found: {meta_path}")

    meta_df = pd.read_csv(meta_path)
    hemi_letter = "L" if hemisphere == "left" else "R"
    allowed_channels = set(
        row["channel_name"]
        for _, row in meta_df.iterrows()
        if row["channel_name"][0].upper() == hemi_letter
    )
    feature_cols = [
        c for c in patient_df_raw.columns
        if c not in ID_COLUMNS and any(c.startswith(ch + "_") for ch in allowed_channels)
    ]
    if not feature_cols:
        raise ValueError(f"No feature columns for hemisphere={hemisphere}, "
                         f"patient={patient_id}")
    print(f"Channels: {len(allowed_channels)}, feature columns: {len(feature_cols)}")
    patient_df_raw = patient_df_raw[ID_COLUMNS + feature_cols]
else:
    n_feat = len([c for c in patient_df_raw.columns if c not in ID_COLUMNS])
    print(f"hemisphere=all: using all {n_feat} feature columns")


# ── decoder ────────────────────────────────────────────────────────────────────
def run_pca_decoder_once(patient_df_input, catdi_override=None):
    patient_df = patient_df_input.copy()
    if catdi_override is not None:
        patient_df["catdi_score"] = catdi_override

    catdi_scores = patient_df["catdi_score"].copy()
    feature_df   = patient_df.drop(columns=ID_COLUMNS)

    sessions   = feature_df.shape[0]
    cutoff     = int(np.ceil(0.70 * sessions))
    sums       = feature_df.notna().sum()
    feature_df = feature_df[sums[sums >= cutoff].index]

    if feature_df.shape[1] == 0:
        return np.array([]), np.array([]), np.inf

    measured_all, predicted_all = [], []

    for train_idx, test_idx in LeaveOneOut().split(feature_df):
        data_train = feature_df.iloc[train_idx].copy()
        data_test  = feature_df.iloc[test_idx].copy()
        y_train    = catdi_scores.iloc[train_idx].values
        y_test     = catdi_scores.iloc[test_idx].values

        train_mean = data_train.mean()
        data_train = data_train.fillna(train_mean)
        data_test  = data_test.fillna(train_mean)

        scaler = StandardScaler()
        scaler.fit(data_train)
        X_train = scaler.transform(data_train)
        X_test  = scaler.transform(data_test)

        pca = PCA(VARIANCE_THRESHOLD)
        pca.fit(X_train)
        pc_train = pca.transform(X_train)
        pc_test  = pca.transform(X_test)

        model = LinearRegression()
        model.fit(pc_train, y_train)
        predicted_all.extend(model.predict(pc_test))
        measured_all.extend(y_test)

    if not measured_all:
        return np.array([]), np.array([]), np.inf

    measured_all  = np.array(measured_all)
    predicted_all = np.array(predicted_all)
    return measured_all, predicted_all, get_nrmse(measured_all, predicted_all)


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
        "patient_id": patient_id,
        "hemisphere": hemisphere,
        "true_nrmse": true_nrmse,
        "r_val":      true_r,
        "pearson_p":  true_pearson_p,
        "mse":        true_mse,
        "r_squared":  true_r2,
        "n_sessions": len(patient_df_raw),
    }]).to_csv(os.path.join(chunks_dir, f"{patient_id}_true_decoder.csv"), index=False)


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

pd.DataFrame(
    [{"patient_id": patient_id, "hemisphere": hemisphere,
      "iteration": idx, "nrmse": nrmse}
     for idx, nrmse in results]
).to_csv(os.path.join(chunks_dir, f"{patient_id}_chunk{chunk_id:03d}.csv"), index=False)
print(f"Saved chunk {chunk_id}")
