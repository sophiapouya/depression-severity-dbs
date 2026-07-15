"""
SEEG PCA permutation test with region + hemisphere support — chunk-based for SLURM array jobs.

Regions are determined from each patient's metadata_greymatter.csv
(channel_name → channel_region). Hemisphere is the first letter of the
channel_name (L = left, R = right).

--region choices : all | acc | amy | dlpfc | vmpfc | ofc
--hemisphere choices : all | left | right

Conditions run:
  acc/all, acc/left, acc/right
  amy/all, amy/left, amy/right
  dlpfc/all, dlpfc/left, dlpfc/right
  vmpfc/all, vmpfc/left, vmpfc/right
  ofc/all, ofc/left, ofc/right
  all/left, all/right

Output:
  {BASE_DIR_SEEG}/seeg_pca_permutation_test_laterality/{region}_{hemisphere}/chunks/
    {patient_id}_true_decoder.csv   (chunk 0 only)
    {patient_id}_chunk{id:03d}.csv

Run from project root:
  python seeg_pca_perm_jobs_laterality/perm_pca_seeg_laterality_chunk.py \\
      --patient_id DBSTRD001 --region acc --hemisphere left \\
      --chunk_id 0 --n_chunks 20
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

REGION_CHOICES = ["all", "acc", "amy", "dlpfc", "vmpfc", "ofc"]
HEMISPHERE_CHOICES = ["all", "left", "right"]

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id",    type=str, required=True)
parser.add_argument("--region",        type=str, required=True, choices=REGION_CHOICES)
parser.add_argument("--hemisphere",    type=str, required=True, choices=HEMISPHERE_CHOICES)
parser.add_argument("--chunk_id",      type=int, required=True)
parser.add_argument("--n_chunks",      type=int, default=20)
parser.add_argument("--n_permutations",type=int, default=1000)
parser.add_argument("--n_jobs",        type=int, default=8)
args = parser.parse_args()

patient_id  = args.patient_id
region      = args.region
hemisphere  = args.hemisphere
chunk_id    = args.chunk_id
n_chunks    = args.n_chunks
n_jobs      = args.n_jobs

condition_label = f"{region}_{hemisphere}"

chunk_size  = int(np.ceil(args.n_permutations / n_chunks))
chunk_start = chunk_id * chunk_size
chunk_end   = min(chunk_start + chunk_size, args.n_permutations)
perm_indices = list(range(chunk_start, chunk_end))

print(f"Patient={patient_id}, condition={condition_label}, chunk={chunk_id}/{n_chunks}, "
      f"perms {chunk_start}-{chunk_end-1}, n_jobs={n_jobs}")

VARIANCE_THRESHOLD = 0.90
ID_COLUMNS = ["patient_id", "session_name", "catdi_score"]
RANDOM_SEED = 42

base_dir   = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "seeg_pca_permutation_test_laterality",
                          condition_label, "chunks")
os.makedirs(chunks_dir, exist_ok=True)

# ── load feature data ──────────────────────────────────────────────────────────
seeg_df = pd.read_csv(str(FEATURES_SEEG_GM))
patient_df_raw = seeg_df[seeg_df["patient_id"] == patient_id].copy()
if patient_df_raw.empty:
    raise ValueError(f"No data found for patient {patient_id}")

# ── load metadata to get channel → region mapping ──────────────────────────────
meta_path = os.path.join(base_dir, patient_id,
                         f"{patient_id}_metadata_greymatter.csv")
if not os.path.exists(meta_path):
    raise FileNotFoundError(f"Metadata not found: {meta_path}")

meta_df = pd.read_csv(meta_path)
# channel_name e.g. "LVPF-OF_2", channel_region e.g. "ofc"

# build set of channel names that pass the region + hemisphere filter
def channel_passes(ch_name, ch_region):
    hemi_letter = ch_name[0].upper()  # 'L' or 'R'
    region_ok = (region == "all") or (ch_region.lower() == region.lower())
    hemi_ok   = (hemisphere == "all") or \
                (hemisphere == "left"  and hemi_letter == "L") or \
                (hemisphere == "right" and hemi_letter == "R")
    return region_ok and hemi_ok

allowed_channels = set(
    row["channel_name"]
    for _, row in meta_df.iterrows()
    if channel_passes(row["channel_name"], row["channel_region"])
)

if not allowed_channels:
    raise ValueError(f"No channels pass filter region={region}, hemisphere={hemisphere} "
                     f"for patient {patient_id}")

print(f"Channels passing filter: {len(allowed_channels)}")
print(f"Allowed channels: {sorted(allowed_channels)}")
non_id_cols = [c for c in patient_df_raw.columns if c not in ID_COLUMNS]
print(f"Total feature cols in CSV for patient: {len(non_id_cols)}")
print(f"Sample feature cols: {non_id_cols[:5]}")

# keep feature columns whose prefix matches an allowed channel name
# feature cols are named "{channel_name}_{feature_name}"
feature_cols = [
    c for c in patient_df_raw.columns
    if c not in ID_COLUMNS and any(c.startswith(ch + "_") for ch in allowed_channels)
]

if not feature_cols:
    raise ValueError(f"No feature columns found for condition {condition_label}, "
                     f"patient {patient_id}")

print(f"Feature columns kept: {len(feature_cols)}")
patient_df_raw = patient_df_raw[ID_COLUMNS + feature_cols]


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
        "patient_id":  patient_id,
        "region":      region,
        "hemisphere":  hemisphere,
        "true_nrmse":  true_nrmse,
        "r_val":       true_r,
        "pearson_p":   true_pearson_p,
        "mse":         true_mse,
        "r_squared":   true_r2,
        "n_sessions":  len(patient_df_raw),
        "n_channels":  len(allowed_channels),
        "n_features":  len(feature_cols),
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
    [{"patient_id": patient_id, "region": region, "hemisphere": hemisphere,
      "iteration": idx, "nrmse": nrmse}
     for idx, nrmse in results]
).to_csv(os.path.join(chunks_dir, f"{patient_id}_chunk{chunk_id:03d}.csv"), index=False)
print(f"Saved chunk {chunk_id}")
