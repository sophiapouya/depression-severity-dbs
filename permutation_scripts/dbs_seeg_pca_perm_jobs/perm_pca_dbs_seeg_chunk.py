"""
DBS + SEEG combined PCA permutation test — chunk-based for SLURM array jobs.

Merges DBS features (filtered by --dbs_region) with SEEG greymatter features.
Region/hemisphere mapping comes from each patient's metadata_greymatter.csv.

Two SEEG selection modes:

  Single mode  (--seeg_right_region not supplied):
    --seeg_region     : all | acc | amy | dlpfc | vmpfc | ofc
    --seeg_hemisphere : all | left | right
    condition label   : dbs_{DBS}_seeg_{region}_{hemisphere}

  Pair mode  (--seeg_right_region supplied):
    Combines LEFT channels from --seeg_region with RIGHT channels from
    --seeg_right_region. --seeg_hemisphere is ignored.
    --seeg_region       : acc | amy | dlpfc | vmpfc | ofc  (left side)
    --seeg_right_region : acc | amy | dlpfc | vmpfc | ofc  (right side)
    condition label     : dbs_{DBS}_seeg_L{left}_R{right}

--dbs_region choices: ALL | VCVS | LVCVS | RVCVS | SCC | LSCC | RSCC

Output:
  {BASE_DIR_SEEG}/dbs_seeg_pca_permutation_test/{condition_label}/chunks/
    {patient_id}_true_decoder.csv   (chunk 0 only)
    {patient_id}_chunk{id:03d}.csv

Run from project root:
  # single mode
  python dbs_seeg_pca_perm_jobs/perm_pca_dbs_seeg_chunk.py \\
      --patient_id DBSTRD001 --dbs_region RVCVS \\
      --seeg_region acc --seeg_hemisphere left \\
      --chunk_id 0 --n_chunks 20

  # pair mode
  python dbs_seeg_pca_perm_jobs/perm_pca_dbs_seeg_chunk.py \\
      --patient_id DBSTRD001 --dbs_region RVCVS \\
      --seeg_region acc --seeg_right_region amy \\
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

from config import BASE_DIR_DBS, BASE_DIR_SEEG, FEATURES_DBS, FEATURES_SEEG_GM
from src.postprocessing_functions import get_nrmse

DBS_REGION_OPTIONS  = ["ALL", "VCVS", "LVCVS", "RVCVS", "SCC", "LSCC", "RSCC"]
SEEG_REGION_OPTIONS = ["all", "acc", "amy", "dlpfc", "vmpfc", "ofc"]
HEMISPHERE_OPTIONS  = ["all", "left", "right"]
ID_COLUMNS          = ["patient_id", "session_name", "catdi_score"]

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id",       type=str, required=True)
parser.add_argument("--dbs_region",       type=str, required=True, choices=DBS_REGION_OPTIONS)
parser.add_argument("--seeg_region",      type=str, required=True, choices=SEEG_REGION_OPTIONS)
parser.add_argument("--seeg_hemisphere",  type=str, default="all",  choices=HEMISPHERE_OPTIONS)
parser.add_argument("--seeg_right_region",type=str, default=None,
                    choices=["acc", "amy", "dlpfc", "vmpfc", "ofc"])
parser.add_argument("--chunk_id",         type=int, required=True)
parser.add_argument("--n_chunks",         type=int, default=20)
parser.add_argument("--n_permutations",   type=int, default=1000)
parser.add_argument("--n_jobs",           type=int, default=8)
args = parser.parse_args()

patient_id        = args.patient_id
dbs_region        = args.dbs_region
seeg_region       = args.seeg_region
seeg_hemisphere   = args.seeg_hemisphere
seeg_right_region = args.seeg_right_region
chunk_id          = args.chunk_id
n_chunks          = args.n_chunks
n_jobs            = args.n_jobs

pair_mode = seeg_right_region is not None
if pair_mode:
    condition_label = f"dbs_{dbs_region}_seeg_L{seeg_region}_R{seeg_right_region}"
else:
    condition_label = f"dbs_{dbs_region}_seeg_{seeg_region}_{seeg_hemisphere}"

chunk_size   = int(np.ceil(args.n_permutations / n_chunks))
chunk_start  = chunk_id * chunk_size
chunk_end    = min(chunk_start + chunk_size, args.n_permutations)
perm_indices = list(range(chunk_start, chunk_end))

print(f"Patient={patient_id}, condition={condition_label}, chunk={chunk_id}/{n_chunks}, "
      f"perms {chunk_start}-{chunk_end-1}, n_jobs={n_jobs}")

VARIANCE_THRESHOLD = 0.90
RANDOM_SEED        = 42

base_dir   = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "dbs_seeg_pca_permutation_test",
                          condition_label, "chunks")
os.makedirs(chunks_dir, exist_ok=True)


# ── DBS column filtering ───────────────────────────────────────────────────────
def subset_dbs_columns(df: pd.DataFrame, region: str) -> pd.DataFrame:
    df = df.drop(columns=["time"], errors="ignore").copy()
    if region == "ALL":
        return df
    if region == "SCC":
        return df.drop(columns=[c for c in df.columns if "VCVS" in c])
    if region == "VCVS":
        return df.drop(columns=[c for c in df.columns if "SCC" in c])
    keep = [c for c in df.columns if region in c] + ID_COLUMNS
    return df[[c for c in keep if c in df.columns]]


# ── load metadata ──────────────────────────────────────────────────────────────
meta_path = os.path.join(base_dir, patient_id,
                         f"{patient_id}_metadata_greymatter.csv")
if not os.path.exists(meta_path):
    raise FileNotFoundError(f"Metadata not found: {meta_path}")

meta_df = pd.read_csv(meta_path)
meta_df["channel_name"]   = meta_df["channel_name"].str.strip()
meta_df["channel_region"] = meta_df["channel_region"].str.strip().str.lower()


# ── build allowed channel set ─────────────────────────────────────────────────
if pair_mode:
    # left channels from seeg_region + right channels from seeg_right_region
    allowed_channels = set(
        row["channel_name"] for _, row in meta_df.iterrows()
        if row["channel_region"] == seeg_region.lower()
        and row["channel_name"][0].upper() == "L"
    ) | set(
        row["channel_name"] for _, row in meta_df.iterrows()
        if row["channel_region"] == seeg_right_region.lower()
        and row["channel_name"][0].upper() == "R"
    )
    if not allowed_channels:
        raise ValueError(f"No channels found for pair L{seeg_region}+R{seeg_right_region}, "
                         f"patient {patient_id}")
else:
    def channel_passes(ch_name, ch_region):
        region_ok = (seeg_region == "all") or (ch_region == seeg_region.lower())
        hemi_ok   = (seeg_hemisphere == "all") or \
                    (seeg_hemisphere == "left"  and ch_name[0].upper() == "L") or \
                    (seeg_hemisphere == "right" and ch_name[0].upper() == "R")
        return region_ok and hemi_ok

    allowed_channels = set(
        row["channel_name"] for _, row in meta_df.iterrows()
        if channel_passes(row["channel_name"], row["channel_region"])
    )
    if not allowed_channels:
        raise ValueError(f"No SEEG channels pass filter seeg_region={seeg_region}, "
                         f"seeg_hemisphere={seeg_hemisphere} for patient {patient_id}")

print(f"SEEG channels passing filter: {len(allowed_channels)}")


# ── load and merge DBS + SEEG data ────────────────────────────────────────────
dbs_df  = pd.read_csv(str(FEATURES_DBS))
seeg_df = pd.read_csv(str(FEATURES_SEEG_GM))

dbs_patient = subset_dbs_columns(dbs_df, dbs_region)
dbs_patient = dbs_patient[dbs_patient["patient_id"] == patient_id].copy()

seeg_patient   = seeg_df[seeg_df["patient_id"] == patient_id].copy()
seeg_feat_cols = [
    c for c in seeg_patient.columns
    if c not in ID_COLUMNS and any(c.startswith(ch + "_") for ch in allowed_channels)
]
if not seeg_feat_cols:
    raise ValueError(f"No SEEG feature columns matched for patient {patient_id}, "
                     f"condition {condition_label}")
seeg_patient = seeg_patient[ID_COLUMNS + seeg_feat_cols]

print(f"DBS feature cols: {len([c for c in dbs_patient.columns if c not in ID_COLUMNS])}, "
      f"SEEG feature cols: {len(seeg_feat_cols)}")

patient_df_raw = pd.merge(dbs_patient, seeg_patient, how="inner", on=ID_COLUMNS)
if patient_df_raw.empty:
    raise ValueError(f"No overlapping sessions after merging DBS and SEEG for {patient_id}")
print(f"Merged sessions: {len(patient_df_raw)}")


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
        "patient_id":       patient_id,
        "dbs_region":       dbs_region,
        "seeg_condition":   condition_label,
        "true_nrmse":       true_nrmse,
        "r_val":            true_r,
        "pearson_p":        true_pearson_p,
        "mse":              true_mse,
        "r_squared":        true_r2,
        "n_sessions":       len(patient_df_raw),
        "n_seeg_channels":  len(allowed_channels),
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
    [{"patient_id": patient_id, "seeg_condition": condition_label,
      "iteration": idx, "nrmse": nrmse}
     for idx, nrmse in results]
).to_csv(os.path.join(chunks_dir, f"{patient_id}_chunk{chunk_id:03d}.csv"), index=False)
print(f"Saved chunk {chunk_id}")
