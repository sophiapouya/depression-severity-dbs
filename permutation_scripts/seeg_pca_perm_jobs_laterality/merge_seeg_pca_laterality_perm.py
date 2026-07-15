"""
Merge chunk outputs into final permutation result files (SEEG PCA, laterality).

Expects chunks at:
  {BASE_DIR_SEEG}/seeg_pca_permutation_test_laterality/{region}_{hemisphere}/chunks/
    {patient_id}_true_decoder.csv
    {patient_id}_chunk000.csv ... chunk{n_chunks-1}.csv

Writes to:
  {BASE_DIR_SEEG}/seeg_pca_permutation_test_laterality/{region}_{hemisphere}/
    {patient_id}_OLS_LOO_seeg_{region}_{hemisphere}_permutation_nrmse.csv
    {patient_id}_OLS_LOO_seeg_{region}_{hemisphere}_permutation_summary.csv

Run from project root:
  python seeg_pca_perm_jobs_laterality/merge_seeg_pca_laterality_perm.py \\
      --patient_id DBSTRD001 --region acc --hemisphere left --n_chunks 20
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from config import BASE_DIR_SEEG

REGION_CHOICES    = ["all", "acc", "amy", "dlpfc", "vmpfc", "ofc"]
HEMISPHERE_CHOICES = ["all", "left", "right"]

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id", type=str, required=True)
parser.add_argument("--region",     type=str, required=True, choices=REGION_CHOICES)
parser.add_argument("--hemisphere", type=str, required=True, choices=HEMISPHERE_CHOICES)
parser.add_argument("--n_chunks",   type=int, default=20)
args = parser.parse_args()

patient_id      = args.patient_id
region          = args.region
hemisphere      = args.hemisphere
n_chunks        = args.n_chunks
condition_label = f"{region}_{hemisphere}"

base_dir   = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "seeg_pca_permutation_test_laterality",
                          condition_label, "chunks")
out_dir    = os.path.join(base_dir, "seeg_pca_permutation_test_laterality",
                          condition_label)
os.makedirs(out_dir, exist_ok=True)

# ── true decoder ───────────────────────────────────────────────────────────────
true_csv = os.path.join(chunks_dir, f"{patient_id}_true_decoder.csv")
if not os.path.exists(true_csv):
    raise FileNotFoundError(f"True decoder file not found: {true_csv}")
true_df   = pd.read_csv(true_csv)
true_nrmse = float(true_df["true_nrmse"].iloc[0])

# ── load and merge chunks ──────────────────────────────────────────────────────
chunk_dfs, missing = [], []
for i in range(n_chunks):
    chunk_csv = os.path.join(chunks_dir, f"{patient_id}_chunk{i:03d}.csv")
    if os.path.exists(chunk_csv):
        chunk_dfs.append(pd.read_csv(chunk_csv))
    else:
        missing.append(i)

if missing:
    print(f"WARNING: missing chunks {missing} — results will be incomplete.")

all_chunks   = pd.concat(chunk_dfs, ignore_index=True)
perm_arr     = all_chunks["nrmse"].values
n_perms      = len(perm_arr)
p_value      = np.mean(perm_arr <= true_nrmse)

print(f"{patient_id} [{condition_label}]: {n_perms} permutations, "
      f"true NRMSE={true_nrmse:.4f}, p={p_value:.4f}")

# ── nrmse CSV ──────────────────────────────────────────────────────────────────
rows = [
    {"patient_id": patient_id, "region": region, "hemisphere": hemisphere,
     "type": "permuted", "iteration": int(r["iteration"]), "nrmse": r["nrmse"]}
    for _, r in all_chunks.iterrows()
]
rows.append({"patient_id": patient_id, "region": region, "hemisphere": hemisphere,
             "type": "true", "iteration": np.nan, "nrmse": true_nrmse})
nrmse_csv = os.path.join(out_dir,
    f"{patient_id}_OLS_LOO_seeg_{condition_label}_permutation_nrmse.csv")
pd.DataFrame(rows).to_csv(nrmse_csv, index=False)

# ── summary CSV ────────────────────────────────────────────────────────────────
summary_csv = os.path.join(out_dir,
    f"{patient_id}_OLS_LOO_seeg_{condition_label}_permutation_summary.csv")
pd.DataFrame([{
    "patient_id":    patient_id,
    "region":        region,
    "hemisphere":    hemisphere,
    "model_choice":  "OLS",
    "cv_choice":     "LOO",
    "r_val":         float(true_df["r_val"].iloc[0]),
    "pearson_p":     float(true_df["pearson_p"].iloc[0]),
    "mse":           float(true_df["mse"].iloc[0]),
    "r_squared":     float(true_df["r_squared"].iloc[0]),
    "true_nrmse":    true_nrmse,
    "perm_p_value":  p_value,
    "n_permutations":n_perms,
    "n_sessions":    int(true_df["n_sessions"].iloc[0]),
    "n_channels":    int(true_df["n_channels"].iloc[0]),
    "n_features":    int(true_df["n_features"].iloc[0]),
}]).to_csv(summary_csv, index=False)

print(f"Saved: {nrmse_csv}")
print(f"Saved: {summary_csv}")
