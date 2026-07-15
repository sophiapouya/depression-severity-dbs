"""
Merge chunk outputs into final permutation result files (DBS + SEEG combined PCA).

Works for both single and pair SEEG modes — the condition label is passed directly.

Expects chunks at:
  {BASE_DIR_SEEG}/dbs_seeg_pca_permutation_test/{condition_label}/chunks/
    {patient_id}_true_decoder.csv
    {patient_id}_chunk000.csv ... chunk{n_chunks-1}.csv

Writes to:
  {BASE_DIR_SEEG}/dbs_seeg_pca_permutation_test/{condition_label}/
    {patient_id}_{condition_label}_permutation_nrmse.csv
    {patient_id}_{condition_label}_permutation_summary.csv

Single mode examples:
  python dbs_seeg_pca_perm_jobs/merge_pca_dbs_seeg_perm.py \\
      --patient_id DBSTRD001 --condition_label dbs_RVCVS_seeg_acc_left

Pair mode example:
  python dbs_seeg_pca_perm_jobs/merge_pca_dbs_seeg_perm.py \\
      --patient_id DBSTRD001 --condition_label dbs_RVCVS_seeg_Lacc_Ramy
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from config import BASE_DIR_SEEG

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id",      type=str, required=True)
parser.add_argument("--condition_label", type=str, required=True)
parser.add_argument("--n_chunks",        type=int, default=20)
args = parser.parse_args()

patient_id      = args.patient_id
condition_label = args.condition_label
n_chunks        = args.n_chunks

base_dir   = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "dbs_seeg_pca_permutation_test",
                          condition_label, "chunks")
out_dir    = os.path.join(base_dir, "dbs_seeg_pca_permutation_test", condition_label)
os.makedirs(out_dir, exist_ok=True)

# ── true decoder ───────────────────────────────────────────────────────────────
true_csv = os.path.join(chunks_dir, f"{patient_id}_true_decoder.csv")
if not os.path.exists(true_csv):
    raise FileNotFoundError(f"True decoder file not found: {true_csv}")
true_df    = pd.read_csv(true_csv)
true_nrmse = float(true_df["true_nrmse"].iloc[0])

# ── merge chunks ───────────────────────────────────────────────────────────────
chunk_dfs, missing = [], []
for i in range(n_chunks):
    chunk_csv = os.path.join(chunks_dir, f"{patient_id}_chunk{i:03d}.csv")
    if os.path.exists(chunk_csv):
        chunk_dfs.append(pd.read_csv(chunk_csv))
    else:
        missing.append(i)

if missing:
    print(f"WARNING: missing chunks {missing} — results will be incomplete.")

all_chunks = pd.concat(chunk_dfs, ignore_index=True)
perm_arr   = all_chunks["nrmse"].values
n_perms    = len(perm_arr)
p_value    = np.mean(perm_arr <= true_nrmse)

print(f"{patient_id} [{condition_label}]: {n_perms} permutations, "
      f"true NRMSE={true_nrmse:.4f}, p={p_value:.4f}")

# ── nrmse CSV ──────────────────────────────────────────────────────────────────
rows = [
    {"patient_id": patient_id, "seeg_condition": condition_label,
     "type": "permuted", "iteration": int(r["iteration"]), "nrmse": r["nrmse"]}
    for _, r in all_chunks.iterrows()
]
rows.append({"patient_id": patient_id, "seeg_condition": condition_label,
             "type": "true", "iteration": np.nan, "nrmse": true_nrmse})
nrmse_csv = os.path.join(out_dir, f"{patient_id}_{condition_label}_permutation_nrmse.csv")
pd.DataFrame(rows).to_csv(nrmse_csv, index=False)

# ── summary CSV ────────────────────────────────────────────────────────────────
summary_csv = os.path.join(out_dir, f"{patient_id}_{condition_label}_permutation_summary.csv")
pd.DataFrame([{
    "patient_id":      patient_id,
    "seeg_condition":  condition_label,
    "model_choice":    "OLS",
    "cv_choice":       "LOO",
    "r_val":           float(true_df["r_val"].iloc[0]),
    "pearson_p":       float(true_df["pearson_p"].iloc[0]),
    "mse":             float(true_df["mse"].iloc[0]),
    "r_squared":       float(true_df["r_squared"].iloc[0]),
    "true_nrmse":      true_nrmse,
    "perm_p_value":    p_value,
    "n_permutations":  n_perms,
    "n_sessions":      int(true_df["n_sessions"].iloc[0]),
    "n_seeg_channels": int(true_df["n_seeg_channels"].iloc[0]),
}]).to_csv(summary_csv, index=False)

print(f"Saved: {nrmse_csv}")
print(f"Saved: {summary_csv}")
