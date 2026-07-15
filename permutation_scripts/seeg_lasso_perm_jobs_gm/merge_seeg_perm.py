"""
Merge chunk outputs into the final permutation result files (greymatter version).

Expects chunks at:
  {BASE_DIR_SEEG}/lasso_permutation_test_gm/chunks/
    {patient_id}_true_decoder.csv
    {patient_id}_chunk000.csv ... chunk{n_chunks-1}.csv

Writes final results to:
  {BASE_DIR_SEEG}/lasso_permutation_test_gm/
    {patient_id}_seeg_permutation_nrmse.csv
    {patient_id}_seeg_permutation_summary.csv

Run from project root:
  python seeg_perm_jobs_gm/merge_seeg_perm.py --patient_id DBSTRD001 --n_chunks 20
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd

from config import BASE_DIR_SEEG

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id", type=str, required=True)
parser.add_argument("--n_chunks", type=int, default=20)
args = parser.parse_args()

patient_id = args.patient_id
n_chunks = args.n_chunks
base_dir = str(BASE_DIR_SEEG)
chunks_dir = os.path.join(base_dir, "lasso_permutation_test_gm", "chunks")
out_dir = os.path.join(base_dir, "lasso_permutation_test_gm")
os.makedirs(out_dir, exist_ok=True)

# ── load true decoder ─────────────────────────────────────────────────────────
true_csv = os.path.join(chunks_dir, f"{patient_id}_true_decoder.csv")
if not os.path.exists(true_csv):
    raise FileNotFoundError(f"True decoder file not found: {true_csv}\n"
                            "Make sure chunk 0 completed successfully.")
true_df = pd.read_csv(true_csv)
true_nrmse = float(true_df["true_nrmse"].iloc[0])

# ── load and merge all chunks ─────────────────────────────────────────────────
chunk_dfs = []
missing = []
for i in range(n_chunks):
    chunk_csv = os.path.join(chunks_dir, f"{patient_id}_chunk{i:03d}.csv")
    if os.path.exists(chunk_csv):
        chunk_dfs.append(pd.read_csv(chunk_csv))
    else:
        missing.append(i)

if missing:
    print(f"WARNING: missing chunks {missing} — results will be incomplete.")

all_chunks = pd.concat(chunk_dfs, ignore_index=True)
perm_arr = all_chunks["nrmse"].values
n_permutations = len(perm_arr)
print(f"{patient_id}: {n_permutations} permutations merged, true NRMSE={true_nrmse:.4f}")

# ── p-value ───────────────────────────────────────────────────────────────────
p_value = np.mean(perm_arr <= true_nrmse)
print(f"p-value = {p_value:.4f}")

# ── write nrmse CSV (matches format used by histogram scripts) ────────────────
rows = [
    {"patient_id": patient_id, "type": "permuted", "iteration": int(row["iteration"]), "nrmse": row["nrmse"]}
    for _, row in all_chunks.iterrows()
]
rows.append({"patient_id": patient_id, "type": "true", "iteration": np.nan, "nrmse": true_nrmse})
nrmse_df = pd.DataFrame(rows)
nrmse_csv = os.path.join(out_dir, f"{patient_id}_seeg_permutation_nrmse.csv")
nrmse_df.to_csv(nrmse_csv, index=False)

# ── write summary CSV ─────────────────────────────────────────────────────────
summary_df = pd.DataFrame([{
    "patient_id": patient_id,
    "r_val": float(true_df["r_val"].iloc[0]),
    "pearson_p": float(true_df["pearson_p"].iloc[0]),
    "mse": float(true_df["mse"].iloc[0]),
    "r_squared": float(true_df["r_squared"].iloc[0]),
    "true_nrmse": true_nrmse,
    "perm_p_value": p_value,
    "n_permutations": n_permutations,
    "n_sessions": int(true_df["n_sessions"].iloc[0]),
}])
summary_csv = os.path.join(out_dir, f"{patient_id}_seeg_permutation_summary.csv")
summary_df.to_csv(summary_csv, index=False)

print(f"Saved: {nrmse_csv}")
print(f"Saved: {summary_csv}")
