"""
SEEG Lasso permutation test (greymatter filtered) — chunk-based for SLURM array jobs.

1000 permutations are split into N_CHUNKS array elements. Each chunk runs
its slice of permutations in parallel via joblib. Chunk 0 also runs and
saves the true decoder result.

Output goes to:
  {BASE_DIR_SEEG}/lasso_permutation_test_gm/chunks/
    {patient_id}_true_decoder.csv          (chunk 0 only)
    {patient_id}_chunk{chunk_id:03d}.csv   (each chunk)

Run from project root:
  python seeg_perm_jobs_gm/perm_lasso_seeg_chunk.py \\
      --patient_id DBSTRD001 --chunk_id 0 --n_chunks 20
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy.stats import pearsonr
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import LeaveOneOut
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from config import BASE_DIR_SEEG, CATDI_SCORES
from src.postprocessing_functions import get_nrmse

# ── args ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--patient_id", type=str, required=True)
parser.add_argument("--chunk_id", type=int, required=True)
parser.add_argument("--n_chunks", type=int, default=20)
parser.add_argument("--n_permutations", type=int, default=1000)
parser.add_argument("--n_jobs", type=int, default=8)
args = parser.parse_args()

patient_id = args.patient_id
chunk_id = args.chunk_id
n_chunks = args.n_chunks
n_permutations = args.n_permutations
n_jobs = args.n_jobs

chunk_size = int(np.ceil(n_permutations / n_chunks))
chunk_start = chunk_id * chunk_size
chunk_end = min(chunk_start + chunk_size, n_permutations)
perm_indices = list(range(chunk_start, chunk_end))

print(f"Patient={patient_id}, chunk={chunk_id}/{n_chunks}, "
      f"perms {chunk_start}-{chunk_end-1}, n_jobs={n_jobs}")

# ── settings ──────────────────────────────────────────────────────────────────
REGIONS = ["amy", "ofc", "vmpfc", "acc", "dlpfc"]
L1_REGS = np.around(np.arange(0.1, 1.1, 0.1), 1)
OUTLIER = 4
RANDOM_SEED = 42

base_dir = str(BASE_DIR_SEEG)
catdi_scores_file = str(CATDI_SCORES)

chunks_dir = os.path.join(base_dir, "lasso_permutation_test_gm", "chunks")
os.makedirs(chunks_dir, exist_ok=True)

# ── load data ─────────────────────────────────────────────────────────────────
patient_power_csv = os.path.join(
    base_dir, patient_id, "bipolar_channels", "power_bipolar",
    f"{patient_id}_bipolar_power_greymatter.csv",
)
patient_df_raw = pd.read_csv(patient_power_csv)

channel_metadata_file = os.path.join(base_dir, patient_id, f"{patient_id}_metadata_greymatter.csv")
channel_region_df = pd.read_csv(channel_metadata_file)
channel_region_df["ch_name"] = channel_region_df["channel_name"]

catdi_excel = pd.read_excel(catdi_scores_file, sheet_name=patient_id)
if patient_id in ["DBSTRD014", "DBSTRD011"]:
    catdi_excel["Name"] = catdi_excel["Name"].astype(str).str.split("task-").str[-1]
true_score_map = dict(zip(catdi_excel["Name"], catdi_excel["Result"]))

valid_sessions = [s for s in patient_df_raw["session"].unique() if s in true_score_map]
true_scores = np.array([true_score_map[s] for s in valid_sessions])


# ── decoder ───────────────────────────────────────────────────────────────────
def run_seeg_decoder_once(patient_df_input, score_map):
    patient_df = patient_df_input.copy()
    patient_df["scores"] = patient_df["session"].map(score_map)
    patient_df = patient_df.dropna(subset=["scores"]).copy()
    patient_df = pd.merge(
        patient_df,
        channel_region_df[["ch_name", "channel_region"]],
        how="left", on="ch_name",
    )

    unique_sessions = patient_df["session"].unique()
    catdi_predict_vals, catdi_measured_vals = [], []

    for train_index, test_index in LeaveOneOut().split(unique_sessions):
        train_sessions = unique_sessions[train_index]
        test_sessions = unique_sessions[test_index]
        data_train = patient_df[patient_df["session"].isin(train_sessions)]
        data_test = patient_df[patient_df["session"].isin(test_sessions)]

        # inner CV: choose best region + alpha
        region_scores, region_alphas = {}, {}
        for region in REGIONS:
            region_df = data_train[data_train["channel_region"] == region]
            pivot = region_df.pivot_table(
                index="session", columns="ch_name",
                values=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"],
            )
            if pivot.shape[0] == 0:
                region_scores[region] = np.inf
                region_alphas[region] = np.nan
                continue

            train_scores_r = region_df.groupby("session")["scores"].first().loc[pivot.index]
            alpha_nrmse = []
            for l1 in L1_REGS:
                inner_pred, inner_true = [], []
                for ii_tr, ii_te in LeaveOneOut().split(pivot):
                    X_tr = pivot.iloc[ii_tr].values
                    X_te = pivot.iloc[ii_te].values
                    y_tr_i = train_scores_r.iloc[ii_tr].values
                    y_te_i = train_scores_r.iloc[ii_te].values
                    mu = np.mean(X_tr, axis=0)
                    sd = np.std(X_tr, axis=0); sd[sd == 0] = 1
                    X_tr = np.where(np.abs(X_tr - mu) > OUTLIER * sd, mu, X_tr)
                    X_te = np.where(np.abs(X_te - mu) > OUTLIER * sd, mu, X_te)
                    m = make_pipeline(StandardScaler(), Lasso(alpha=l1, random_state=0, max_iter=2000))
                    m.fit(X_tr, y_tr_i)
                    inner_pred.extend(m.predict(X_te))
                    inner_true.extend(y_te_i)
                alpha_nrmse.append(get_nrmse(inner_true, inner_pred))
            best_idx = int(np.argmin(alpha_nrmse))
            region_scores[region] = alpha_nrmse[best_idx]
            region_alphas[region] = L1_REGS[best_idx]

        best_region = min(region_scores, key=region_scores.get)
        best_alpha = region_alphas[best_region]

        # outer fit
        train_r = data_train[data_train["channel_region"] == best_region]
        pivot_tr = train_r.pivot_table(
            index="session", columns="ch_name",
            values=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"],
        )
        test_r = data_test[data_test["channel_region"] == best_region]
        pivot_te = test_r.pivot_table(
            index="session", columns="ch_name",
            values=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"],
        ).reindex(columns=pivot_tr.columns)

        y_tr = train_r.groupby("session")["scores"].first().loc[pivot_tr.index]
        y_te = test_r.groupby("session")["scores"].first().loc[pivot_te.index]

        X_tr = pivot_tr.values; X_te = pivot_te.values
        mu = np.mean(X_tr, axis=0); sd = np.std(X_tr, axis=0); sd[sd == 0] = 1
        X_tr = np.where(np.abs(X_tr - mu) > OUTLIER * sd, mu, X_tr)
        X_te = np.where(np.abs(X_te - mu) > OUTLIER * sd, mu, X_te)

        final_model = make_pipeline(
            StandardScaler(), Lasso(alpha=best_alpha, random_state=42, max_iter=10000)
        )
        final_model.fit(X_tr, y_tr)
        catdi_predict_vals.extend(final_model.predict(X_te))
        catdi_measured_vals.extend(y_te.values)

    return (
        np.array(catdi_measured_vals),
        np.array(catdi_predict_vals),
        get_nrmse(catdi_measured_vals, catdi_predict_vals),
    )


# ── true decoder (chunk 0 only) ───────────────────────────────────────────────
if chunk_id == 0:
    print("Running true decoder...")
    true_measured, true_predicted, true_nrmse = run_seeg_decoder_once(
        patient_df_raw, true_score_map
    )
    true_r, true_pearson_p = pearsonr(true_measured, true_predicted)
    true_mse = mean_squared_error(true_measured, true_predicted)
    true_r2 = r2_score(true_measured, true_predicted)
    print(f"True NRMSE={true_nrmse:.4f}, R={true_r:.4f}")

    true_df = pd.DataFrame([{
        "patient_id": patient_id,
        "true_nrmse": true_nrmse,
        "r_val": true_r,
        "pearson_p": true_pearson_p,
        "mse": true_mse,
        "r_squared": true_r2,
        "n_sessions": len(valid_sessions),
    }])
    true_df.to_csv(
        os.path.join(chunks_dir, f"{patient_id}_true_decoder.csv"), index=False
    )


# ── permutation chunk ─────────────────────────────────────────────────────────
def run_one_permutation(perm_idx):
    rng = np.random.default_rng(RANDOM_SEED * 10000 + perm_idx)
    shuffled = rng.permutation(true_scores)
    perm_score_map = dict(zip(valid_sessions, shuffled))
    _, _, perm_nrmse = run_seeg_decoder_once(patient_df_raw, perm_score_map)
    if (perm_idx + 1) % 10 == 0:
        print(f"  Done permutation {perm_idx}")
    return perm_idx, perm_nrmse


print(f"Running {len(perm_indices)} permutations with {n_jobs} workers...")
results = Parallel(n_jobs=n_jobs, backend="loky")(
    delayed(run_one_permutation)(i) for i in perm_indices
)

chunk_df = pd.DataFrame(
    [{"patient_id": patient_id, "iteration": idx, "nrmse": nrmse}
     for idx, nrmse in results]
)
chunk_csv = os.path.join(chunks_dir, f"{patient_id}_chunk{chunk_id:03d}.csv")
chunk_df.to_csv(chunk_csv, index=False)
print(f"Saved chunk: {chunk_csv}")
