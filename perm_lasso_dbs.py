import argparse
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import LeaveOneOut
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from scipy.stats import pearsonr

from src.postprocessing_functions import get_nrmse
from config import CATDI_SCORES, BASE_DIR_DBS

# -----------------------------
# Settings
# -----------------------------
cv_choice = "LOO"
catdi_scores_file = str(CATDI_SCORES)
base_dir = str(BASE_DIR_DBS)

parser = argparse.ArgumentParser()
parser.add_argument("--patient_id", type=str, required=True)
args = parser.parse_args()
patient_id = args.patient_id

l1_reg_list = np.around(np.arange(0.1, 1.1, 0.1), 1)
outlier = 4

n_permutations = 100
random_seed = 42

output_dir = os.path.join(base_dir, "lasso_permutation_test")
os.makedirs(output_dir, exist_ok=True)


# -----------------------------
# Load patient data
# -----------------------------
patient_power_csv = os.path.join(
    base_dir,
    patient_id,
    "bipolar_alternating_channels",
    "power_bipolar_alternating",
    f"{str(patient_id)}_bipolar_alternating_power.csv"
)
patient_df_raw = pd.read_csv(patient_power_csv)

catdi_excel = pd.read_excel(catdi_scores_file, sheet_name=patient_id)
if patient_id in ["DBSTRD014", "DBSTRD011"]:
    catdi_excel["Name"] = catdi_excel["Name"].astype(str).str.split("task-").str[-1]

true_score_map = dict(zip(catdi_excel["Name"], catdi_excel["Result"]))


# -----------------------------
# Helper: run your full nested decoder once
# -----------------------------
def run_decoder_for_one_patient(
    patient_df_input: pd.DataFrame,
    score_map: dict[str, float],
    l1_reg_list: np.ndarray,
    outlier: float,
):
    patient_df = patient_df_input.copy()
    patient_df["scores"] = patient_df["session"].map(score_map)
    patient_df = patient_df.dropna(subset=["scores"]).copy()

    folds = LeaveOneOut()
    unique_sessions = patient_df["session"].unique()

    catdi_predict_vals = []
    catdi_measured_vals = []

    for train_index, test_index in folds.split(unique_sessions):
        train_sessions = unique_sessions[train_index]
        test_sessions = unique_sessions[test_index]

        data_train = patient_df[patient_df["session"].isin(train_sessions)]
        data_test = patient_df[patient_df["session"].isin(test_sessions)]

        data_train_pivot = data_train.pivot_table(
            index="session",
            columns="ch_name",
            values=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
        )

        # fill NaNs using training means only
        train_feature_means = data_train_pivot.mean(axis=0).fillna(0)
        data_train_pivot = data_train_pivot.fillna(train_feature_means)

        train_catdi_scores = (
            data_train.groupby("session")["scores"].first().loc[data_train_pivot.index]
        )

        alpha_nrmse = []
        for l1_reg in l1_reg_list:
            inner_predictions = []
            inner_test_vals = []
            inner_loo = LeaveOneOut()

            for train_idx, test_idx in inner_loo.split(data_train_pivot):
                inner_train_data = data_train_pivot.iloc[train_idx].values
                inner_test_data = data_train_pivot.iloc[test_idx].values
                inner_train_scores = train_catdi_scores.iloc[train_idx].values
                inner_test_scores = train_catdi_scores.iloc[test_idx].values

                mu = np.mean(inner_train_data, axis=0)
                std = np.std(inner_train_data, axis=0)
                std[std == 0] = 1

                inner_train_data = np.where(
                    np.abs(inner_train_data - mu) > outlier * std,
                    mu,
                    inner_train_data
                )
                inner_test_data = np.where(
                    np.abs(inner_test_data - mu) > outlier * std,
                    mu,
                    inner_test_data
                )

                inner_model = make_pipeline(
                    StandardScaler(),
                    Lasso(alpha=l1_reg, random_state=0, max_iter=2000)
                )
                inner_model.fit(inner_train_data, inner_train_scores)

                inner_prediction = inner_model.predict(inner_test_data)
                inner_predictions.extend(inner_prediction)
                inner_test_vals.extend(inner_test_scores)

            current_alpha_nrmse = get_nrmse(inner_test_vals, inner_predictions)
            alpha_nrmse.append(current_alpha_nrmse)

        best_alpha_idx = np.argmin(alpha_nrmse)
        best_alpha = l1_reg_list[best_alpha_idx]

        data_test_pivot = data_test.pivot_table(
            index="session",
            columns="ch_name",
            values=["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
        ).reindex(columns=data_train_pivot.columns)

        data_test_pivot = data_test_pivot.fillna(train_feature_means)

        catdi_scores_outer_train = (
            data_train.groupby("session")["scores"].first().loc[data_train_pivot.index]
        )
        catdi_scores_outer_test = (
            data_test.groupby("session")["scores"].first().loc[data_test_pivot.index]
        )

        data_train_pivot = data_train_pivot.values
        data_test_pivot = data_test_pivot.values

        mu_outer = np.mean(data_train_pivot, axis=0)
        std_outer = np.std(data_train_pivot, axis=0)
        std_outer[std_outer == 0] = 1

        outer_train_clamped = np.where(
            np.abs(data_train_pivot - mu_outer) > outlier * std_outer,
            mu_outer,
            data_train_pivot
        )
        outer_test_clamped = np.where(
            np.abs(data_test_pivot - mu_outer) > outlier * std_outer,
            mu_outer,
            data_test_pivot
        )

        outer_model = make_pipeline(
            StandardScaler(),
            Lasso(alpha=best_alpha, random_state=42, max_iter=100000)
        )
        outer_model.fit(outer_train_clamped, catdi_scores_outer_train)
        catdi_predicted = outer_model.predict(outer_test_clamped)

        catdi_predict_vals.extend(catdi_predicted)
        catdi_measured_vals.extend(catdi_scores_outer_test)

    true_nrmse = get_nrmse(catdi_measured_vals, catdi_predict_vals)

    return (
        np.array(catdi_measured_vals),
        np.array(catdi_predict_vals),
        true_nrmse,
    )


# -----------------------------
# Run true decoder
# -----------------------------
true_measured, true_predicted, true_nrmse = run_decoder_for_one_patient(
    patient_df_input=patient_df_raw,
    score_map=true_score_map,
    l1_reg_list=l1_reg_list,
    outlier=outlier,
)

true_r, true_pearson_p = pearsonr(true_measured, true_predicted)
true_mse = mean_squared_error(true_measured, true_predicted)
true_r2 = r2_score(true_measured, true_predicted)

# -----------------------------
# Permutation test
# -----------------------------
rng = np.random.default_rng(random_seed)

valid_sessions = patient_df_raw["session"].unique()
valid_sessions = [s for s in valid_sessions if s in true_score_map]
true_scores = np.array([true_score_map[s] for s in valid_sessions])

perm_nrmse_list = []

for i in range(n_permutations):
    shuffled_scores = rng.permutation(true_scores)
    perm_score_map = dict(zip(valid_sessions, shuffled_scores))

    _, _, perm_nrmse = run_decoder_for_one_patient(
        patient_df_input=patient_df_raw,
        score_map=perm_score_map,
        l1_reg_list=l1_reg_list,
        outlier=outlier
    )
    perm_nrmse_list.append(perm_nrmse)

    if (i + 1) % 10 == 0:
        print(f"Finished permutation {i + 1}/{n_permutations}")

perm_nrmse_array = np.array(perm_nrmse_list)

# lower NRMSE = better
perm_p_value = np.mean(perm_nrmse_array <= true_nrmse)

# format p-value: scientific if very small, else 4 decimal places
if perm_p_value < 1e-4:
    p_str = f"{perm_p_value:.4e}"
else:
    p_str = f"{perm_p_value:.4f}"


# -----------------------------
# Save summary CSV
# -----------------------------
summary_df = pd.DataFrame([{
    "patient_id": patient_id,
    "r_val": true_r,
    "pearson_p": true_pearson_p,
    "mse": true_mse,
    "r_squared": true_r2,
    "true_nrmse": true_nrmse,
    "perm_p_value": perm_p_value,
    "n_permutations": n_permutations,
    "n_sessions": len(valid_sessions),
}])

summary_csv = os.path.join(output_dir, f"{patient_id}_permutation_summary.csv")
summary_df.to_csv(summary_csv, index=False)


# -----------------------------
# Save NRMSE distribution CSV for plotting
# -----------------------------
rows = []

for i, val in enumerate(perm_nrmse_array):
    rows.append({
        "patient_id": patient_id,
        "type": "permuted",
        "iteration": i,
        "nrmse": val
    })

rows.append({
    "patient_id": patient_id,
    "type": "true",
    "iteration": np.nan,
    "nrmse": true_nrmse
})

perm_csv_df = pd.DataFrame(rows)
perm_csv_path = os.path.join(output_dir, f"{patient_id}_permutation_nrmse.csv")
perm_csv_df.to_csv(perm_csv_path, index=False)


# -----------------------------
# Scatter plot: measured vs predicted
# -----------------------------
scatter_path = os.path.join(output_dir, f"{patient_id}_true_decoder_scatter.png")

plt.figure(figsize=(5, 5))
min_lim = min(np.min(true_measured), np.min(true_predicted))
max_lim = max(np.max(true_measured), np.max(true_predicted))

plt.plot([min_lim, max_lim], [min_lim, max_lim], linestyle="--", color="black")
plt.scatter(true_measured, true_predicted, color="blue", s=50)
plt.xlim(min_lim, max_lim)
plt.ylim(min_lim, max_lim)
plt.xlabel("Measured CATDI")
plt.ylabel("Predicted CATDI")
plt.title(patient_id)
plt.gca().set_box_aspect(1)
plt.text(
    0.05, 0.85,
    f"R={true_r:.2f}\nperm p={p_str}",
    transform=plt.gca().transAxes,
    fontsize=12
)
plt.tight_layout()
plt.savefig(scatter_path, dpi=300)
plt.close()


# -----------------------------
# Histogram: permuted NRMSEs + true NRMSE
# -----------------------------
hist_path = os.path.join(output_dir, f"{patient_id}_permutation_histogram.png")

plt.figure(figsize=(6, 4))
plt.hist(perm_nrmse_array, bins=30, color="lightgray", edgecolor="black")
plt.axvline(true_nrmse, color="red", linestyle="--", linewidth=2)
y_max = plt.ylim()[1]
plt.text(
    true_nrmse,
    y_max * 0.9,
    f"p = {p_str}",
    color="black",
    ha="left",   # keeps it just left of the line
    va="top",
    fontsize=12,
    bbox=dict(facecolor="white", edgecolor="none", alpha=0.7)
)
plt.xlabel("Permuted NRMSE")
plt.ylabel("Count")
plt.title(f"{patient_id} permutation test")
plt.tight_layout()
plt.savefig(hist_path, dpi=300)
plt.close()

