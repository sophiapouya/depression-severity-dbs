import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# paths
base_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData"
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")
output_dir = os.path.join(base_dir, "pca", "scree_plots")
os.makedirs(output_dir, exist_ok=True)

# load csv
all_patient_df = pd.read_csv(csv_path)

# remove non-feature column
clean_df = all_patient_df.drop(columns=["time"], errors="ignore")

# get all patient IDs
patient_ids = clean_df["patient_id"].unique()

# make subplot grid
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
axes = axes.flatten()

for i, patient in enumerate(patient_ids):
    # get one patient's data
    patient_df = clean_df[clean_df["patient_id"] == patient].copy()

    # keep only feature columns
    feature_df = patient_df.drop(columns=["patient_id", "catdi_score", "session_name"], errors="ignore")

    # apply 70% rule
    total_sessions = feature_df.shape[0]
    cutoff = int(np.ceil(0.70 * total_sessions))
    columns_to_keep = feature_df.columns[feature_df.notna().sum() >= cutoff]
    feature_df = feature_df[columns_to_keep]

    # fill remaining NaNs with column means
    feature_df = feature_df.fillna(feature_df.mean())

    # standardize
    scaler = StandardScaler()
    feature_matrix = scaler.fit_transform(feature_df)

    # fit PCA
    pca = PCA()
    pca.fit(feature_matrix)

    # cumulative explained variance
    cumulative_variance = np.cumsum(pca.explained_variance_ratio_)

    # x values = PC numbers
    pcs = np.arange(1, len(cumulative_variance) + 1)

    # find number of PCs for 90% variance
    pcs_90 = np.argmax(cumulative_variance >= 0.90) + 1

    # plot
    axes[i].plot(pcs, cumulative_variance, marker="o")
    axes[i].axhline(0.90, linestyle="--", color="red")
    axes[i].axvline(pcs_90, linestyle="--",color="black")
    axes[i].text(pcs_90 + 0.2, 0.2, f"{pcs_90} PCs")
    axes[i].set_title(patient)
    axes[i].set_xlabel("Number of Principal Components")
    axes[i].set_ylabel("Cumulative Explained Variance")
    axes[i].set_ylim(0, 1.05)

# tidy layout
fig.tight_layout()

# save figure
fig_path = os.path.join(output_dir, "all_patients_scree_plot_ALL.png")
fig.savefig(fig_path, dpi=300)
plt.close()

print(f"Saved figure to: {fig_path}")