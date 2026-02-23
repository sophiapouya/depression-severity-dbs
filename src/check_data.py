import pandas as pd
import numpy as np
import os

# Setup paths
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData'
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")

# Load data
master_df = pd.read_csv(csv_path)

# Separate metadata
meta_cols = ['patient_id', 'session_name', 'catdi_score', 'time']
feature_cols = [c for c in master_df.columns if c not in meta_cols]

# Per-patient feature counts
feature_counts = (
    master_df
    .groupby('patient_id')[feature_cols]
    .apply(lambda x: x.notnull().sum())
)

# Example: look at patient 1
patient1_counts = feature_counts.loc["DBSTRD001"].sort_values()

pid = "DBSTRD001"
n_sessions = master_df.loc[master_df["patient_id"] == pid].shape[0]
print("DBSTRD001 total sessions:", n_sessions)

print(patient1_counts.head(20))
print(patient1_counts.tail(20))