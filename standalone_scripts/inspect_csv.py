import pandas as pd

input_csv ="/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/final_mast_feat.csv"
csv_df = pd.read_csv(input_csv) 

patients = csv_df["patient_id"].unique()

for patient in patients:
    empties = csv_df[csv_df["patient_id"] == patient].isna().sum().sum()
    non_empties = csv_df[csv_df["patient_id"] == patient].notna().sum().sum()
    print(f"{patient} has {empties} empty cells and {non_empties} cells")
