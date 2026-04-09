import pandas as pd
import os
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# read in features csv
features_csv_path = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/CATDI_master_features.csv'
feature_pd = pd.read_csv(features_csv_path)

# output directory
save_dir = "/Users/sophiapouya/workspace/bcm/CATDI/depression-severity-dbs/feature_corr_plots"
os.makedirs(save_dir,exist_ok=True)

features_dict= { "Hjorth":
                    ["RawHjorth_Mobility", "RawHjorth_Complexity", "RawHjorth_Activity"],
                "Fft":
                    ["fft_delta_mean", "fft_theta_mean", "fft_alpha_mean",
                    "fft_beta_mean","fft_low_gamma_mean","fft_high_gamma_mean"],
                "Sharpwave":
                    ["Sharpwave_Mean_prominence_range_12_30", "Sharpwave_Max_prominence_range_12_30", 
                    "Sharpwave_Mean_interval_range_12_30", "Sharpwave_Mean_sharpness_range_12_30",
                    "Sharpwave_Max_sharpness_range_12_30"],
                "Fooof": 
                    ["fooof_a_exp","fooof_a_offset"]
}

for patient in feature_pd["patient_id"].unique():
    patient_df = feature_pd[feature_pd["patient_id"] == patient]
    patient_df_data= patient_df.drop(columns=["time", "patient_id", "session_name"])
    patient_df_filled = patient_df_data.fillna(patient_df_data.mean(axis=0))
    
    # patient specific output directory
    patient_save_dir = os.path.join(save_dir,patient)
    os.makedirs(patient_save_dir, exist_ok=True)
    
    for feature, values in features_dict.items():
        correlation_vals = []
        if feature == "Coherence":
            # handle this a bit differently
            continue
        else: 
            for value in values:
               for col in patient_df_filled.columns:
                   if value in col:
                       col_parts = col.split("_")
                       contact_name = col_parts[0]+"_"+col_parts[1]
                       single_feature_vals = patient_df_filled[col]
                       catdi_scores = patient_df_filled["catdi_score"]
                       # perform the correlation
                       r, p = pearsonr(single_feature_vals, catdi_scores)
                       correlation_row = {
                           "contact": contact_name,
                           "feature": value,
                           "r": r,
                           "p":p
                       }
                       correlation_vals.append(correlation_row)

        correlation_df = pd.DataFrame(correlation_vals)
        corr_matrix = correlation_df.pivot(index="feature", columns="contact", values="r")
        pval_matrix = correlation_df.pivot(index="feature", columns="contact", values="p")
        if feature == "Fft":
            # drawn from top to bottom
            desired_order = ["fft_delta_mean", "fft_theta_mean","fft_alpha_mean", 
                             "fft_beta_mean","fft_low_gamma_mean", "fft_high_gamma_mean"]
            corr_matrix=corr_matrix.loc[desired_order]
            pval_matrix=pval_matrix.loc[desired_order]
        annotations = pval_matrix.map(lambda p: '*' if p<0.05 else "")

        plt.figure(figsize=(12,10))
        plt.title(f"{patient}_{feature}_Heatmap")
        plt.xlabel("Contacts")
        plt.ylabel("Features")
        sns.heatmap(corr_matrix, cmap='RdBu_r',annot=annotations,fmt="", vmin=-1, vmax=1)
        plt.tight_layout()
        fig_path = os.path.join(patient_save_dir,f"{patient}_{feature}_correlation.png")
        plt.savefig(fig_path)
        # plt.show()
        plt.close()

    