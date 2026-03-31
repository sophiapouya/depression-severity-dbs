import pandas as pd
import os
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt


# subj_name = "DBSTRD001"
# ref_type = "bipolar_alternating"

catdi_scores_excel = "/Users/sophiapouya/workspace/bcm/CATDI/CATDI_scores.xlsx"
all_subjs = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD014"]
ref_types = ["bipolar","car","esr","bipolar_alternating"]
base_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
plot_dir = "/Users/sophiapouya/workspace/bcm/CATDI/depression-severity-dbs/corr_plots"

for subj_name in all_subjs: 
    for ref_type in ref_types: 
        sbj_dir = os.path.join(base_dir, subj_name)
        if ref_type == "bipolar":
            bipolar_power_dir = os.path.join(sbj_dir, "bipolar_channels","power_bipolar")
            file = os.path.join(bipolar_power_dir, f"{subj_name}_bipolar_power.csv")
        elif ref_type == "car":
            car_power_dir = os.path.join(sbj_dir, "car_channels","power_car")
            file = os.path.join(car_power_dir, f"{subj_name}_car_power.csv")
        elif ref_type == "esr":
            esr_power_dir = os.path.join(sbj_dir, "esr_channels","power_esr")
            file = os.path.join(esr_power_dir, f"{subj_name}_esr_power.csv")
        elif ref_type == "bipolar_alternating":
            bipolar_alternating_power_dir = os.path.join(sbj_dir, "bipolar_alternating_channels","power_bipolar_alternating")
            file = os.path.join(bipolar_alternating_power_dir, f"{subj_name}_bipolar_alternating_power.csv")

        # add the catdi scores to the df saved
        subj_data = pd.read_csv(file)

        # bring in the catdi scores
        catdi_excel = pd.read_excel(catdi_scores_excel, sheet_name=subj_name)

        # clean up the names if they have text before CATDI
        if subj_name in ["DBSTRD011", "DBSTRD014"]:
            catdi_excel["session"] = catdi_excel["Name"].str.split("_task-").str[1]

        else:
            catdi_excel["session"] = catdi_excel["Name"]

        # Force both 'session' columns to be strings so they match perfectly
        subj_data["session"] = subj_data["session"].astype(str)
        catdi_excel["session"] = catdi_excel["session"].astype(str)

        subj_data = pd.merge(subj_data, catdi_excel[["session", "Result"]], how="left", on="session")

        # rename the column 
        subj_data = subj_data.rename(columns={"Result":"score"})

        # bands (make sure these match your CSV column names exactly)
        bands = ['delta', 'theta', 'alpha', 'beta', 'low_gamma', 'high_gamma']

        # Use LaTeX labels for the plot
        band_labels = ['$\delta$', '$\\theta$', '$\\alpha$', '$\\beta$', '$\gamma$', 'h$\gamma$']

        plot_data = {}
        p_val_data = {}

        # loop through each unique contact (e.g., 'LSCC-1', 'LSCC-2')
        for contact in subj_data['ch_name'].unique():

            # filter data for just this one contact
            contact_df = subj_data[subj_data['ch_name'] == contact].dropna(subset=['score'])
            r_values = []
            p_values = []

            for band in bands:

                # correlation between the power band and the merged score
                r, p = pearsonr(contact_df[band], contact_df['score'])
                r_values.append(r)
                p_values.append(p)
            
            plot_data[contact] = r_values
            p_val_data[contact] = p_values

        # create the dataframes for plotting
        plot_df = pd.DataFrame(plot_data, index=band_labels).sort_index(axis=1)
        p_val_df = pd.DataFrame(p_val_data, index=band_labels).sort_index(axis=1)

        # annotations for significant p values
        annotations = p_val_df.map(lambda p: '*' if p < 0.05 else '')

        # plot
        plt.figure(figsize=(12, 8))
        ax = sns.heatmap(data=plot_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1, annot=annotations, fmt='')
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=20)

        plt.title(f"Correlation: {subj_name} ({ref_type})",fontsize=20)
        plt.xlabel("Contact", fontsize=16)
        plt.ylabel("Frequency Band", fontsize=16)
        plt.tight_layout()

        # save and show
        save_path = os.path.join(plot_dir, f"{subj_name}_{ref_type}_correlation_segmented_sections.png")
        plt.savefig(save_path)
        plt.close()

