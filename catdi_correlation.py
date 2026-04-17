import pandas as pd
import os
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import fdrcorrection
from config import BASE_DIR_DBS, BASE_DIR_SEEG, CATDI_SCORES, CATDI_ELECTRODES

ELEC_TYPE = "DBS"  # either DBS or SEEG

catdi_scores_excel = str(CATDI_ELECTRODES)
all_subjs = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
#all_subjs = ["DBSTRD011", "DBSTRD014"]

if ELEC_TYPE == "SEEG":
    ref_types = ["bipolar"]
    base_dir = str(BASE_DIR_SEEG)
    plot_dir = os.path.join(base_dir, "corr_plots_seeg")
    os.makedirs(plot_dir,exist_ok=True)
else: 
    ref_types = ["bipolar_alternating"]
    base_dir = str(BASE_DIR_DBS)
    plot_dir = os.path.join(base_dir, "corr_plots_dbs")
    os.makedirs(plot_dir, exist_ok=True)

for subj_name in all_subjs: 

    if ELEC_TYPE == "SEEG":
        metadata_file = os.path.join(base_dir, subj_name, f"{subj_name}_metadata.csv")
        # csv with columns channel_name and channel_region
        metadata_df = pd.read_csv(metadata_file)
       
        # ensure the name is consistent with power data file name
        metadata_df["ch_name"] = metadata_df["channel_name"]

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
        
        if ELEC_TYPE == "SEEG":
            subj_data = pd.merge(subj_data, metadata_df, how="left", on="ch_name")
            # create the channel mapping 
            channel_region_map = dict(zip(subj_data["ch_name"], subj_data["channel_region"]))
        
        # rename the column 
        subj_data = subj_data.rename(columns={"Result":"score"})

        # bands (make sure these match your CSV column names exactly)
        bands = ['delta', 'theta', 'alpha', 'beta', 'low_gamma', 'high_gamma']

        # Use LaTeX labels for the plot
        band_labels = [r'$\delta$', r'$\theta$', r'$\alpha$', r'$\beta$', r'$\gamma$', r'h$\gamma$']

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

        if ELEC_TYPE == "SEEG":
            sorted_chans = sorted(
                plot_df.columns, 
                key= lambda ch: channel_region_map[ch]
            )
            plot_df = plot_df[sorted_chans]
            p_val_df = p_val_df[sorted_chans]

        # fdr correction for the pvals
        p_vals_flat = p_val_df.values.flatten()
        rejected, correct_p_vals=fdrcorrection(p_vals_flat)
        p_val_df_corrected = pd.DataFrame(
            correct_p_vals.reshape(p_val_df.shape),
            index=p_val_df.index,
            columns=p_val_df.columns
        )

        # annotations for significant p values
        annotations = p_val_df_corrected.map(lambda p: '*' if p < 0.05 else '')

        # plot
        plt.figure(figsize=(12, 8))


        if ELEC_TYPE == "SEEG":
            region_labels = [channel_region_map[ch].upper() for ch in plot_df.columns]
            ax = sns.heatmap(data=plot_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1, annot=annotations, annot_kws={"color": "black", "fontsize":14}, fmt='',  xticklabels=region_labels)
        else:
            ax = sns.heatmap(data=plot_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1, annot=annotations, annot_kws={"color": "black", "fontsize":14}, fmt='')
            ax.set_xticklabels(ax.get_xticklabels(), fontsize=16)

        ax.set_yticklabels(ax.get_yticklabels(), fontsize=20)

        plt.title(f"Correlation: {subj_name} ({ref_type})",fontsize=20)
        plt.xlabel("Contact", fontsize=16)
        plt.ylabel("Frequency Band", fontsize=16)
        plt.tight_layout()

        # save and show
        save_path = os.path.join(plot_dir, f"{subj_name}_{ref_type}.png")
        plt.savefig(save_path)
        plt.close()

