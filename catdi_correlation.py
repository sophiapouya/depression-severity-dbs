import pandas as pd
import os
import numpy as np
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt


subj_name = "DBSTRD014"
ref_type = "bipolar_alternating"

catdi_scores_excel = "/Users/sophiapouya/workspace/bcm/CATDI/CATDI_scores.xlsx"
base_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
plot_dir = "/Users/sophiapouya/workspace/bcm/CATDI/depression-severity-dbs/corr_plots"

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
sns.heatmap(data=plot_df, cmap='RdBu_r', center=0, vmin=-1, vmax=1, annot=annotations, fmt='')
plt.title(f"Correlation: {subj_name} ({ref_type})")
plt.xlabel("Contact")
plt.ylabel("Frequency Band")

# save and show
save_path = os.path.join(plot_dir, f"{subj_name}_{ref_type}_correlation.png")
plt.savefig(save_path)
plt.show()

