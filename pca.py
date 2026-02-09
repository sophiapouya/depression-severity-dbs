import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# define input csv
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData'
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")
all_patient_df = pd.read_csv(csv_path)
output_dir = os.path.join(base_dir,"pca")

# remove columns that aren't features 
columns_to_exclude = ['session_name', 'time']
clean_df =all_patient_df.drop(columns_to_exclude, axis=1)

# apply the 70% rule
# go through each patient
patient_pca_list = {}
for patient in clean_df['patient_id'].unique():
    patient_df = clean_df[clean_df['patient_id'] == patient]
    # remove catdi and patient id columns 
    patient_df = patient_df.drop(['patient_id','catdi_score'], axis=1)
    
    # implement 70% cutoff-> 70% of sessions need a value for a feature to be analyzed
    sessions= patient_df.shape[0]
    cutoff = int(np.ceil(.70*sessions))
    boolean_df = patient_df.notna()
    sums = boolean_df.sum()
    columns_to_keep = sums[sums >= cutoff].index
    patient_df = patient_df[columns_to_keep]
    # fill in the remaining gaps with average values from columns since pca and lasso cannot deal with blank values 
    patient_df = patient_df.fillna(patient_df.mean())

    # perform pca analysis 
    pca = PCA()
    pca.fit_transform(patient_df)
    eigenvalues = pca.explained_variance_
    patient_pca_list[patient] = eigenvalues

patient_list = clean_df['patient_id'].unique()

# figure with 3 rows and 2 columns
fig, axes = plt.subplots(3,2, figsize=(12,8))
axes = axes.flatten()

for index, patient in enumerate(patient_list):

    #find the point where 90% of the variance is explained
    total_variance = np.sum(patient_pca_list[patient])
    variance_ratios = patient_pca_list[patient]/total_variance

    # 90% variance explained cutoff
    sum_total= 0
    cutoff_index= 0
    for idx, ratio in enumerate(variance_ratios):
        sum_total = sum_total + ratio
        if sum_total >= 0.90:
            cutoff_index = idx + 1
            break

    x_axis= np.arange(1, len(patient_pca_list[patient][:cutoff_index])+1)
    # bar graph for components and eigenvalues 
    axes[index].bar(x_axis, patient_pca_list[patient][:cutoff_index])
    # scatter to see the points above the bars and connect them with a dashed line
    axes[index].scatter(x_axis, patient_pca_list[patient][:cutoff_index],color='black')
    # add a line connecting the dots above the plots for 90% variance explained
    axes[index].plot(x_axis, patient_pca_list[patient][:cutoff_index], linestyle="--", color='black')
    # label the points with the values above for 90% variance explained
    for idx, value in enumerate(patient_pca_list[patient][:cutoff_index]):
        axes[index].text(idx+1, 1.05*value, f"{value:.2f}", ha="center", va="bottom",color='black')
    
    axes[index].set_xticks(x_axis)
    axes[index].set_xlabel("Dimensions")
    axes[index].set_ylabel("Eigenvalues")
    axes[index].spines["top"].set_visible(False)
    axes[index].set_title(f"{patient} Scree Plot")
    
plt.tight_layout()
output_file = os.path.join(output_dir,"scree_plots.png")
plt.savefig(output_file)
plt.show()




