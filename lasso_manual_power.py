

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.linear_model import LassoCV
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr 
import matplotlib.pyplot as plt


# cv choice
cv_choice="LOO" # choices -> LOO or KFOLD

# define input files
catdi_scores_file = '/Users/sophiapouya/workspace/bcm/CATDI/CATDI_scores.xlsx'

base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/'
all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD014"]


# list for keeping track of performance for output csv
performance_metrics = []

# list for patient catdi predicted vs catdi measured for plotting
decoding_plot_info = []

# go through each patient
for sbj in all_subjs:
    patient_power_csv = os.path.join(base_dir,sbj,"bipolar_alternating_channels","power_bipolar_alternating",f"{sbj}_bipolar_alternating_power.csv")
    patient_df = pd.read_csv(patient_power_csv)

    # take an average for the power values across probes and contacts so that it's just session x power values
    patient_df_avg=patient_df.groupby("session", as_index=False).mean(numeric_only=True)
    
    # add a column for the catdi scores
    catdi_excel = pd.read_excel(catdi_scores_file, sheet_name=sbj)
    if sbj == "DBSTRD014":
        catdi_excel["Name"]= catdi_excel["Name"].astype(str).str.split("task-").str[-1]
    # create a column of scores catdi_excel["Result"] that matches the patient_power_csv["session"]
    catdi_session_and_scores = dict(zip(catdi_excel["Name"],catdi_excel["Result"]))
    
    catdi_scores = patient_df_avg["session"].map(catdi_session_and_scores)

    # remove columns that aren't features 
    patient_df_avg =patient_df_avg.drop(['session'], axis=1)

    # drop any session that don't have catdi scores
    valid_score_indices = catdi_scores.notna()
    patient_df_avg=patient_df_avg.loc[valid_score_indices].reset_index(drop=True)
    catdi_scores = catdi_scores.loc[valid_score_indices].reset_index(drop=True)
    
    # initialize kfold 
    if cv_choice == "KFOLD":
        folds = KFold(n_splits=5, shuffle=True, random_state=42)
    else:
        folds = LeaveOneOut()
    
    inner_cv = KFold(n_splits=5, shuffle=True, random_state=42)

    fold_num = 0
    catdi_predict_vals, catdi_measured_vals = [],[]
    for train_index, test_index in folds.split(patient_df_avg):
        # keep track of the current fold
        fold_num = fold_num+1
        data_train, data_test = patient_df_avg.iloc[train_index,:], patient_df_avg.iloc[test_index,:]
        catdi_train, catdi_test = catdi_scores.iloc[train_index], catdi_scores.iloc[test_index]

        # standardize the training data
        scaler = StandardScaler()
        # get the mean and std dev of the training data
        scaler.fit(data_train)
        # z score the training data
        data_train_transformed = scaler.transform(data_train)
        # z score the testing data w/ mean and std dev from the training data
        data_test_transformed= scaler.transform(data_test)

        # perform lasso cv regression
        model = LassoCV(cv=inner_cv, random_state=42,max_iter=100000)
        # fit the model 
        model.fit(data_train_transformed, catdi_train)

        # predict on the test set
        catdi_predicted = model.predict(data_test_transformed)

        if cv_choice == "KFOLD":
            # pearson correlation
            r_value, p_value = pearsonr(catdi_test, catdi_predicted)
            # mean squared error
            mse = mean_squared_error(catdi_test, catdi_predicted)
            # r-squared 
            r2 = r2_score(catdi_test, catdi_predicted)
            fold_data = { "patient_id":sbj,
                                "fold_number": fold_num,
                                "r_value": r_value,
                                "p_value": p_value,
                                "alpha": model.alpha_,
                                "mse": mse,
                                "r_squared": r2
            }
            performance_metrics.append(fold_data)

        catdi_predict_vals.extend(catdi_predicted)
        catdi_measured_vals.extend(catdi_test)

    decoding_plot_dict = {
        'patient_id': sbj,
        'measured_catdi': catdi_measured_vals,
        'predicted_catdi': catdi_predict_vals
    }
    decoding_plot_info.append(decoding_plot_dict)

    if cv_choice == "LOO":
        r_val, p_val = pearsonr(catdi_measured_vals,catdi_predict_vals)
        mse_val = mean_squared_error(catdi_measured_vals,catdi_predict_vals)
        r2 = r2_score(catdi_measured_vals, catdi_predict_vals)
        performance_metrics.append({"patient_id": sbj,
                                "mse": mse_val,
                                "r_val": r_val,
                                "p_val":p_val,
                                "r_squared": r2,
                                "total_sessions": len(catdi_scores)})
# output stats to csv
output_dir = os.path.join(base_dir,"lasso")
os.makedirs(output_dir, exist_ok=True)
csv = os.path.join(output_dir,f"lasso_manual_power_{cv_choice}.csv")
performance_df = pd.DataFrame(performance_metrics)
performance_df.to_csv(csv, index=False)

# output figure for all patients 
fig_path = os.path.join(output_dir, f"lasso_decoding_results_manual_power_{cv_choice}.png")

fig, axes = plt.subplots(2,3,figsize=(10,6))
axes = axes.flatten()   # 1 through 6 instead of the grid

for index, dict in enumerate(decoding_plot_info):
    # define x and y points
    measured = dict['measured_catdi']
    predicted = dict['predicted_catdi']
    # find the r and p value with pearson correlation
    r_val, p_val = pearsonr(measured, predicted)
    # show y =x line
    min_lim= min((min(measured), min(predicted)))
    max_lim =max((max(measured), max(predicted)))
    axes[index].plot((min_lim, max_lim),(min_lim,max_lim), linestyle='--', color='black')
    # show scatter plot of measured and predicted catdi values
    axes[index].scatter(measured, predicted, color='blue')
    # label R value and P value from pearson correlation
    axes[index].text(.10,.80,f"R={r_val:.2f} \nP={p_val:.4}",transform=axes[index].transAxes)
    axes[index].set_xlabel('Measured CATDI')
    axes[index].set_ylabel('Predicted CATDI')
    axes[index].set_title(f"{dict['patient_id']}")

fig.tight_layout()
fig.savefig(fig_path)
plt.show()


        


    



