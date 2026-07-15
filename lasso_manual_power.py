

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from src.postprocessing_functions import get_nrmse
from scipy.stats import pearsonr 
import matplotlib.pyplot as plt
import math
from config import CATDI_SCORES, BASE_DIR_DBS


# cv choice
cv_choice="LOO" # choices -> LOO or KFOLD

# define input files
catdi_scores_file = str(CATDI_SCORES)

base_dir = str(BASE_DIR_DBS)
all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]


# list for keeping track of performance for output csv
performance_metrics = []

# list for patient catdi predicted vs catdi measured for plotting
decoding_plot_info = []

l1_reg_list = np.around(np.arange(0.1, 1.1, 0.1), 1)
outlier = 4

# go through each patient
for sbj in all_subjs:
    patient_power_csv = os.path.join(base_dir,sbj,"bipolar_alternating_channels","power_bipolar_alternating",f"{sbj}_bipolar_alternating_power.csv")
    patient_df = pd.read_csv(patient_power_csv)

    # add a column for the catdi scores
    catdi_excel = pd.read_excel(catdi_scores_file, sheet_name=sbj)
    if sbj in ["DBSTRD014", "DBSTRD011"]:
        catdi_excel["Name"]= catdi_excel["Name"].astype(str).str.split("task-").str[-1]
    # create a column of scores catdi_excel["Result"] that matches the patient_power_csv["session"]
    catdi_session_and_scores = dict(zip(catdi_excel["Name"],catdi_excel["Result"]))
    
    patient_df["scores"] = patient_df["session"].map(catdi_session_and_scores)
    # drop any rows without scores
    patient_df=patient_df.dropna(subset=["scores"])

    
    # initialize kfold 
    if cv_choice == "KFOLD":
        folds = KFold(n_splits=5, shuffle=True, random_state=42)
    else:
        folds = LeaveOneOut()
    
    fold_num = 0
    catdi_predict_vals, catdi_measured_vals = [],[]
    unique_sessions = patient_df["session"].unique()
    for train_index, test_index in folds.split(unique_sessions):
        # keep track of the current fold
        fold_num = fold_num+1
        train_sessions = unique_sessions[train_index]
        test_sessions = unique_sessions[test_index]

        data_train = patient_df[patient_df["session"].isin(train_sessions)]
        data_test = patient_df[patient_df["session"].isin(test_sessions)]
        
        data_train_pivot = data_train.pivot_table(
            index="session",
            columns="ch_name",
            values = ["delta","theta","alpha","beta","low_gamma","high_gamma"]
        )
        # use the means of contacts to fill in nans
        train_feature_means =data_train_pivot.mean(axis=0)
        data_train_pivot = data_train_pivot.fillna(train_feature_means)
        train_catdi_scores = data_train.groupby("session")["scores"].first().loc[data_train_pivot.index]

        alpha_nrmse = []
        for l1_reg in l1_reg_list:

            inner_predictions, inner_test_vals =[],[]
            inner_loo = LeaveOneOut()

            for train_idx, test_idx in inner_loo.split(data_train_pivot):
                inner_train_data, inner_test_data = data_train_pivot.iloc[train_idx].values, data_train_pivot.iloc[test_idx].values
                inner_train_scores, inner_test_scores = train_catdi_scores.iloc[train_idx].values, train_catdi_scores.iloc[test_idx].values

                mu, std = np.mean(inner_train_data, axis=0), np.std(inner_train_data, axis=0)
                std[std==0] = 1

                # clamp the training data
                inner_train_data =np.where(np.abs(inner_train_data -mu)> outlier * std, mu, inner_train_data)
                inner_test_data = np.where(np.abs(inner_test_data-mu)>outlier * std, mu, inner_test_data)
                inner_model = make_pipeline(StandardScaler(), Lasso(alpha=l1_reg, random_state=0, max_iter=2000))
                inner_model.fit(inner_train_data, inner_train_scores)

                inner_prediction = inner_model.predict(inner_test_data)
                inner_predictions.extend(inner_prediction)
                inner_test_vals.extend(inner_test_scores)
            
            current_alpha_nrmse = get_nrmse(inner_test_vals, inner_predictions)
            alpha_nrmse.append(current_alpha_nrmse)
        
        best_alpha_idx = np.argmin(alpha_nrmse)
        best_alpha = l1_reg_list[best_alpha_idx]

        data_test_pivot=data_test.pivot_table(
            index="session",
            columns="ch_name",
            values = ["delta","theta","alpha","beta","low_gamma","high_gamma"]
        ).reindex(columns=data_train_pivot.columns)
        
        #fill any nans in the test session with averages from the training data
        data_test_pivot = data_test_pivot.fillna(train_feature_means)
        catdi_scores_outer_train = data_train.groupby("session")["scores"].first().loc[data_train_pivot.index]
        catdi_scores_outer_test = data_test.groupby("session")["scores"].first().loc[data_test_pivot.index]

        data_train_pivot = data_train_pivot.values
        data_test_pivot = data_test_pivot.values

        # outer loop clamping
        mu_outer, std_outer = np.mean(data_train_pivot, axis=0), np.std(data_train_pivot, axis=0)
        std_outer[std_outer==0] = 1

        outer_train_clamped = np.where(np.abs(data_train_pivot - mu_outer)> outlier * std_outer, mu_outer, data_train_pivot)
        outer_test_clamped = np.where(np.abs(data_test_pivot - mu_outer) > outlier * std_outer, mu_outer, data_test_pivot)

        # perform lasso cv regression
        outer_model = make_pipeline(StandardScaler(), Lasso(alpha=best_alpha, random_state=42,max_iter=100000))
        # fit the model 
        outer_model.fit(outer_train_clamped, catdi_scores_outer_train)
        # predict on the test set
        catdi_predicted = outer_model.predict(outer_test_clamped)

        catdi_predict_vals.extend(catdi_predicted)
        catdi_measured_vals.extend(catdi_scores_outer_test)

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
                                "total_sessions": len(unique_sessions)})
# output stats to csv
output_dir = os.path.join(base_dir,"lasso")
os.makedirs(output_dir, exist_ok=True)
csv = os.path.join(output_dir,f"lasso_manual_power_{cv_choice}.csv")
performance_df = pd.DataFrame(performance_metrics)
performance_df.to_csv(csv, index=False)

# output figure for all patients 
fig_path = os.path.join(output_dir, f"lasso_{cv_choice}.png")
n_rows = 2
n_cols = math.ceil(len(all_subjs)/n_rows)

fig, axes = plt.subplots(n_rows, n_cols,figsize=(4 * n_cols, 4 * n_rows),
    constrained_layout=True)
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
    axes[index].text(.10,.80,f"R={r_val:.2f}",fontsize=14, transform=axes[index].transAxes)
    axes[index].set_xlabel('Measured CATDI')
    axes[index].set_ylabel('Predicted CATDI')
    axes[index].set_title(f"{dict['patient_id']}")
    axes[index].set_box_aspect(1)

# remove unused subplots
for i in range(len(decoding_plot_info), len(axes)):
    fig.delaxes(axes[i])

fig.tight_layout()
fig.savefig(fig_path)
plt.show()


        


    



