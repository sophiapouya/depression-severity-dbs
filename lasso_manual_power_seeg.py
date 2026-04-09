

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.linear_model import Lasso
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.pipeline import make_pipeline
from scipy.stats import pearsonr
from src.postprocessing_functions import get_nrmse
import matplotlib.pyplot as plt
import math

# cv choice
cv_choice="LOO" # choices -> LOO or KFOLD

# define input files
label = "new_labels_all_best_region_including_none_contacts"
catdi_scores_file = '/Users/sophiapouya/workspace/bcm/CATDI/CATDI_scores.xlsx'
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/'
all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
#all_subjs= ["DBSTRD001"]
#regions = ["acc"]
regions= ["acc","vmpfc","dlpfc","ofc"]

# list for keeping track of performance for output csv
performance_metrics = []

# list for patient catdi predicted vs catdi measured for plotting
decoding_plot_info = []

l1_reg_list = np.around(np.arange(0.1, 1.1, 0.1), 1)
outlier = 4

# go through each patient
for sbj in all_subjs:
    patient_power_csv = os.path.join(base_dir,sbj,"bipolar_channels","power_bipolar",f"{sbj}_bipolar_power.csv")
    patient_df = pd.read_csv(patient_power_csv)

    # add a column for the catdi scores
    catdi_excel = pd.read_excel(catdi_scores_file, sheet_name=sbj)
    if sbj in ["DBSTRD014","DBSTRD011"]:
        catdi_excel["Name"]= catdi_excel["Name"].astype(str).str.split("task-").str[-1]
    # create a column of scores catdi_excel["Result"] that matches the patient_power_csv["session"]
    catdi_session_and_scores = dict(zip(catdi_excel["Name"],catdi_excel["Result"]))
    patient_df["scores"] = patient_df["session"].map(catdi_session_and_scores)
    # drop any session that don't have catdi scores
    patient_df = patient_df.dropna(subset=["scores"])
    # # residualize the catdi scores
    # m, b, timedelta_name_map = get_residual_line(score_csv=catdi_scores_file, patient_id=sbj)

    # # add a timedelta columns to the patient_df
    # patient_df["timedeltas"] = patient_df["session"].map(timedelta_name_map)

    # # residualize the scores
    # patient_df["scores"] = patient_df["scores"] - (m * patient_df["timedeltas"] + b)

    # add the region metadata to the patient's power value dataframe
    channel_metadata_file = os.path.join(base_dir, sbj, f"{sbj}_metadata.csv")
    channel_region_df = pd.read_csv(channel_metadata_file)
    channel_region_df["ch_name"] = channel_region_df["channel_name"]
    patient_df = pd.merge(patient_df, channel_region_df[["ch_name","channel_region"]], how= "left", on="ch_name")
    
    # initialize kfold 
    if cv_choice == "KFOLD":
        folds = KFold(n_splits=5, shuffle=True, random_state=42)
    else:
        folds = LeaveOneOut()

    fold_num = 0
    catdi_predict_vals, catdi_measured_vals = [],[]

    unique_sessions = patient_df["session"].unique()
    print(f"num of unique sessions = {len(unique_sessions)}")

    for train_index, test_index in folds.split(unique_sessions):
        # keep track of the current fold
        fold_num = fold_num+1
        # get the session names for training and testing
        train_sessions = unique_sessions[train_index]
        test_sessions = unique_sessions[test_index]

        data_train, data_test = patient_df[patient_df["session"].isin(train_sessions)], patient_df[patient_df["session"].isin(test_sessions)]
        
        # inner CV to select best region with training data only
        region_scores = {}
        region_alphas = {}
        for region in regions:
            inner_folds = LeaveOneOut()
            region_df = data_train[data_train["channel_region"]==region]
            region_df_pivot = region_df.pivot_table(
                index="session",
                columns="ch_name",
                values=["delta","theta","alpha","beta","low_gamma","high_gamma"]
            )
            train_scores = region_df.groupby("session")["scores"].first().loc[region_df_pivot.index]
            
            alpha_nrmse = []
            for l1_reg in l1_reg_list:
                inner_loo = LeaveOneOut()
                inner_predictions,inner_test_vals = [], []

                for train_idx, test_idx in inner_folds.split(region_df_pivot):
                    inner_train_data, inner_test_data = region_df_pivot.iloc[train_idx].values, region_df_pivot.iloc[test_idx].values
                    inner_train_scores, inner_test_scores = train_scores.iloc[train_idx].values, train_scores.iloc[test_idx].values 
                    
                    mu = np.mean(inner_train_data, axis=0)
                    std = np.std(inner_train_data, axis=0)

                    std[std==0]=1

                    # clamp the training data
                    inner_train_data = np.where(np.abs(inner_train_data - mu)> outlier * std, mu, inner_train_data)
                    inner_test_data = np.where(np.abs(inner_test_data - mu) > outlier * std, mu, inner_test_data)

                    inner_model = make_pipeline(StandardScaler(), Lasso(alpha=l1_reg, random_state=0, max_iter=2000))
                    inner_model.fit(inner_train_data, inner_train_scores)

                    # predict with the model
                    inner_prediction = inner_model.predict(inner_test_data)
                    inner_predictions.extend(inner_prediction)
                    inner_test_vals.extend(inner_test_scores)

                current_alpha_nrmse = get_nrmse(inner_test_vals, inner_predictions)
                alpha_nrmse.append(current_alpha_nrmse)
                
            # get the region's nrmse
            best_alpha_idx = np.argmin(alpha_nrmse)
            region_scores[region] = alpha_nrmse[best_alpha_idx]
            region_alphas[region] = l1_reg_list[best_alpha_idx]
        
        # select the best region
        best_region = min(region_scores, key=region_scores.get)
        best_alpha = region_alphas[best_region]

        print(f"\n{sbj} fold {fold_num}")
        print("test_sessions:", list(test_sessions))
        print("best_region:", best_region)

        print(
            "test rows overall:",
            len(data_test),
            "test rows in best region:",
            len(data_test[data_test['channel_region'] == best_region])
        )

        print(
            "test sessions overall:",
            data_test['session'].nunique(),
            "test sessions in best region:",
            data_test[data_test['channel_region'] == best_region]['session'].nunique()
        )
        
        # training data
        train_best_region = data_train[data_train["channel_region"]==best_region]

        pivot_train_best_region = train_best_region.pivot_table(
            index="session",
            columns="ch_name",
            values =["delta","theta","alpha","beta","low_gamma","high_gamma"])
        
        data_test_best_region = data_test[data_test["channel_region"] == best_region]
        data_test_pivot = data_test_best_region.pivot_table(
            index="session",
            columns="ch_name",
            values=["delta","theta","alpha","beta","low_gamma","high_gamma"]
        ).reindex(columns=pivot_train_best_region.columns)


        catdi_scores_train = train_best_region.groupby("session")["scores"].first()
        catdi_scores_best_region_train = catdi_scores_train.loc[pivot_train_best_region.index]
        catdi_scores_test = data_test_best_region.groupby("session")["scores"].first()
        catdi_scores_test = catdi_scores_test.loc[data_test_pivot.index]

        # final outlier clamping
        pivot_train_best_region = pivot_train_best_region.values
        data_test_pivot = data_test_pivot.values

        mu_outer = np.mean(pivot_train_best_region, axis=0)
        std_outer = np.std(pivot_train_best_region, axis=0)
        std_outer[std_outer==0] = 1

        outer_train_clamped = np.where(np.abs(pivot_train_best_region - mu_outer)> outlier * std_outer, mu_outer, pivot_train_best_region)
        outer_test_clamped = np.where(np.abs(data_test_pivot - mu_outer) > outlier * std_outer, mu_outer, data_test_pivot)

        final_model = make_pipeline(StandardScaler(), Lasso(alpha=best_alpha, random_state=42, max_iter=10000))
        final_model.fit(outer_train_clamped, catdi_scores_best_region_train)
        catdi_predicted = final_model.predict(outer_test_clamped)

        catdi_predict_vals.extend(catdi_predicted)
        catdi_measured_vals.extend(catdi_scores_test)

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
                                "nrmse":  region_scores[best_region],
                                "r_val": r_val,
                                "p_val":p_val,
                                "r_squared": r2, 
                                "alpha": best_alpha,
                                "total_sessions": len(unique_sessions),
                                "best_region": best_region})
# output stats to csv
output_dir = os.path.join(base_dir,"lasso_test")
os.makedirs(output_dir, exist_ok=True)
csv = os.path.join(output_dir,f"{label}_lasso_manual_power_{cv_choice}.csv")
performance_df = pd.DataFrame(performance_metrics)
performance_df.to_csv(csv, index=False)

# output figure for all patients 
fig_path = os.path.join(output_dir, f"{label}_lasso_decoding_results_manual_power_{cv_choice}.png")

n_rows = 2
n_cols = math.ceil(len(all_subjs)/n_rows)

fig, axes = plt.subplots(n_rows, n_cols, figsize=(10,6))
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

# remove unused subplots
for i in range(len(decoding_plot_info), len(axes)):
    fig.delaxes(axes[i])

fig.tight_layout()
fig.savefig(fig_path)
plt.show()


        


    



