

import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.linear_model import LinearRegression, RidgeCV
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr 
import matplotlib.pyplot as plt
import math


probe_type = "DBS"      # choices: "DBS", "SEEG"

# regression model
model_choice = "OLS"    # choices: "RIDGE", "OLS"

# cv method
cv_choice = "LOO"     # choices: "KFOLD", "LOO"

# # exclude probes
# included_probes = "ALL"   # choices: "LEFT", "RIGHT", "ALL", "LSCC", "RSCC", "RVCVS", "LVCVS", "SCC", "VCVS"
if probe_type == "DBS":
    base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData'
    #probe_options = ["LEFT", "RIGHT", "ALL", "LSCC", "RSCC", "RVCVS", "LVCVS", "SCC", "VCVS"]
    probe_options= ["ALL"]
else:
    base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData'
    probe_options = ["ALL"]

# define input csv
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")
all_patient_df = pd.read_csv(csv_path)


for included_probes in probe_options:

    # remove columns that aren't features 
    columns_to_exclude = ['time']
    clean_df =all_patient_df.drop(columns_to_exclude, axis=1)

    # decode on specific regions
    if included_probes == "ALL":
        clean_df=clean_df
    elif included_probes == "LEFT":
        columns_to_drop = [col for col in clean_df.columns if ("RS" in col) or ("RV" in col)]
        clean_df = clean_df.drop(columns=columns_to_drop)
        print(f"number of columns for {included_probes} = {len(clean_df.columns)}")
    elif included_probes == "RIGHT":
        columns_to_drop = [col for col in clean_df.columns if ("LS" in col) or ("LV" in col)]
        clean_df = clean_df.drop(columns=columns_to_drop)
        print(f"number of columns for {included_probes} = {len(clean_df.columns)}")
    elif included_probes == "SCC":
        columns_to_drop = [col for col in clean_df.columns if ("VCVS" in col)]
        clean_df = clean_df.drop(columns=columns_to_drop)
        print(f"number of columns for {included_probes} = {len(clean_df.columns)}")
    elif included_probes == "VCVS":
        columns_to_drop = [col for col in clean_df.columns if ("SCC" in col)]
        clean_df = clean_df.drop(columns=columns_to_drop)
        print(f"number of columns for {included_probes} = {len(clean_df.columns)}")
    else:
        columns_to_keep = [col for col in clean_df.columns if included_probes in col]
        columns_to_keep.append("catdi_score")
        columns_to_keep.append("patient_id")
        columns_to_keep.append("session_name")
        clean_df = clean_df[columns_to_keep]
        print(f"number of columns for {included_probes} = {len(clean_df.columns)}")

    # list for keeping track of performance for output csv
    performance_metrics = []

    # list for patient catdi predicted vs catdi measured for plotting
    decoding_plot_info = []

    # apply the 70% rule
    # go through each patient
    for patient in clean_df['patient_id'].unique():
        patient_df = clean_df[clean_df['patient_id'] == patient]
        session_names = patient_df["session_name"]
        catdi_scores = patient_df['catdi_score']

        # remove catdi and patient id columns 
        patient_df = patient_df.drop(['patient_id','catdi_score','session_name'], axis=1)
        
        # implement 70% cutoff-> 70% of sessions need a value for a feature to be analyzed
        sessions= patient_df.shape[0]
        cutoff = int(np.ceil(.70*sessions))
        boolean_df = patient_df.notna()
        sums = boolean_df.sum()
        columns_to_keep = sums[sums >= cutoff].index
        patient_df = patient_df[columns_to_keep]

        # initialize kfold 
        if cv_choice == "KFOLD":
            folds = KFold(n_splits=5, shuffle=True, random_state=42)
        else:
            folds = LeaveOneOut()

        fold_num = 0
        catdi_predict_vals, catdi_measured_vals = [],[]
        for train_index, test_index in folds.split(patient_df):
            # keep track of the current fold
            fold_num = fold_num+1
            data_train, data_test = patient_df.iloc[train_index,:], patient_df.iloc[test_index,:]
            catdi_train, catdi_test = catdi_scores.iloc[train_index], catdi_scores.iloc[test_index]
            
            # fill in the remaining gaps with average values from columns of training data ONLY since pca and ols/ridge cannot deal with blank values 
            train_mean = data_train.mean()
            data_train = data_train.fillna(train_mean)
            data_test = data_test.fillna(train_mean)

            # standardize the training data
            scaler = StandardScaler()
            # get the mean and std dev of the training data
            scaler.fit(data_train)
            # z score the training data
            data_train_transformed = scaler.transform(data_train)
            # z score the testing data w/ mean and std dev from the training data
            data_test_transformed= scaler.transform(data_test)

            # perform pca -> as many pcs to explain 90% variance in the data
            pca = PCA(0.90)
            # fit the pca to the training data only
            pca.fit(data_train_transformed)
            pc_scores_train = pca.transform(data_train_transformed)
            pc_scores_test = pca.transform(data_test_transformed)
            
            if model_choice == "OLS":
                model = LinearRegression()
            elif model_choice == "RIDGE":
                model = RidgeCV(alphas=np.logspace(-4, 4, 50), cv=5)
            else:
                print("need to specify a model")
                break
            # fit the model 
            model.fit(pc_scores_train, catdi_train)

            # predict on the test set
            catdi_predicted = model.predict(pc_scores_test)

            if cv_choice == "KFOLD":
                # pearson correlation
                r_value, p_value = pearsonr(catdi_test, catdi_predicted)
                # mean squared error
                mse = mean_squared_error(catdi_test, catdi_predicted)
                # r-squared 
                r2 = r2_score(catdi_test, catdi_predicted)
                fold_data = { "patient_id":patient,
                                    "fold_number": fold_num,
                                    "number of pcs": pca.n_components_,
                                    "r_value": r_value,
                                    "p_value": p_value,
                                    "alpha": max(abs(model.coef_)) if model_choice == "RIDGE" else None,
                                    "mse": mse,
                                    "r_squared": r2,
                                    "total_sessions": len(catdi_scores)
                }
                performance_metrics.append(fold_data)
            
            catdi_predict_vals.extend(catdi_predicted)
            catdi_measured_vals.extend(catdi_test)

        decoding_plot_dict = {
            'patient_id': patient,
            'measured_catdi': catdi_measured_vals,
            'predicted_catdi': catdi_predict_vals
        }
        decoding_plot_info.append(decoding_plot_dict)

        if cv_choice == "LOO":
            r_val, p_val = pearsonr(catdi_measured_vals,catdi_predict_vals)
            mse_val = mean_squared_error(catdi_measured_vals,catdi_predict_vals)
            r2 = r2_score(catdi_measured_vals, catdi_predict_vals)
            performance_metrics.append({"patient_id": patient,
                                    "mse": mse_val,
                                    "r_val": r_val,
                                    "p_val":p_val,
                                    "r2": r2,
                                    "pcs": pca.n_components_,
                                    "total_sessions": len(catdi_scores)})

    # output stats to csv
    output_dir = os.path.join(base_dir,"pca")
    os.makedirs(output_dir, exist_ok=True)
    csv = os.path.join(output_dir,f"{model_choice}_{cv_choice}_{included_probes}_probes.csv")
    performance_df = pd.DataFrame(performance_metrics)
    performance_df.to_csv(csv, index=False)

    # output figure for all patients 
    fig_path = os.path.join(output_dir, f"{model_choice}_{cv_choice}_{included_probes}_probes_decoding_results.png")

    n_rows = 2
    n_cols = math.ceil(len(clean_df['patient_id'].unique())/2)

    fig, axes = plt.subplots(n_rows,n_cols,figsize=(4 * n_cols, 4 * n_rows),constrained_layout=True)
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
        axes[index].text(.10,.80,f"R={r_val:.2f}",fontsize=14,transform=axes[index].transAxes)
        axes[index].set_xlabel('Measured CATDI')
        axes[index].set_ylabel('Predicted CATDI')
        axes[index].set_title(f'{dict['patient_id']}')


    # remove unused subplots
    for i in range(len(decoding_plot_info), len(axes)):
        fig.delaxes(axes[i])

    fig.tight_layout()
    fig.savefig(fig_path)
    plt.close()


            


        



