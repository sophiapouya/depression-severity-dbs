

import pandas as pd
import numpy as np
import os
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.linear_model import LassoCV
from sklearn.metrics import mean_squared_error, r2_score
from scipy.stats import pearsonr 
import matplotlib.pyplot as plt


# define input csv
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData'
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")
all_patient_df = pd.read_csv(csv_path)

# remove columns that aren't features 
columns_to_exclude = ['session_name', 'time']
clean_df =all_patient_df.drop(columns_to_exclude, axis=1)

# list for keeping track of performance for output csv
performance_metrics = []

# list for patient catdi predicted vs catdi measured for plotting
decoding_plot_info = []

# apply the 70% rule
# go through each patient
for patient in clean_df['patient_id'].unique():
    patient_df = clean_df[clean_df['patient_id'] == patient]
    catdi_scores = patient_df['catdi_score']

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

    # initialize kfold 
    folds = KFold(n_splits=5, shuffle=True, random_state=42)

    fold_num = 0
    catdi_predict_vals, catdi_measured_vals = [],[]
    for train_index, test_index in folds.split(patient_df):
        # keep track of the current fold
        fold_num = fold_num+1
        data_train, data_test = patient_df.iloc[train_index,:], patient_df.iloc[test_index,:]
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
        model = LassoCV(cv=5, random_state=42,max_iter=100000)
        # fit the model 
        model.fit(data_train_transformed, catdi_train)

        # predict on the test set
        catdi_predicted = model.predict(data_test_transformed)

        # pearson correlation
        r_value, p_value = pearsonr(catdi_test, catdi_predicted)

        # mean squared error
        mse = mean_squared_error(catdi_test, catdi_predicted)

        # r-squared 
        r2 = r2_score(catdi_test, catdi_predicted)

        fold_data = { "patient_id":patient,
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
        'patient_id': patient,
        'measured_catdi': catdi_measured_vals,
        'predicted_catdi': catdi_predict_vals
    }
    decoding_plot_info.append(decoding_plot_dict)

# output stats to csv
output_dir = os.path.join(base_dir,"lasso")
os.makedirs(output_dir, exist_ok=True)
csv = os.path.join(output_dir,"lasso.csv")
performance_df = pd.DataFrame(performance_metrics)
performance_df.to_csv(csv, index=False)

# output figure for all patients 
fig_path = os.path.join(output_dir, "lasso_decoding_results.png")

fig, axes = plt.subplots(3,2,figsize=(15,10))
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
    axes[index].set_title(f'Lasso on Features for {dict['patient_id']}')

fig.tight_layout()
fig.savefig(fig_path)
plt.show()


        


    



