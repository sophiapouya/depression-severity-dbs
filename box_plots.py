import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import os

input_folder = Path("/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/pca/")
output_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/pca/"

pca_list = []
for file in input_folder.glob("*.csv"):
    # files named like -> OLS_LOO_LEFT_probes.csv
    name_parts = file.stem.split("_")
    regression_method= name_parts[0]
    cv_method = name_parts[1]
    probe = name_parts[2]

    # read in the csv
    csv_df = pd.read_csv(file)

    if (cv_method != "LOO") or (regression_method != "OLS"):
        continue
    # append the data to the pca_list as a dictionary
    for index, row in csv_df.iterrows():
        csv_dict = {
            'patient': row['patient_id'],
            'probe': probe.lower(),
            'r':  row['r_val'],
            'p': row['p_val']
        }
        pca_list.append(csv_dict)

pca_df = pd.DataFrame(pca_list)
patient_r_df = pca_df.pivot(index="patient", columns="probe", values="r" )
patient_p_df = pca_df.pivot(index="patient", columns="probe", values="p" )
desired_patient_order = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008","DBSTRD010", "DBSTRD014"]
patient_r_df = patient_r_df.loc[desired_patient_order,:]
patient_p_df = patient_p_df.loc[desired_patient_order,:]

# # create comparative box plot w/ r-values connected b/w patients 
# def box_plot_compare(df, columns, output_dir):
#     cat_1 = columns[0]
#     cat_2 = columns[1]
#     output_file = os.path.join(output_dir, f"boxplot_{cat_1}_vs_{cat_2}.png")
#     df_col1 = df[cat_1]
#     df_col2 = df[cat_2]

#     plt.boxplot([df_col1,df_col2] ,labels=[columns[0], columns[1]])

#     # connect paired points
#     for patient in df.index:
#         plt.plot([1,2], [df_col1.loc[patient], df_col2[patient]])
#         plt.text(1.05, df_col1.loc[patient], patient )

#     plt.tight_layout()
#     plt.savefig(output_file)
#     plt.show()


def box_plot_all(r_df, p_df, output_dir):
    output_file = os.path.join(output_dir, "ALL_boxplots.png")
    expected_col_order = ["all","vcvs","lvcvs","rvcvs","scc", "lscc", "rscc"]
    r_df = r_df[expected_col_order]
    p_df = p_df[expected_col_order]

    patient_colormap={
        "DBSTRD001": 'red',
        "DBSTRD002":'gray',
        "DBSTRD006": 'orange',
        "DBSTRD008": 'pink',
        "DBSTRD010": 'green',
        "DBSTRD014": 'blue',
    }
    
    plt.figure(figsize=(15,10))

    ordered_cols, ordered_p_vals = [], []
    for feature_num, col in enumerate(r_df.columns):
        ordered_col = r_df[col].sort_values()
        ordered_p_val = p_df[col].loc[ordered_col.index]
        count = 1
        for index, value in enumerate(ordered_col):
            patient_label = ordered_col.index[index]
            if count % 2 == 0:
                text_offset= 0.55
            else:
                text_offset = 1.05
            # plot the patient labels
            plt.text(feature_num+text_offset, value, patient_label)
            if ordered_p_val.iloc[index] < 0.05:
                plt.plot(feature_num+1, value, marker="*", markersize=10, color=patient_colormap[patient_label])
            else: 
                plt.plot(feature_num+1, value, marker="o", markersize=4, color="black")
            count = count+1
        ordered_cols.append(ordered_col)
        ordered_p_vals.append(ordered_p_val)

    plt.boxplot(ordered_cols,labels=r_df.columns)
    plt.xlabel("Region features used to decode", fontsize=16)
    plt.ylabel("R-values", fontsize=16)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()
    plt.close()
                
        
#box_plot_all(patient_r_df, patient_p_df, output_dir)

# create the box plots
# box_plot_compare(r_df=patient_r_df, p_df=patient_p_df, columns=["rscc", "rvcvs"], output_dir=output_dir)






from scipy.stats import ttest_rel

def run_paired_tests(df):
    results = []

    cols = df.columns

    for i in range(len(cols)):
        for j in range(i+1, len(cols)):
            col1 = cols[i]
            col2 = cols[j]

            # drop NaNs pairwise
            paired = df[[col1, col2]].dropna()

            if len(paired) < 3:
                continue

            t_stat, p_val = ttest_rel(paired[col1], paired[col2])

            results.append({
                "comparison": f"{col1} vs {col2}",
                "mean_diff": (paired[col1] - paired[col2]).mean(),
                "t_stat": t_stat,
                "p_val": p_val,
                "n": len(paired)
            })

    return pd.DataFrame(results)

ttest_results = run_paired_tests(patient_r_df)
print(ttest_results)