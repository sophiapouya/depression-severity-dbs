import py_neuromodulation as nm
import os
import mne
import pandas as pd
import pprint
from config import BASE_DIR_SEEG, CATDI_SCORES

# define the list of sessions to exclude 
EXCLUDED_SESSIONS = {
    "DBSTRD001": ["CATDI_run-08_blk-04"],
    "DBSTRD002": ["CATDI_run-Day5_blk-04", "CATDI_run-Day7_blk-05", "CATDI_run-Day3_blk-02","CATDI_run-Day3_blk-03"],
    "DBSTRD006": ["CATDI_date-02-08-2022_time-12-42-20"],
    "DBSTRD008": ["CATDI_date-10-25-2022_time-20-20-50","CATDI_date-10-25-2022_time-13-42-28","CATDI_date-10-25-2022_time-11-27-18",
                  "CATDI_date-10-26-2022_time-07-36-57","CATDI_date-10-26-2022_time-16-16-07","CATDI_date-10-26-2022_time-14-58-43",
                  "CATDI_date-10-25-2022_time-16-51-50","CATDI_date-10-25-2022_time-14-50-44","CATDI_date-10-25-2022_time-08-23-49", 
                  "CATDI_date-10-26-2022_time-18-24-14"],
    "DBSTRD010": ["CATDI_date-05-11-2023_time-16-20-04"],
    "DBSTRD011": ["CATDI_date-20240724_time-145011", "CATDI_date-20240719_time-213319","CATDI_date-20240724_time-111946"],
    "DBSTRD014": ["CATDI_date-20250307_time-130552","CATDI_date-20250312_time-183153","CATDI_date-20250307_time-195211",
                  "CATDI_date-20250308_time-084929"]
}

# default template setttings
settings = nm.get_default_settings()

settings["frequency_ranges_hz"] = {
    'delta': [1, 4], 'theta': [4, 8], 'alpha': [8, 12], 
    'beta': [12, 30], 'low_gamma': [35, 50], 'high_gamma': [70, 150]
}

settings.sampling_rate_features_hz = 1.0 # Match your stream call

settings.coherence_settings.nperseg = 512 # Provides better frequency resolution at 1000Hz
settings.coherence_settings.noverlap = 256  # 50% overlap
settings["segment_length_features_ms"] = 2000  # 2 second windows

# disable preprocessing
settings["preprocessing"] = []
settings["features"]["linelength"] = False
settings["features"]["welch"] = False
settings["features"]["adaptive_filter"] = False
settings.postprocessing.feature_normalization = False

# enable specific feature modules
settings["features"]["fft"] = True
settings["features"]["hjorth"] = True
settings["features"]["sharpwave_analysis"] = True
settings["features"]["fooof"] = True 

# FFT settings
settings["fft_settings"]["windowlength_ms"] = 2000

# FOOOF settings
settings["fooof_settings"]["knee"] = False  # This is being ignored due to bug
settings["fooof_settings"]["aperiodic"]["knee"] = False  # Try this too
settings["fooof_settings"]["freq_range_hz"] = [2, 45]
settings["fooof_settings"]["windowlength_ms"] = 2000
settings["fooof_settings"]["max_n_peaks"] = 4

# Burst settings
settings["features"]["bursts"] = False  # leaving false for now b/c don't know how to specify it

# sharpwave settings
settings["sharpwave_analysis_settings"]["filter_ranges_hz"] = [[12, 30], [70,150]]
settings["sharpwave_analysis_settings"]["estimator"]["mean"] = ['interval', 'prominence', 'sharpness']

# coherence 
settings["features"]["coherence"] = False 

# print settings to verify everything is okay
pprint.pprint(settings)

all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
#all_subjs = ["DBSTRD006"]
base_dir = str(BASE_DIR_SEEG)
catdi_scores_excel =str(CATDI_SCORES)
all_session_results = []

for subj in all_subjs:
    sbj_dir = os.path.join(base_dir,subj)
    fif_dir = os.path.join(sbj_dir, "bipolar_channels_greymatter")

    # bring in the catdi scores
    catdi_excel = pd.read_excel(catdi_scores_excel, sheet_name=subj)
    # clean up the names if they have text before CATDI
    if subj in ["DBSTRD014","DBSTRD011"]:
        catdi_excel["session"] = catdi_excel["Name"].str.split("_task-").str[1]
    else:
        catdi_excel["session"] = catdi_excel["Name"]

    for file in os.listdir(fif_dir):
        if file.endswith(".fif"):
            # average session data
            sesh = os.path.splitext(file)[0]
            
            # skip over the session if it's excluded 
            if sesh in EXCLUDED_SESSIONS[subj]:
                continue
            
            score_row = catdi_excel.loc[catdi_excel["session"] == sesh]

            # skip over file if it doesn't have a catdi score
            if score_row.empty:
                print(f"skipping over {file} because no catdi score found in excel")
                continue

            file_path = os.path.join(fif_dir, file)
            raw = mne.io.read_raw_fif(file_path, preload=True, verbose=False)
            data = raw.get_data()
  
            # channel df
            channels_df = pd.DataFrame({
                "name": raw.ch_names, 
                "new_name": raw.ch_names,
                "type": "seeg", 
                "used": 1, 
                "target": 0, 
                "status": "good",
                "rereference": "None"
            })

            # setup stream
            stream = nm.Stream(
                sfreq=raw.info['sfreq'],
                data=data,
                settings=settings,
                sampling_rate_features_hz= 1.0,  # slow down step size to 1hz (since patient average -> one calculation every second
                channels=channels_df,
                verbose=False
            )

            # run extraction 
            df = stream.run()
            session_avg = df.mean(numeric_only=True).to_frame().T

            session_avg['session_name'] = sesh
            session_avg['patient_id'] = subj
            session_avg['catdi_score'] = score_row["Result"].values[0]
            all_session_results.append(session_avg)

master_df = pd.concat(all_session_results, ignore_index=True)

# Check which features have the most NaNs
nan_counts = master_df.isna().sum()
print(nan_counts[nan_counts > 0].sort_values(ascending=False))

# Check which patients have the most NaNs
print(master_df.groupby('patient_id').apply(lambda x: x.isna().sum().sum()))

# reorder columns so the metadata is at the front (easier to read)
id_cols = ['patient_id', 'session_name', 'catdi_score']
other_cols = [c for c in master_df.columns if c not in id_cols]
master_df = master_df[id_cols + other_cols]

# save as pkl and csv
csv_path = os.path.join(base_dir, "CATDI_master_features_gm_bursting.csv")
master_df.to_csv(csv_path, index=False)

pkl_path = os.path.join(base_dir, "CATDI_master_features_gm_bursting.pkl")
master_df.to_pickle(pkl_path)