import py_neuromodulation as nm
import os
import mne
import pandas as pd
import pprint
import numpy as np


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

# enable specific feature modules
settings["features"]["fft"] = True
settings["features"]["hjorth"] = True
settings["features"]["sharpwave_analysis"] = True
settings["features"]["coherence"] = True
settings["features"]["fooof"] = True 

# FOOOF settings
settings["fooof_settings"]["knee"] = False  # This is being ignored due to bug
settings["fooof_settings"]["aperiodic"]["knee"] = False  # Try this too
settings["fooof_settings"]["freq_range_hz"] = [2, 45]
settings["fooof_settings"]["windowlength_ms"] = 2000
settings["fooof_settings"]["max_n_peaks"] = 4

# Burst settings
settings["features"]["bursts"] = False  # leaving false for now b/c don't know how to specify it
# settings["bursts_settings"]["threshold"] = 75  # Amplitude threshold (percentile)
# settings["bursts_settings"]["time_duration_s"] = [0.1, 1.0]  # Min and max burst duration
# settings["burst_settings"]["frequency_bands"] = ['theta', 'alpha', 'beta']  # Which bands to analyze

# sharpwave settings
settings["sharpwave_analysis_settings"]["filter_ranges_hz"] = [[12, 30]]
settings["sharpwave_analysis_settings"]["estimator"]["mean"] = ['interval', 'prominence', 'sharpness']

# coherence   
settings["coherence_settings"]["frequency_bands"] = ['theta', 'alpha', 'beta']
settings["coherence_settings"]["channels"] = [
    ["LSCC_1", "LVCVS_1"], # left hemisphere, bp contact 1
    ["LSCC_2", "LVCVS_2"], # left hemisphere, bp contact 2
    ["LSCC_3", "LVCVS_3"], # left hemisphere, bp contact 3
    ["RSCC_1", "RVCVS_1"], # right hemisphere, bp contact 1
    ["RSCC_2", "RVCVS_2"], # right hemisphere, bp contact 2
    ["RSCC_3", "RVCVS_3"], # right hemisphere, bp contact 3
]

# print settings to verify everything is okay
pprint.pprint(settings)

all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD014"]
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData'
catdi_scores_excel = "/Users/sophiapouya/workspace/bcm/CATDI/CATDI_scores.xlsx"
all_session_results = []

for subj in all_subjs:
    sbj_dir = os.path.join(base_dir,subj)
    fif_dir = os.path.join(sbj_dir, "bipolar_alternating_channels")
    
    # all features for one patient
    subj_features = []

    # bring in the catdi scores
    catdi_excel = pd.read_excel(catdi_scores_excel, sheet_name=subj)
    # clean up the names if they have text before CATDI
    if subj =="DBSTRD014":
        catdi_excel["session"] = catdi_excel["Name"].str.split("_task-").str[1]
    else:
        catdi_excel["session"] = catdi_excel["Name"]

    for file in os.listdir(fif_dir):
        if file.endswith(".fif"):
            # average session data
            sesh = os.path.splitext(file)[0]
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
                "type": "dbs", 
                "used": 1, 
                "target": 0, 
                "status": "good",
                "rereference": "None"
            })
            settings.postprocessing.feature_normalization = False
            # build file-specific coherence pairs + safe nperseg ---
            file_settings = settings.model_copy(deep=True)
            requested_pairs = file_settings.coherence_settings.channels
            valid_pairs = [[a, b] for a, b in requested_pairs if a in raw.ch_names and b in raw.ch_names]

            if len(valid_pairs) == 0:
                file_settings.features.coherence = False
            else:
                file_settings.coherence_settings.channels = valid_pairs

            # run the analysis separately on segments if the data has artifacts
            annotations = raw.annotations
            if len(annotations.onset) > 0:
                # get the bad intervals
                bads = []
                for onset, duration in zip(annotations.onset, annotations.duration):
                    bads.append((onset,onset+duration))
                
                # sort pairs chronologically by onset time
                bads = sorted(bads, key=lambda x: x[0])
                # good intervals of time to do analysis with
                goods = []
                total_time = raw.times[-1]
                current_time = 0.0
                for bad_start, bad_end in bads:
                    if bad_start > current_time:
                        goods.append((current_time, bad_start))
                    current_time = bad_end
                if current_time < total_time:
                    goods.append((current_time, total_time))  
                
                segment_dfs = []
                for start_time, end_time in goods:

                    # crop the actual data
                    raw_segment = raw.copy().crop(tmin=start_time, tmax=end_time)
                    
                    recording_duration_s = raw_segment.n_times/raw_segment.info['sfreq']
                    # check if the segment meets the requirements
                    if (recording_duration_s < 2.5):
                        continue
                    
                    data_segment=raw_segment.get_data()

                    # setup stream
                    stream = nm.Stream(
                        sfreq=raw_segment.info['sfreq'],
                        data=data_segment,
                        settings=file_settings,
                        sampling_rate_features_hz= 1.0,  # slow down step size to 1hz (since patient average -> one calculation every second
                        channels=channels_df,
                        verbose=False
                    )

                    # run extraction 
                    seg_df = stream.run()
                    segment_dfs.append(seg_df)
                    # testing
                    print("seg_df shape:", seg_df.shape)
                    print("nan frac:", seg_df.isna().mean().mean())
                all_windows_df = pd.concat(segment_dfs, ignore_index=True)
                session_avg=all_windows_df.mean(numeric_only=True).to_frame().T

                #testing
                print("all_windows_df shape:", all_windows_df.shape)
                print("session nan frac:", all_windows_df.isna().mean().sort_values(ascending=False).head(10))


            else:       # run everything normally if there are no artifacts that segment the time series
                recording_duration_s = data.shape[1] / raw.info['sfreq']
                # setup stream
                stream = nm.Stream(
                    sfreq=raw.info['sfreq'],
                    data=data,
                    settings=file_settings,
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
csv_path = os.path.join(base_dir, "CATDI_master_features.csv")
master_df.to_csv(csv_path, index=False)

pkl_path = os.path.join(base_dir, "CATDI_master_features.pkl")
master_df.to_pickle(pkl_path)