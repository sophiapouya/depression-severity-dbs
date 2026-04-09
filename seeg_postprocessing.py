from src.postprocessing_functions import save_power_data_from_fif_seeg
import os
from scipy.signal import butter

base_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/"
all_subjs = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
#all_subjs = ["DBSTRD011"]
REF_TYPE = "bipolar"

# define the list of sessions to exclude 
EXCLUDED_SESSIONS = {
    "DBSTRD001": ["CATDI_run-08_blk-04"],
    
    #"DBSTRD002": ["CATDI_run-Day6_blk-02"],
    "DBSTRD002": ["CATDI_run-Day5_blk-04", "CATDI_run-Day7_blk-05", "CATDI_run-Day3_blk-02","CATDI_run-Day3_blk-03"],
    "DBSTRD006": ["CATDI_date-02-08-2022_time-12-42-20"],
    "DBSTRD008": ["CATDI_date-10-25-2022_time-20-20-50","CATDI_date-10-25-2022_time-13-42-28","CATDI_date-10-25-2022_time-11-27-18",
                  "CATDI_date-10-26-2022_time-07-36-57","CATDI_date-10-26-2022_time-16-16-07","CATDI_date-10-26-2022_time-14-58-43",
                  "CATDI_date-10-25-2022_time-16-51-50","CATDI_date-10-25-2022_time-14-50-44","CATDI_date-10-25-2022_time-08-23-49", 
                  "CATDI_date-10-26-2022_time-18-24-14"],
    "DBSTRD010": ["CATDI_date-05-11-2023_time-16-20-04"],
    "DBSTRD011": ["CATDI_date-20240724_time-145011", "CATDI_date-20240719_time-213319","CATDI_date-20240724_time-111946"],
    "DBSTRD014": ["CATDI_date-20250307_time-130552","CATDI_date-20250312_time-183153","CATDI_date-20250311_time-100834"]
}
# perform hilbert power series in freq bands
FEATURE_BANDS = {
    'delta': [1, 4], 'theta': [4, 8], 'alpha': [8, 12], 
    'beta': [12, 30], 'low_gamma': [35, 50], 'high_gamma': [70, 150]
}
SFREQ = 1000

# precompute butterworth filter coefficients
band_coefficients = {}
for band, values in FEATURE_BANDS.items():
    band_coefficients[band] = butter(N=2, Wn=values, btype='bandpass', fs=SFREQ)

# iterate over all subjects to calculate the power in each frequency band
for sbj in all_subjs:
    sbj_dir = os.path.join(base_dir, sbj)
    save_power_data_from_fif_seeg(band_coefficients=band_coefficients, FEATURE_BANDS=FEATURE_BANDS, EXCLUDED_SESSIONS=EXCLUDED_SESSIONS, subj_name=sbj, sbj_dir=sbj_dir)


