# Save off the sessions with specific referencing in a way that py neurmodulation expects
import os
from utils.preprocessing_functions import save_py_neuromod_chans

base_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
subj_name = "DBSTRD001"
all_subjs = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
ref_type = "bipolar_alternating"
ref_types = ["bipolar","car","esr","bipolar_alternating"]

EXCLUDED_SESSIONS = {
    "DBSTRD001": ["CATDI_run-08_blk-04", "CATDI_run-09_blk-01", "CATDI_run-08_blk-01"],
    "DBSTRD002": ["CATDI_run-Day6_blk-02", "CATDI_run-Day7_blk-02"],
    "DBSTRD006": ["CATDI_date-02-14-2022_time-08-34-02"],
    "DBSTRD008": ["CATDI_date-10-25-2022_time-14-50-44", "CATDI_date-10-26-2022_time-14-58-43", "CATDI_date-10-25-2022_time-08-23-49", 
                  "CATDI_date-10-24-2022_time-14-42-53", "CATDI_date-10-24-2022_time-12-27-04", "CATDI_date-10-25-2022_time-20-20-50", 
                  "CATDI_date-10-25-2022_time-16-51-50", "CATDI_date-10-24-2022_time-17-01-51", "CATDI_date-10-25-2022_time-13-42-28",
                  "CATDI_date-10-25-2022_time-11-27-18", "CATDI_date-10-24-2022_time-13-50-41", "CATDI_date-10-24-2022_time-10-57-12" ],
    "DBSTRD010": [],
    "DBSTRD011": ["CATDI_date-20240717_time-135720", "CATDI_date-20240720_time-121427", "CATDI_date-20240723_time-183759"],
    "DBSTRD014": []
}

for subj in all_subjs:
    sbj_dir = os.path.join(base_dir, subj)
    save_py_neuromod_chans(ref_type=ref_type, sbj_dir=sbj_dir, subj_name=subj, EXCLUDED_SESSIONS=EXCLUDED_SESSIONS)