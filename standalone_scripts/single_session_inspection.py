import os
from src.preprocessing_functions import *
# 191
# 194
# 395

# assumes preprocessed .fif file (bandpass and notch filtered)

input_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
sbj_name = "DBSTRD006"
session_name = "CATDI_date-02-14-2022_time-15-59-03"
fif_file = os.path.join(input_dir,sbj_name,"rerefData/",session_name,f"fiEEG_dbs_{session_name}.fif")

raw_file = mne.io.read_raw_fif(fif_file)
raw_file.plot(block=True)