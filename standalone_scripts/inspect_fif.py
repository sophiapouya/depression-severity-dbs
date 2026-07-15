# inspect the fif file
import mne
import os

PLOTTING_SCALE = 200e-6


# single instances
# fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/DBSTRD014/raw_fif_files/"
# sessions = ["CATDI_date-20250312_time-110624.fif","CATDI_date-20250310_time-174228.fif"]

fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/DBSTRD011/raw_fif_files/"
sessions = ["CATDI_date-20240717_time-193928.fif"]
for session in sessions:
    fif_file = os.path.join(fif_dir, session)
    raw = mne.io.read_raw_fif(fif_file, preload=True)
    raw.plot(block=True,scalings=PLOTTING_SCALE)

# patient level sifting, define patient 
electrode_type = "SEEG"  # choices: DBS or SEEG
sbj = "DBSTRD014"
mode = "raw"    # choices: raw or reref

if mode == "reref":
    if electrode_type == "DBS":
        fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
        fif_patient_dir = os.path.join(fif_dir, sbj,"bipolar_alternating_channels")
    else:
        fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/"
        fif_patient_dir = os.path.join(fif_dir, sbj,"bipolar_channels")
else:
    if electrode_type == "DBS":
        fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
        fif_patient_dir = os.path.join(fif_dir, sbj,"raw_fif_files")
    else:
        fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/"
        fif_patient_dir = os.path.join(fif_dir, sbj,"raw_fif_files")

for file in os.listdir(fif_patient_dir):
    if file.endswith(".fif"):
        fif_file = os.path.join(fif_patient_dir, file)
        raw = mne.io.read_raw_fif(fif_file, preload=True)
        raw.plot(block=True, scalings=PLOTTING_SCALE)
    
    else:
        continue