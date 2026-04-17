# inspect the fif file
import mne
import os

PLOTTING_SCALE = 200e-6


# single instance
# fif_file = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/DBSTRD001/py_neuro_files/CATDI_run-03_blk-01_pynm.fif"
# raw = mne.io.read_raw_fif(fif_file, preload=True)
# raw.plot(block=True,scalings=PLOTTING_SCALE)

# patient level sifting, define patient 
electrode_type = "DBS"  # choices: DBS or SEEG
sbj = "DBSTRD002"

if electrode_type == "DBS":
    fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/"
    fif_patient_dir = os.path.join(fif_dir, sbj,"bipolar_alternating_channels")
else:
    fif_dir = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData/"
    fif_patient_dir = os.path.join(fif_dir, sbj,"bipolar_channels")

for file in os.listdir(fif_patient_dir):
    if file.endswith(".fif"):
        fif_file = os.path.join(fif_patient_dir, file)
        raw = mne.io.read_raw_fif(fif_file, preload=True)
        raw.plot(block=True, scalings=PLOTTING_SCALE)
    
    else:
        continue