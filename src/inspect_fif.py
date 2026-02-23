# inspect the fif file
import mne

fif_file = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/dbsData/DBSTRD001/py_neuro_files/CATDI_run-03_blk-01_pynm.fif"
raw = mne.io.read_raw_fif(fif_file, preload=True)
print("finished")