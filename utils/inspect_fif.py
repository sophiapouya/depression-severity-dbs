# inspect the fif file
import mne

fif_file = "/Users/sophiapouya/workspace/bcm/CATDI/preprocessing_reworked/preprocessing_results/DBSTRD006/EMU-194_subj-DBSTRD006_task-CATDI_date-02-11-2022_time-08-41-17_clean_annotated_ieeg.fif"
raw = mne.io.read_raw_fif(fif_file, preload=True)
print("finished")