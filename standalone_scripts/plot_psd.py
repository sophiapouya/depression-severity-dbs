import mne
import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go

all_subjs = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
#all_subjs= ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD014"]
base_dir = '/Users/sophiapouya/workspace/bcm/CATDI/neuralData/seegData'

# psd settings
fmin = 1
fmax = 150
n_fft = 2048

for subj in all_subjs:
    fif_dir = os.path.join(base_dir,subj,"bipolar_channels")
    
    # plotly object
    fig = go.Figure()

    for file in os.listdir(fif_dir):
        if not file.endswith(".fif"):
            continue
        fif_path = os.path.join(fif_dir,file)
        raw_obj = mne.io.read_raw_fif(fif_path, preload=True, verbose=False)

        # compute the psd
        psd = raw_obj.compute_psd(method="welch",fmin=1, fmax=150, n_fft=2048, verbose=False)
        psd_data = psd.get_data()   # n_channels x freqs
        # average the frequencies across channels so there's only one curve per session
        mean_psd = psd_data.mean(axis=0)
        # convert to dB for better understanding
        mean_psd_db = 10 * np.log10(mean_psd)
        
        # plot the data
        # plt.plot(psd.freqs, mean_psd_db, alpha=0.6)

        fig.add_trace(
            go.Scatter(x=psd.freqs, y=mean_psd_db, mode="lines", name=file.replace(".fif","").split("TDI_")[-1])
        )
    
    fig.update_layout(
        title=f"PSD for {subj}", xaxis_title="Frequency (Hz)", yaxis_title="Power (dB)")
    
    fig.show()

    # plot for all sessions with one patient
    # plt.title(f"PSD for {subj}")
    # plt.xlabel("Frequency (Hz)")
    # plt.ylabel("Power (dB)")
    # plt.show(block=True)
