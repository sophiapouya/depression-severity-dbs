# blackrock files -> .nev and .ns3 (2kHz sampling)
from neo import rawio
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, iirnotch
   
base_directory = '/Users/sophiapouya/workspace/bcm/depression_research/CATDI/neuralData/originalData/'
sfreq = 2000    # sampling frequency
harmonics = np.arange(60, 501, 60)

# generic 4th order butterworth bandpass filter
def bandpass_filter(raw_signal, order= 4, sfreq=sfreq, low_cutoff=0.3, high_cutoff=500):
    nyquist = sfreq / 2
    wn= [low_cutoff/nyquist, high_cutoff/nyquist]
    b, a = butter(N=order, Wn=wn, btype='bandpass')
    bandpassed_signal = filtfilt(b, a, raw_signal)
    return bandpassed_signal

# generic notch filter
def notch_filter(raw_signal, harmonic):
    w0 = harmonic / (sfreq/2)
    b, a = iirnotch(w0=w0, Q=30)
    notch_filtered = filtfilt(b,a,raw_signal)
    return notch_filtered

# perform bandpass and notch filtering (60hz harmonics) all at once
def preprocess_signal(raw_signal, harmonics=harmonics):
    filtered_signal = bandpass_filter(raw_signal=raw_signal)
    for harmonic in harmonics:
        filtered_signal = notch_filter(filtered_signal, harmonic)
    return filtered_signal
    

# for all dbs leads (3)
for dbs_directory in os.listdir(base_directory):
    dbs_path = os.path.join(base_directory,dbs_directory)

    # for all events in the event directory
    for event_directory in os.listdir(dbs_path):
        event_path = os.path.join(dbs_path,event_directory)

        # booleans to keep track of .ns3 and .nev files
        nev_file=False
        ns3_file = False
        file_basename = None

        for file in os.listdir(event_path):
            if file.endswith(".nev"):
                nev_file = True 
                file_basename = os.path.splitext(file)[0]
            elif file.endswith(".ns3"):
                ns3_file = True

        # only process the folder if it has both necessary files
        if nev_file and ns3_file: 
            full_path = os.path.join(event_path,file_basename)
            reader = rawio.BlackrockRawIO(filename=full_path, nsx_to_load=3)
            reader.parse_header()

            nsx_key = list(reader.nsx_datas.keys())[0]
            block_key = list(reader.nsx_datas[nsx_key].keys())[0]

            # nsx_datas shape -> (time points, channels)
            data_length = np.shape(reader.nsx_datas[nsx_key][block_key])[0]
            channels = np.shape(reader.nsx_datas[nsx_key][block_key])[1]

            times = np.arange(data_length) / sfreq
            
            # for all channels in file
            for channel in range(channels):
                
                # for all time points in a given channel, bandpass and notch filter the data
                raw_signal = reader.nsx_datas[nsx_key][block_key][:,channel]                
                filtered_data = preprocess_signal(raw_signal=raw_signal)
                
                # make a plot for each event and channel
                plt.plot(times, filtered_data)
                plt.xlabel("Seconds (s)")
                plt.ylabel("Microvolts (uV)")
                
                plot_dir = os.path.join("/Users/sophiapouya/workspace/bcm/depression_research/filtered_channel_plots", dbs_directory, str(channel+1))
                os.makedirs(plot_dir, exist_ok=True)
                
                plot_name = f"{dbs_directory}_{event_directory}_channel_{channel+1}"
                plot_path = os.path.join(plot_dir,plot_name+".png")
                plt.title(plot_name)
                plt.savefig(plot_path)
                plt.close()

#####################
# Single File Example
#####################

# #contains a .nev and .ns3 file
# test_file_directory = 'CATDI/neuralData/originalData/DBSTRD001/EMU-027_task-CATDI_run-03_blk-02/'

# filename_base = 'EMU-027_subj-DBSTRD001_task-CATDI_run-03_blk-02'
# full_path = os.path.join(test_file_directory,filename_base)
# reader = rawio.BlackrockRawIO(filename=full_path, nsx_to_load=3)
# reader.parse_header()

# #nsx_datas -> 427442 time points
# data_length = np.shape(reader.nsx_datas[3][0])
# data_vec = np.arange(0, data_length[0])
# times = np.arange(0,data_length[0]*.0005, .0005)
# uV_data = []

# #channel 1 data ->
# for data_pt in data_vec:
#     data = reader.nsx_datas[3][0][data_pt][0]
#     uV_data.append(data)

# # make a plot of time vs voltage for this channel during this session
# plt.plot(times,uV_data)
# plt.ylabel("Voltage (uV)")
# plt.xlabel("Seconds (s)")
# plt.title("channel_1")




# required preprocessing: amplification and bandpass filter from .3 to 500hz w/ fourth order butterworth filter
    # additionally: notch filtering (60hz and harmonics)
    # bipolar referencing -> subtracting voltage of neighboring contact on each electrode 
    # Hilbert transform to estimate spectral power features
# spectral power from 6 frequency bands for each channel
# bands-> delta (1-4hz), theta (4-8hz), alpha (8-12hz), beta (12-30hz), gamma (35-50hz), high gamma (70-150hz)
