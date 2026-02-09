from utils.preprocessing_functions import *
import json

# paths
DATA_ROOT = "/Users/sophiapouya/workspace/bcm"
PROJECT_NAME = "CATDI"
SBJ_NAME = "DBSTRD001" # Subject to process

ORIGINAL_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'originalData', SBJ_NAME)
DBS_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'dbsData', SBJ_NAME)
REREF_DATA_ROOT = os.path.join(DBS_DATA_ROOT, 'rerefData')

# channel patterns
DBS_CHANNEL_PATTERNS = ['*scc*', '*vcvs*'] 
SEEG_CHANNEL_PATTERNS = ['*of*', '*vm*',' *acc*', '*amy*', '*pf*']

# params
PLOTTING_SCALE = 200e-6
TARGET_SFREQ = 2000  

# flags
OVERWRITE = False
 
# main loop
if __name__ == "__main__":

    file_list = create_nsx_file_list(data_path=ORIGINAL_DATA_ROOT)

    for nsx_path in file_list: 

        block_folder_name = os.path.basename(os.path.dirname(nsx_path))
        simplified_block_name = block_folder_name.split('task-')[-1]

        # load nsx, dbs channels
        X_counts, ch_names_all, fs, ext_headers_all, original_elec_ids = load_blackrock_data(nsx_path)
        
        # remove all sessions less than 60 seconds total
        if (X_counts.shape[-1])/fs < 60:
            continue    

        # use just the dbs leads for now
        dbs_chans, dbs_indices = find_channels(channel_names=ch_names_all, patterns=DBS_CHANNEL_PATTERNS)
        
        ext_headers_dbs = [ext_headers_all[i] for i in dbs_indices]
        X_counts_dbs = X_counts[dbs_indices,:]

        raw_voltage_dbs = scale_to_volts(X_counts_dbs, ext_headers_dbs)
        info_dbs = mne.create_info(ch_names=dbs_chans, sfreq=fs, ch_types='dbs')
        raw_dbs = mne.io.RawArray(raw_voltage_dbs, info_dbs, verbose=False, preload=True)
        
        # downsample if necessary (should this go after filtering?)
        if fs > TARGET_SFREQ:
            raw_dbs.resample(TARGET_SFREQ)
        final_fs = raw_dbs.info['sfreq']
        
        # bandpass filter 
        raw_dbs.filter(l_freq=0.3, h_freq=500.0, method="iir", iir_params=dict(order=4, ftype="butter"))
        
        # detect line noise
        freqs_list = detect_line_noise_peaks(data=raw_dbs._data, fs=final_fs)
        
        if len(freqs_list) > 0:
            freqs_list_int = [int(freq) for freq in freqs_list]
            # notch filter
            raw_dbs.notch_filter(freqs=freqs_list_int, notch_widths=4)
        
        # make this step dependent on whether the file exists or overwrite is selected
        json_file = os.path.join(DBS_DATA_ROOT+"/annotations", f"{block_folder_name}_annotations.json")

        if not os.path.exists(json_file) or OVERWRITE:
            # create annotation category to store the artifact segments
            annotation = mne.Annotations(
                onset=[0.0],
                duration=[0.0],
                description=["BAD_artifact"]
            )
            raw_dbs.set_annotations(raw_dbs.annotations + annotation)

            # visual graph to eliminate bad chans or bad time segments
            raw_dbs.plot(block=True, scalings=PLOTTING_SCALE, title=f"{block_folder_name}" )
        
            # save off removed chans/time segments to json
            artifacts = []
            for onset, duration, description in zip(raw_dbs.annotations.onset, raw_dbs.annotations.duration, raw_dbs.annotations.description):
                if description.startswith('BAD'):
                    artifact_dict = {
                        "onset": float(onset),
                        "duration": float(duration)
                    }
                    artifacts.append(artifact_dict)

            channel_metadata = {
                "bad_chans": raw_dbs.info['bads'],
                "artifacts": artifacts
            }
            os.makedirs(DBS_DATA_ROOT+"/annotations", exist_ok=True)
        
            with open(json_file, "w") as file:
                json.dump(channel_metadata, file, indent=4)

        # don't repeat the visual plotting and manual annotations if an annotation file already exists
        else:
            with open(json_file, "r") as file:
               metadata = json.load(file)

            onsets, durations, descriptions = [],[], []

            # grab the bad channels
            raw_dbs.info['bads'] = metadata["bad_chans"]

            # remake the annotations
            for item in metadata["artifacts"]:
                if not (item["onset"] == 0.0 and item["duration"] == 0.0):  # don't include filler 
                    onsets.append(float(item["onset"]))
                    durations.append(float(item["duration"]))
                    descriptions.append("BAD_artifact")
            annotations = mne.Annotations(
                onset = onsets,
                duration= durations,
                description = descriptions
            )

            raw_dbs.set_annotations(annotations)

        # remove the bad chans if there are any
        if raw_dbs.info['bads']:
            raw_dbs.drop_channels(raw_dbs.info['bads'])
                
        reref_dir = os.path.join(REREF_DATA_ROOT, simplified_block_name)
        
        # save off cleaned data
        # if not os.path.exists(reref_dir):
        #     os.makedirs(reref_dir, exist_ok=True)
        #     fiData_path_dbs = os.path.join(reref_dir, f"fiEEG_dbs_{simplified_block_name}.fif")
        #     raw_dbs.save(fiData_path_dbs, overwrite=True)

        # os.makedirs(reref_dir, exist_ok=True)
        # fiData_path_dbs = os.path.join(reref_dir, f"fiEEG_dbs_{simplified_block_name}.fif")
        # raw_dbs.save(fiData_path_dbs, overwrite=True)

        # common average reference the data
        probes = create_dbs_probes(raw_file=raw_dbs)

        # # bipolar reference the data
        # bipolar_dir = os.path.join(DBS_DATA_ROOT, "bipolar_channels")
        # os.makedirs(bipolar_dir, exist_ok=True)
        # save_bipolar_chans(probes=probes, raw_data=raw_dbs, block_name= simplified_block_name, save_dir=bipolar_dir, mode= "regular")

        # bipolar alternating referencing for the data
        bipolar_alternating_dir = os.path.join(DBS_DATA_ROOT, "bipolar_alternating_channels")
        os.makedirs(bipolar_alternating_dir, exist_ok = True)
        save_bipolar_chans(probes=probes, raw_data=raw_dbs, block_name= simplified_block_name, save_dir=bipolar_alternating_dir, mode="alternating")

        # # common average reference
        # car_dir = os.path.join(DBS_DATA_ROOT, "car_channels")
        # os.makedirs(car_dir, exist_ok=True)
        # save_car_chans(probes=probes, raw_data=raw_dbs, block_name= simplified_block_name, save_dir=car_dir)

        # # esr -> averaging 
        # esr_dir = os.path.join(DBS_DATA_ROOT, "esr_channels")
        # os.makedirs(esr_dir, exist_ok = True)
        # save_esr_chans(probes=probes, raw_data=raw_dbs, block_name= simplified_block_name, save_dir=esr_dir)





                





