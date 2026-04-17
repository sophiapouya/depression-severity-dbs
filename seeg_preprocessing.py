from src.preprocessing_functions import create_nsx_file_list, load_blackrock_data, scale_to_volts, find_seeg_channels, save_seeg_chans, create_seeg_probes, get_seeg_metadata, output_metadata
import mne
import os
import re
from config import ROOT_DIR, CATDI_ELECTRODES

# paths
DATA_ROOT = str(ROOT_DIR)
PROJECT_NAME = "CATDI"
ALL_SUBJ = ["DBSTRD001","DBSTRD002","DBSTRD006","DBSTRD008","DBSTRD010","DBSTRD011","DBSTRD014"]
# ALL_SUBJ = ["DBSTRD011"]
ELECTRODE_INFO_EXCEL= str(CATDI_ELECTRODES)

# params
PLOTTING_SCALE = 200e-6
TARGET_SFREQ = 2000  

# seeg probe regex pattern
SEEG_PATTERN = re.compile(r"[A-Za-z]+-[A-Za-z]+\d+-\d{3}")

# flags
OVERWRITE = True

# bad channels for seeg
BAD_CHANNELS = {
    "DBSTRD001": ['RdPF-AMY05-197','RdPF-AMY06-198','RdPF-AMY07-199''LdPF-ACC07-053','LdPF-ACC08-054', 'LdPF-mOF11-041', 'LdPF-mOF12-042', 'LdPF-mPF11-025', 'LdPF-mPF12-026', 'LdPF-mPF13-027'], 
    "DBSTRD002": ['lVPF-OF13-013', 'LdPF-mOF14-030','LSTG-Amy03-047','LdPF-mPF09-073','LdPF-mPF10-074','LdPF-mPF12-076','RdPF-mPF08-200','RdPF-mPF15-207'],
    "DBSTRD006": ['CLpDLP-AC08-024','ZLpDLP-AC09-025','LpDLP-AC13-029', 'LpDLP-AC14-030', 'LMTG-Amy01-065', 'LMTG-Amy02-066', 'RaDLP-vmP03-131', 'RVLPF-OF16-192'],
    "DBSTRD008": ['ZLdPF-ACC08-024','CLdPF-ACC09-025','RSTG-Amy14-206','RSTG-Amy09-201','RSTG-Amy07-199', 'RdPF-mPF03-131', 'LdPF-ACC10-026'],
    "DBSTRD010": ['ZLdPF-ACC04-020','CLdPF-ACC05-021','LdPF-ACC06-022'],
    "DBSTRD011": ['RMTG-Amy02-066','RMTG-Amy03-067','RMTG-Amy04-068','RMTG-Amy05-069'],
    "DBSTRD014": ['CSub-Gale08-100','Sub-Gale01-093', 'Sub-Gale02-094', 'Sub-Gale03-095', 'Sub-Gale04-096', 'Sub-Gale05-097', 'Sub-Gale06-098', 'Sub-Gale07-099']
}

# preprocessing code
for SBJ_NAME in ALL_SUBJ:

    ORIGINAL_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'originalData', SBJ_NAME)
    SEEG_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'seegData', SBJ_NAME)
    # get the metadata from the electrodes excel file
    seeg_metadata = get_seeg_metadata(file_path = ELECTRODE_INFO_EXCEL, patient=SBJ_NAME)

    file_list = create_nsx_file_list(data_path=ORIGINAL_DATA_ROOT)
    skipped_sessions = []
    for nsx_path in file_list: 

        block_folder_name = os.path.basename(os.path.dirname(nsx_path))
        simplified_block_name = block_folder_name.split('task-')[-1]

        # load nsx, dbs channels
        X_counts, ch_names_all, fs, ext_headers_all, original_elec_ids = load_blackrock_data(nsx_path)
        
        # remove all sessions less than 45 seconds total
        if (X_counts.shape[-1])/fs < 45:
            skipped_sessions.append((simplified_block_name,(X_counts.shape[-1])/fs))
            continue    

        # use just the dbs leads for now
        seeg_chans, seeg_indices = find_seeg_channels(channel_names=ch_names_all, pattern=SEEG_PATTERN)
        
        ext_headers_dbs = [ext_headers_all[i] for i in seeg_indices]
        X_counts_dbs = X_counts[seeg_indices,:]

        raw_voltage_dbs = scale_to_volts(X_counts_dbs, ext_headers_dbs)
        info_dbs = mne.create_info(ch_names=seeg_chans, sfreq=fs, ch_types='seeg')
        raw_seeg = mne.io.RawArray(raw_voltage_dbs, info_dbs, verbose=False)
        
        # downsample if necessary (should this go after filtering?)
        if fs > TARGET_SFREQ:
            raw_seeg.resample(TARGET_SFREQ)
        final_fs = raw_seeg.info['sfreq']
        
        # bandpass filter 
        raw_seeg.filter(l_freq=0.3, h_freq=500.0, method="iir", iir_params=dict(order=4, ftype="butter"))
        
        # detect line noise
        freqs_list = [60, 120, 180]
        
        if len(freqs_list) > 0:
            # notch filter
            raw_seeg.notch_filter(freqs=freqs_list, notch_widths=4)

        # grab the bad channels from the definition above
        bads = [chan for chan in BAD_CHANNELS[SBJ_NAME] if chan in seeg_chans]
        raw_seeg.info['bads'] = bads

        # remove the bad chans if there are any
        if raw_seeg.info['bads']:
            raw_seeg.drop_channels(raw_seeg.info['bads'])

        # common average reference the data
        probes = create_seeg_probes(raw_file=raw_seeg)

        # bipolar reference the data
        bipolar_dir = os.path.join(SEEG_DATA_ROOT, "bipolar_channels")
        os.makedirs(bipolar_dir, exist_ok=True)
        output_data = save_seeg_chans(probes=probes, raw_data=raw_seeg, block_name= simplified_block_name, save_dir=bipolar_dir, seeg_metadata=seeg_metadata)

    for session, length in skipped_sessions:
        print(f"session skipped: {session}")
        print(f"length of session: {length}")
    # save the metadata 
    output_file = os.path.join(SEEG_DATA_ROOT,f"{SBJ_NAME}_metadata.csv")
    output_metadata(output_file, output_data)







                





