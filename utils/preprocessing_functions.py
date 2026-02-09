import os
import numpy as np
import mne
import scipy.io as sio
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from .brpylib import NsxFile
from scipy.signal import filtfilt, welch, firwin, convolve, find_peaks, decimate, hilbert
import fnmatch
import pandas as pd


def _get(d, keys, default=None):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return default

# create list of nsx files (ns3 if it exists, ns5 if there is no ns3 file)
def create_nsx_file_list(data_path):    #data path assumed format: original_data_root/session_folder/files
    file_list = []  # create a final list of files to process

    for directory, sessions, files in os.walk(data_path): 
        ns3found=False    
        # first pass through files to see if ns3 exists
        for file in files: 
            if file.endswith(".ns3"):
                file_list.append(os.path.join(directory, file))
                ns3found = True
                break 

        if ns3found == False:
            # second pass to see append ns5 in the case that ns3 doesn't exist
            for file in files:
                if file.endswith(".ns5"):
                    file_list.append(os.path.join(directory,file))
                    break

    return file_list

def load_blackrock_data(nsx_path):
    nsx = NsxFile(nsx_path)
    basic_header = getattr(nsx, "basic_header", {}) or {}
    period = basic_header.get("Period", basic_header.get("period"))
    fs = 30000.0 / float(period) if period else 2000.0
    data_dict = nsx.getdata()
    ext_headers = getattr(nsx, "extended_headers", [])
    original_elec_ids = [eh.get("ElectrodeID", i+1) for i, eh in enumerate(ext_headers)]
    ch_names_all = [eh.get('ElectrodeLabel', f'ch{i+1}') for i, eh in enumerate(ext_headers)]
    arr = data_dict.get("data")
    X_counts = np.concatenate([np.asarray(a, dtype=np.float32) for a in arr], axis=1) if isinstance(arr, list) else np.asarray(arr, dtype=np.float32)
    nsx.close()
    return X_counts, ch_names_all, fs, ext_headers, original_elec_ids

def scale_to_volts(X_counts, ext_headers):
    n_ch = X_counts.shape[0]
    X = X_counts.astype(np.float32, copy=True)
    if not ext_headers or len(ext_headers) < n_ch:
        return X * 1e-6  # assume µV
    scales = np.empty(n_ch, dtype=np.float64)
    offsets = np.empty(n_ch, dtype=np.float64)
    post_mul = np.ones(n_ch, dtype=np.float64)
    ok = True
    for i, eh in enumerate(ext_headers[:n_ch]):
        min_d = _get(eh, ["MinDigiValue","MinDigitalValue"])
        max_d = _get(eh, ["MaxDigiValue","MaxDigitalValue"])
        min_a = _get(eh, ["MinAnalogValue"])
        max_a = _get(eh, ["MaxAnalogValue"])
        units = _get(eh, ["Units","units"])
        if None in (min_d, max_d, min_a, max_a) or float(max_d) == float(min_d):
            ok = False; break
        s = (float(max_a) - float(min_a)) / (float(max_d) - float(min_d))
        o = float(min_a) - float(min_d) * s
        if units and str(units).lower() in {"uv","µv","microvolts","microvolt"}:
            post_mul[i] = 1e-6
        scales[i] = s; offsets[i] = o
    if not ok:
        return X * 1e-6
    for i in range(n_ch):
        X[i] = X[i] * scales[i] + offsets[i]
        X[i] = X[i] * post_mul[i]
    return X
       
def find_channels(channel_names, patterns):
    dbs_chans = []
    dbs_indices = []
    for index, ch in enumerate(channel_names):
        # clean up the names
        ch_cleaned=ch.upper()
        ch_cleaned=ch.lower()
        for pat in patterns:
            if fnmatch.fnmatch(ch_cleaned, pat):
                dbs_chans.append(ch)
                dbs_indices.append(index)
    
    return dbs_chans, dbs_indices

def save_bipolar_chans(probes, raw_data, block_name, save_dir, mode):
    
    ch_names = []
    bp_chans = []
    for probe in probes:
        avg_1, avg_2 = [], []
        ch1, ch8 = None, None
        for chan in probes[probe]:
            prefix = chan.split("-")[0]
            # average 2, 3, and 4
            if prefix[-1] in ("2","3","4"):
                avg_1.append(chan)
            elif prefix[-1] in ("5","6","7"):
                avg_2.append(chan)
            elif prefix[-1] == "1":
                ch1 = chan
            elif prefix[-1] == "8": 
                ch8 = chan

        avg_1_total = raw_data.get_data(picks =avg_1).mean(axis=0)   # 1 channel, time series data
        avg_2_total = raw_data.get_data(picks=avg_2).mean(axis=0)    # 1 channel, time series data

        if mode == "regular":
            bp1 = raw_data.get_data(picks=ch1).flatten() - avg_1_total
            bp2 = avg_1_total - avg_2_total
            bp3 = avg_2_total - raw_data.get_data(picks=ch8).flatten()
        elif mode == "alternating":
            bp1 = raw_data.get_data(picks=ch1).flatten() - raw_data.get_data(picks=ch8).flatten()
            bp2 = raw_data.get_data(picks=ch1).flatten() - avg_2_total
            bp3 = avg_1_total - raw_data.get_data(picks=ch8).flatten()
        ch_names.append(f"{probe}_1")
        ch_names.append(f"{probe}_2")
        ch_names.append(f"{probe}_3")
        bp_chans.append(bp1)
        bp_chans.append(bp2)
        bp_chans.append(bp3)
    bp_chans_arr = np.stack(bp_chans,axis=0)

    # save the session level fif file
    info = mne.create_info(ch_names=ch_names, sfreq=raw_data.info['sfreq'], ch_types='dbs')
    raw_bp=mne.io.RawArray(bp_chans_arr, info, verbose=False)

    # deal with the annotations
    annotations = raw_data.annotations
    # remove the onset 0.0, duration 0.0 used when visually inspecting
    if annotations is not None and len(annotations)>0:
        new_onset, new_duration, new_description = [],[],[]
        for onset, duration, description in zip(annotations.onset, annotations.duration, annotations.description):
            if not (float(onset) == 0.0 and float(duration) == 0.0):
                new_onset.append(float(onset))
                new_duration.append(float(duration))
                new_description.append(description)
        if len(new_onset)>0:
            raw_bp.set_annotations(mne.Annotations(onset=new_onset, duration=new_duration,description=new_description))
        else:
            raw_bp.set_annotations(mne.Annotations(onset=[],duration=[], description=[]))
    
    # downsample to 1000hz 
    raw_bp.resample(1000)
    
    # save as fif file
    session_file_name = os.path.join(save_dir, f"{block_name}.fif")
    raw_bp.save(session_file_name, overwrite=True)


def save_car_chans(probes, raw_data, block_name, save_dir):
    avg_ref = raw_data.get_data(picks=raw_data.ch_names).mean(axis=0)

    for probe in probes:
        for chan in probes[probe]:
            car_ch = raw_data.get_data(picks=chan).flatten() - avg_ref
            num = chan.split("-")[0]
            # save the car channel
            car_file_name = os.path.join(save_dir, f"{block_name}_{probe}_carCh_{num[-1]}.npy")
            np.save(car_file_name, car_ch)

def save_esr_chans(probes, raw_data, block_name, save_dir):

    for probe in probes:
        probe_chans = probes[probe]
        probe_avg = raw_data.get_data(picks=probe_chans).mean(axis=0)
        for chan in probes[probe]:
            num = chan.split("-")[0]
            esr_chan = raw_data.get_data(picks=chan).flatten() - probe_avg
            # save the esr chan
            esr_file_name = os.path.join(save_dir, f"{block_name}_{probe}_esrCh_{num[-1]}.npy")
            np.save(esr_file_name, esr_chan)

def detect_line_noise_peaks(data, fs,
                             fmin=50.0,
                             fmax=500.0,
                             prominence_db=10.0,
                             min_distance_hz=20.0):
    """
    Detect narrow, high-power peaks (line noise / harmonics) in the PSD
    of band-passed data.

    data: ndarray, shape (n_channels, n_samples)
    fs:   sampling frequency
    returns: 1D numpy array of peak frequencies in Hz
    """

    # Use first channel as representative (line noise is global)
    x = data[0, :]

    # Compute PSD
    f, Pxx = welch(
        x,
        fs=fs,
        window='hamming',
        nperseg=int(2 * fs),
        noverlap=int(fs),
        nfft=int(2 * fs),
    )
    Pxx_db = 10 * np.log10(Pxx + 1e-20)

    # Limit to desired frequency range
    fmax = min(fmax, fs / 2.0 - 1.0)
    band_mask = (f >= fmin) & (f <= fmax)
    f_band = f[band_mask]
    Pxx_band = Pxx_db[band_mask]

    if len(f_band) < 3:
        return np.array([])

    # Convert min_distance_hz to bins
    df = f_band[1] - f_band[0]
    min_distance_bins = max(1, int(min_distance_hz / df))

    # Find peaks with specified prominence
    peak_inds, _ = find_peaks(Pxx_band,
                              prominence=prominence_db,
                              distance=min_distance_bins)

    peak_freqs = f_band[peak_inds]

    # Optional: round a bit so they're nice numbers
    peak_freqs = np.round(peak_freqs, 1)

    return peak_freqs

def create_dbs_probes(raw_file):
    probes = {}
    for ch in raw_file.ch_names:
        prefix = ch.split("-")[0]
        if prefix[0] not in ("L","R"):
            prefix = prefix[1:]
        probe_name = prefix[:-2]
        
        if probe_name not in probes:
            probes[probe_name] = []
        probes[probe_name].append(ch)
    return probes

def save_png(raw, filename, scaling):
    with mne.viz.use_browser_backend('matplotlib'):
        fig = raw.plot(
            show=False,
            n_channels=len(raw.ch_names),
            duration=1000,
            scalings=scaling,
        )
    fig.savefig(filename)  # e.g., "…/DBS_block01.png"

def save_raw_voltages(master_block_dir, block_name, dbs_voltage, seeg_voltage):
    os.makedirs(master_block_dir, exist_ok=True)
    raw_vol_path_dbs = os.path.join(master_block_dir, f"RawVoltage_dbs_{block_name}.mat")
    raw_vol_path_seeg = os.path.join(master_block_dir, f"RawVoltage_seeg_{block_name}.mat")
    sio.savemat(raw_vol_path_dbs, {'RawVoltage': dbs_voltage}, do_compression=True)
    sio.savemat(raw_vol_path_seeg, {'RawVoltage': seeg_voltage}, do_compression=True)

def save_master_vars(master_block_dir, block_name, master_vars_dbs, master_vars_seeg):
    os.makedirs(master_block_dir, exist_ok=True)
    master_vars_path_dbs = os.path.join(master_block_dir, f"master_dbs_{block_name}.mat")
    master_vars_path_seeg = os.path.join(master_block_dir, f"master_seeg_{block_name}.mat")
    sio.savemat(master_vars_path_dbs, {'master_vars': master_vars_dbs}, oned_as='column')
    sio.savemat(master_vars_path_seeg, {'master_vars': master_vars_seeg}, oned_as='column')

def save_psd_session(path, raw_file, block_name, fs):
    
    psds = []
    for ch in raw_file.ch_names:
        voltage = raw_file.get_data(picks=ch)[0]
        f, Pxx = welch(voltage, fs=fs, 
                    window='hamming', 
                    nperseg=int(2*fs),   # 2 second window
                    noverlap=int(fs),   # 50% overlap
                    nfft = int(2*fs))
        Pxx_db = 10 *np.log10(Pxx)
        psds.append(Pxx_db)
    
    psd_arr = np.array(psds)
    avg_psd = np.average(psd_arr, axis=0)
    fig, ax = plt.subplots(1, 1, figsize=(8,6))
    ax.plot(f, avg_psd)
    ax.set_title(f"Average PSD for {block_name}")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_xlim(0, 500)
    ax.set_ylabel("Power (dB)")
    fig.savefig(path)
    plt.close(fig)
    

def save_psd_plots(pdf_path, raw_file, block_name, fs, rows=8, cols=4):
    with PdfPages(pdf_path) as pdf:

        fig, axes = plt.subplots(rows, cols, figsize=(20,20))
        axes = axes.ravel()

        for dbs_ch, ax in zip(raw_file.ch_names, axes):
            
            voltage = raw_file.get_data(picks=dbs_ch)[0]
            f, Pxx = welch(voltage, fs=fs, 
                        window='hamming', 
                        nperseg=int(2*fs),   # 2 second window
                        noverlap=int(fs),   # 50% overlap
                        nfft = int(2*fs))
            Pxx_db = 10*np.log10(Pxx)
            ax.plot(f, Pxx_db)
            ax.set_title(dbs_ch, fontsize=8)
            ax.set_xlabel("Hz", fontsize=6)
            ax.set_xlim(0,500)
            ax.set_ylabel("Power (dB)", fontsize=6)

        fig.suptitle(f"PSD_{block_name}", fontsize=12)
        fig.tight_layout(pad=1.0)
        pdf.savefig(fig)

    
def save_power_data(band_coefficients, FEATURE_BANDS, EXCLUDED_SESSIONS, ref_type, subj_name, sbj_dir):
    if ref_type == "bipolar":
        working_dir = os.path.join(sbj_dir, "bipolar_channels")
        power_dir = os.path.join(working_dir, "power_bipolar")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "car":
        working_dir = os.path.join(sbj_dir, "car_channels")
        power_dir = os.path.join(working_dir, "power_car")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "esr":
        working_dir = os.path.join(sbj_dir, "esr_channels")
        power_dir = os.path.join(working_dir, "power_esr")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "bipolar_alternating":
        working_dir = os.path.join(sbj_dir, "bipolar_alternating_channels")
        power_dir = os.path.join(working_dir, "power_bipolar_alternating")
        os.makedirs(power_dir, exist_ok=True)
    

    master_list = []
    with os.scandir(working_dir) as files:
        for file in files:
            session_dict = {}
            # make sure it's actually a file
            if not file.name.endswith(".npy"):
                continue

            name, ext = os.path.splitext(file.name)
            name_parts = name.split("_")
            ch_num= name_parts[-1]
            probe_name = name_parts[-3]
            session_name = "_".join(name_parts[:-3])
            probe_ch_num = probe_name + "_" + ch_num
            
            # check if it's an excluded session
            if session_name in EXCLUDED_SESSIONS[subj_name]:
                continue

            session_dict["session"] = session_name
            session_dict["probe"] = probe_name
            session_dict["ch_name"] = probe_ch_num

            bipolar_data = np.load(file.path)
            
            #downsample the data
            bipolar_data=decimate(bipolar_data,2)
            
            for band, values in FEATURE_BANDS.items():
        
                b,a = band_coefficients[band]
                band_data = filtfilt(b,a,bipolar_data)
                analytic_signal = hilbert(band_data)
                power = (np.abs(analytic_signal))**2
                log_power=np.log10(power)
                avg_log_power= np.mean(log_power)
                session_dict[band] = avg_log_power

            #append session dictionary to master list
            master_list.append(session_dict)
            
    df = pd.DataFrame(master_list)
    df.to_csv(os.path.join(power_dir,f"{subj_name}_{ref_type}_power.csv"))

def save_power_data_from_fif(band_coefficients, FEATURE_BANDS, EXCLUDED_SESSIONS, ref_type, subj_name, sbj_dir):
    if ref_type == "bipolar":
        working_dir = os.path.join(sbj_dir, "bipolar_channels")
        power_dir = os.path.join(working_dir, "power_bipolar")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "car":
        working_dir = os.path.join(sbj_dir, "car_channels")
        power_dir = os.path.join(working_dir, "power_car")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "esr":
        working_dir = os.path.join(sbj_dir, "esr_channels")
        power_dir = os.path.join(working_dir, "power_esr")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "bipolar_alternating":
        working_dir = os.path.join(sbj_dir, "bipolar_alternating_channels")
        power_dir = os.path.join(working_dir, "power_bipolar_alternating")
        os.makedirs(power_dir, exist_ok=True)
    
    # window settings
    SFREQ = 1000
    WINDOW_S = 2.0
    STEP_S = 1.0
    WINDOW_SAMPLES = int(WINDOW_S * SFREQ)
    STEP_SAMPLES = int(STEP_S * SFREQ)

    master_list = []
    with os.scandir(working_dir) as files:
        for file in files:
            session_dict = {}
            # make sure it's a fif file
            if not file.name.endswith(".fif"):
                continue

            session_name = os.path.splitext(file.name)[0]

            if session_name in EXCLUDED_SESSIONS[subj_name]:
                continue

            file_path = file.path

            raw = mne.io.read_raw_fif(file_path, preload=True, verbose=False)

            data = raw.get_data()
            ch_names = raw.ch_names

            total_time_s = raw.times[-1]
            n_samples = data.shape[1]

            # bad intervals list
            bads = []
            for onset, duration, desc in zip(raw.annotations.onset,
                                             raw.annotations.duration,
                                             raw.annotations.description):
                if str(desc) == "BAD_artifact":
                    bads.append((float(onset), float(onset + duration)))
            
            bads = sorted(bads, key=lambda x: x[0])

            # good intervals list
            goods = []
            current = 0.0

            for bad_start, bad_end in bads:
                if bad_start > current:
                    goods.append((current, bad_start))
                current = bad_end

            if current < total_time_s:
                goods.append((current, total_time_s))

            # If there were no BAD annotations, goods = whole session
            if len(bads) == 0:
                goods = [(0.0, total_time_s)]

            # Keep only good intervals that can fit at least ONE full window
            usable_goods = []
            for g_start, g_end in goods:
                if (g_end - g_start) >= WINDOW_S:
                    usable_goods.append((g_start, g_end))

            if len(usable_goods) == 0:
                # no usable data after artifact removal
                continue

            # Per channel power calc
            for ch_idx, ch in enumerate(ch_names):

                probe_name = ch.split("_")[0]

                session_dict = {
                    "session": session_name,
                    "probe": probe_name,
                    "ch_name": ch
                }
                x = data[ch_idx, :]

                for band in FEATURE_BANDS.keys():

                    b, a = band_coefficients[band]
                    window_means = []

                    # loop over each good interval
                    for g_start_s, g_end_s in usable_goods:

                        start_sample = int(np.round(g_start_s * SFREQ))
                        end_sample = int(np.round(g_end_s * SFREQ))

                        # make sure bounds are valid
                        start_sample = max(0, start_sample)
                        end_sample = min(n_samples, end_sample)

                        interval_len = end_sample - start_sample
                        n_windows = (interval_len - WINDOW_SAMPLES) // STEP_SAMPLES + 1

                        if n_windows == 0:
                            continue
                        
                        x_interval = x[start_sample:end_sample]
                        band_data = filtfilt(b, a, x_interval)

                        for w in range(n_windows):
                            s = w * STEP_SAMPLES
                            e = s + WINDOW_SAMPLES

                            segment = band_data[s:e]

                            # Safety: skip if window isn't full length
                            if segment.shape[0] != WINDOW_SAMPLES:
                                continue

                            analytic = hilbert(segment)
                            power = np.abs(analytic) ** 2
                            log_power = np.log10(power)
                            window_means.append(np.mean(log_power))

                    # If a channel/band has no windows (should be rare), store NaN
                    if len(window_means) == 0:
                        session_dict[band] = np.nan
                    else:
                        session_dict[band] = float(np.mean(window_means))

                master_list.append(session_dict)

    df = pd.DataFrame(master_list)
    df.to_csv(os.path.join(power_dir, f"{subj_name}_{ref_type}_power.csv"), index=False)


def save_fif_chans(ref_type, sbj_dir, subj_name, EXCLUDED_SESSIONS):
    if ref_type == "bipolar":
        working_dir = os.path.join(sbj_dir, "bipolar_channels")
        power_dir = os.path.join(working_dir, "power_bipolar")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "car":
        working_dir = os.path.join(sbj_dir, "car_channels")
        power_dir = os.path.join(working_dir, "power_car")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "esr":
        working_dir = os.path.join(sbj_dir, "esr_channels")
        power_dir = os.path.join(working_dir, "power_esr")
        os.makedirs(power_dir, exist_ok=True)
    elif ref_type == "bipolar_alternating":
        working_dir = os.path.join(sbj_dir, "bipolar_alternating_channels")
        power_dir = os.path.join(working_dir, "power_bipolar_alternating")
        os.makedirs(power_dir, exist_ok=True)
    
    # create output directory for py_neuromod files
    py_neuro_dir = os.path.join(sbj_dir,"py_neuro_files")
    os.makedirs(py_neuro_dir, exist_ok=True)

    with os.scandir(working_dir) as files:
        session_info = {}
        for file in files:
            # make sure it's actually a file
            if not file.name.endswith(".npy"):
                continue

            name, ext = os.path.splitext(file.name)
            name_parts = name.split("_")
            ch_num= name_parts[-1]
            probe_name = name_parts[-3]
            session_name = "_".join(name_parts[:-3])
            probe_ch_num = probe_name + "_" + ch_num
            
            # check if it's an excluded session
            if session_name in EXCLUDED_SESSIONS[subj_name]:
                continue

            # load the file
            bipolar_data = np.load(file.path)
            
            #downsample the data
            bipolar_data=decimate(bipolar_data,2)

            # if the category doesn't exist, create it:
            if session_name not in session_info:
                session_info[session_name] = {'data': [], 'ch_names': []}
            
            # load the time series for the data category
            session_info[session_name]['data'].append(bipolar_data)
            session_info[session_name]['ch_names'].append(probe_ch_num)

    
    for session_name, category_data in session_info.items():
        # create the time series for the session (channels x time points)
        session_data = np.array(category_data['data'])

        # create the session info (the channel names)
        info = mne.create_info(ch_names=category_data['ch_names'], sfreq=1000)

        # create the pynm mne object
        pynm_raw_obj = mne.io.RawArray(session_data, info)

        # save as a fif file
        save_path = os.path.join(py_neuro_dir, f'{session_name}_pynm.fif')
        pynm_raw_obj.save(save_path, overwrite=True)

