import os
import numpy as np
import mne
import scipy.io as sio
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from brpylib import NsxFile
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

def save_bipolar_chans(probes, raw_data, block_name, save_dir):
    for probe in probes:
        avg_1, avg_2 = [], []
        ch1, ch8 = [], []
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

        # save off each bipolar channel
        bp_chans = []
        bp1 = raw_data.get_data(picks=ch1).flatten() - avg_1_total
        bp_chans.append(bp1)
        bp2 = avg_1_total - avg_2_total
        bp_chans.append(bp2)
        bp3 = avg_2_total - raw_data.get_data(picks=ch8).flatten()
        bp_chans.append(bp3)

        for i, bp in enumerate(bp_chans):
            bipolar_file_name = os.path.join(save_dir, f"{block_name}_{probe}_bipolarCh_{i+1}.npy")
            np.save(bipolar_file_name, bp)

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

def matlab_notch_filter(data, fs, notch_freqs, N=300):
    # data: shape (n_channels, n_samples)
    # fs: sampling rate
    # notch_freqs: list like [60] or [60, 120, 180]

    if notch_freqs is None or len(notch_freqs)==0:
        return data

    freqconv = 2 / fs
    filters = []

    for f0 in notch_freqs:
        wn = f0 * freqconv
        h = firwin(N+1, wn, window='hann')     # FIR lowpass
        nf = 2*h - np.concatenate([
            np.zeros(N//2),
            np.array([1.0]),
            np.zeros(N//2)
        ])
        filters.append(nf)

    # combine notches if needed
    if len(filters) > 1:
        nf = filters[0]
        for f in filters[1:]:
            nf = convolve(nf, f)
    else:
        nf = filters[0]

    # apply filtfilt to each channel
    out = np.zeros_like(data)
    for i in range(data.shape[0]):
        out[i] = filtfilt(nf, [1.0], data[i])

    return out

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