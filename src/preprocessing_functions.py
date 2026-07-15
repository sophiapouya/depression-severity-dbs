import os
import numpy as np
import mne
from .brpy.brpylib import NsxFile
from scipy.signal import welch, find_peaks
import fnmatch
from typing import Any, Literal
import re
import pandas as pd
import json

# create list of nsx files (ns3 if it exists, ns5 if there is no ns3 file)
def create_nsx_file_list(data_path: str) -> list[str]:    #data path assumed format: original_data_root/session_folder/files
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

# helper function for loading blackrock data
def _get(d: dict[str, Any], keys: list[str]) -> Any | None:
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return None

# load time series blackrock data
def load_blackrock_data(nsx_path: str) -> tuple[np.ndarray, list[str], float, list[dict[str, Any]], list[int]]:
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

# scale the time series correctly to volts
def scale_to_volts(X_counts: np.ndarray, ext_headers: list[dict[str, Any]]) -> np.ndarray:
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

# find the dbs channel names and indices       
def find_dbs_channels(
    channel_names: list[str], 
    patterns: list[str]
) -> tuple[list[str], list[int]]:
    
    dbs_chans = []
    dbs_indices = []
    for index, ch in enumerate(channel_names):
        # clean up the names
        ch_cleaned=ch.upper()
        ch_cleaned=ch_cleaned.lower()
        found = False
        for pat in patterns:
            if fnmatch.fnmatch(ch_cleaned, pat):
                dbs_chans.append(ch)
                dbs_indices.append(index)
                found = True
            if found:
                break

    return dbs_chans, dbs_indices

# find the seeg channel names and indices       
def find_seeg_channels(
    channel_names: list[str],
    pattern: re.Pattern,
) -> tuple[list[str], list[int]]:
    
    seeg_chans = []
    seeg_indices = []

    for index, ch in enumerate(channel_names):
        if pattern.search(ch):
            seeg_chans.append(ch)
            seeg_indices.append(index)

    return seeg_chans, seeg_indices

# reference and save the dbs channels as a fif file
def save_dbs_chans(
    patient: str,
    probes: dict[str, list[str]], 
    raw_data: mne.io.BaseRaw, 
    block_name: str, 
    save_dir: str, 
    mode: Literal["bipolar_alternating", "bipolar_regular", "esr", "car"]
) -> None:  
    ch_names, chans = [], []
    expected_pairs = [(1,2), (2,3), (3,4), (4,5), (5,6), (6,7), (7,8)]  #relevant only to dbstrd011 and dbstrd014 vcvs probe
    avg_ref = raw_data.get_data(picks=raw_data.ch_names).mean(axis=0)
    for probe in probes:

        if mode in ["bipolar_alternating","bipolar_regular"]:
            probe_exception = False
            if patient in ["DBSTRD011", "DBSTRD014"] and "VCVS" in probe:
                print(f"{probe} for {patient}, doing bipolar referencing")
                probe_exception = True
                contact_map = {}
                for chan in probes[probe]:
                    contact_num = int((chan.split("-")[0])[-2:])
                    contact_map[contact_num] = chan
                
                for c1, c2 in expected_pairs:
                    chan_name = f"{probe}_{c2}-{c1}"

                    if (c1 in contact_map) and (c2 in contact_map):
                        data_1 = raw_data.get_data(picks=contact_map[c1]).flatten() 
                        data_2 = raw_data.get_data(picks=contact_map[c2]).flatten()
                        data = data_2 - data_1
                        ch_names.append(chan_name)
                        chans.append(data)

                    else:
                        print(f"Skipping {chan_name} in {block_name}")

            if not probe_exception: 
                avg_1, avg_2 = [], []
                ch1, ch8 = None, None

                for chan in probes[probe]:
                    prefix = chan.split("-")[0]
                    contact = int(prefix[-2:])
                    # average 2, 3, and 4
                    if contact in (2,3,4):
                        avg_1.append(chan)
                    elif contact in (5,6,7):
                        avg_2.append(chan)
                    elif contact ==1:
                        ch1 = chan
                    elif contact == 8: 
                        ch8 = chan
                
                # if any of the channels are not present, skip this session b/c there isn't enough data
                if (ch1 is None) and (ch8 is None) and (len(avg_1) == 0) and (len(avg_2) == 0):
                    print(f"Skipping probe {probe} in {block_name} for bipolar channel since at least all channels are missing or excluded")
                    continue
                
                # compute averages only if the contributing contacts exist
                avg_1_total = None
                avg_2_total = None

                if len(avg_1) > 0:
                    avg_1_total = raw_data.get_data(picks=avg_1).mean(axis=0)

                if len(avg_2) > 0:
                    avg_2_total = raw_data.get_data(picks=avg_2).mean(axis=0)

                ch1_data = None
                ch8_data = None

                if ch1 is not None:
                    ch1_data = raw_data.get_data(picks=ch1).flatten()

                if ch8 is not None:
                    ch8_data = raw_data.get_data(picks=ch8).flatten()

                if mode == "bipolar_regular":
                    # probe_1 = ch1 - avg(2,3,4)
                    if (ch1_data is not None) and (avg_1_total is not None):
                        ch_names.append(f"{probe}_1")
                        chans.append(ch1_data - avg_1_total)
                    else:
                        print(f"Skipping {probe}_1 in {block_name}")

                    # probe_2 = avg(2,3,4) - avg(5,6,7)
                    if (avg_1_total is not None) and (avg_2_total is not None):
                        ch_names.append(f"{probe}_2")
                        chans.append(avg_1_total - avg_2_total)
                    else:
                        print(f"Skipping {probe}_2 in {block_name}")

                    # probe_3 = avg(5,6,7) - ch8
                    if (avg_2_total is not None) and (ch8_data is not None):
                        ch_names.append(f"{probe}_3")
                        chans.append(avg_2_total - ch8_data)
                    else:
                        print(f"Skipping {probe}_3 in {block_name}")

                elif mode == "bipolar_alternating":
                    # probe_1 = ch1 - ch8
                    if (ch1_data is not None) and (ch8_data is not None):
                        ch_names.append(f"{probe}_1")
                        chans.append(ch1_data - ch8_data)
                    else:
                        print(f"Skipping {probe}_1 in {block_name}")

                    # probe_2 = ch1 - avg(5,6,7)
                    if (ch1_data is not None) and (avg_2_total is not None):
                        ch_names.append(f"{probe}_2")
                        chans.append(ch1_data - avg_2_total)
                    else:
                        print(f"Skipping {probe}_2 in {block_name}")

                    # probe_3 = avg(2,3,4) - ch8
                    if (avg_1_total is not None) and (ch8_data is not None):
                        ch_names.append(f"{probe}_3")
                        chans.append(avg_1_total - ch8_data)
                    else:
                        print(f"Skipping {probe}_3 in {block_name}")

        elif mode == "esr":
            probe_chans = probes[probe]
            probe_avg = raw_data.get_data(picks=probe_chans).mean(axis=0)
            
            for chan in probes[probe]:
                prefix = chan.split("-")[0]
                num = int(prefix[-2:])
                esr_chan = raw_data.get_data(picks=chan).flatten() - probe_avg
                ch_names.append(f"{probe}_{num}")
                chans.append(esr_chan)

        elif mode == "car":
            for chan in probes[probe]:
                car_ch = raw_data.get_data(picks=chan).flatten() - avg_ref
                prefix = chan.split("-")[0]
                num = int(prefix[-2:])
                ch_names.append(f"{probe}_{num}")
                chans.append(car_ch)

    chans_arr = np.stack(chans,axis=0)

    # save the session level fif file
    info = mne.create_info(ch_names=ch_names, sfreq=raw_data.info['sfreq'], ch_types='dbs')
    raw_bp=mne.io.RawArray(chans_arr, info, verbose=False)

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

# reference and save the seeg channels as a fif file
def save_seeg_chans(
    probes: dict[str, list[str]], 
    raw_data: mne.io.BaseRaw, 
    block_name: str, 
    save_dir: str,
    seeg_metadata: dict[str,str],
) -> dict[str,list[str]]:  
   
    ch_names, chans_arr, output_dict = bipolar_seeg(raw_data=raw_data, probes=probes, seeg_metadata=seeg_metadata)

    # save the session level fif file
    info = mne.create_info(ch_names=ch_names, sfreq=raw_data.info['sfreq'], ch_types='seeg')
    raw_seeg=mne.io.RawArray(chans_arr, info, verbose=False)

    # downsample to 1000hz 
    raw_seeg.resample(1000)
    
    # save as fif file
    session_file_name = os.path.join(save_dir, f"{block_name}.fif")
    raw_seeg.save(session_file_name, overwrite=True)

    return output_dict

def get_matched_name(
    chan_name: str
) -> str:
    chan_name_parts = chan_name.split("-")
    name = str(chan_name_parts[0]+"-"+chan_name_parts[1])
    final_name = name.lower()
    return final_name

# bipolar referencing for seeg contacts
def bipolar_seeg(
    *, 
    probes: dict[str, list[str]], 
    raw_data: mne.io.BaseRaw,
    seeg_metadata: dict[str, str],
) -> tuple[list, np.ndarray]:
    
    chan_names, chans_data, chans_regions = [], [], []

    # pre pull the data
    data = raw_data.get_data()
    ch_to_idx = {ch:index for index, ch in enumerate(raw_data.ch_names)}

    for probe_name, probe_channels in probes.items():
        sorted_chans = sorted(probe_channels, key = lambda ch: int(ch.split("-")[-1]))
        for i in range(len(sorted_chans) -1):  
            # get the matching channel name to compare against the seeg metadata
            chan_1_name = get_matched_name(sorted_chans[i])
            chan_2_name = get_matched_name(sorted_chans[i+1])
            chan_1_region = seeg_metadata.get(chan_1_name, "none")
            chan_2_region = seeg_metadata.get(chan_2_name, "none")

            # look up the region for the contact
            # only include if both regions are the same OR one region is none and the other contact has a region
            # continue on if both regions are None
            if (chan_1_region == "none" and chan_2_region == "none"):
                continue
            
            # same region 
            elif ((chan_1_region == chan_2_region) or ((chan_1_region == "none") or (chan_2_region == "none"))):
            
                chan_1_data = data[ch_to_idx[sorted_chans[i]]]
                chan_2_data = data[ch_to_idx[sorted_chans[i+1]]]
                chan_name = f"{probe_name.upper()}_{i+1}"
                
                chan_data = chan_2_data - chan_1_data
                chan_names.append(chan_name)
                chans_data.append(chan_data)
                
                if chan_1_region == "none":
                    region = chan_2_region
                else:
                    region = chan_1_region

                chans_regions.append(region)
            
            # regions are different, don't include
            else:
                continue

    output_dict = {
        "channel_name": chan_names,
        "channel_region": chans_regions
    }

    chans_data_arr = np.stack(chans_data, axis=0)
    return chan_names, chans_data_arr, output_dict


# detect line noise for further filtering for each session
def detect_line_noise_peaks(data: np.ndarray, fs: float,
                             fmin=50.0,
                             fmax=300.0,
                             prominence_db=10.0,
                             min_distance_hz=20.0) -> np.ndarray:
    # Compute PSD
    f, Pxx = welch(
        data,
        fs=fs,
        window='hamming',
        nperseg=int(2 * fs),
        noverlap=int(fs),
        nfft=int(2 * fs),
    )

    # collapse all channels to one "typical" channel for the session
    Pxx_median = np.median(Pxx, axis=0)
    Pxx_db = 10 * np.log10(Pxx_median + 1e-20)

    # Limit to desired frequency range
    fmax = min(fmax, fs / 2.0)
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

    # round a bit so they're nice numbers
    peak_freqs = np.round(peak_freqs, 1)

    return peak_freqs

# create a dictionary of dbs probes with clean names
def create_dbs_probes(raw_file: mne.io.BaseRaw) -> dict[str, list[str]]:
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

# create a dictionary of seeg probes with clean names
def create_seeg_probes(
    raw_file: mne.io.BaseRaw,
) -> dict[str, list[str]]:
    probes = {}
    for channel in raw_file.ch_names:
        match = re.match(r"^[^\d]+", channel)
        probe_prefix = match.group(0)
        if probe_prefix not in probes:
            probes[probe_prefix] = []
        
        probes[probe_prefix].append(channel)
    return probes

def get_seeg_metadata(
    file_path: str,
    patient: str,
) -> dict[str, str]:
    
    final_dict = {}
    excel_df = pd.read_excel(file_path, sheet_name=patient)
    
    # Drop all the DBS contacts
    clean_excel = excel_df[excel_df['Type'] != "DBS"].copy()  # .copy() fixes the SettingWithCopyWarning too

    # Fill missing areas and drop rows with no label
    clean_excel["area"] = clean_excel["area"].fillna("none")
    clean_excel = clean_excel.dropna(subset=["Label"])

    for i, row in clean_excel.iterrows():
        label = str(row["Label"]).strip().lower()
        region = str(row["area"]).strip().lower()
        # exclude white matter contacts regardless of region label
        if str(row.get("Grey v White", "Grey")).strip() == "White":
            region = "none"
        final_dict[label] = region

    return final_dict

def output_metadata(
    file_path: str,
    output_dict: dict[str, list[str]]
) -> None: 
    
    df = pd.DataFrame(output_dict)
    df.to_csv(file_path, index=False)

    
    
        
