#!/usr/bin/env python

"""
Batch processing script for DBS data with segregated output directories.
- FIX: Replaced MNE's set_eeg_reference with robust NumPy Manual CAR due to internal errors.
"""

import os
import re
import fnmatch
import numpy as np
import mne
import scipy.io as sio
from brpylib import NsxFile
from typing import Dict, Any, List

# ==============================
# --- 1. USER CONFIGURATION ---
# ==============================

# --- Paths (Set these to match your environment) ---
DATA_ROOT = "/Users/sophiapouya/workspace/bcm"
PROJECT_NAME = "CATDI"
SBJ_NAME = "DBSTRD002" # Subject to process

# --- Dynamic Path Assembly ---
ORIGINAL_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'originalData', SBJ_NAME)
DBS_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'dbsData', SBJ_NAME)
REREF_DATA_ROOT = os.path.join(DBS_DATA_ROOT, 'rerefData')
POWER_DATA_ROOT = os.path.join(DBS_DATA_ROOT, 'powerData')
MASTER_DATA_ROOT = os.path.join(DBS_DATA_ROOT, 'originalData')

# --- Channel Definitions (Case-insensitive matching) ---
DBS_CHANNEL_PATTERNS = ['*scc*', '*vcvs*'] 
NONNEURAL_PATTERNS = ['empty*', 'ref*']

# --- Processing Parameters ---
TARGET_SFREQ = 1000  
NOTCH_FREQS = [60, 120, 180]
APPLY_LEGACY_0_25_SCALING = True 

# --- Feature Parameters (Matching MATLAB's Frequencies) ---
FEATURE_BANDS = {
    'delta': [1, 4], 'theta': [4, 8], 'alpha': [8, 12], 
    'beta': [12, 30], 'low_gamma': [35, 70], 'high_gamma': [70, 150]
}

# ==============================
# --- 2. HELPER FUNCTIONS ---
# ==============================

def _load_blackrock_data(nsx_path):
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

def _scale_to_volts(X_counts, ext_headers):
    X_volts = X_counts.astype(np.float32, copy=True)
    if not ext_headers or len(ext_headers) < X_counts.shape[0]:
        return X_volts * 1e-6
    for i, eh in enumerate(ext_headers):
        min_d, max_d = eh.get("MinDigiValue"), eh.get("MaxDigiValue")
        min_a, max_a = eh.get("MinAnalogValue"), eh.get("MaxAnalogValue")
        units = eh.get("Units", "µV")
        if any(v is None for v in [min_d, max_d, min_a, max_a]) or float(max_d) == float(min_d): continue
        s = (float(max_a) - float(min_a)) / (float(max_d) - float(min_d))
        o = float(min_a) - float(min_d) * s
        X_volts[i] = X_volts[i] * s + o
        if str(units).lower() in {"uv", "µv", "microvolts"}:
            X_volts[i] = X_volts[i] * 1e-6
    return X_volts

def _create_probes(ch_names: List[str], patterns: List[str]) -> Dict[str, List[str]]:
    """
    Groups channels by physical lead name (e.g., 'LSCC', 'RVCVS') for CAR.
    The final output names are upper-cased and stripped to match the MNE object.
    """
    probes = {}
    
    strip_pattern = re.compile(r'[-_.]*\d+[-_.]*\d*$|[-_.]*\d*$')
    target_stems = [pat.strip('*') for pat in patterns] 

    for ch in ch_names:
        ch_clean = ch.strip().upper() 
        ch_lower = ch_clean.lower()
        
        # Strip contact number and unique ID
        group_name = strip_pattern.sub('', ch_lower)

        # Sanity check: ensure the group name contains one of the targets
        if group_name and (any(target in group_name for target in target_stems)):
            if group_name not in probes:
                probes[group_name] = []
            
            # Append the cleaned, uppercase name (MUST match raw.ch_names)
            probes[group_name].append(ch_clean) 

    return probes

# ==============================
# --- 3. EXECUTION FUNCTION ---
# ==============================

def process_dbs_pipeline(nsx_path: str, block_folder_name: str) -> None:
    
    block_name = block_folder_name 
    simplified_block_name = block_name.split('task-')[-1] 
    print(f"\n--- Processing Block: **{block_name}** ---")
    
    # --- STEP 1: Load NSx, Select DBS Channels, Scale ---
    X_counts, ch_names_all, fs, ext_headers_all, original_elec_ids = _load_blackrock_data(nsx_path)

    dbs_indices_to_keep: List[int] = []
    ch_names_dbs: List[str] = []
    
    for i, ch in enumerate(ch_names_all):
        ch_clean = ch.strip().upper() # Cleaned and upper-cased name
        ch_lower = ch_clean.lower()
        
        is_dbs = any(fnmatch.fnmatch(ch_lower, pat) for pat in DBS_CHANNEL_PATTERNS) 
        is_nonneural = any(fnmatch.fnmatch(ch_lower, pat) for pat in NONNEURAL_PATTERNS)
        
        if is_dbs and not is_nonneural:
            dbs_indices_to_keep.append(i)
            ch_names_dbs.append(ch_clean) # Use the cleaned, UPPERCASE name

    if not ch_names_dbs:
        print("  ❌ Skipping: No DBS channels found matching patterns.")
        return

    X_counts_dbs = X_counts[dbs_indices_to_keep, :]
    ext_headers_dbs = [ext_headers_all[i] for i in dbs_indices_to_keep]
    original_elecs = np.array([original_elec_ids[i] for i in dbs_indices_to_keep], dtype=np.uint16).reshape(-1, 1)

    RawVoltage = _scale_to_volts(X_counts_dbs, ext_headers_dbs) 
    
    if APPLY_LEGACY_0_25_SCALING:
        RawVoltage = RawVoltage * 0.25
        print("  - Applied legacy 0.25 voltage scaling factor.")
        
    info = mne.create_info(ch_names=ch_names_dbs, sfreq=fs, ch_types='dbs')
    raw = mne.io.RawArray(RawVoltage, info, verbose=False)
    
    print(f"  - Loaded **{raw.info['nchan']} DBS channels** at {raw.info['sfreq']} Hz.")

    # --- STEP 2: Save master/RawVoltage (omitted details) ---
    master_block_dir = os.path.join(MASTER_DATA_ROOT, block_name)
    os.makedirs(master_block_dir, exist_ok=True)
    raw_vol_path = os.path.join(master_block_dir, f"RawVoltage_{block_name}.mat")
    sio.savemat(raw_vol_path, {'RawVoltage': RawVoltage}, do_compression=True)
    master_elecs = np.arange(1, len(ch_names_dbs) + 1, dtype=np.uint16).reshape(-1, 1)
    bad_channels_matlab = np.array([], dtype=np.uint16).reshape(1, -1) 
    master_vars_data: Dict[str, Any] = {'ecog_srate': np.array([[int(fs)]]), 'pdiode_srate': np.array([[30000]]), 'compress': np.array([[2]]), 'badchan': bad_channels_matlab, 'originalelecs': original_elecs, 'nchan': np.array([[len(ch_names_dbs)]]), 'elecs': master_elecs, 'onset': np.array([[]]), 'comments': np.array([[]]),}
    master_vars_path = os.path.join(master_block_dir, f"master_{block_name}.mat")
    sio.savemat(master_vars_path, {'master_vars': master_vars_data}, oned_as='column')
    print(f"  - Saved master/RawVoltage to: {master_block_dir.split(SBJ_NAME)[-1]}")

    # --- STEP 3: Preprocessing (Notch -> CAR -> Downsample) ---
    
    # 3a. Notch Filter (MNE handles filtering efficiently)
    raw.notch_filter(NOTCH_FREQS, verbose=False)
    
    # Save fiData.mat (Post-Notch)
    reref_dir = os.path.join(REREF_DATA_ROOT, simplified_block_name)
    os.makedirs(reref_dir, exist_ok=True)
    fiData_path = os.path.join(reref_dir, f"fiEEG{simplified_block_name}.mat")
    sio.savemat(fiData_path, {'fiData': raw.get_data()}, do_compression=True)
    
    # =================================================================
    # 3b. Manual Common Average Rereference (CAR) per Probe (NUMPY FIX)
    # =================================================================
    
    # 1. Get the channel groups
    probes_dict = _create_probes(raw.ch_names, DBS_CHANNEL_PATTERNS) 
    car_groups = [ch_list for ch_list in probes_dict.values() if len(ch_list) > 0]

    # 2. Extract the NumPy array from the MNE object (Post-Notch)
    data = raw.get_data() 

    # 3. Perform the CAR subtraction for each lead group
    for group in car_groups:
        # Get the indices of the channels in this group
        ch_indices = [raw.ch_names.index(ch) for ch in group]
        
        # Calculate the average reference for this lead (mean across channels, axis=0)
        avg_ref = np.mean(data[ch_indices, :], axis=0)
        
        # Subtract the average from every channel in the group
        data[ch_indices, :] = data[ch_indices, :] - avg_ref

    # 4. Update the MNE Raw object with the Rereferenced data
    raw._data = data 
    print(f"  - Applied Manual Common Average Reference (CAR) for **{len(car_groups)}** leads.")

    # 3c. Downsample (Decimation)
    if raw.info['sfreq'] > TARGET_SFREQ:
        raw.resample(TARGET_SFREQ, verbose=False)
        print(f"  - Downsampled to {TARGET_SFREQ} Hz.")

    # --- STEP 4: Feature Extraction and Save Power.mat ---
    
    power_matrix = np.empty((len(FEATURE_BANDS), raw.info['nchan']))
    
    for i, (band_name, (fmin, fmax)) in enumerate(FEATURE_BANDS.items()):
        raw_band = raw.copy().filter(fmin, fmax, fir_design='firwin', verbose=False)
        raw_band.apply_hilbert(envelope=True, verbose=False)
        band_power_data = raw_band.get_data()**2
        avg_power = np.mean(band_power_data, axis=1)
        power_matrix[i, :] = avg_power

    power_matrix = power_matrix.T 
    
    power_file_name = f"CAR_power_sub-{SBJ_NAME}_task-{simplified_block_name}.mat"
    power_path = os.path.join(POWER_DATA_ROOT, power_file_name)
    
    sio.savemat(power_path, {'power': power_matrix}, do_compression=True)
    print(f"  - Saved final **CAR power data** to: {power_path.split(SBJ_NAME)[-1]}")
    print("---------------------------------------")


# ==============================
# --- 4. BATCH RUN SCRIPT ---
# ==============================
if __name__ == "__main__":
    
    print(f"Starting batch processing for subject: **{SBJ_NAME}**")
    
    # 1. Setup output directories
    os.makedirs(REREF_DATA_ROOT, exist_ok=True)
    os.makedirs(POWER_DATA_ROOT, exist_ok=True)
    os.makedirs(MASTER_DATA_ROOT, exist_ok=True)
    
    file_list = []
    # 2. Iterate through all subdirectories and find NSx files
    for root, dirs, files in os.walk(ORIGINAL_DATA_ROOT):
        for file in files:
            if file.endswith('.ns3') or file.endswith('.ns5'):
                file_list.append(os.path.join(root, file))

    if not file_list:
        print(f"🚫 Error: No .ns3 or .ns5 files found under: {ORIGINAL_DATA_ROOT}")
    
    # 3. Process each file
    for nsx_path in file_list:
        block_folder_name = os.path.basename(os.path.dirname(nsx_path)) 
        process_dbs_pipeline(nsx_path, block_folder_name)

    print(f"\n✅ **Batch Processing Complete** for subject **{SBJ_NAME}**.")