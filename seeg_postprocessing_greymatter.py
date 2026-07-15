"""
Extract power from greymatter bipolar channels.
Reads from bipolar_channels_greymatter/ and saves to the same dir.

Run from project root: python seeg_postprocessing_greymatter.py
"""
import os
import numpy as np
import pandas as pd
from scipy.signal import butter, filtfilt, hilbert
import mne
from config import BASE_DIR_SEEG

base_dir = str(BASE_DIR_SEEG)
#all_subjs = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008", "DBSTRD010", "DBSTRD011", "DBSTRD014"]
all_subjs = ["DBSTRD006"]

EXCLUDED_SESSIONS = {
    "DBSTRD001": ["CATDI_run-08_blk-04"],
    "DBSTRD002": ["CATDI_run-Day5_blk-04", "CATDI_run-Day7_blk-05",
                  "CATDI_run-Day3_blk-02", "CATDI_run-Day3_blk-03"],
    "DBSTRD006": ["CATDI_date-02-08-2022_time-12-42-20"],
    "DBSTRD008": ["CATDI_date-10-25-2022_time-20-20-50","CATDI_date-10-25-2022_time-13-42-28",
                  "CATDI_date-10-25-2022_time-11-27-18","CATDI_date-10-26-2022_time-07-36-57",
                  "CATDI_date-10-26-2022_time-16-16-07","CATDI_date-10-26-2022_time-14-58-43",
                  "CATDI_date-10-25-2022_time-16-51-50","CATDI_date-10-25-2022_time-14-50-44",
                  "CATDI_date-10-25-2022_time-08-23-49","CATDI_date-10-26-2022_time-18-24-14"],
    "DBSTRD010": ["CATDI_date-05-11-2023_time-16-20-04"],
    "DBSTRD011": ["CATDI_date-20240724_time-145011", "CATDI_date-20240719_time-213319",
                  "CATDI_date-20240724_time-111946"],
    "DBSTRD014": ["CATDI_date-20250307_time-130552","CATDI_date-20250312_time-183153",
                  "CATDI_date-20250307_time-195211","CATDI_date-20250308_time-084929"]
}

FEATURE_BANDS = {
    'delta': [1, 4], 'theta': [4, 8], 'alpha': [8, 12],
    'beta': [12, 30], 'low_gamma': [35, 50], 'high_gamma': [70, 150]
}
SFREQ = 1000

band_coefficients = {
    band: butter(N=2, Wn=vals, btype='bandpass', fs=SFREQ)
    for band, vals in FEATURE_BANDS.items()
}

for sbj in all_subjs:
    print(f"\n{sbj}:")
    sbj_dir = os.path.join(base_dir, sbj)
    working_dir = os.path.join(sbj_dir, "bipolar_channels_greymatter")
    power_dir = os.path.join(working_dir, "power_bipolar_greymatter")

    if not os.path.exists(working_dir):
        print(f"  Skipped (no bipolar_channels_greymatter)")
        continue

    os.makedirs(power_dir, exist_ok=True)

    master_list = []
    fif_count = 0
    with os.scandir(working_dir) as files:
        for file in sorted(files, key=lambda f: f.name):
            if not file.name.endswith(".fif"):
                continue

            session_name = os.path.splitext(file.name)[0]

            if session_name in EXCLUDED_SESSIONS.get(sbj, []):
                continue

            fif_count += 1
            raw = mne.io.read_raw_fif(file.path, preload=True, verbose=False)
            data = raw.get_data()
            ch_names = raw.ch_names

            for ch_idx, ch in enumerate(ch_names):
                session_dict = {
                    "session": session_name,
                    "probe": ch.split("_")[0],
                    "ch_name": ch
                }
                x = data[ch_idx, :]

                for band in FEATURE_BANDS:
                    b, a = band_coefficients[band]
                    band_data = filtfilt(b, a, x)
                    analytic = hilbert(band_data)
                    power = np.abs(analytic) ** 2
                    session_dict[band] = np.log10(np.mean(power))

                master_list.append(session_dict)

    if master_list:
        df = pd.DataFrame(master_list)
        out_csv = os.path.join(power_dir, f"{sbj}_bipolar_power_greymatter.csv")
        df.to_csv(out_csv, index=False)
        print(f"  Saved: {fif_count} files, {len(df)} channels")
    else:
        print(f"  No sessions processed")
