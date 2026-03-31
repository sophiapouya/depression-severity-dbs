import os
import mne
import numpy as np
import pandas as pd
from scipy.signal import filtfilt, butter, hilbert

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
                            mean_power = np.mean(power)
                            window_means.append(mean_power)

                    # If a channel/band has no windows (should be rare), store NaN
                    if len(window_means) == 0:
                        session_dict[band] = np.nan
                    else:
                        session_mean_power = np.mean(window_means)
                        session_dict[band] = np.log10(session_mean_power)

                master_list.append(session_dict)

    df = pd.DataFrame(master_list)
    df.to_csv(os.path.join(power_dir, f"{subj_name}_{ref_type}_power.csv"), index=False)


