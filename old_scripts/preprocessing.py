# store all path variables here


# save master variables
    # master variables must include:
        # .ecog_srate
        # .piode_srate
        # .compress
        # .bad_chans for each patient
        # montageInfo -> from blackrockMontageInfo_modified
        # .originalelecs
        # .probes
        # .nchan
        # .elecs

        # open NS3 file and process it (with NMPK library)
        # save rawvoltage by multiplying by 0.25 after reading in NS3.Data(master_vars.originalelecs)

# save timestamps
    
# visual inspection for artifacts
# channels excluded 
# blocks excluded
# each channel notch filtered
# each channel rereferenced through common average referencing (CAR)
# down sampling to 1000hz
# hilbert transform to estimate spectral power features in 6 bands -> 
# log transformed average power during catdi test for each freq band


#!/usr/bin/env python

## for large files can run in terminal with:
## (base) katkab@MacBook-Pro-(2) py_neuromodulation % /Users/katkab/Documents/GitHub/py_neuromodulation/.venv/bin/python \ /Users/katkab/Documents/GitHub/py_neuromodulation/clean_raw_neural_data_preprocess_replot_extract_features_memory_efficient_10_13_25.py 2>&1 | tee -a run.log
"""
NSX (.ns3) → Clean (pass-1) → Preview (pass-2; CAR+filters VIEW-ONLY) → Merge → Features (chunked, incremental)

This script can:
1) Start from raw .ns3 when no clean FIF exists (first-pass cleaning + save).
2) Resume from an existing clean, unfiltered multi-part FIF (second-pass preview/QC).
In both cases, it proceeds to feature extraction with py_neuromodulation.

Outputs (in OUTPUT_DIR):
- *_clean_annotated_ieeg.fif (+ -1/-2/-3 parts if large)  [UNFILTERED; bads+annotations only]
- PSD PNGs (initialQC / preprocPreview)
- *_bad_channels_snapshot.csv, *_annotations_snapshot.csv, *_annotations_summary.csv (after each pass)
- sub_FEATURES.csv  (incrementally appended; memory-safe)
- sub_channels.csv, sub_SETTINGS.yaml, sub_SIDECAR.json
"""

import os
import re
import gc
import json
import fnmatch
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne
from brpylib import NsxFile, NevFile
import py_neuromodulation as nm

# ============================== USER CONFIG ==================================


NONNEURAL_CHANS = ['Empty*','empty*','Fz*','Cz*','Pz*','Ref*', 'x*','CL*','ZL*','E*']
DBS_CHANS = ['dLVCVS*','dLSCC*','dRVCVS*','dRSCC*']

ECOG_SRATE = 2000
PDIODE_SRATE = 30000
COMPRESS = 2

# PSD preview slice for speed; None = use full (of the preview window)
PSD_SEGMENT_SEC = None


# ============================== HELPERS ======================================
def _get(d, keys, default=None):
    for k in keys:
        if k in d and d[k] is not None:
            return d[k]
    return default

def _compute_fs(basic_header):
    period = _get(basic_header, ["Period", "period"])
    if period and float(period) != 0:
        return float(30000.0 / float(period))  # 30 kHz base clock
    return int(_get(basic_header, ["Fs","SamplingFreq","SamplingFrequency","SamplingRate"], 2000))

_name_ok_pat = re.compile(r"[^A-Za-z0-9_\-\.]+")
def _sanitize_label(x: str) -> str:
    if x is None:
        return ""
    if isinstance(x, (bytes, bytearray)):
        x = x.decode("utf-8", "ignore")
    x = x.replace("\x00", "").strip()
    x = _name_ok_pat.sub("", x)
    return x if x else ""

def _unique(names, fallback_ids):
    seen = {}
    out = []
    for i, nm in enumerate(names):
        base = nm if nm else (str(fallback_ids[i]) if fallback_ids is not None else f"ch{i+1:03d}")
        cand = base; c = 1
        while cand in seen:
            c += 1; cand = f"{base}_{c}"
        seen[cand] = True; out.append(cand)
    return out

def _scale_to_volts(X_counts, ext_headers):
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

def _labels_from_ext_headers(ext_headers, n_ch, elec_ids):
    labels = []
    for i in range(n_ch):
        eh = ext_headers[i] if (ext_headers and i < len(ext_headers)) else {}
        lab = _get(eh, ["ElectrodeLabel","Label","ChanLabel","ChannelLabel","Name","name"])
        labels.append(_sanitize_label(lab))
    return _unique(labels, elec_ids)

def _annotations_to_df(ann: mne.Annotations) -> pd.DataFrame:
    df = pd.DataFrame({"onset_s": ann.onset, "duration_s": ann.duration, "description": ann.description})
    if ann.orig_time is not None:
        if isinstance(ann.orig_time, tuple):
            base_posix = ann.orig_time[0] + ann.orig_time[1] * 1e-6
        elif isinstance(ann.orig_time, datetime):
            base_posix = ann.orig_time.timestamp()
        else:
            base_posix = float(ann.orig_time)
        df["onset_abs_utc"] = [
            datetime.fromtimestamp(base_posix + o, tz=timezone.utc).isoformat()
            for o in df["onset_s"].to_list()
        ]
    return df

def _save_snapshots(raw_obj, base_name, out_dir):
    bad_csv = os.path.join(out_dir, f"{base_name}_bad_channels_snapshot.csv")
    ann_csv = os.path.join(out_dir, f"{base_name}_annotations_snapshot.csv")
    ann_summary_csv = os.path.join(out_dir, f"{base_name}_annotations_summary.csv")

    df_channels = pd.DataFrame({
        "index": range(len(raw_obj.ch_names)),
        "name": raw_obj.ch_names,
        "type": raw_obj.get_channel_types(),
        "is_bad": [nm in set(raw_obj.info["bads"]) for nm in raw_obj.ch_names],
    })
    df_channels.to_csv(bad_csv, index=False)
    print(f"[SNAPSHOT] Bad-channel list → {bad_csv}")

    ann_df = _annotations_to_df(raw_obj.annotations)
    ann_df.to_csv(ann_csv, index=False)
    print(f"[SNAPSHOT] Annotations → {ann_csv}")
    if not ann_df.empty:
        summary = (ann_df.groupby("description", as_index=False)
                         .agg(n_events=("description", "size"),
                              total_duration_s=("duration_s", "sum"))
                         .sort_values(["n_events", "total_duration_s"], ascending=False))
        summary.to_csv(ann_summary_csv, index=False)
        print(f"[SNAPSHOT] Annotations summary → {ann_summary_csv}")
    else:
        print("[SNAPSHOT] No annotations present.")

# ========================== INTERACTIVE PSD (READ-ONLY) ======================
def clickable_psd_show_label(raw_obj, label, out_dir, base_name, segment_sec=None, fmin=1.0, fmax_cap=200.0, follow_cursor=True):
    """
    Interactive PSD for QC (read-only):
      • Click a line → shows label "channel_name (index)" INSIDE the figure.
      • Does NOT print, highlight, or modify bads.
    """
    # Optional crop for PSD
    if segment_sec is not None:
        duration = raw_obj.n_times / raw_obj.info['sfreq']
        if segment_sec < duration:
            raw_for_psd = raw_obj.copy().crop(tmin=0, tmax=segment_sec); seg = f"{int(segment_sec)}s"
        else:
            raw_for_psd = raw_obj; seg = "full"
    else:
        raw_for_psd, seg = raw_obj, "full"

    fmax = min(fmax_cap, raw_for_psd.info["sfreq"]/2 - 1.0)
    n_fft = int(2 ** np.ceil(np.log2(raw_for_psd.info["sfreq"] * 2.0)))

    psd = raw_for_psd.compute_psd(method="welch", fmin=fmin, fmax=fmax, n_fft=n_fft, picks="seeg")
    freqs = psd.freqs
    spectra = psd.get_data()
    ch_names = list(psd.ch_names) if hasattr(psd, "ch_names") else raw_for_psd.ch_names

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title(f"PSD ({label}) — click a line to show channel name/index; close to continue", fontsize=12)
    ax.set_xlabel("Frequency (Hz)"); ax.set_ylabel("Power Spectral Density (dB)")
    ax.set_xlim([fmin, fmax])

    line_to_idx = {}
    eps = np.finfo(float).tiny
    for i, spec in enumerate(spectra):
        y = 10.0 * np.log10(np.maximum(spec, eps))
        ln, = ax.plot(freqs, y, lw=0.6, alpha=0.9, picker=5)
        line_to_idx[ln] = i

    psd_png = os.path.join(out_dir, f"{base_name}_PSD_{label}_{seg}.png")
    fig.savefig(psd_png, dpi=150, bbox_inches="tight")
    print(f"[PSD] Saved (snapshot) → {psd_png}")

    ann = ax.annotate("", xy=(0,0), xytext=(12,12), textcoords="offset points",
                      bbox=dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9),
                      arrowprops=dict(arrowstyle="->", alpha=0.6), fontsize=10)
    ann.set_visible(False)

    def on_pick(event):
        ln = event.artist
        i = line_to_idx.get(ln, None)
        if i is None:
            return
        if follow_cursor and event.mouseevent.xdata is not None and event.mouseevent.ydata is not None:
            mx, my = float(event.mouseevent.xdata), float(event.mouseevent.ydata)
        else:
            xdata = ln.get_xdata(); ydata = ln.get_ydata()
            mid = len(xdata) // 2
            mx, my = float(xdata[mid]), float(ydata[mid])
        ann.xy = (mx, my)
        ann.set_text(f"{ch_names[i]}  (index {i})")
        ann.set_visible(True)
        fig.canvas.draw_idle()

    fig.canvas.mpl_connect('pick_event', on_pick)
    plt.show(block=True)
    plt.close(fig)

# ====================== PREVIEW PREPROCESSING (COPY ONLY) ====================


# ========================== STAGE 1: FIRST-PASS CLEAN ========================
def stage1_clean_from_nsx():
    print("Reading NSX via brpylib…")
    nsx = NsxFile(NSX_PATH)
    # nev = NevFile(NEV_PATH)

    try:
        # # get the comments from the nev file to find the onset sample
        # nev_data = nev.getdata()
        # comments = nev_data.get("comments")

        # timestamps_list = comments.get("TimeStamps")

        # if comments["TimeStampsStarted"] is not None:
        #     if NS3 == True:
        #         ONSET_SAMPLE = int(timestamps_list[0]/(PDIODE_SRATE/ECOG_SRATE)) - 1
        #     else:
        #         ONSET_SAMPLE = int(timestamps_list[0]) -1
        # else:
        #     ONSET_SAMPLE = 0
        
        # session_metadata = {
        #     "onset_sample": ONSET_SAMPLE,
        #     "compress": COMPRESS,
        #     "ecog_srate": ECOG_SRATE,
        #     "pdiode_srate": PDIODE_SRATE
        # }

        bh = getattr(nsx, "basic_header", {}) or {}
        fs = _compute_fs(bh)
        data_dict = nsx.getdata()
        arr = data_dict.get("data")
        # Concat segments if presencdt
        if isinstance(arr, list):
            X_counts = np.concatenate([np.asarray(a, dtype=np.float32) for a in arr], axis=1)
        else:
            X_counts = np.asarray(arr, dtype=np.float32)
        elec_ids = data_dict.get("elec_ids") or data_dict.get("electrode_ids")
        if elec_ids is not None:
            elec_ids = list(elec_ids)
        ext_headers = getattr(nsx, "extended_headers", None)
        if ext_headers and len(ext_headers) >= X_counts.shape[0]:
            ch_names = _labels_from_ext_headers(ext_headers, X_counts.shape[0], elec_ids)
        else:
            ch_names = _unique([str(e) for e in elec_ids], elec_ids) if elec_ids is not None \
                       else [f"ch{i+1:03d}" for i in range(X_counts.shape[0])]
    finally:
        nsx.close()

    
    # filter channels to exclude nonneural channels like empty
    ch_names_cleaned, ch_types = [], []
    indices_to_exclude, indices_to_keep = [], []
    dbs_indices, seeg_indices = [], []
    dbs_ch_names, seeg_ch_names = [], []

    for i, ch in enumerate(ch_names):
        match = False
        match_dbs = False
        for ch_nn in NONNEURAL_CHANS:
            if fnmatch.fnmatch(ch, ch_nn): 
                match = True
                indices_to_exclude.append(i)
                break
        if match == True: 
            continue
        else:
            #global list for seeg and dbs
            ch_names_cleaned.append(ch) 
            indices_to_keep.append(i)

            for dbs_chan in DBS_CHANS:
                if fnmatch.fnmatch(ch, dbs_chan):
                    dbs_ch_names.append(ch)
                    dbs_indices.append(i)
                    ch_types.append("dbs")
                    match_dbs = True
                    break
            if match_dbs == True:
                continue
            else:
                seeg_ch_names.append(ch)
                seeg_indices.append(i)
                ch_types.append("seeg")
    
    X_counts = X_counts[indices_to_keep,:]
    
    duration_s = X_counts.shape[1] / fs
    print(f"Loaded NSX: channels={X_counts.shape[0]}, samples={X_counts.shape[1]}, fs≈{fs:.2f} Hz, duration≈{duration_s/60:.2f} min")

    ext_headers_filtered = []
    for i in indices_to_keep:
        ext_headers_filtered.append(ext_headers[i])

    X_volts = _scale_to_volts(X_counts, ext_headers_filtered)
    # do I multiple by .25 here???
    # X_volts = X_volts * 0.25

    info = mne.create_info(ch_names=ch_names_cleaned, sfreq=fs, ch_types=ch_types)
    raw = mne.io.RawArray(X_volts, info, verbose=True)

    print("[QC] Light high-pass at 1 Hz for initial visual QC (view only)…")
    raw_qc = raw.copy().filter(l_freq=1.0, h_freq=None, fir_design="firwin", filter_length="auto", verbose=True)

    print("Interactive PSD (initial QC): click a line to show channel name/index; close to continue…")
    clickable_psd_show_label(raw_qc, label="initialQC", out_dir=OUTPUT_DIR, base_name=BASE_NAME, segment_sec=PSD_SEGMENT_SEC)

    print("Open browser: mark bad channels & add annotations; close when done…")
    raw_qc.plot(block=True, duration=10.0, n_channels=min(32, raw_qc.info["nchan"]),
                show_scrollbars=True, scalings="auto")

    # Save CLEAN (unfiltered) file with bads + annotations
    raw_clean = raw.copy()
    raw_clean.set_annotations(raw_qc.annotations, emit_warning=False)
    raw_clean.info["bads"] = list(set(raw_qc.info["bads"]))

    # #serialize the metadata with the onset sample
    # session_metadata_string = json.dumps(session_metadata)
    # raw_clean.info["description"] = session_metadata_string

    print(f"Saving CLEAN (unfiltered) annotated FIF → {FIF_CLEAN_PATH}")
    raw_clean.save(FIF_CLEAN_PATH, overwrite=True, fmt="single", buffer_size_sec=300.0)

    # Snapshots
    base1 = os.path.splitext(os.path.basename(FIF_CLEAN_PATH))[0]
    _save_snapshots(raw_clean, base1, OUTPUT_DIR)
 
# # ===================== STAGE 2 =====================
# notch filtering for 60, 120, and 180 hz
# FIR filter with hanning windows
# rereference probe data
# pull in onset (if onset is not just 1, onset = onset/compress)
# define 6 canonical freq bands
# decimate the signal (by compress -> 2) with FIR filter
# ecog decomposition with hilbert
# data only chopped when taking the average power -> 


# load in the fif data

# # ====================== STAGE 3 ====================


# ================================== MAIN =====================================
if __name__ == "__main__":

    # -- Paths
    # ## SINGLE USE ##
    # NSX_PATH = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/originalData/DBSTRD006/EMU-001_task-CATDI_date-02-08-2022_time-12-42-20/EMU-001_subj-DBSTRD006_task-CATDI_date-02-08-2022_time-12-42-20.ns3"
    # NEV_PATH = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/originalData/DBSTRD006/EMU-001_task-CATDI_date-02-08-2022_time-12-42-20/EMU-001_subj-DBSTRD006_task-CATDI_date-02-08-2022_time-12-42-20.nev"
    # OUTPUT_DIR = "/Users/sophiapouya/workspace/bcm/CATDI/preprocessing_reworked/preprocessing_results/DBSTRD006"
    # os.makedirs(OUTPUT_DIR, exist_ok=True)
    # BASE_NAME = os.path.splitext(os.path.basename(NSX_PATH))[0]  # e.g., EMU-0104_Convo_NSP-1

    # # Base of the multi-part clean FIF (MNE will auto-load -1.fif, -2.fif, -3.fif)
    # FIF_CLEAN_PATH = os.path.join(OUTPUT_DIR, f"{BASE_NAME}_clean_annotated_ieeg.fif")

    # # SCENARIO HANDLER:

    # ## SINGLE USE ##
    # # If a clean (unfiltered) annotated FIF already exists → skip Stage 1.
    # # Else → build it from NSX (Stage 1).
    # if os.path.exists(FIF_CLEAN_PATH):
    #     print(f"[RESUME] Using existing CLEAN (unfiltered) FIF → {FIF_CLEAN_PATH}")
    # else:
    #     if not os.path.exists(NSX_PATH):
    #         raise FileNotFoundError(f"No clean FIF found and NSX missing.\n  Expected clean: {FIF_CLEAN_PATH}\n  NSX_PATH: {NSX_PATH}")
    #     stage1_clean_from_nsx()

    ## DIRECTORY USE ##
    ####################################################################################################
    PATIENT = "DBSTRD006"
    NEURAL_BASE_DIR = "/Users/sophiapouya/workspace/bcm/CATDI/neuralData/originalData"
    PATIENT_BASE_DIR = os.path.join(NEURAL_BASE_DIR, PATIENT)
    OUTPUT_BASE = "/Users/sophiapouya/workspace/bcm/CATDI/preprocessing_reworked/preprocessing_results/"
    OUTPUT_DIR = os.path.join(OUTPUT_BASE, PATIENT)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for exp_dir in os.listdir(PATIENT_BASE_DIR):
        
        exp_path = os.path.join(PATIENT_BASE_DIR,exp_dir)

        for file in os.listdir(exp_path):
            # find the nev file first
            if file.endswith(".nev") == True:
                BASE_NAME = file.split('.')[0]
                NEV_PATH = os.path.join(exp_path,f"{BASE_NAME}.nev")
                NSX_PATH = os.path.join(exp_path,f"{BASE_NAME}.ns3")
                NS3 = True
                FIF_CLEAN_PATH = os.path.join(OUTPUT_DIR, f"{BASE_NAME}_clean_annotated_ieeg.fif")

                if not os.path.exists(NSX_PATH):
                    NSX_PATH = os.path.join(exp_path,f"{BASE_NAME}.ns5")
                    NS3= False

                if os.path.exists(FIF_CLEAN_PATH):
                    print(f"[RESUME] Using existing CLEAN (unfiltered) FIF → {FIF_CLEAN_PATH}")
                    break
                else:
                    if not os.path.exists(NSX_PATH):
                        raise FileNotFoundError(f"No clean FIF found and NSX missing.\n  Expected clean: {FIF_CLEAN_PATH}\n  NSX_PATH: {NSX_PATH}")
                    stage1_clean_from_nsx()
                    break
            else:
                continue
    

    # Stage 2: preview (second pass) + merge updates into clean FIF
    #stage2_preview_and_merge()

    # Stage 3: features (chunked; incremental CSV)
    #stage3_features_incremental()

    # print("\nDone. CLEAN (unfiltered) FIF, PSDs, snapshots, sidecars, and sub_FEATURES.csv saved in:")
    # print(f"  {OUTPUT_DIR}")

