"""
Create grey-matter-only versions of SEEG metadata, power, and features CSVs.
Uses same preprocessing but with grey-matter-only metadata.

Run from project root: python create_greymatter_seeg_files.py
"""
import os
import pandas as pd
from src.preprocessing_functions import (
    save_seeg_chans, create_seeg_probes, output_metadata
)
import mne
import re
from config import ROOT_DIR, CATDI_ELECTRODES, BASE_DIR_SEEG, FEATURES_SEEG

DATA_ROOT = str(ROOT_DIR)
PROJECT_NAME = "CATDI"
#ALL_SUBJ = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008", "DBSTRD010", "DBSTRD011", "DBSTRD014"]
ALL_SUBJ = ["DBSTRD006"]
ELECTRODE_FILE = str(CATDI_ELECTRODES)
SEEG_PATTERN = re.compile(r"[A-Za-z]+-[A-Za-z]+\d+-\d{3}")
BASE_DIR_SEEG_STR = str(BASE_DIR_SEEG)

def get_seeg_metadata_greymatter(file_path: str, patient: str) -> dict:
    """Get metadata, marking white matter contacts as 'none'."""
    final_dict = {}
    excel_df = pd.read_excel(file_path, sheet_name=patient)

    clean_excel = excel_df[excel_df['Type'] != "DBS"].copy()
    clean_excel["area"] = clean_excel["area"].fillna("none")
    clean_excel = clean_excel.dropna(subset=["Label"])

    for _, row in clean_excel.iterrows():
        label = str(row["Label"]).strip().lower()
        region = str(row["area"]).strip().lower()

        # Mark white matter contacts as "none"
        if str(row.get("Grey v White", "Grey")).strip() != "Grey":
            region = "none"

        final_dict[label] = region

    return final_dict


for SBJ_NAME in ALL_SUBJ:
    print(f"\n{SBJ_NAME}:")

    SEEG_DATA_ROOT = os.path.join(DATA_ROOT, PROJECT_NAME, 'neuralData', 'seegData', SBJ_NAME)

    # Get grey matter metadata
    seeg_metadata_gm = get_seeg_metadata_greymatter(ELECTRODE_FILE, SBJ_NAME)
    gm_count = sum(1 for v in seeg_metadata_gm.values() if v != "none")
    print(f"  Grey matter contacts: {gm_count}")

    # Load existing raw .fif files and rerun bipolar creation with greymatter metadata
    raw_fif_dir = os.path.join(SEEG_DATA_ROOT, "raw_fif_files")
    fif_files = [f for f in os.listdir(raw_fif_dir) if f.endswith(".fif")]

    output_data_gm = None
    for fif_file in fif_files:
        raw_seeg = mne.io.read_raw_fif(os.path.join(raw_fif_dir, fif_file),
                                       preload=False, verbose=False)
        simplified_block_name = os.path.splitext(fif_file)[0]

        probes = create_seeg_probes(raw_file=raw_seeg)
        bipolar_dir_gm = os.path.join(SEEG_DATA_ROOT, "bipolar_channels_greymatter")
        os.makedirs(bipolar_dir_gm, exist_ok=True)

        output_data_gm = save_seeg_chans(
            probes=probes, raw_data=raw_seeg,
            block_name=simplified_block_name,
            save_dir=bipolar_dir_gm,
            seeg_metadata=seeg_metadata_gm
        )

    # Save metadata
    if output_data_gm:
        meta_out = os.path.join(SEEG_DATA_ROOT, f"{SBJ_NAME}_metadata_greymatter.csv")
        output_metadata(meta_out, output_data_gm)

        # Count ACC
        meta_gm = pd.read_csv(meta_out)
        acc_all = [ch for ch in meta_gm["channel_name"] if "ACC" in ch]
        print(f"  ACC channels: {len(acc_all)}")
        print(f"  Total grey matter channels: {len(meta_gm)}")


        # Filter features
        features = pd.read_csv(str(FEATURES_SEEG))
        features_subj = features[features["patient_id"] == SBJ_NAME]

        gm_feature_cols = ["patient_id", "session_name", "catdi_score"]
        for col in features_subj.columns:
            if col not in gm_feature_cols:
                for ch in meta_gm["channel_name"]:
                    if col.startswith(ch + "_"):
                        gm_feature_cols.append(col)
                        break

        if len(gm_feature_cols) > 3:
            print(f"  Features: {len(features_subj.columns)-3} → {len(gm_feature_cols)-3}")



