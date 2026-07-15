"""
All-patient CATDI correlation heatmap for paper.
X-axis: one centered region label per region group, patient IDs below dividers.
Run from project root: python paper_figures/catdi_correlation_paper.py
"""
import os
import sys

# add project root to path so config and src imports work from any subfolder
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import pearsonr
from statsmodels.stats.multitest import fdrcorrection

from config import BASE_DIR_SEEG, CATDI_SCORES

# ── style ──────────────────────────────────────────────────────────────────────
FS = 5
matplotlib.rcParams.update({
    "font.size": FS,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "svg.fonttype": "none",
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": FS,
    "ytick.labelsize": FS,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
})

# ── constants ─────────────────────────────────────────────────────────────────
BASE_DIR = str(BASE_DIR_SEEG)
CATDI_FILE = str(CATDI_SCORES)
ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
BANDS = ["delta", "theta", "alpha", "beta", "low_gamma", "high_gamma"]
BAND_LABELS = ["δ", "θ", "α", "β", "γ", "hγ"]

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "catdi_correlation_seeg")
os.makedirs(OUT_DIR, exist_ok=True)


# ── data loading ───────────────────────────────────────────────────────────────
def load_subject_data(subj):
    meta = pd.read_csv(os.path.join(BASE_DIR, subj, f"{subj}_metadata_greymatter.csv"))
    meta["ch_name"] = meta["channel_name"]

    power = pd.read_csv(os.path.join(
        BASE_DIR, subj, "bipolar_channels_greymatter", "power_bipolar_greymatter",
        f"{subj}_bipolar_power_greymatter.csv",
    ))

    catdi = pd.read_excel(CATDI_FILE, sheet_name=subj)
    if subj in ["DBSTRD011", "DBSTRD014"]:
        catdi["session"] = catdi["Name"].str.split("_task-").str[1]
    else:
        catdi["session"] = catdi["Name"]

    power["session"] = power["session"].astype(str)
    catdi["session"] = catdi["session"].astype(str)

    power = pd.merge(power, catdi[["session", "Result"]], how="left", on="session")
    power = pd.merge(power, meta[["ch_name", "channel_region"]], how="left", on="ch_name")
    power = power.rename(columns={"Result": "score"})
    return power


def compute_correlations(subj_data):
    region_map = dict(zip(subj_data["ch_name"], subj_data["channel_region"]))
    # sort by region, then hemisphere (L before R), then channel name
    contacts = sorted(subj_data["ch_name"].unique(),
                      key=lambda ch: (region_map.get(ch, ""), ch[0].upper() != "L", ch))

    r_data, p_data = {}, {}
    for ch in contacts:
        ch_df = subj_data[subj_data["ch_name"] == ch].dropna(subset=["score"])
        r_vals, p_vals = [], []
        for band in BANDS:
            sub = ch_df[[band, "score"]].dropna()
            if len(sub) < 2 or sub[band].nunique() < 2 or sub["score"].nunique() < 2:
                r_vals.append(float("nan"))
                p_vals.append(1.0)
                continue
            r, p = pearsonr(sub[band], sub["score"])
            r_vals.append(r)
            p_vals.append(p)
        r_data[ch] = r_vals
        p_data[ch] = p_vals

    r_df = pd.DataFrame(r_data, index=BAND_LABELS)
    p_df = pd.DataFrame(p_data, index=BAND_LABELS)
    if r_df.empty:
        return r_df, p_df, region_map

    _, corrected = fdrcorrection(p_df.values.flatten())
    corrected_df = pd.DataFrame(
        corrected.reshape(p_df.shape), index=p_df.index, columns=p_df.columns
    )
    return r_df, corrected_df.map(lambda v: "*" if v < 0.05 else ""), region_map


# ── collect data ───────────────────────────────────────────────────────────────
all_r, all_ann = [], []
patient_labels, patient_centers, patient_boundaries = [], [], []
hemi_label_info   = []  # (center_col, 'L'/'R') for each hemisphere sub-group
region_label_info = []  # (left_col, right_col, region) for each region group
region_inner_bounds = []  # thick-ish lines between regions
hemi_inner_bounds   = []  # thin dotted lines between L/R within a region

running = 0
for subj in ALL_SUBJS:
    try:
        subj_data = load_subject_data(subj)
    except FileNotFoundError:
        print(f"Skipping {subj}: data not found")
        continue

    r_df, ann_df, region_map = compute_correlations(subj_data)
    if r_df.empty:
        continue

    all_r.append(r_df)
    all_ann.append(ann_df)

    contacts = r_df.columns.tolist()
    n = len(contacts)
    abs_offset = running

    patient_centers.append(running + n / 2)
    running += n
    patient_boundaries.append(running)
    patient_labels.append(subj.replace('DBSTRD', 'TRD'))

    prev_region, prev_hemi = None, None
    region_start = abs_offset
    hemi_start   = abs_offset

    for i, ch in enumerate(contacts):
        region = region_map.get(ch, "unknown").upper()
        hemi   = ch[0].upper()  # 'L' or 'R'

        if region != prev_region:
            if prev_hemi is not None:
                hemi_label_info.append(((hemi_start + abs_offset + i) / 2, prev_hemi))
            if prev_region is not None:
                region_label_info.append((region_start, abs_offset + i, prev_region))
                region_inner_bounds.append(abs_offset + i)
            prev_region  = region
            region_start = abs_offset + i
            prev_hemi    = hemi
            hemi_start   = abs_offset + i
        elif hemi != prev_hemi:
            hemi_label_info.append(((hemi_start + abs_offset + i) / 2, prev_hemi))
            hemi_inner_bounds.append(abs_offset + i)
            prev_hemi  = hemi
            hemi_start = abs_offset + i

    # close last groups
    if prev_hemi is not None:
        hemi_label_info.append(((hemi_start + running) / 2, prev_hemi))
    if prev_region is not None:
        region_label_info.append((region_start, running, prev_region))

combined_r = pd.concat(all_r, axis=1)
combined_ann = pd.concat(all_ann, axis=1)

# ── plot ───────────────────────────────────────────────────────────────────────
n_contacts = combined_r.shape[1]
fig, ax = plt.subplots(figsize=(160 / 25.4 * (96 / 72), 27 / 25.4 * (96 / 72)))

sns.heatmap(
    data=combined_r,
    ax=ax,
    cmap="RdBu_r",
    center=0,
    vmin=-1,
    vmax=1,
    annot=combined_ann,
    annot_kws={"color": "black", "fontsize": 5, "fontweight": "bold"},
    fmt="",
    xticklabels=False,
    cbar=False,
)

# thick patient dividers
for b in patient_boundaries[:-1]:
    ax.axvline(b, color="black", linewidth=0.8)

# thin dashed region dividers within patients
for b in region_inner_bounds:
    if b not in patient_boundaries:
        ax.axvline(b, color="black", linewidth=0.3, linestyle="--")

# very thin dotted L/R dividers within regions
hemi_only_bounds = set(hemi_inner_bounds) - set(region_inner_bounds) - set(patient_boundaries)
for b in hemi_only_bounds:
    ax.axvline(b, color="grey", linewidth=0.3, linestyle=":")

# L/R labels (first level, close to axis)
for col_center, hemi in hemi_label_info:
    ax.text(
        col_center, -0.04, hemi,
        transform=ax.get_xaxis_transform(),
        ha="center", va="top",
        fontsize=3.5, clip_on=False,
    )

# region labels — tick from axis then diagonal label centered under L+R group
for left_col, right_col, region in region_label_info:
    center = (left_col + right_col) / 2
    ax.plot([center, center], [-0.09, -0.12],
            transform=ax.get_xaxis_transform(),
            color="black", linewidth=0.5, clip_on=False)
    ax.text(
        center, -0.13, region.lower(),
        transform=ax.get_xaxis_transform(),
        ha="center", va="top",
        fontsize=5, rotation=-45, clip_on=False,
    )

# patient labels above heatmap
for col_center, label in zip(patient_centers, patient_labels):
    ax.text(
        col_center, 1.02, label,
        transform=ax.get_xaxis_transform(),
        ha="center", va="bottom",
        fontsize=7, weight="bold", clip_on=False,
    )

ax.set_yticklabels(ax.get_yticklabels(), fontsize=7, rotation=0)
ax.tick_params(axis="x", length=0)
ax.set_ylabel("Frequency Band", fontsize=7)
fig.subplots_adjust(top=0.88)

fig.savefig(os.path.join(OUT_DIR, "catdi_correlation_seeg_all_patients_greymatter.png"),
            dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT_DIR, "catdi_correlation_seeg_all_patients_greymatter.svg"),
            format="svg", bbox_inches="tight")
print(f"Saved to {OUT_DIR}/catdi_correlation_seeg_all_patients.png")
plt.show()
