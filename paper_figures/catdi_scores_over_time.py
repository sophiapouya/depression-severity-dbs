"""
CATDI scores over time — one scatter subplot per patient, single row.
X-axis: calendar day rank (1 = first day of recording).
Run from project root: python paper_figures/catdi_scores_over_time.py
"""
import os
import sys
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from config import CATDI_SCORES

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
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})

# ── constants ─────────────────────────────────────────────────────────────────
CATDI_FILE = str(CATDI_SCORES)
ALL_SUBJS = [
    "DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
    "DBSTRD010", "DBSTRD011", "DBSTRD014",
]
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "catdi_scores_over_time")
os.makedirs(OUT_DIR, exist_ok=True)


def _parse_dates(df, subj):
    """Return a list of datetime.date objects, one per row of df."""
    if subj in ["DBSTRD001", "DBSTRD002"]:
        dates = []
        for val in df["EMUdate"]:
            if isinstance(val, str):
                dates.append(datetime.strptime(val.strip(), "%m/%d/%y").date())
            else:
                dates.append(pd.Timestamp(val).date())
        return dates

    names = df["Name"].astype(str)
    if subj in ["DBSTRD011", "DBSTRD014"]:
        names = names.str.split("task-").str[-1]

    dates = []
    for name in names:
        date_part = name.split("_date-")[1].split("_time-")[0]
        if "-" in date_part:
            dates.append(datetime.strptime(date_part.replace("-", ""), "%m%d%Y").date())
        else:
            dates.append(datetime.strptime(date_part, "%Y%m%d").date())
    return dates


def load_scores_with_days(subj):
    df = pd.read_excel(CATDI_FILE, sheet_name=subj).dropna(subset=["Result"])
    dates = _parse_dates(df, subj)
    df["_date"] = dates

    unique_dates = sorted(set(dates))
    date_to_day = {d: i + 1 for i, d in enumerate(unique_dates)}
    df["day_num"] = df["_date"].map(date_to_day)
    return df[["day_num", "Result"]].dropna().sort_values("day_num")


# ── patient colors (matches all other paper figures) ─────────────────────────
PATIENT_COLORS = {
    subj: plt.cm.tab10(np.linspace(0, 1, len(ALL_SUBJS)))[i]
    for i, subj in enumerate(ALL_SUBJS)
}

# ── load all data first to get shared x range ─────────────────────────────────
rng = np.random.default_rng(42)
patient_data = {}
for subj in ALL_SUBJS:
    try:
        patient_data[subj] = load_scores_with_days(subj)
    except Exception as e:
        print(f"Skipping {subj}: {e}")

max_day = max(df["day_num"].max() for df in patient_data.values())
day_ticks = list(range(1, max_day + 1))

# ── plot ───────────────────────────────────────────────────────────────────────
n = len(ALL_SUBJS)
fig, axes = plt.subplots(
    1, n,
    figsize=(160 / 25.4 * (96 / 72), 1.6),
    gridspec_kw={"wspace": 0.25},
)

for i, subj in enumerate(ALL_SUBJS):
    ax = axes[i]
    if subj not in patient_data:
        ax.set_visible(False)
        continue

    df = patient_data[subj]

    # alternating gray/white column backgrounds
    for d in day_ticks:
        if d % 2 == 0:
            ax.axvspan(d - 0.5, d + 0.5, color="lightgray", alpha=0.4, zorder=0)

    # jittered scatter within each day bucket
    jitter = rng.uniform(-0.25, 0.25, size=len(df))
    ax.scatter(df["day_num"] + jitter, df["Result"],
               s=4, color=PATIENT_COLORS[subj], zorder=3)

    ax.set_xlim(0.5, max_day + 0.5)
    ax.set_xticks(day_ticks)
    ax.set_xticklabels(day_ticks, fontsize=5)
    ax.set_title(subj.replace('DBSTRD', 'TRD'), fontsize=6, fontweight="bold", pad=4)
    ax.set_xlabel("Time (Day)", fontsize=6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    if i == 0:
        ax.set_ylabel("Depression Severity", fontsize=6)
    else:
        ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0)

fig.savefig(os.path.join(OUT_DIR, "catdi_scores_over_time.png"),
            dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(OUT_DIR, "catdi_scores_over_time.svg"),
            format="svg", bbox_inches="tight")
print(f"Saved: {OUT_DIR}/catdi_scores_over_time.png")
plt.show()
