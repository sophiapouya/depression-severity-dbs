"""
Pairwise paired t-test for all DBS × SEEG combined PCA conditions.

Each condition is one (DBS region, SEEG combo) pair — 7 × 35 = 245 total.
Per-patient permutation p-value is the observation for each condition.
Prints all significant pairs (p < 0.05) sorted by p-value, and saves to CSV.

Reads from the cache built by dbs_seeg_perm_heatmaps_lateral.py:
  paper_figures/pca_perm_heatmaps_dbs_seeg/_all_results_cache.csv

Run from project root: python paper_figures/ttest_grid_dbs_seeg_perm_pvalue.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
from scipy.stats import ttest_rel

# ── condition ordering (matches dbs_seeg_perm_heatmaps_lateral.py) ────────────
DBS_REGIONS = ["ALL", "LSCC", "RSCC", "LVCVS", "RVCVS", "SCC", "VCVS"]
DBS_LABELS  = {"ALL": "ALL", "LSCC": "L-SCC", "RSCC": "R-SCC",
               "LVCVS": "L-VCVS", "RVCVS": "R-VCVS", "SCC": "SCC", "VCVS": "VCVS"}

SEEG_BASE = ["acc", "amy", "dlpfc", "ofc", "vmpfc"]
SEEG_UNI  = [f"{r}_{s}" for r in SEEG_BASE for s in ("left", "right")]
SEEG_BI   = [f"L{r1}_R{r2}" for r1 in SEEG_BASE for r2 in SEEG_BASE]
SEEG_ALL  = SEEG_UNI + SEEG_BI   # 35

SEEG_LABEL = {}
for s in SEEG_UNI:
    SEEG_LABEL[s] = ("L-" if "_left" in s else "R-") + s.split("_")[0]
for s in SEEG_BI:
    parts = s.split("_")
    SEEG_LABEL[s] = f"L-{parts[0][1:]}/R-{parts[1][1:]}"

CONDITIONS = [(dbs, seeg) for dbs in DBS_REGIONS for seeg in SEEG_ALL]
N = len(CONDITIONS)   # 245

SIG_THRESHOLD = 0.05

CACHE_CSV = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "pca_perm_heatmaps_dbs_seeg", "_all_results_cache.csv",
)
OUT_DIR = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "ttest_grid_dbs_seeg_perm_pvalue"
)
os.makedirs(OUT_DIR, exist_ok=True)

# ── load data ──────────────────────────────────────────────────────────────────
if not os.path.exists(CACHE_CSV):
    raise FileNotFoundError(
        f"Cache not found:\n  {CACHE_CSV}\n"
        "Run paper_figures/dbs_seeg_perm_heatmaps_lateral.py first."
    )

df = pd.read_csv(CACHE_CSV)
df["condition"] = list(zip(df["dbs_region"], df["seeg_combo"]))

wide = df.pivot_table(index="patient_id", columns="condition", values="perm_p_value")
wide = wide.reindex(columns=CONDITIONS)

print(f"{N} conditions × {len(wide)} patients")

# ── pairwise paired t-tests ────────────────────────────────────────────────────
results = []

for i in range(N):
    for j in range(i + 1, N):
        paired = wide[[CONDITIONS[i], CONDITIONS[j]]].dropna()
        if len(paired) >= 2:
            _, p_val = ttest_rel(
                paired[CONDITIONS[i]].values,
                paired[CONDITIONS[j]].values,
            )
            results.append({
                "cond_a":     f"{DBS_LABELS[CONDITIONS[i][0]]} / {SEEG_LABEL[CONDITIONS[i][1]]}",
                "cond_b":     f"{DBS_LABELS[CONDITIONS[j][0]]} / {SEEG_LABEL[CONDITIONS[j][1]]}",
                "p_value":    p_val,
                "n_patients": len(paired),
            })

results_df = pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)
sig_df = results_df[results_df["p_value"] < SIG_THRESHOLD]

n_total = len(results_df)
n_sig   = len(sig_df)
print(f"Significant pairs (p < {SIG_THRESHOLD}): {n_sig} / {n_total}  ({100*n_sig/n_total:.1f}%)\n")

# ── print significant pairs ────────────────────────────────────────────────────
if sig_df.empty:
    print("No significant pairs found.")
else:
    print(f"{'Condition A':<35}  {'Condition B':<35}  {'p-value':>10}  n")
    print("-" * 90)
    for _, row in sig_df.iterrows():
        print(f"{row['cond_a']:<35}  {row['cond_b']:<35}  {row['p_value']:>10.4f}  {int(row['n_patients'])}")

# ── save full results to CSV ───────────────────────────────────────────────────
out_csv = os.path.join(OUT_DIR, "ttest_all_pairs.csv")
results_df.to_csv(out_csv, index=False)
print(f"\nFull results saved: {out_csv}")

out_sig_csv = os.path.join(OUT_DIR, "ttest_significant_pairs.csv")
sig_df.to_csv(out_sig_csv, index=False)
print(f"Significant pairs saved: {out_sig_csv}")
