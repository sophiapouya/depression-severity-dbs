"""
Analyze PCA loadings per patient, matching the exact method from pca_regression_seeg.py.

Feature type importance is computed as:
  importance[i] = sum_k( |loading[i,k]| * |regression_coef[k]| )
averaged across LOO folds, then normalized by feature count per method.
This reflects each feature type's actual contribution to PREDICTING depression severity.

Run from project root: python paper_figures/pca_loadings_analysis.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

plt.rcParams.update({
    "svg.fonttype": "none",
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
})
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler

from config import FEATURES_SEEG_GM, BASE_DIR_SEEG

# ── constants ──────────────────────────────────────────────────────
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_loadings")
os.makedirs(OUT_DIR, exist_ok=True)

VARIANCE_THRESHOLD = 0.90
PRESENCE_THRESHOLD = 0.70
NON_FEATURE_COLS = ["patient_id", "catdi_score", "session_name", "time"]

PATIENTS = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
            "DBSTRD010", "DBSTRD011", "DBSTRD014"]

METHOD_MAP = {
    'FFT':       lambda c: 'fft' in c,
    'FOOOF':     lambda c: 'fooof' in c,
    'Hjorth':    lambda c: 'RawHjorth' in c,
    'Sharpwave': lambda c: 'Sharpwave' in c,
}
METHOD_COLORS = {
    'FFT':       '#8ECFC9',
    'FOOOF':     '#7BAFD4',
    'Hjorth':    '#D9D9D9',
    'Sharpwave': '#FFE680',
}

REGION_BASE  = ['ofc', 'vmpfc', 'acc', 'dlpfc', 'amy']
REGION_ORDER = [f"{r}_{s}" for r in REGION_BASE for s in ('left', 'right')]

REGION_COLORS = {
    'ofc_left':    '#AAAAAA',
    'ofc_right':   '#AAAAAA',
    'vmpfc_left':  '#AAAAAA',
    'vmpfc_right': '#AAAAAA',
    'acc_left':    '#AAAAAA',
    'acc_right':   '#AAAAAA',
    'dlpfc_left':  '#AAAAAA',
    'dlpfc_right': '#AAAAAA',
    'amy_left':    '#AAAAAA',
    'amy_right':   '#AAAAAA',
    'unknown':     '#AAAAAA',
}

def region_label(region_lat):
    """'ofc_left' → 'L-ofc', 'ofc_right' → 'R-ofc'"""
    if '_left' in region_lat:
        return 'L-' + region_lat.replace('_left', '')
    if '_right' in region_lat:
        return 'R-' + region_lat.replace('_right', '')
    return region_lat

# 21 subcategories: FFT(6) + FOOOF(2) + Sharpwave(10, both bands) + Hjorth(3)
SUBCATEGORIES = [
    ('FFT delta',                 lambda c: 'fft_delta'      in c),
    ('FFT theta',                 lambda c: 'fft_theta'      in c),
    ('FFT alpha',                 lambda c: 'fft_alpha'      in c),
    ('FFT beta',                  lambda c: 'fft_beta'       in c),
    ('FFT low gamma',             lambda c: 'fft_low_gamma'  in c),
    ('FFT high gamma',            lambda c: 'fft_high_gamma' in c),
    ('FOOOF exponent',            lambda c: 'fooof_a_exp'    in c),
    ('FOOOF offset',              lambda c: 'fooof_a_offset' in c),
    ('SW max prominence (beta)',  lambda c: 'Sharpwave' in c and 'Max_prominence'  in c and 'range_12_30'  in c),
    ('SW max sharpness (beta)',   lambda c: 'Sharpwave' in c and 'Max_sharpness'   in c and 'range_12_30'  in c),
    ('SW mean interval (beta)',   lambda c: 'Sharpwave' in c and 'Mean_interval'   in c and 'range_12_30'  in c),
    ('SW mean prominence (beta)', lambda c: 'Sharpwave' in c and 'Mean_prominence' in c and 'range_12_30'  in c),
    ('SW mean sharpness (beta)',  lambda c: 'Sharpwave' in c and 'Mean_sharpness'  in c and 'range_12_30'  in c),
    ('SW max prominence (hi-γ)',  lambda c: 'Sharpwave' in c and 'Max_prominence'  in c and 'range_70_150' in c),
    ('SW max sharpness (hi-γ)',   lambda c: 'Sharpwave' in c and 'Max_sharpness'   in c and 'range_70_150' in c),
    ('SW mean interval (hi-γ)',   lambda c: 'Sharpwave' in c and 'Mean_interval'   in c and 'range_70_150' in c),
    ('SW mean prominence (hi-γ)', lambda c: 'Sharpwave' in c and 'Mean_prominence' in c and 'range_70_150' in c),
    ('SW mean sharpness (hi-γ)',  lambda c: 'Sharpwave' in c and 'Mean_sharpness'  in c and 'range_70_150' in c),
    ('Hjorth activity',           lambda c: 'RawHjorth_Activity'   in c),
    ('Hjorth mobility',           lambda c: 'RawHjorth_Mobility'   in c),
    ('Hjorth complexity',         lambda c: 'RawHjorth_Complexity' in c),
]


def _lighten(color, factor=0.72):
    """Blend color toward white. factor=1 → white, factor=0 → original."""
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


def draw_violin(ax, data_list, positions, color_list, rng, dot_s=12, ylabel="|loading|×|coef|",
               fontsize=7, dot_color_list=None, body_color_list=None, median_color='#222222'):
    """Violin plots with individual dots on top. Groups with < 2 points show only dots.
    color_list: used for edge color and as fallback for dots/body.
    body_color_list: if provided, used directly as violin fill (no lightening).
    dot_color_list: if provided, overrides dot color per position.
    """
    for idx, (pos, vals, color) in enumerate(zip(positions, data_list, color_list)):
        dot_color  = dot_color_list[idx]  if dot_color_list  is not None else color
        body_color = body_color_list[idx] if body_color_list is not None else _lighten(color)
        if len(vals) >= 2:
            vp = ax.violinplot([vals], positions=[pos], showmedians=True, showextrema=True)
            for body in vp['bodies']:
                body.set_facecolor(body_color)
                body.set_alpha(0.9)
                body.set_edgecolor(color)
                body.set_linewidth(0.4)
            for part in ['cmins', 'cmaxes', 'cbars']:
                if part in vp:
                    vp[part].set_color(color)
                    vp[part].set_linewidth(0.4)
            if 'cmedians' in vp:
                vp['cmedians'].set_color(median_color)
                vp['cmedians'].set_linewidth(0.6)
                vp['cmedians'].set_zorder(5)
        jitter = rng.uniform(-0.1, 0.1, size=len(vals))
        ax.scatter(np.array([pos] * len(vals)) + jitter, vals,
                   color=dot_color, s=dot_s, alpha=0.9, linewidths=0, zorder=3)
    ax.set_ylabel(ylabel, fontsize=fontsize)
    ax.tick_params(labelsize=fontsize)
    ax.grid(axis='y', alpha=0.2, linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# ── load data ──────────────────────────────────────────────────────
print("Loading SEEG greymatter features...")
df = pd.read_csv(FEATURES_SEEG_GM)
df = df.drop(columns=["time"], errors="ignore")
print(f"Loaded {len(df)} sessions, {len(df.columns)} columns\n")

pooled_records = []

# ── process each patient ───────────────────────────────────────────
for patient in PATIENTS:
    print(f"\n{'='*70}")
    print(f"Patient: {patient}")
    print(f"{'='*70}")

    patient_df = df[df["patient_id"] == patient].copy()
    if len(patient_df) == 0:
        print(f"  No data found")
        continue

    print(f"  Sessions: {len(patient_df)}")

    catdi_scores = patient_df["catdi_score"].reset_index(drop=True)
    feature_df = patient_df.drop(columns=[c for c in NON_FEATURE_COLS if c in patient_df.columns])
    feature_df = feature_df.reset_index(drop=True)

    # 70% presence rule
    cutoff = int(np.ceil(PRESENCE_THRESHOLD * feature_df.shape[0]))
    feature_df = feature_df[feature_df.columns[feature_df.notna().sum() >= cutoff]]
    feat_cols = list(feature_df.columns)

    print(f"  Features after 70% presence rule: {len(feat_cols)}")

    # ── LOO CV: compute weighted feature importance per fold ──────────
    fold_importances = []

    for tr_idx, te_idx in LeaveOneOut().split(feature_df):
        data_tr = feature_df.iloc[tr_idx]
        y_tr = catdi_scores.iloc[tr_idx]

        train_mean = data_tr.mean()
        data_tr = data_tr.fillna(train_mean)

        scaler = StandardScaler()
        X_tr = scaler.fit_transform(data_tr)

        pca = PCA(n_components=VARIANCE_THRESHOLD)
        pc_tr = pca.fit_transform(X_tr)

        model = LinearRegression()
        model.fit(pc_tr, y_tr)

        loadings = pca.components_.T
        coef = np.abs(model.coef_)
        importance = np.sum(np.abs(loadings) * coef[np.newaxis, :], axis=1)
        fold_importances.append(importance)

    mean_importance = np.mean(fold_importances, axis=0)

    # ── aggregate by method, normalized by feature count ──────────────
    method_importance = {m: 0.0 for m in METHOD_MAP}
    method_counts = {m: 0 for m in METHOD_MAP}

    for i, feat in enumerate(feat_cols):
        for method, check in METHOD_MAP.items():
            if check(feat):
                method_importance[method] += mean_importance[i]
                method_counts[method] += 1
                break

    method_importance_norm = {
        m: method_importance[m] / method_counts[m] if method_counts[m] > 0 else 0.0
        for m in METHOD_MAP
    }

    print(f"  Feature counts per method: { {m: method_counts[m] for m in METHOD_MAP} }")
    for m, v in method_importance_norm.items():
        print(f"    {m}: {v:.4f}")

    # ── also fit PCA on all data (for n_components reporting only) ─────
    all_mean = feature_df.mean()
    feature_filled = feature_df.fillna(all_mean)
    scaler_all = StandardScaler()
    X_all = scaler_all.fit_transform(feature_filled)
    pca_all = PCA(n_components=VARIANCE_THRESHOLD)
    pca_all.fit(X_all)
    print(f"  PCs (all data): {pca_all.n_components_}")

    # ── load metadata: map channel → region_laterality ───────────────
    meta_path = os.path.join(str(BASE_DIR_SEEG), patient, f"{patient}_metadata_greymatter.csv")
    contact_to_region = {}
    if os.path.exists(meta_path):
        meta_df = pd.read_csv(meta_path)
        for _, row in meta_df.iterrows():
            ch   = row['channel_name']
            reg  = row.get('channel_region', 'unknown')
            side = 'left' if str(ch).startswith('L') else ('right' if str(ch).startswith('R') else None)
            contact_to_region[ch] = f"{reg}_{side}" if (reg != 'unknown' and side) else 'unknown'

    # ── aggregate mean_importance by contact, normalize by feature count ─
    contact_names = []
    for feat in feat_cols:
        parts = feat.split("_")
        if len(parts) >= 2 and parts[1].isdigit():
            contact = f"{parts[0]}_{parts[1]}"
        else:
            contact = parts[0]
        contact_names.append(contact)

    contact_importance = {}
    contact_feat_counts = {}
    for i, contact in enumerate(contact_names):
        contact_importance[contact] = contact_importance.get(contact, 0.0) + mean_importance[i]
        contact_feat_counts[contact] = contact_feat_counts.get(contact, 0) + 1

    contact_importance_norm = {
        c: contact_importance[c] / contact_feat_counts[c]
        for c in contact_importance
    }

    print(f"\n  Top 15 CONTACTS by weighted importance:")
    for i, (contact, imp) in enumerate(sorted(contact_importance_norm.items(), key=lambda x: x[1], reverse=True)[:15]):
        reg_lat = contact_to_region.get(contact, 'unmapped')
        print(f"    {i+1:2d}. {contact:20s} {imp:.4f}  ({reg_lat})")

    # ── save contact importance CSV ───────────────────────────────────
    patient_out_dir = os.path.join(OUT_DIR, patient)
    os.makedirs(patient_out_dir, exist_ok=True)

    contact_imp_df = pd.DataFrame([
        {"contact": c, "importance_norm": contact_importance_norm[c],
         "region": contact_to_region.get(c, "unknown"), "n_features": contact_feat_counts[c]}
        for c in contact_importance_norm
    ]).sort_values("importance_norm", ascending=False)
    contact_imp_df.to_csv(os.path.join(patient_out_dir, "contact_importance.csv"), index=False)

    # ── build feature → region_laterality lookup ─────────────────────
    feat_region = []
    for feat in feat_cols:
        parts   = feat.split("_")
        contact = f"{parts[0]}_{parts[1]}" if len(parts) >= 2 and parts[1].isdigit() else parts[0]
        feat_region.append(contact_to_region.get(contact, 'unknown'))

    # ── collect pooled data (excluding DBSTRD011) ─────────────────────
    if patient != 'DBSTRD011':
        for i, feat in enumerate(feat_cols):
            subcat_name = 'other'
            for label, fn in SUBCATEGORIES:
                if fn(feat):
                    subcat_name = label
                    break
            pooled_records.append({
                'importance': mean_importance[i],
                'region':     feat_region[i],
                'subcat':     subcat_name,
            })

    rng = np.random.default_rng(42)

    # ── plot: per-region violin with all individual feature points ────
    region_vals_pat = {r: [] for r in REGION_ORDER}
    for i in range(len(feat_cols)):
        r = feat_region[i]
        if r in region_vals_pat:
            region_vals_pat[r].append(mean_importance[i])

    present_pat   = [r for r in REGION_ORDER if len(region_vals_pat[r]) > 0]
    positions_pat = list(range(1, len(present_pat) + 1))

    fig, ax = plt.subplots(figsize=(max(7, len(present_pat) * 0.8), 4))
    draw_violin(ax,
                [region_vals_pat[r] for r in present_pat],
                positions_pat,
                [REGION_COLORS.get(r, REGION_COLORS['unknown']) for r in present_pat],
                rng, dot_s=10, fontsize=9)
    ax.set_xticks(positions_pat)
    ax.set_xticklabels([region_label(r) for r in present_pat], fontsize=9)
    ax.set_ylabel("|loading|×|coef|", fontsize=9)
    plt.tight_layout()
    fig.savefig(os.path.join(patient_out_dir, "region_importance_violin.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(patient_out_dir, "region_importance_violin.svg"), format="svg", bbox_inches="tight")
    plt.close()

    # ── plot: 3×7 grid — one violin per feature subcategory ──────────
    nrow, ncol = 3, 7
    fig, axes = plt.subplots(nrow, ncol, figsize=(28, 12), sharey=True)
    axes = axes.flatten()

    for ax_idx, (subcat_label, subcat_fn) in enumerate(SUBCATEGORIES):
        ax = axes[ax_idx]

        region_vals = {r: [] for r in REGION_ORDER}
        for i, feat in enumerate(feat_cols):
            if subcat_fn(feat) and feat_region[i] in region_vals:
                region_vals[feat_region[i]].append(mean_importance[i])

        present = [r for r in REGION_ORDER if len(region_vals[r]) > 0]
        if not present:
            ax.set_visible(False)
            continue

        positions        = list(range(1, len(present) + 1))
        data_by_region   = [region_vals[r] for r in present]
        colors_by_region = [REGION_COLORS.get(r, REGION_COLORS['unknown']) for r in present]

        draw_violin(ax, data_by_region, positions, colors_by_region, rng)
        ax.set_xticks(positions)
        ax.set_xticklabels([region_label(r) for r in present], fontsize=7)
        ax.set_title(subcat_label, fontsize=9)

    fig.suptitle("", fontsize=13, y=1.01)
    plt.tight_layout()
    fig.savefig(os.path.join(patient_out_dir, "subcategory_region_grid.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(patient_out_dir, "subcategory_region_grid.svg"), format="svg", bbox_inches="tight")
    plt.close()

    # ── plot: stacked panels, one per region+laterality, 21 subcategory violins ─
    subcat_labels = [s[0] for s in SUBCATEGORIES]
    n_subcat      = len(SUBCATEGORIES)
    x_positions   = list(range(1, n_subcat + 1))

    present_regions = [r for r in REGION_ORDER
                       if any(feat_region[i] == r for i in range(len(feat_cols)))]

    fig, axes = plt.subplots(len(present_regions), 1,
                              figsize=(24, 5 * len(present_regions)),
                              sharey=True)
    if len(present_regions) == 1:
        axes = [axes]

    for ri, region in enumerate(present_regions):
        ax = axes[ri]

        region_subcat_data = []
        for _, subcat_fn in SUBCATEGORIES:
            vals = [mean_importance[i] for i, feat in enumerate(feat_cols)
                    if subcat_fn(feat) and feat_region[i] == region]
            region_subcat_data.append(vals)

        non_empty_positions = [x_positions[k] for k, v in enumerate(region_subcat_data) if len(v) > 0]
        non_empty_data      = [v for v in region_subcat_data if len(v) > 0]

        if not non_empty_positions:
            ax.set_visible(False)
            continue

        region_color = REGION_COLORS.get(region, REGION_COLORS['unknown'])
        draw_violin(ax, non_empty_data, non_empty_positions,
                    [region_color] * len(non_empty_data), rng, fontsize=8)
        ax.set_xticks(x_positions)
        ax.set_xticklabels(subcat_labels, rotation=45, ha='right', fontsize=7)
        ax.set_xlim(0.5, n_subcat + 0.5)
        ax.set_title(region_label(region), fontsize=11,
                     color=region_color)

    fig.suptitle("", fontsize=13)
    plt.tight_layout()
    fig.savefig(os.path.join(patient_out_dir, "region_feature_stacked.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(patient_out_dir, "region_feature_stacked.svg"), format="svg", bbox_inches="tight")
    plt.close()

    print(f"\n  Saved to: {patient_out_dir}/")


# ── pooled plots (all patients except DBSTRD011) ──────────────────────────────
print(f"\n{'='*70}")
print("Generating pooled plots (all patients except DBSTRD011)...")
print(f"{'='*70}")

pooled_out_dir = os.path.join(OUT_DIR, "pooled")
os.makedirs(pooled_out_dir, exist_ok=True)

pooled_df = pd.DataFrame(pooled_records)
rng = np.random.default_rng(42)

# sort regions by descending median importance across all pooled features
_region_medians = {
    r: pooled_df[pooled_df['region'] == r]['importance'].median()
    for r in REGION_ORDER
    if (pooled_df['region'] == r).any()
}
pooled_region_order = sorted(_region_medians, key=_region_medians.get, reverse=True)

# shared y-axis range for pooled region + top10 plots
_pooled_ymin = pooled_df['importance'].min()
_pooled_ymax = pooled_df['importance'].max()
_pooled_ypad = (_pooled_ymax - _pooled_ymin) * 0.05
POOLED_YLIM  = (_pooled_ymin - _pooled_ypad, _pooled_ymax + _pooled_ypad)

# ── pooled plot 1: per-region violin with all points ─────────────────────────
PINK_BODY    = '#fff2f8'   # inner violin fill
PINK_EDGE    = '#ff80aa'   # outer violin edge and dots
MEDIAN_BLACK = 'black'

fig, ax = plt.subplots(figsize=(7, 3))
region_vals_p = {r: pooled_df[pooled_df['region'] == r]['importance'].tolist() for r in pooled_region_order}
present_p     = [r for r in pooled_region_order if len(region_vals_p[r]) > 0]
positions_p   = list(range(1, len(present_p) + 1))
draw_violin(ax,
            [region_vals_p[r] for r in present_p],
            positions_p,
            [PINK_EDGE] * len(present_p),
            rng, dot_s=6, fontsize=10,
            dot_color_list=[PINK_EDGE] * len(present_p),
            body_color_list=[PINK_BODY] * len(present_p),
            median_color=MEDIAN_BLACK)
ax.set_xticks(positions_p)
ax.set_xticklabels([region_label(r) for r in present_p], fontsize=10)
ax.set_xlim(0.3, len(present_p) + 0.7)
ax.set_ylim(*POOLED_YLIM)
ax.set_ylabel("|loading|×|coef|", fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(pooled_out_dir, "region_importance_violin.png"), dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(pooled_out_dir, "region_importance_violin.svg"), format="svg", bbox_inches="tight")
plt.close()

# ── pooled plot 2: 3×7 grid — one violin per subcategory ─────────────────────
nrow, ncol = 3, 7
fig, axes = plt.subplots(nrow, ncol, figsize=(28, 12), sharey=True)
axes = axes.flatten()

for ax_idx, (subcat_label, _) in enumerate(SUBCATEGORIES):
    ax = axes[ax_idx]

    sub_df = pooled_df[pooled_df['subcat'] == subcat_label]
    region_vals = {r: sub_df[sub_df['region'] == r]['importance'].tolist() for r in REGION_ORDER}

    present = [r for r in REGION_ORDER if len(region_vals[r]) > 0]
    if not present:
        ax.set_visible(False)
        continue

    positions        = list(range(1, len(present) + 1))
    data_by_region   = [region_vals[r] for r in present]
    colors_by_region = [REGION_COLORS.get(r, REGION_COLORS['unknown']) for r in present]

    draw_violin(ax, data_by_region, positions, colors_by_region, rng)
    ax.set_xticks(positions)
    ax.set_xticklabels([region_label(r) for r in present], fontsize=7)
    ax.set_title(subcat_label, fontsize=9)

fig.suptitle("", fontsize=13, y=1.01)
plt.tight_layout()
fig.savefig(os.path.join(pooled_out_dir, "subcategory_region_grid.png"), dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(pooled_out_dir, "subcategory_region_grid.svg"), format="svg", bbox_inches="tight")
plt.close()

# ── pooled plot 3: stacked panels, one per region+laterality ─────────────────
fig, axes = plt.subplots(len(pooled_region_order), 1,
                          figsize=(24, 5 * len(pooled_region_order)),
                          sharey=True)
if len(pooled_region_order) == 1:
    axes = [axes]

for ri, region in enumerate(pooled_region_order):
    ax = axes[ri]

    subcat_items = []
    for subcat_label, _ in SUBCATEGORIES:
        mask = (pooled_df['region'] == region) & (pooled_df['subcat'] == subcat_label)
        vals = pooled_df[mask]['importance'].tolist()
        if vals:
            subcat_items.append((subcat_label, vals, np.mean(vals)))

    if not subcat_items:
        ax.set_visible(False)
        continue

    subcat_items.sort(key=lambda x: x[2], reverse=True)
    sorted_labels = [s[0] for s in subcat_items]
    sorted_data   = [s[1] for s in subcat_items]
    positions     = list(range(1, len(subcat_items) + 1))

    region_color = REGION_COLORS.get(region, REGION_COLORS['unknown'])
    draw_violin(ax, sorted_data, positions,
                [region_color] * len(sorted_data), rng, fontsize=8)
    ax.set_xticks(positions)
    ax.set_xticklabels(sorted_labels, rotation=45, ha='right', fontsize=7)
    ax.set_xlim(0.5, len(subcat_items) + 0.5)
    ax.set_title(region_label(region), fontsize=11,
                 color=region_color)

plt.tight_layout()
fig.savefig(os.path.join(pooled_out_dir, "region_feature_stacked.png"), dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(pooled_out_dir, "region_feature_stacked.svg"), format="svg", bbox_inches="tight")
plt.close()

# ── pooled plot 4: top-10 subcategories across ALL regions combined ───────────
subcat_items_all = []
for subcat_label, _ in SUBCATEGORIES:
    vals = pooled_df[pooled_df['subcat'] == subcat_label]['importance'].tolist()
    if vals:
        subcat_items_all.append((subcat_label, vals, np.median(vals)))

subcat_items_all.sort(key=lambda x: x[2], reverse=True)
top10 = subcat_items_all[:10]

def _wrap_label(label, max_len=13):
    if len(label) <= max_len:
        return label
    if ' (' in label:
        return label.replace(' (', '\n(')
    mid = len(label) // 2
    left  = label.rfind(' ', 0, mid + 1)
    right = label.find(' ', mid)
    if left == -1 and right == -1:
        return label
    if left == -1:   idx = right
    elif right == -1: idx = left
    else: idx = left if (mid - left) <= (right - mid) else right
    return label[:idx] + '\n' + label[idx + 1:]

top10_labels = [s[0] for s in top10]
top10_data   = [s[1] for s in top10]

fig, ax = plt.subplots(figsize=(7, 3))
positions = list(range(1, len(top10) + 1))
draw_violin(ax, top10_data, positions, [PINK_EDGE] * len(top10), rng, dot_s=6, fontsize=9,
            dot_color_list=[PINK_EDGE] * len(top10),
            body_color_list=[PINK_BODY] * len(top10),
            median_color=MEDIAN_BLACK)
ax.set_xticks(positions)
ax.set_xticklabels([_wrap_label(l) for l in top10_labels], rotation=0, ha='center', fontsize=8)
ax.set_xlim(0.3, len(top10) + 0.7)
ax.set_ylim(*POOLED_YLIM)
ax.set_ylabel("|loading|×|coef|", fontsize=10)
plt.tight_layout()
fig.savefig(os.path.join(pooled_out_dir, "top10_all_regions.png"), dpi=200, bbox_inches="tight")
fig.savefig(os.path.join(pooled_out_dir, "top10_all_regions.svg"), format="svg", bbox_inches="tight")
plt.close()

print(f"\nPooled plots saved to: {pooled_out_dir}/")

print(f"\n{'='*70}")
print("Done! All results saved to: paper_figures/pca_loadings/")
print(f"{'='*70}")
