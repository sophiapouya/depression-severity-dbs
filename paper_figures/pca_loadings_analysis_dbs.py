"""
Analyze PCA loadings per patient for DBS contacts (ALL and RVCVS probe types).

Feature type importance is computed as:
  importance[i] = sum_k( |loading[i,k]| * |regression_coef[k]| )
averaged across LOO folds, then normalized by feature count per method/contact.
This reflects each feature type's actual contribution to PREDICTING depression severity.

Contact region is derived directly from the contact name prefix (LSCC, RSCC, LVCVS, RVCVS).

Run from project root: python paper_figures/pca_loadings_analysis_dbs.py
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

from config import FEATURES_DBS

# ── constants ──────────────────────────────────────────────────────
BASE_OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pca_loadings_dbs")

VARIANCE_THRESHOLD = 0.90
PRESENCE_THRESHOLD = 0.70
NON_FEATURE_COLS   = ["patient_id", "catdi_score", "session_name", "time"]

PATIENTS = ["DBSTRD001", "DBSTRD002", "DBSTRD006", "DBSTRD008",
            "DBSTRD010", "DBSTRD011", "DBSTRD014"]

PROBE_TYPES = ["ALL"]

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

# individual contact prefix colors (used for top-contacts bar chart)
REGION_COLORS = {
    'LSCC':    '#D4B896',
    'RSCC':    '#D4B896',
    'LVCVS':   '#4A9E68',
    'RVCVS':   '#4A9E68',
    'unknown': '#CCCCCC',
}

# DBS electrode groupings: each feature belongs to all groups whose fn(prefix) is True
DBS_GROUP_ORDER = ['all', 'left', 'right', 'scc', 'vcvs', 'lscc', 'rscc', 'lvcvs', 'rvcvs']
DBS_GROUP_FNS = {
    'all':   lambda p: True,
    'left':  lambda p: p in ('LSCC', 'LVCVS'),
    'right': lambda p: p in ('RSCC', 'RVCVS'),
    'scc':   lambda p: p in ('LSCC', 'RSCC'),
    'vcvs':  lambda p: p in ('LVCVS', 'RVCVS'),
    'lscc':  lambda p: p == 'LSCC',
    'rscc':  lambda p: p == 'RSCC',
    'lvcvs': lambda p: p == 'LVCVS',
    'rvcvs': lambda p: p == 'RVCVS',
}
DBS_GROUP_COLORS = {
    'all':   '#AAAAAA',
    'left':  '#AAAAAA',
    'right': '#AAAAAA',
    'scc':   '#D4B896',
    'vcvs':  '#4A9E68',
    'lscc':  '#D4B896',
    'rscc':  '#D4B896',
    'lvcvs': '#4A9E68',
    'rvcvs': '#4A9E68',
}

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


def get_contact_and_region(feat):
    """'LVCVS_1_RawHjorth_Activity' → contact='LVCVS_1', region='LVCVS'"""
    parts = feat.split("_")
    if len(parts) >= 2 and parts[1].isdigit():
        contact = f"{parts[0]}_{parts[1]}"
        region  = parts[0]
    else:
        contact = parts[0]
        region  = parts[0]
    return contact, region


def _lighten(color, factor=0.72):
    r, g, b = mcolors.to_rgb(color)
    return (r + (1 - r) * factor, g + (1 - g) * factor, b + (1 - b) * factor)


def _darken(color, factor=0.5):
    r, g, b = mcolors.to_rgb(color)
    return (r * factor, g * factor, b * factor)


PINK_BODY    = '#fff2f8'   # inner violin fill
PINK_EDGE    = '#ff80aa'   # outer violin edge and dots
MEDIAN_BLACK = 'black'


def draw_violin(ax, data_list, positions, color_list, rng,
                dot_s=12, ylabel="|loading|×|coef|", fontsize=7,
                dot_color_list=None, body_color_list=None, median_color='#222222'):
    """Violin plots with individual dots on top. Groups with < 2 points show only dots.
    color_list: edge color and fallback for dots/body.
    body_color_list: direct fill color (no lightening applied).
    dot_color_list: per-position dot color override.
    """
    for idx, (pos, vals, color) in enumerate(zip(positions, data_list, color_list)):
        dot_color  = dot_color_list[idx]  if dot_color_list  is not None else _darken(color)
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
print("Loading DBS features...")
df = pd.read_csv(str(FEATURES_DBS))
df = df.drop(columns=["time"], errors="ignore")
print(f"Loaded {len(df)} sessions, {len(df.columns)} columns\n")

# ── loop over probe types ──────────────────────────────────────────
for probe_type in PROBE_TYPES:
    print(f"\n{'#'*70}")
    print(f"PROBE TYPE: {probe_type}")
    print(f"{'#'*70}")

    out_dir = os.path.join(BASE_OUT_DIR, probe_type)
    os.makedirs(out_dir, exist_ok=True)

    pooled_records = []

    for patient in PATIENTS:
        print(f"\n{'='*70}")
        print(f"Patient: {patient}  |  Probe: {probe_type}")
        print(f"{'='*70}")

        patient_df = df[df["patient_id"] == patient].copy()
        if len(patient_df) == 0:
            print(f"  No data found")
            continue

        # Apply probe filter
        if probe_type != "ALL":
            feat_cols_all = [c for c in patient_df.columns if c not in NON_FEATURE_COLS]
            keep_feats = [c for c in feat_cols_all if probe_type in c]
            keep_cols  = [c for c in NON_FEATURE_COLS if c in patient_df.columns] + keep_feats
            patient_df = patient_df[keep_cols]

        print(f"  Sessions: {len(patient_df)}")

        catdi_scores = patient_df["catdi_score"].reset_index(drop=True)
        feature_df   = patient_df.drop(columns=[c for c in NON_FEATURE_COLS if c in patient_df.columns])
        feature_df   = feature_df.reset_index(drop=True)

        # 70% presence rule
        cutoff = int(np.ceil(PRESENCE_THRESHOLD * feature_df.shape[0]))
        feature_df = feature_df[feature_df.columns[feature_df.notna().sum() >= cutoff]]
        feat_cols  = list(feature_df.columns)

        if len(feat_cols) == 0:
            print(f"  No features after 70% presence rule — skipping")
            continue

        print(f"  Features after 70% presence rule: {len(feat_cols)}")

        # ── LOO CV: compute weighted feature importance per fold ──────
        fold_importances = []

        for tr_idx, te_idx in LeaveOneOut().split(feature_df):
            data_tr = feature_df.iloc[tr_idx]
            y_tr    = catdi_scores.iloc[tr_idx]

            train_mean = data_tr.mean()
            data_tr    = data_tr.fillna(train_mean)

            scaler = StandardScaler()
            X_tr   = scaler.fit_transform(data_tr)

            pca    = PCA(n_components=VARIANCE_THRESHOLD)
            pc_tr  = pca.fit_transform(X_tr)

            model  = LinearRegression()
            model.fit(pc_tr, y_tr)

            loadings   = pca.components_.T
            coef       = np.abs(model.coef_)
            importance = np.sum(np.abs(loadings) * coef[np.newaxis, :], axis=1)
            fold_importances.append(importance)

        mean_importance = np.mean(fold_importances, axis=0)

        # ── aggregate by method ────────────────────────────────────────
        method_importance = {m: 0.0 for m in METHOD_MAP}
        method_counts     = {m: 0   for m in METHOD_MAP}

        for i, feat in enumerate(feat_cols):
            for method, check in METHOD_MAP.items():
                if check(feat):
                    method_importance[method] += mean_importance[i]
                    method_counts[method]     += 1
                    break

        method_importance_norm = {
            m: method_importance[m] / method_counts[m] if method_counts[m] > 0 else 0.0
            for m in METHOD_MAP
        }

        print(f"  Feature counts per method: { {m: method_counts[m] for m in METHOD_MAP} }")
        for m, v in method_importance_norm.items():
            print(f"    {m}: {v:.4f}")

        # ── aggregate by contact ───────────────────────────────────────
        contact_importance  = {}
        contact_feat_counts = {}
        contact_to_region   = {}

        for i, feat in enumerate(feat_cols):
            contact, region = get_contact_and_region(feat)
            contact_importance[contact]  = contact_importance.get(contact, 0.0) + mean_importance[i]
            contact_feat_counts[contact] = contact_feat_counts.get(contact, 0) + 1
            contact_to_region[contact]   = region

        contact_importance_norm = {
            c: contact_importance[c] / contact_feat_counts[c]
            for c in contact_importance
        }

        print(f"\n  Top 5 CONTACTS by weighted importance:")
        for i, (contact, imp) in enumerate(
            sorted(contact_importance_norm.items(), key=lambda x: x[1], reverse=True)[:5]
        ):
            region = contact_to_region.get(contact, 'unknown')
            print(f"    {i+1:2d}. {contact:20s} {imp:.4f}  ({region})")

        # ── save ───────────────────────────────────────────────────────
        patient_out_dir = os.path.join(out_dir, patient)
        os.makedirs(patient_out_dir, exist_ok=True)

        contact_imp_df = pd.DataFrame([
            {"contact": c, "importance_norm": contact_importance_norm[c],
             "region": contact_to_region.get(c, "unknown"), "n_features": contact_feat_counts[c]}
            for c in contact_importance_norm
        ]).sort_values("importance_norm", ascending=False)
        contact_imp_df.to_csv(os.path.join(patient_out_dir, "contact_importance.csv"), index=False)

        # ── plot: Feature type contributions ───────────────────────────
        fig, ax = plt.subplots(figsize=(6, 5))
        methods = list(METHOD_MAP.keys())
        values  = [method_importance_norm[m] for m in methods]
        colors  = [METHOD_COLORS[m] for m in methods]
        ax.bar(methods, values, color=colors, edgecolor='black', linewidth=0.5)
        ax.set_ylabel("Mean Weighted Importance per Feature", fontsize=11)
        ax.set_xlabel("Analysis Method", fontsize=11)
        ax.grid(axis='y', alpha=0.2)
        plt.tight_layout()
        fig.savefig(os.path.join(patient_out_dir, "feature_type_contributions.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(patient_out_dir, "feature_type_contributions.svg"), format="svg", bbox_inches="tight")
        plt.close()

        # ── plot: Top 5 contacts ────────────────────────────────────────
        top_n = sorted(contact_importance_norm.items(), key=lambda x: x[1], reverse=True)[:5]
        labels, vals = zip(*top_n)
        colors_bar, labels_with_region = [], []
        for label in labels:
            region = contact_to_region.get(label, 'unknown')
            colors_bar.append(REGION_COLORS.get(region, REGION_COLORS['unknown']))
            labels_with_region.append(f"{label} ({region})")

        fig, ax = plt.subplots(figsize=(6, 3))
        ax.barh(range(len(labels_with_region)), vals, color=colors_bar, edgecolor="black", linewidth=0.3)
        ax.set_yticks(range(len(labels_with_region)))
        ax.set_yticklabels(labels_with_region, fontsize=9)
        ax.set_xlabel("Mean Weighted Importance per Feature\n(|loading| × |coef|, LOO avg, normalized)", fontsize=9)
        ax.invert_yaxis()
        ax.grid(axis='x', alpha=0.2)
        plt.tight_layout()
        fig.savefig(os.path.join(patient_out_dir, "top_contacts.png"), dpi=300, bbox_inches="tight")
        fig.savefig(os.path.join(patient_out_dir, "top_contacts.svg"), format="svg", bbox_inches="tight")
        plt.close()

        # ── extract prefix for each feature ───────────────────────────
        feat_prefix = [get_contact_and_region(f)[1] for f in feat_cols]

        # ── collect pooled data (excluding DBSTRD011) ─────────────────
        if patient != 'DBSTRD011':
            for i, feat in enumerate(feat_cols):
                subcat_name = 'other'
                for label, fn in SUBCATEGORIES:
                    if fn(feat):
                        subcat_name = label
                        break
                pooled_records.append({
                    'importance': mean_importance[i],
                    'prefix':     feat_prefix[i],
                    'subcat':     subcat_name,
                })

        rng = np.random.default_rng(42)

        # ── group importance by DBS grouping ──────────────────────────
        group_vals = {g: [mean_importance[i] for i in range(len(feat_cols))
                          if DBS_GROUP_FNS[g](feat_prefix[i])]
                      for g in DBS_GROUP_ORDER}
        present = [g for g in DBS_GROUP_ORDER if len(group_vals[g]) > 0]

        # ── plot: region_importance_violin ────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 4))
        positions = list(range(1, len(present) + 1))
        draw_violin(ax, [group_vals[g] for g in present], positions,
                    [DBS_GROUP_COLORS.get(g, '#CCCCCC') for g in present],
                    rng, dot_s=10, fontsize=9)
        ax.set_xticks(positions)
        ax.set_xticklabels([g.upper() for g in present], fontsize=9)
        ax.set_ylabel("|loading|×|coef|", fontsize=9)
        plt.tight_layout()
        fig.savefig(os.path.join(patient_out_dir, "region_importance_violin.png"), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(patient_out_dir, "region_importance_violin.svg"), format="svg", bbox_inches="tight")
        plt.close()

        # ── plot: subcategory_region_grid (3×7) ───────────────────────
        nrow, ncol = 3, 7
        fig, axes = plt.subplots(nrow, ncol, figsize=(28, 12), sharey=True)
        axes = axes.flatten()

        for ax_idx, (subcat_label, subcat_fn) in enumerate(SUBCATEGORIES):
            ax = axes[ax_idx]

            grp_vals_sub = {g: [mean_importance[i] for i, feat in enumerate(feat_cols)
                                if subcat_fn(feat) and DBS_GROUP_FNS[g](feat_prefix[i])]
                            for g in DBS_GROUP_ORDER}
            present_sub = [g for g in DBS_GROUP_ORDER if len(grp_vals_sub[g]) > 0]
            if not present_sub:
                ax.set_visible(False)
                continue

            pos_sub = list(range(1, len(present_sub) + 1))
            draw_violin(ax, [grp_vals_sub[g] for g in present_sub], pos_sub,
                        [DBS_GROUP_COLORS.get(g, '#CCCCCC') for g in present_sub], rng)
            ax.set_xticks(pos_sub)
            ax.set_xticklabels([g.upper() for g in present_sub], fontsize=7)
            ax.set_title(subcat_label, fontsize=9)

        plt.tight_layout()
        fig.savefig(os.path.join(patient_out_dir, "subcategory_region_grid.png"), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(patient_out_dir, "subcategory_region_grid.svg"), format="svg", bbox_inches="tight")
        plt.close()

        # ── plot: region_feature_stacked ──────────────────────────────
        n_grps = len(present)
        fig, axes_st = plt.subplots(n_grps, 1, figsize=(24, 5 * n_grps), sharey=True)
        if n_grps == 1:
            axes_st = [axes_st]

        for ri, grp in enumerate(present):
            ax = axes_st[ri]

            subcat_items = []
            for subcat_label, subcat_fn in SUBCATEGORIES:
                vals = [mean_importance[i] for i, feat in enumerate(feat_cols)
                        if subcat_fn(feat) and DBS_GROUP_FNS[grp](feat_prefix[i])]
                if vals:
                    subcat_items.append((subcat_label, vals, np.mean(vals)))

            if not subcat_items:
                ax.set_visible(False)
                continue

            subcat_items.sort(key=lambda x: x[2], reverse=True)
            sorted_labels = [s[0] for s in subcat_items]
            sorted_data   = [s[1] for s in subcat_items]
            positions_st  = list(range(1, len(subcat_items) + 1))

            grp_color = DBS_GROUP_COLORS.get(grp, '#CCCCCC')
            draw_violin(ax, sorted_data, positions_st,
                        [grp_color] * len(sorted_data), rng, fontsize=8)
            ax.set_xticks(positions_st)
            ax.set_xticklabels(sorted_labels, rotation=45, ha='right', fontsize=7)
            ax.set_xlim(0.5, len(subcat_items) + 0.5)
            ax.set_title(grp.upper(), fontsize=11, color=grp_color)

        plt.tight_layout()
        fig.savefig(os.path.join(patient_out_dir, "region_feature_stacked.png"), dpi=200, bbox_inches="tight")
        fig.savefig(os.path.join(patient_out_dir, "region_feature_stacked.svg"), format="svg", bbox_inches="tight")
        plt.close()

        print(f"\n  Saved to: {patient_out_dir}/")

    # ── pooled plots (all patients except DBSTRD011) ──────────────────
    print(f"\n{'='*70}")
    print(f"Generating pooled plots for probe type: {probe_type}...")
    print(f"{'='*70}")

    pooled_out_dir = os.path.join(out_dir, "pooled")
    os.makedirs(pooled_out_dir, exist_ok=True)

    pooled_df = pd.DataFrame(pooled_records)
    rng = np.random.default_rng(42)

    # sort groups by descending median importance across all pooled features
    _region_medians = {
        g: pooled_df[pooled_df['prefix'].apply(DBS_GROUP_FNS[g])]['importance'].median()
        for g in DBS_GROUP_ORDER
        if pooled_df['prefix'].apply(DBS_GROUP_FNS[g]).any()
    }
    pooled_group_order = sorted(_region_medians, key=_region_medians.get, reverse=True)

    # ── pooled plot 1: region_importance_violin ────────────────────────
    group_vals_p = {
        g: pooled_df[pooled_df['prefix'].apply(DBS_GROUP_FNS[g])]['importance'].tolist()
        for g in pooled_group_order
    }
    present_p   = [g for g in pooled_group_order if len(group_vals_p[g]) > 0]
    positions_p = list(range(1, len(present_p) + 1))
    n_p = len(present_p)

    fig, ax = plt.subplots(figsize=(12, 4))
    draw_violin(ax, [group_vals_p[g] for g in present_p], positions_p,
                [PINK_EDGE] * n_p,
                rng, dot_s=6, fontsize=10,
                dot_color_list=[PINK_EDGE] * n_p,
                body_color_list=[PINK_BODY] * n_p,
                median_color=MEDIAN_BLACK)
    ax.set_xticks(positions_p)
    ax.set_xticklabels([g.upper() for g in present_p], fontsize=10)
    ax.set_ylabel("|loading|×|coef|", fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(pooled_out_dir, "region_importance_violin.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(pooled_out_dir, "region_importance_violin.svg"), format="svg", bbox_inches="tight")
    plt.close()

    # ── pooled plot 2: subcategory_region_grid (3×7) ──────────────────
    nrow, ncol = 3, 7
    fig, axes = plt.subplots(nrow, ncol, figsize=(28, 12), sharey=True)
    axes = axes.flatten()

    for ax_idx, (subcat_label, _) in enumerate(SUBCATEGORIES):
        ax = axes[ax_idx]

        sub_df = pooled_df[pooled_df['subcat'] == subcat_label]
        grp_vals_sub = {
            g: sub_df[sub_df['prefix'].apply(DBS_GROUP_FNS[g])]['importance'].tolist()
            for g in DBS_GROUP_ORDER
        }
        present_sub = [g for g in DBS_GROUP_ORDER if len(grp_vals_sub[g]) > 0]
        if not present_sub:
            ax.set_visible(False)
            continue

        pos_sub = list(range(1, len(present_sub) + 1))
        draw_violin(ax, [grp_vals_sub[g] for g in present_sub], pos_sub,
                    [DBS_GROUP_COLORS.get(g, '#CCCCCC') for g in present_sub], rng)
        ax.set_xticks(pos_sub)
        ax.set_xticklabels([g.upper() for g in present_sub], fontsize=7)
        ax.set_title(subcat_label, fontsize=9)

    plt.tight_layout()
    fig.savefig(os.path.join(pooled_out_dir, "subcategory_region_grid.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(pooled_out_dir, "subcategory_region_grid.svg"), format="svg", bbox_inches="tight")
    plt.close()

    # ── pooled plot 3: region_feature_stacked ─────────────────────────
    n_grps_p = len(pooled_group_order)
    fig, axes_st = plt.subplots(n_grps_p, 1, figsize=(24, 5 * n_grps_p), sharey=True)
    if n_grps_p == 1:
        axes_st = [axes_st]

    for ri, grp in enumerate(pooled_group_order):
        ax = axes_st[ri]
        grp_mask = pooled_df['prefix'].apply(DBS_GROUP_FNS[grp])

        subcat_items = []
        for subcat_label, _ in SUBCATEGORIES:
            vals = pooled_df[grp_mask & (pooled_df['subcat'] == subcat_label)]['importance'].tolist()
            if vals:
                subcat_items.append((subcat_label, vals, np.mean(vals)))

        if not subcat_items:
            ax.set_visible(False)
            continue

        subcat_items.sort(key=lambda x: x[2], reverse=True)
        sorted_labels = [s[0] for s in subcat_items]
        sorted_data   = [s[1] for s in subcat_items]
        positions_st  = list(range(1, len(subcat_items) + 1))

        grp_color = DBS_GROUP_COLORS.get(grp, '#CCCCCC')
        draw_violin(ax, sorted_data, positions_st,
                    [grp_color] * len(sorted_data), rng, fontsize=8)
        ax.set_xticks(positions_st)
        ax.set_xticklabels(sorted_labels, rotation=45, ha='right', fontsize=7)
        ax.set_xlim(0.5, len(subcat_items) + 0.5)
        ax.set_title(grp.upper(), fontsize=11, color=grp_color)

    plt.tight_layout()
    fig.savefig(os.path.join(pooled_out_dir, "region_feature_stacked.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(pooled_out_dir, "region_feature_stacked.svg"), format="svg", bbox_inches="tight")
    plt.close()

    # ── pooled plot 4: top-10 subcategories across ALL groups combined ────────
    subcat_items_all = []
    for subcat_label, _ in SUBCATEGORIES:
        vals = pooled_df[pooled_df['subcat'] == subcat_label]['importance'].tolist()
        if vals:
            subcat_items_all.append((subcat_label, vals, np.median(vals)))

    subcat_items_all.sort(key=lambda x: x[2], reverse=True)
    top10 = subcat_items_all[:10]

    top10_labels = [s[0] for s in top10]
    top10_data   = [s[1] for s in top10]
    n_top10 = len(top10_labels)

    fig, ax = plt.subplots(figsize=(11, 6))
    positions = list(range(1, len(top10) + 1))
    draw_violin(ax, top10_data, positions, [PINK_EDGE] * n_top10, rng, dot_s=6, fontsize=9,
                dot_color_list=[PINK_EDGE] * n_top10,
                body_color_list=[PINK_BODY] * n_top10,
                median_color=MEDIAN_BLACK)
    ax.set_xticks(positions)
    ax.set_xticklabels(top10_labels, rotation=45, ha='right', fontsize=8)
    ax.set_xlim(0.3, len(top10) + 0.7)
    ax.set_ylabel("|loading|×|coef|", fontsize=10)
    plt.tight_layout()
    fig.savefig(os.path.join(pooled_out_dir, "top10_all_groups.png"), dpi=200, bbox_inches="tight")
    fig.savefig(os.path.join(pooled_out_dir, "top10_all_groups.svg"), format="svg", bbox_inches="tight")
    plt.close()

    print(f"\nPooled plots saved to: {pooled_out_dir}/")

print(f"\n{'='*70}")
print(f"Done! Results saved to: paper_figures/pca_loadings_dbs/")
print(f"{'='*70}")
