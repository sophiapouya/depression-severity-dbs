import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import itertools
import re
import scipy.io
import glob
import warnings
from datetime import datetime, date

# ====================================================================
# --- 1. CONFIGURATION ---
# ====================================================================

project_name = 'CATDI'
temp_sbj_name = 'DBSTRD002'
path = '/Users/sophiapouya/workspace/bcm' # Path to the 'bcm' directory

# Dynamic Path Assembly
project_root = os.path.join(path, project_name)
data_root = os.path.join(project_root, 'neuralData')
result_root = os.path.join(project_root, 'Results')

POWER_DATA_ROOT = os.path.join(data_root, 'dbsData', temp_sbj_name, 'powerData')

log_tran = True
bandlist = np.array(['delta', 'theta', 'alpha', 'beta', 'gamma', 'highgamma'])
warnings.filterwarnings('ignore', category=RuntimeWarning)

# ====================================================================
# --- 2. HELPER FUNCTIONS (Unchanged) ---
# ====================================================================

def find_original_folder(path, project_name, sbj_name, block_name):
    """Finds the originalData folder corresponding to a ScoreList block name."""
    subject_dir = os.path.join(path, project_name, 'neuralData', 'originalData', sbj_name)
    matches = []

    # 1. NEW LOGIC: Direct Match for DBSTRD011 / DBSTRD014 (ScoreList name IS the folder name)
    if sbj_name in ['DBSTRD011', 'DBSTRD014']:
        # Create a flexible pattern to handle potential inconsistencies (e.g., 'TRD014')
        name_pattern = block_name.replace('subj-TRD014', 'subj-*').replace('subj-DBSTRD014', 'subj-*')
        pattern = os.path.join(subject_dir, name_pattern)
        matches.extend(glob.glob(pattern))
        if matches:
            return sorted(matches)[0]
        
    # 2. Run/Block Style (DBSTRD001 / DBSTRD002)
    m_run = re.search(r'run-(?:Day)?(\d+)', block_name)
    m_blk = re.search(r'blk-(\d+)', block_name)
    if m_run and m_blk:
        rnum = int(m_run.group(1)); bnum = int(m_blk.group(1))
        # Account for common zero-padding and 'Day' variations
        run_tokens = [f"run-Day{rnum}", f"run-{rnum:02d}", f"run-{rnum}"]
        blk_tokens = [f"blk-{bnum:02d}", f"blk-{bnum}"]
        for cat in ["CATDI", "CAT-DI"]:
            for r in run_tokens:
                for b in blk_tokens:
                    pattern = os.path.join(subject_dir, f"*{cat}*{r}*{b}*")
                    matches.extend(glob.glob(pattern))
        if matches:
            return sorted(matches)[0]
    
    if ("date-" in block_name and "time-" in block_name):
        date_part = re.search(r'(date-\d{2}-\d{2}-\d{4})', block_name)
        time_part = re.search(r'(time-\d{2}-\d{2}-\d{2})', block_name)
        if date_part and time_part:
            for cat in ["CATDI", "CAT-DI"]:
                pattern = os.path.join(subject_dir, f"*{cat}*{date_part.group(1)}*{time_part.group(1)}*")
                matches.extend(glob.glob(pattern))
    if not matches:
        raise FileNotFoundError(f"No originalData folder for '{block_name}' found.")
    return sorted(matches)[0]

def find_power_file_by_folder(data_root, sbj_name, folder_basename):
    """Maps originalData folder name -> powerData file, handling all token formats."""
    base = os.path.join(data_root, 'dbsData', sbj_name, 'powerData')
    files = glob.glob(os.path.join(base, "*.mat"))
    if not files:
        raise FileNotFoundError(f"No .mat files under {base}")

    m = re.search(r'task-(.+)$', folder_basename)
    task_token = m.group(1) if m else folder_basename
    
    # 1. Specialized Token Creation for DBSTRD011 / DBSTRD014
    if sbj_name in ['DBSTRD011', 'DBSTRD014']:
        date_m = re.search(r'(date-\d+)', folder_basename)
        time_m = re.search(r'(time-\d+)', folder_basename)
        if date_m and time_m:
            task_token_fixed = f"CATDI_{date_m.group(0)}_{time_m.group(0)}"
        else:
            task_token_fixed = task_token
    else:
        task_token_fixed = task_token

    # A. Exact Match Search (Best attempt)
    for f in files:
        fname = os.path.basename(f)
        if f"_task-{task_token_fixed}.mat" in fname:
            return f
    
    # B. Fallback: Run/Block Style Match (DBSTRD001 / DBSTRD002)
    run_m = re.search(r'run-(?:Day)?(\d+)', task_token)
    blk_m = re.search(r'blk-(\d+)', task_token)
    if run_m and blk_m:
        rnum = int(run_m.group(1)); bnum = int(blk_m.group(1))
        run_alts = [f"run-Day{rnum}", f"run-{rnum:02d}", f"run-{rnum}"]
        blk_alts = [f"blk-{bnum:02d}", f"blk-{bnum}"]
        for f in files:
            fname = os.path.basename(f)
            if "_task-CATDI" in fname and any(r in fname for r in run_alts) and any(b in fname for b in blk_alts):
                return f

    raise FileNotFoundError(f"No powerData file for task '{task_token}' in {base}")

# ====================================================================
# --- 2.5. SCORE LIST LOADING & FILTERING (CRITICAL MISSING BLOCK) ---
# ====================================================================

print(f"Loading behavioral data for subject: **{temp_sbj_name}**")

# Load and filter ScoreList (y)
Scoredata = os.path.join(project_root, 'CATDI_scores.xlsx')
ScoreList = pd.read_excel(Scoredata, sheet_name=temp_sbj_name)
# CRITICAL: Filter out incomplete results
ScoreList = ScoreList[ScoreList.Result != 'incomplete'] 


# --- Subject-Specific Score Parsing and Filtering ---
if temp_sbj_name == 'DBSTRD001':
    # CRITICAL FIX for '93.3 (S)': Extract numeric part before any parentheses/text.
    ScoreList['Score'] = ScoreList['Result'].astype(str).str.extract(r'([\d\.]+)').astype(float)
    # Exclusions
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_run-08_blk-04']
    
elif temp_sbj_name == 'DBSTRD002':
    # Regular float conversion is fine for 002.
    ScoreList['Score'] = ScoreList['Result'].astype(float) 
    # Exclusions
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_run-05_blk-04']
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_run-07_blk-05']
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_run-03_blk-02']
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_run-03_blk-03']
    
elif temp_sbj_name == 'DBSTRD006':
    # Regular float conversion is fine for 006.
    ScoreList['Score'] = ScoreList['Result'].astype(float) 
    # Exclusions
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_date-02-13-2022_time-15-50-56']
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_date-02-09-2022_time-12-45-13']
    ScoreList = ScoreList[ScoreList['Name'] != 'CATDI_date-02-09-2022_time-14-01-20']

else:
    # Default for all other subjects (DBSTRD011, DBSTRD014, etc.)
    ScoreList['Score'] = ScoreList['Result'].astype(float) 

# These lines define the variables 'y' and 'bn' which the loop needs!
y = ScoreList['Score'].values
bn = ScoreList['Name'].values

print(f"Initial score list size (y, bn): {len(y)} blocks.")


# --- 3. DATA LOADING & AGGREGATION (Continued from Score Filtering) ---
# Assuming ScoreList, y, and bn have been created and filtered for exclusions (e.g., 34 blocks)

# --- Synchronized Power Data Aggregation ---
power_data_list = []
y_valid = []
bn_valid = []
processed_blocks_count = 0

for ii in range(0, len(bn)):
    block_name = bn[ii]
    try:
        # 1. Map Block Name to Folder Name
        orig_dir = find_original_folder(path, project_name, temp_sbj_name, block_name)
        folder_basename = os.path.basename(orig_dir)
        
        # 2. Map Folder Name to Power File Path
        powerpath = find_power_file_by_folder(data_root, temp_sbj_name, folder_basename) 
        
        # 3. Load Power Data
        power = scipy.io.loadmat(powerpath)
        allband_task = power['power'] # (channels x bands)
        
        # CRITICAL: Transpose and expand for consistent stacking (bands x channels x 1)
        allband_task = allband_task.T 
        allband_task = np.expand_dims(allband_task, axis=2) 
        
        # 4. Success: Append both power data and score to the *valid* lists
        power_data_list.append(allband_task)
        y_valid.append(y[ii])
        bn_valid.append(block_name)
        processed_blocks_count += 1
            
    except FileNotFoundError as e:
        # This means the score exists, but the power file is missing. We skip it.
        print(f"  - WARNING: Skipping block {block_name} (Index {ii}) - Power file missing or failed mapping. Error: {e}")
    except Exception as e:
        # Catch any other unexpected loading errors
        print(f"  - WARNING: Skipping block {block_name} (Index {ii}) due to unexpected error: {e}")


# 5. Final Aggregation and Consistency Check
if not power_data_list:
    print("FATAL ERROR: No power data files were loaded successfully after filtering.")
    sys.exit(1)

allblk_task = np.concatenate(power_data_list, axis=2)

# Overwrite the original y and bn variables with the validated subset
y = np.array(y_valid)
bn = np.array(bn_valid)

print(f"Aggregated power data shape: {allblk_task.shape} (Bands x Channels x Blocks)")
print(f"Validated scores (y) count: {len(y)}")

if allblk_task.shape[2] != len(y):
    print("FATAL ERROR: Array shape inconsistency after aggregation. Debug is needed.")
    sys.exit(1)
    
print(f"Aggregated power data shape: {allblk_task.shape} (Bands x Channels x Blocks)")

# ====================================================================
# --- 3.5. BEHAVIORAL SCORE DETRENDING (Time Regression) ---
# ====================================================================

# 1. Generate Datetime objects from the validated ScoreList
datetimes = []

# CRITICAL: We need to use the subset of ScoreList that corresponds to the VALIDATED 'bn'
ScoreList_valid = ScoreList[ScoreList['Name'].isin(bn)].copy()

if temp_sbj_name == 'DBSTRD001':
    for idx in range(len(ScoreList_valid)):
        date_str = str(ScoreList_valid.iloc[idx]['EMUdate'])
        time_str = str(ScoreList_valid.iloc[idx]['Timestamp'])

        if '/' in date_str:
            parts = date_str.split('/', 1)
            date_str = parts[1].strip()

        if date_str.startswith("Thur"):
            date_str = date_str.replace("Thur", "Thu", 1)

        combined = f"{date_str} {time_str}"

        # Try AM/PM first, then fallback to 24-hour with seconds
        try:
            dt_obj = datetime.strptime(combined, "%a %b %d, %Y %I:%M %p")
        except ValueError:
            dt_obj = datetime.strptime(combined, "%a %b %d, %Y %H:%M:%S")
        datetimes.append(dt_obj)
        
# For other subjects, the logic for datetimes would need to be re-added here...
# For DBSTRD001, this block is sufficient based on your original code's logic.

elif temp_sbj_name == 'DBSTRD002': # <--- NEW/CORRECT LOCATION FOR DBSTRD002 LOGIC
    for idx in range(len(ScoreList_valid)):
        raw_date = ScoreList_valid.iloc[idx]['EMUdate']
        raw_time = str(ScoreList_valid.iloc[idx]['Timestamp']).strip()

        # --- date: robustly get a date() from Excel datetime or string ---
        if isinstance(raw_date, (pd.Timestamp, datetime, date)):
            d = pd.to_datetime(raw_date).date()
        else:
            d = pd.to_datetime(str(raw_date), errors='coerce')
            if pd.isna(d):
                d = datetime.strptime(str(raw_date).strip(), "%Y-%m-%d %H:%M:%S")
            d = d.date()

        # --- time: accept various formats ---
        t = None
        raw_time_nospace = raw_time.replace(' ', '').upper()
        for fmt, candidate in [
            ("%I:%M%p", raw_time_nospace),     # "3:55PM"
            ("%I:%M %p", raw_time.upper()),    # "3:55 PM"
            ("%H:%M:%S", raw_time),            # "15:55:00"
            ("%H:%M", raw_time),               # "15:55"
        ]:
            try:
                t = datetime.strptime(candidate, fmt).time()
                break
            except ValueError:
                continue
        if t is None:
            # pandas fallback for any weird cases
            t = pd.to_datetime(raw_time, errors='raise').time()

        dt_obj = datetime.combine(d, t)
        datetimes.append(dt_obj)

elif temp_sbj_name in ['DBSTRD011', 'DBSTRD014']:
    for name in ScoreList_valid.Name:
        date_match = re.search(r'date-(\d{8})', name)
        time_match = re.search(r'time-(\d{6})', name)
        
        if date_match and time_match:
            dt_str = f"{date_match.group(1)} {time_match.group(1)}"
            dt_obj = datetime.strptime(dt_str, "%Y%m%d %H%M%S")
            datetimes.append(dt_obj)
        else:
            raise ValueError(f"Could not parse date/time from block name: {name}")

else: # Fallback logic (for DBSTRD006)
    for name in ScoreList_valid.Name:
        date_part = name[-24:-14]
        time_part = name[-8:]
        d = datetime.strptime(str(date_part), "%m-%d-%Y").date()
        t = datetime.strptime(str(time_part), "%H-%M-%S").time()
        datetimes.append(datetime.combine(d, t))


# 2. Calculate time in seconds elapsed since the first block
times = [(mydatetimes - datetimes[0]).total_seconds() for mydatetimes in datetimes]

# 3. Fit linear regression (Score ~ Time)
m, b = np.polyfit(times, y, 1)  # m = slope, b = intercept

# 4. Calculate residuals
y_residual = y - (m * np.array(times) + b)
print(f"Calculated y_residual by detrending scores (y) against session time.")

# Apply Log Transformation to Power Data (x)
if log_tran:
    eps = np.finfo(np.float32).tiny
    x = 10 * np.log10(np.maximum(allblk_task, eps))
    print("Applied 10*log10 transformation to power data.")
else:
    x = allblk_task
# ====================================================================

# ====================================================================
# --- 4. CHANNEL LIST FIX: Using Data Structure for Labels ---
# ====================================================================

MAX_CHANNELS = allblk_task.shape[1] # This is 32

# Manually define the channel groups based on your debugging screenshots
GROUP_NAMES = ['DLVCVS', 'DLSCC', 'DRVCVS', 'DRSCC']
CHANNELS_PER_GROUP = 8

# Create the sequential list of labels (e.g., DLVCVS01, DLVCVS02, ..., DRSCC08)
labels = []
for group_name in GROUP_NAMES:
    for i in range(1, CHANNELS_PER_GROUP + 1):
        labels.append(f'{group_name}{i:02d}')

# Create the DataFrame used for plotting
df = pd.DataFrame({
    'Electrode': np.arange(1, MAX_CHANNELS + 1), # 1 to 32
    'Label': labels 
})

if len(df) != MAX_CHANNELS:
    print("FATAL ERROR: Programmatic DF creation failed.")
    sys.exit(1)

print(f"Programmatic channel list created with {len(df)} channels and descriptive labels.")


# Apply Log Transformation
if log_tran:
    eps = np.finfo(np.float32).tiny
    x = 10 * np.log10(np.maximum(allblk_task, eps))
    print("Applied 10*log10 transformation to power data.")
else:
    x = allblk_task


# ====================================================================
# --- 4.5. CRITICAL DEBUGGING CHECK ---
# ====================================================================
print("\n--- CRITICAL DATA INTEGRITY CHECK ---")
print(f"Score Data (y) Min: {np.min(y)}, Max: {np.max(y)}")

x_finite = x[np.isfinite(x)]
if x_finite.size > 0:
    print(f"Power Data (x) Shape: {x.shape}")
    print(f"Power Data (x) Min: {np.min(x_finite):.2f}, Max: {np.max(x_finite):.2f}, Mean: {np.mean(x_finite):.2f}")
    # Check for variance: if variance is zero, correlation will be zero.
    variance_check = np.var(x_finite)
    print(f"Power Data (x) Variance (Should be > 0): {variance_check:.2f}")
else:
    print("Power Data (x) contains no finite values (all NaN or Inf).")

if np.var(y) < eps:
    print("WARNING: Score data (y) has zero variance (all scores are the same). Correlation will be zero.")
    
print("-------------------------------------\n")
# ====================================================================


# ====================================================================
# --- 5. PLOTTING FUNCTION & EXECUTION (Unchanged) ---
# ====================================================================

def plot_r_allchannel(allblk_task_log, y, df, bandlist, ax, patient_label=None):
    """Plots the correlation (r-value) between log-power and scores."""
    colorbarlabel = 'Correlation Coefficient'

    chan_indices = df['Electrode'].values - 1  
    
    r_all = np.empty([len(bandlist), len(df)]) 
    p_all = np.empty([len(bandlist), len(df)])

    # Calculate Pearson Correlation
    for bandindex in range(len(bandlist)):
        for i, original_elec_index in enumerate(chan_indices):
            x_band_chan = allblk_task_log[bandindex, original_elec_index, :] 

            finite_mask = np.isfinite(x_band_chan) & np.isfinite(y)
            x_f = x_band_chan[finite_mask]
            y_f = y[finite_mask]

            if x_f.size < 3 or np.allclose(np.var(x_f), 0) or np.allclose(np.var(y_f), 0):
                rvalue, pvalue = 0.0, 1.0 
            else:
                # The core correlation function: should produce non-zero 'rvalue' if data varies
                rvalue, pvalue = stats.pearsonr(y_f, x_f)

            r_all[bandindex, i] = rvalue
            p_all[bandindex, i] = pvalue

    # FDR correction
    _, p_corr = fdrcorrection(np.reshape(p_all, [-1]))
    p_corr = p_corr.reshape(p_all.shape) # CORRECT: Reshape 1D array back to original shape (6, 32)

    # Labels
    y_axis_labels = ['$\delta$', '$\\theta$', '$\\alpha$', '$\\beta$', '$\gamma$', 'h$\gamma$']
    x_axis_labels = df['Label'].values 

    # Determine significance for annotation
    sign_matrix = np.select([p_corr >= 0.05, p_corr < 0.05], ['', '*'], default='')
    
    matplotlib_rc = {'font.size': 12, 'font.family': "Arial", 'axes.labelsize': 12, 'legend.fontsize': 12, 
                     'axes.titlesize': 14, 'xtick.labelsize': 8, 'ytick.labelsize': 12}
    sns.set(style="white", rc=matplotlib_rc)

    # Plot Heatmap
    ax = sns.heatmap(
        r_all, ax=ax, annot_kws={"size": 10, "color": 'black', "weight": 'bold'}, 
        annot=sign_matrix, fmt="", yticklabels=y_axis_labels,
        cmap=sns.diverging_palette(220, 20, as_cmap=True), center=0, vmin=-1, vmax=1,
        xticklabels=x_axis_labels, cbar_kws={'label': colorbarlabel, 'ticks':[-1, 0, 1]}, cbar=True
    )
    
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha='center')
    plt.yticks(rotation=0)
    ax.set_title(patient_label or "Correlation (All Channels)")
    ax.tick_params(axis='x', length=0)
    
    print(f"\nFDR corrected sig count (q<0.05): {np.sum(p_corr < 0.05)}")
    return ax

# --- Execution ---
fig, ax = plt.subplots(1, 1, figsize=(14, 7)) 

plot_r_allchannel(
    allblk_task_log=x, 
    y=y_residual,  # <-- CHANGE THIS
    df=df, 
    bandlist=bandlist, 
    ax=ax, 
    patient_label=f'Patient {temp_sbj_name[-3:]}: Power-Score Correlation (FDR-Corrected)'
)

plt.tight_layout()
plt.savefig(f"dbs_corr_plots/{temp_sbj_name}.svg")