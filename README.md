# Predicting Depression Severity From Intracranial Neural Data

Python pipeline for preprocessing intracranial neural recordings, extracting signal features, and building statistical and machine learning models to study depression severity.

This project works with two neural recording modalities:

- **DBS** (deep brain stimulation leads)
- **sEEG** (stereoelectroencephalography leads)

It was built as an end-to-end workflow for real biomedical time-series data: from raw signal ingestion, to quality control, to feature engineering, to predictive modeling.

Built with Python- key libraries include MNE, NumPy, pandas, scikit-learn, SciPy

---

## Overview

This repository converts raw intracranial recordings into structured, model-ready data.

The pipeline includes:

- raw Blackrock data loading
- signal preprocessing and rereferencing
- artifact and bad-channel handling
- feature extraction from DBS and sEEG recordings
- correlation analysis with depression severity scores
- regression modeling with LASSO and PCA-based workflows
- permutation testing for model validation

Although this project is grounded in a neuroscience application, the technical work is broadly relevant to AI/ML and data science:

- time-series preprocessing
- feature engineering
- statistical modeling
- reproducible Python workflows
- analysis of noisy real-world data

---

## Why this project is interesting

Most clinical symptom measures are subjective and sampled infrequently. This project explores whether depression severity can be linked to measurable patterns in brain activity using intracranial recordings.

From an engineering perspective, the main challenge is not just modeling. It is building a reliable pipeline that can take raw biomedical signals, clean them, transform them into usable features, and support downstream analysis in a reproducible way.

---

## What the pipeline does

### 1. Preprocessing

The preprocessing scripts load raw Blackrock neural recordings and convert them into analysis-ready files.

Implemented steps include:

- loading raw `.ns3` and `.ns5` files
- bandpass filtering from **0.3 Hz to 500 Hz**
- notch filtering to reduce line noise
- visual inspection for noisy channels and artifacts
- saving quality-control decisions to JSON
- removing bad channels
- rereferencing signals
- saving cleaned session-level `.fif` outputs

Supported rereferencing methods:

- bipolar
- alternating bipolar
- ESR
- CAR (common average reference)

Main scripts:

- `dbs_preprocessing.py`
- `seeg_preprocessing.py`

### 2. Feature extraction

After preprocessing, the pipeline extracts numerical features that summarize neural activity for downstream analysis and modeling.

Implemented workflows include:

- feature extraction for DBS data
- feature extraction for sEEG data
- chunk-level computation
- session-level aggregation
- export to structured CSV files

Main scripts:

- `feature_extraction_dbs.py`
- `feature_extraction_seeg.py`

### 3. Postprocessing and analysis

The repository includes analysis scripts for relating extracted neural features to depression severity scores.

Implemented workflows include:

- CAT-DI correlation analysis
- modality-specific feature correlation
- visualization and summary plotting
- downstream preparation of outputs for modeling

Main scripts:

- `dbs_postprocessing.py`
- `seeg_postprocessing.py`
- `catdi_correlation.py`
- `feature_correlation_dbs.py`
- `feature_correlation_seeg.py`
- `box_plots.py`

### 4. Modeling and validation

The repository also supports predictive modeling workflows.

Implemented workflows include:

- LASSO regression
- PCA + regression
- modeling from manually prepared power features
- permutation testing to assess model significance

Main scripts:

- `lasso.py`
- `lasso_manual_power.py`
- `lasso_manual_power_seeg.py`
- `pca_regression.py`
- `perm_lasso_dbs.py`
- `perm_lasso_seeg.py`

---

### 5. Results

The final models acheived statistically significant depression severity decoding in 5 of 6 patients (83% of cohort) with PCA reducing ~300 features to 7-20 prinicpal component explaining ~90% of variance. Permuation testing confirmed model significance above chance.

---

## Repository structure

```text
.
├── src/
├── standalone_scripts/
├── README.md
├── box_plots.py
├── catdi_correlation.py
├── config.py
├── config.json (must create this file for your paths)
├── dbs_postprocessing.py
├── dbs_preprocessing.py
├── environment.yml
├── feature_correlation_dbs.py
├── feature_correlation_seeg.py
├── feature_extraction_dbs.py
├── feature_extraction_seeg.py
├── lasso.py
├── lasso_manual_power.py
├── lasso_manual_power_seeg.py
├── pca_regression.py
├── perm_lasso_dbs.py
├── perm_lasso_seeg.py
├── seeg_postprocessing.py
└── seeg_preprocessing.py