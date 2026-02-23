# Predicting depression severity from DBS probe data 

## Data preprocessing 
dbs_preprocessing.py
   - loads raw blackrock data (.ns3 or .ns5 files)
   - band pass filters raw data from 0.3hz to 500hz (mne filter function)
   - notch filters raw data (find_peaks from scipy.signal)
   - visual data inspection to remove bad artifacts or bad channels
   - saves json of visual inspection (bad channels and/or bad artifacts)
   - removes the channels marked as bad
   - performs the following referencing methods
      - bipolar
      - alternating bipolar
      - esr
      - car
   - saves the output of the referencing methods as session level fif which contains (if any) the artifact times in the annotations 
   

## Feature Extraction
dbs_power.py
  - calculates the power within each band using a hilbert transform
  - saves a csv per patient with averaged log power values in each frequency band

    frequency bands defined as:
       
       - delta: 1-4 hz
       - theta: 4-8 hz
       - alpha: 8-12 hz
       - beta: 12-30 hz
       - low gamma: 35-50 hz
       - high gamma: 70-150 hz

feature_exactraction.py
  - calculates a range of features using the pyneuromodulation package
  - calculates features on chunks of data and averages over sessions
  - returns a master CSV of all session, patients, and features

## Regression

## Correlation
catdi_correlation.py 
  - performs pearson correlation b/w band averaged log power and CATDI session score
  - plots correlation per patient 


## Source
src/
  - brpy/
    - brMiscFxns.py & brpylib.py -> Blackrock python files (edits made to brpy library, required to function correctly) 
  - preprocessing_functions.py -> functions used in dbs_preprocessing.py stored here 
  - inspect_npy.py -> inspecting a single npy file
  - inspect_fif.py -> inspecting a single fif file

## Software requirements
- see environment.yml for a full list of versions and packages used

