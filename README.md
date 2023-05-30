# analysis for updates to Deutsch et al 2019
```
data collection : Yiqin Gao
analysis: Kyle Thieringer
```

## order of operations
1. exportslp.py
   - exports .slp files to .h5 files
2. createfeatures.py
   - calculates kinematic features and saves them to .h5 files
3. createMetatable.ipynb
   - uses exported googledrive csv file and creates a metadata sheet 
   to be used by the createAllData.m script
4. exptSync.ipynb
   - uses the vfaas (video frame at audio sample) to interpolate the audio sample 
   at each video frame which is needed by createAllData.m
5. checkFtrs.ipynb
   - visualizes basic features as a sanity check to assess quality of data
6. createAllData.m
   - creates a mat file containing kinematic features, audio data, and sync data for each fly
   - too large for github upload
7. Analyze_Behavior_Jans.m
   - uses output of createAllData to bin behavior and audio data for later analysis
8. rankCorrelation.ipynb
   - make plots using the output of Analyze_Behavior_Jans.m

---
**some notable distinctions in processing the data**
- I smoothed data during the createfeatures.py step using a kalman filter.
  - in the original analysis, behavior features were median filtered
- female speed was zscored for each fly then binned rather than using the raw speed