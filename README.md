# Transient and steady-state acoustic unit discovery


## Steps for AUD:
1. Extract features inside featExtract directory
2. Segment the audio using the segmentation binary from syllable_segmentation_gd repo
3. Cluster the syllable-like segments inside clusterDir
4. Iterative self-training using the scripts in main_expt (requires kaldi)

## Requirements
The repo uses Kaldi for feature extraction and acoustic modelling;
