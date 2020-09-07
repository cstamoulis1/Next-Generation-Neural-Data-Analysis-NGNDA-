# Next-Generation-Neural-Data-Analysis (NGNDA)-fMRI

1. construct_raw_connectome_main.m %this is the code submitted for each brain to preprocess th
e fMRI data
Calls functions construct_raw_connectome_helper.m, parcellation_generator.m and decomposition 
codes ceemdan and emd.m  
function construct_raw_connectome_helper containing all preprocessing steps (registration, seg
mentation, motion-correction, parcellation, filtering and signal decomposition

function parcellation_generator.m generates a parcellation atlas used to reduce the spatial di
mensionality of the fMRI data to set of brain regions (parcels)
functions ceemdan (implements the complete Ensemble Empirical Mode Decomposition to address th
e issue of mode mixing) and eemd and emd (implement the ensemble EMD and original EMD, respect
ively) are used for the fMRI decomposition into dominant components but were written by others
 (see citations in the functions themselves)

2. calculate_mode_characteristics.m %this code calculates the characteristic frequency and amp
litude of each narrowband component (mode) estimated from each fMRI (node)
signal
Calls function save_modeInfo.m to save these estimates for each brain

3. modes_classification_and_connectome_reduction.m %this code denoises fMRI data by removing s
ignal components (modes) that are outside of the
expected frequency range and those with amplitudes outside a statistically estimated range.
Based on their characteristic (dominant) frequencies, modes are classified in distinct mode/fr
equency bands that are consistent across brains, to enable comparison between brains
The output is both the denoised broadband data, obtained as a linear superposition of the rema
ining narrowband components of each fMRI time series (synthesis) and the individual narrowband
 components with their characteristics for further analysis.

Calls functions calculate_mode_thresh.m to estimate the relevant frequency and amplitude thres
holds for denoising and save_reduced_connectomes.m to save the denoised data

4. process_connectivity_matrix.m %this code computes connectivity matrices for denoised (mode-
reduced and resynthesized) raw data and individual mode band data,
using multiple approaches (signal cross-correlation, cross-coherence and Mutual Information (MI)

5. calculate_network_properties.m % this code calculates multiple network properties (includin
g those in the Brain Connectivity Toolbox as well as an additional measure of network robustne
ss (natural connectivity) using a user-defined set of thresholds to threshold connectivity mat
rics and create corresponding adjacency matrices (binary and weighted) from this network prope
rties are estimated.
Calls function calculate_network_properties_helper.m which contains the implementation of or c
all to these measures.

Connectivity thresholds both for the broadband and narrowband (mode) connectivity matrices are
 estimated using auxiliary codes
compute_raw_connectivity_thresholds.m 
compute_mode_connectivity_thresholds.m (which calls another auxiliary function most_common_mod
e_bands.m)





