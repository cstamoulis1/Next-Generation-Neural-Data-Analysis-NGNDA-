# Next-Generation-Neural-Data-Analysis-NGNDA-EEG
This part of the repository contains codes for high-dimensional analysis of human electrophysiology data
The main codes are:
1. sc2_analysis_allbands: this code performs signal decomposition and connectivity analysis 
for a particular file and segment of continuous EEG data

It decomposes broadband electrophysiological signals using a modified Empirical Mode Decomposition
approach, and unsupervised approach for narrowband signal decomposition/filtering. This analysis uses various modifications
of the original EMD approach.
References:
[1] N. E. Huang et al., "The empirical mode decomposition and the
Hilbert spectrum for non-linear and non stationary time series analysis",
Proc. Royal Soc. London A, Vol. 454, pp. 903-995, 1998
[2] Stamoulis, C, Betensky, RA, A novel signal processing approach for the detection of copy-number variations 
in the human genome, Bionformatics, 27(17):2338-2345, 2011, PMID: 21752800
[3] Torres, M. E., Colominas, M. A., Schlotthauer, G., and Flandrin, P. (2011, May):
A complete ensemble empirical mode decomposition with adaptive noise
(ICASSP), 2011 (pp. 4144-4147). DOI: 10.1109/ICASSP.2011.594726
[4] Stamoulis, C., Betensky, RA, Optimization of Signal Decomposition Matched Filtering for improved detection 
of copy-number variations in genomic data, IEEE Trans Comput Biol Bioinform, 13(3): 584-591, 2016, PMID: 27295643

Signals analysis using this code represent segments of continuous (typically multi-hour) recordings at multiple electodes 
(often >100 in invasive studies). Even in the case of high-resolution probes, the temporal dimension will always be the highest
requiring parallel processing.
The code uses both cross-correlation and mutual information as measures of connectivity, estimated for each narrowband mode,
and saves these multiband connectivity matrices, as well as mode parameters as part of a single structure.

Given that the highest dimension in this analysis is in time, true parallelization is done in that dimension.
So, sc2_analysis_allbands takes as arguments the file of interest and segment of interest and 
is called through a parfor loop by sc2_masterjob.

2. master_submit_sc2.m sets the cluster, queue and sc2_master parameters for processing and su
bmits one batch job (sc2_masterjob).

