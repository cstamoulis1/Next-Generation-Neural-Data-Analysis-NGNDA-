# Next-Generation-Neural-Data-Analysis (NGNDA)
Tools for massive parallel analysis of brain data
The NGNDA repository includes Matlab codes for multiscale (multiband) connectivity analysis of human electrophysiological (EEG) and fMRI data. Version 1.0 of the pipeline includes codes for preprocessing large-scale fMRI datasets such as the Adolescent Brain Cognitive Development (ABCD) dataset and estimating resting-state connectomes and their network properties, offering options for different connectivity methods and comparisons. Version 1.0 also includes codes for connectivity analysis of very high-dimensional time series (human EEG) using high-resolution processing in time. 

The repository is organized into 2 folders, NGNDA-EEG and NGNDA-fMRI, each containing modality-specific codes and a list of the codes (README) with details on use. All codes include parallelizations in different dimensions (always the highest dimension), depending on the data modality and the implemented task. They have been designed for deployment in a High Performance Cluster (HPC) with shared resources. All fMRI related codes require modest resources (<20 CPUs). In contrast, although EEG codes can be deployed using modest resources as well (and each time-domain computation is relatively quick), given the very high number of time segments processed, extensive simulations using different allocation of resources have shown that 40-100 CPUs is optimal for these codes.  

Version 1.0 was released September 5, 2020
Use of the NGNDA codes is entirely unrestricted. However, please acknowledge the Computational Neuroscience Laboratory at Boston Childrens Hospital/Harvard Medical School in any research products using the codes. 
