function calculate_mode_characteristics(brains2process, cortical_atlas, frames_threshold)
% Each broadband fMRI time series is decomposed into its dominant narrowband components with 
% distinct characteristic frequencies and contributions (amplitude) to the original signal
% purpose: to calculate characteristic frequency and amplitude for all
% modes of all nodes; some modes will later be excluded as part of denoising using this frequency and amplitude information
% Also calls crossing.m to compute mode frequency based on the zero crossings of the data and save_modeInfo.m to save mode information to the 
% correct folder
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: September 5, 2020
% INPUT: brains2process = cell array containing a subject id string in each
%                         individual cell according to the naming convention 'sub-NDARINVXXXXXXXX'
%                         These consist of all brains that need mode info
%                         to be calculated for them
%        cortical_atlas = parcellation atlas/template used to parcellate cortex
%        frames_threshold = threshold of % frames censored above which a data run is excluded
% OUTPUT: for a subject, saves file of mode characteristic properties of
%         characteristic frequency and characteristic amplitude for each
%         mode of each node (region) in the parcellated raw connectome data

% paths to where connectome data are located and where to save info on
% modes
    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    cd(main_directory);
    modeinfo_loc = ''; % path to directory where you want to save mode info (ie '/pathtomodeinfodirectory/modeinfodirectory')
    connectome_loc = ''; % path to directory where connectomes are located (ie '/pathtoconnectomedirectory/connectomedirectory')
    
    if cortical_atlas == 1
        identifier = '_Gordon';
    elseif cortical_atlas == 2
        identifier = '_Schaefer_1000';
    elseif cortical_atlas == 3
        identifier = '_Schaefer_300';
    elseif cortical_atlas == 4
        identifier = '_Schaefer_800';
    elseif cortical_atlas == 5
        identifier = '_Schaefer_100';
    end
    
% parfor loop through each brain that needs mode info calculated; True parallelization occurs at the brain level 
%In the ABCD and other large-scale fMRI studies, the highest dimension is that of the cohort once the fMRI data are parcellated
%In the case of voxel-based analysis (not performed in this version of the pipeline) the parallelization can be changed
    parfor s = 1:length(brains2process)
        cd(connectome_loc); % cds into location where mode info is located
        S = dir(strcat(brains2process{s},'*')); 
        RS_modeInfo = cell(1,4); 
        filename = S(1).name;
        RS_net = load(filename, 'RS_net'); % loads subject specific connectome file
        cd(main_directory); % cd to main folder where code is located
        if isfield(RS_net,'RS_net') == 1 
            RS_net = RS_net.RS_net;
        end
        norm_exists = 0; % whether raw data is normalized
        for j = 1:length(RS_net)
            if ~isempty(RS_net{j})
                RS_modeInfo{j}.percentCensored = RS_net{j}.percentCensored;
                if RS_net{j}.percentCensored < frames_threshold % if percent censored for run is below frames threshold
                    if isfield(RS_net{j},'rawData_normalized')
                        norm_exists = 1; 
                        if ~isempty(RS_net{j}.rawData_normalized)
                            TR = RS_net{j}.TR;
                            fs = 1/TR;
                            numNodes = size(RS_net{j}.rawData_normalized,1); % number of regions
                            npts = size(RS_net{j}.rawData_normalized,2); % number of time points
                            decomposedNodes = RS_net{j}.decomposedNodes; % gets modal data
                            maxModes = max(cellfun('size',decomposedNodes,1)); % max number of modes out of any region
                            modeAmp = NaN(numNodes,maxModes);
                            modeFreq = NaN(numNodes,maxModes);
                            
                            for node = 1:numNodes
                                modes = decomposedNodes{node};
                                for m = 1:size(modes,1)
                                    modeAmp(node,m) = nanmedian(abs(modes(m,:))); % median amplitude of mode (characteristic amplitude)
                                    modeFreq(node,m) = (fs/2)*length(crossing(modes(m,:)))/npts; % rough estimate for characteristic frequency of mode
                                end
                            end
                            RS_modeInfo{j}.amp = modeAmp; % stores mode info for each run
                            RS_modeInfo{j}.freq = modeFreq;
                        end
                    end
                end
            end
        end
        if norm_exists == 1
            k = regexp(filename, '_', 'once')-1;
            name = filename(1:k);
            str = strcat(name,'_RS_modeInfo',identifier,'.mat'); % specify file name
            save_modeInfo(str,RS_modeInfo,modeinfo_loc); % uses helper function to save file within parfor loop
        end
    end
end
