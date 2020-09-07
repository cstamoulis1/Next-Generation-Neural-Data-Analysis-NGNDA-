function [common_mode_num, average_freqs] = most_common_mode_bands(frames_threshold)
    % Purpose: to compute most common number of mode bands across all
    % brains and runs, and the average frequencies of these mode bands needed for mode classification
    % and comparison between brains
  % Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
  % Last version: September 5, 2020
  % INPUT:     frames_threshold = Threshold of % frames censored for excessive motion above which a data run is excluded from analysis

    %
    % OUTPUTS:   common_mode_num = most common number of mode bands for a given run out of all
    %                           runs of all subjects
    %            average_freqs = median frequency for each mode band
    %            calculated as the cohort median of all run-specific median (over nodes) frequencies across

    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    connectivity_loc = '';  % enter path to directory where connectivity files are (ie '/path_to_connectivity_directory/connectivity_directory')
    cd(connectivity_loc); % cd into directory where connectivities are located
    S = dir('sub*'); % get list of all subjects with connectivities in directory
    
	freqs_nums = NaN(length(S),4); % array to get number of mode bands present for each run of each subject
    freqs_all = cell(length(S),4); % each cell gets average (median) frequencies for each mode band for particular run of particular subject
    parfor i = 1:length(S)
        filename = S(i).name;
        RS_connectivity = load(filename); % loads RS_connectivity file
        if isfield(RS_connectivity,'RS_connectivity') == 1
            RS_connectivity = RS_connectivity.RS_connectivity;
        end   
        
        for r = 1:4 % goes through each run
            if length(RS_connectivity) >= r
                if ~isempty(RS_connectivity{r})
                    if RS_connectivity{r}.percentCensored < frames_threshold % if run is not empty and meets percent censored threshold
                        if isfield(RS_connectivity{r},'modes') % if modes is a field
                    	    freqs = RS_connectivity{r}.modes.freqs; % extracts mode frequencies for run
                            
                            freqs_i = nanmedian(freqs,1); % gets median frequency for each mode band 
                            freqs_all(i,r) = {freqs_i}; % stores these median values in cell
                    	    freqs_nums(i,r) = sum(~isnan(freqs_i)); % stores number of mode bands for this run
                            
                        end
                    end
                
                end
            end
        end

    end
    
	common_mode_num = mode(freqs_nums,'all'); % most common mode number (number of mode bands) is the mode of this array
    freqs_all = freqs_all(:); % linearizes cell array
    freqs_undone = []; % unpacks each cell and stores in array
    for i = 1:length(freqs_all)
        freqs_undone(i,1:length(freqs_all{i})) = freqs_all{i};
    end
    
    % freqs undone is an array where each row represents a run, and the
    % columns represent the mode bands, all values are median frequencies
    freqs_undone(freqs_undone == 0) = NaN; % any 0s are set to nan
    average_freqs = nanmedian(freqs_undone,1); % gets median of median frequencies for each mode band
end
