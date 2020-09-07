function calculate_network_properties(subID, measure, cortical_atlas, raw_thresholds, threshold_labels, mode_thresholds, frames_threshold, common_mode_num, average_freqs)
    % Purpose: to compute network properties of a resting connectome
    % using one or multiple population level statistics to threshold the broadband connectivity and the mode-specific connectivity matrices
    % and obtain corresponding adjacency matrices
    % Also calls calculate_network_properties_helper.m  
    % Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
    % Last version: September 5, 2020
    % INPUTS: measure = 1 is peak cross-correlation
    %         measure = 2 is peak coherence
    %         measure = 3 is mutual information
    %
    %         subID = 'sub-NDARINVXXXXXXXX'
    %
    %         cortical_atlas = 1 is Gordon
    %         cortical_atlas = 2 is Schaefer 1000
    %         cortical_atlas = 3 is Schaefer 300
    %         cortical_atlas = 4 is Schaefer 800
    %         cortical_atlas = 5 is Schaefer 100
    %         
  %         raw_thresholds = vector of thresholds for raw connectivity; can be a single threshold
    %         These thresholds are necessary in order to estimate the weighted (or binary) adjacency matrix from which network properties 
    %         will be computed
    %
    %         threshold_labels = cell array containing labels for each
    %                            threshold used
    %         mode_thresholds = matrix of thresholds for modal connectivity
    %                           data. Each row is a different mode band and
    %                           each column is a different threshold; can be a vector of thresholds in the case one only one threshold per mode
    %         frames_threshold = threshold of frames censored above which a run is excluded
    %
    %
    %         The next two parameters are obtained by running 
    %         compute_mode_connectivity_thresholds.m; these values will
    %         be contained in the file saved in the main directory,
    %         or by running the function most_common_mode_bands.m (called by compute_mode_connectivity_thresholds)
    %         and using its output (common_mode_num and average_freqs)
    %
    %         common_mode_num = most common number of mode bands for a given run out of all
    %                           runs of all subjects
    %         average_freqs = median frequency for each mode band
    %                         calculated as the median of all median frequencies across
    %                         all subjects and runs for that mode band
    %         
    %
    
    % OUTPUTS: file with network properties variable called network_props which contains the network properties for all desired thresholds for
    %          modes and raw data for subject
    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    connectivity_loc = '';  % enter path to directory where connectivity files are (ie '/path_to_connectivity_directory/connectivity_directory')
    network_loc = ''; % enter path to directory where you want to save network property files (ie '/path_to_networkprops_directory/networkprops_directory')
    
    if measure == 1
        measStr = 'xcorr/';
    elseif measure == 2
        measStr = 'coherence/';
    elseif measure == 3
        measStr = 'mi/';
    end
    
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
    
cd(connectivity_loc);   % cd into directory where connectivity files are located

S = dir(strcat(subID,'*')); % gets filename of subject specific connectivity file

pathToSpornsCode = ''; % path to sporns code in order to use network properties functions
% path is of form '/pathtoSpornsCode/SpornsCodeDirectory'
addpath(pathToSpornsCode); % adds path to sporns network property code

filename = S(1).name; 

RS_connectivity = load(filename); % loads connectivity file for subject
if isfield(RS_connectivity,'RS_connectivity') == 1
   RS_connectivity = RS_connectivity.RS_connectivity;
end   

network_props = cell(1,4); 

for r = 1:4
    if length(RS_connectivity) >= r
        if ~isempty(RS_connectivity{r})
            if RS_connectivity{r}.percentCensored < frames_threshold % if run meets percent censored criteria and isn't empty
                if isfield(RS_connectivity{r},'modes') % if modes is a field
                    
                    % extracts connectivity matrix for modes and raw data dependent on measure
                   
                    if measure == 1
                        cmat = RS_connectivity{r}.raw.rmat;
                        modes_temp = RS_connectivity{r}.modes.rmat;
                    elseif measure == 2
                        cmat = RS_connectivity{r}.raw.cxymat;
                        modes_temp = RS_connectivity{r}.modes.cxymat;    
                    elseif measure == 3
                        cmat = RS_connectivity{r}.raw.mimat;
                        modes_temp = RS_connectivity{r}.modes.mimat;
                    end
                    
                    modes_cmat = cell(1,length(average_freqs)); % cell array which will contain the connectivity matrices for each mode band
                    freqs_i = NaN(1,length(average_freqs)); % empty vector size of max number of mode bands across all subjects
                    freqs = RS_connectivity{r}.modes.freqs; % extracts frequency data for all modes
                    network_props{r}.percentCensored = RS_connectivity{r}.percentCensored; % stores percentCensored value
                    props = calculate_network_properties_helper(1,cmat,raw_thresholds,threshold_labels,measure); % gets network properties for raw connectivity matrix for desired thresholds
                    network_props{r}.raw.props = props; % stores in network_props structure
                    
                    modes_cmat(1:length(modes_temp)) = modes_temp; % stores modal connectivity matrices in cell array
                    freqs_i(1:size(freqs,2)) = nanmedian(freqs,1); % populates vector with average frequency for available mode bands
                    
                    if sum(~isnan(freqs_i)) < common_mode_num %if less mode bands than number of most common number of mode bands
                                % need to check to see that mode bands for
                                % run of this particular subject are in the
                                % right location, or whether they need to
                                % be shifted to better match in frequency
                        num_bands = sum(~isnan(freqs_i)); % number of mode bands for subject run
                        temp_freqs_i = NaN(1,length(average_freqs)); % temporary variable
                        shift_band = 0; % how much to shift the mode bands over
                        diff_vals = Inf; % absolute difference between subject mode band frequencies and average frequencies across all subjects
                                % want to minimize this absolute difference
                                
                                % goes through each possible shift in mode
                                % bands and calculates absolute difference
                                % if absolute difference is less than
                                % diff_vals, this new shift value is
                                % considered best
                        for band_iter = 1:(length(average_freqs) - num_bands + 1)
                            temp_freqs_i(band_iter:(band_iter+num_bands-1)) = freqs_i(1:num_bands);
                            temp_diff = nansum(abs(temp_freqs_i - average_freqs).*(1./average_freqs));
                            if temp_diff < diff_vals
                                diff_vals = temp_diff;
                                shift_band = band_iter - 1;
                            end
                            temp_freqs_i = NaN(1,length(average_freqs));
                        end   
                        
                        % shifts all values for mode bands over by
                                % shift_band which was the shift that
                                % minimized the absolute difference between
                                % average frequency of mode band across
                                % subjects and subject specific median frequency
                        freqs_temp = freqs_i;
                        freqs_i(shift_band+1:shift_band+num_bands) = freqs_temp(1:num_bands);
                        modes_cmat(shift_band+1:shift_band+num_bands) = modes_temp(1:num_bands);
                        if shift_band ~= 0
                            freqs_i(1:shift_band) = NaN;
                            modes_cmat(1:shift_band) = {};
                        end
                    end
                    
                    network_props{r}.modes = cell(1,length(average_freqs)); % empty cell array for network properties for each mode band
                    for m = 1:length(modes_cmat)
                        if ~isempty(modes_cmat{m})
                            cmat = modes_cmat{m}; % connectivity matrix for individual mode band
                            thresholds_modes_m = mode_thresholds(m,:); % thresholds for individual mode band
                            props = calculate_network_properties_helper(0,cmat,thresholds_modes_m,threshold_labels, measure); % gets network properties for mode band connectivity matrix for desired thresholds
                            network_props{r}.modes{m}.props = props;
                            network_props{r}.modes{m}.freq = freqs_i(m); % stores in network_props structure as well as average frequency of that mode band for this particular run and subject
                        end
                    end
                end
            end
        end
    end
end

savename = strcat(subID,'_network_properties_',measStr(1:end-1),identifier,'.mat');
cd(network_loc);
save(savename,'network_props'); % saves network property file in desired location for network properties
end
