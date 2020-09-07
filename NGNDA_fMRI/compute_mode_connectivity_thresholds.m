function compute_mode_connectivity_thresholds(measure, cortical_atlas, frames_threshold)
  % purpose: To compute statistical thresholds for median, third quartile, moderate outlying and extreme outlying connectivity from narrowband mode connectomes
% Each statistic is calculated from the available sample of runs and bootstrapping is used to estimate corresponding confidence intervals
% The upper confidence interval for each statistic is used as the threshold         
% Also calls most_common_mode_bands.m
  % When estimating connectivity thresholds for individual modes, the classification of these modes into the same frequency bands needs to be ensured prior to
  % the estimation; the function most_common_modes_bands calculates both the most common (across the sample) number of modes and their individual frequencies
  % any misclassified mode is moved to the correct frequency band prior to the connectivity threshold estimation (a very rare occurrence)
  % Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
  % Last version: September 5, 2020 

  % INPUTS: 
    %         measure = 1 is peak cross-correlation
    %         measure = 2 is peak coherence
    %         measure = 3 is mutual information
    %
    %
    %         cortical_atlas = 1 is Gordon
    %         cortical_atlas = 2 is Schaefer 1000
    %         cortical_atlas = 3 is Schaefer 300
    %         cortical_atlas = 4 is Schaefer 800
    %         cortical_atlas = 5 is Schaefer 100
    
    %         frames_threshold = Threshold  of % frames with excessive motion above which a data run is excluded

    
    % OUTPUTS: medianThresh = upper confidence interval of median for each mode band of all median connectivity
    %                        value for each subject
    %         thirdThresh = upper confidence interval of median for each mode band of all third quartile of connectivity
    %                        value for each subject
    %         modThresh = upper confidence interval of median for each mode band of all upper moderate outlier connectivity
    %                        value for each subject
    %         extThresh = upper confidence interval of median for each mode band of all upper extreme outlier connectivity
    %                        value for each subject
    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    connectivity_loc = '';  % enter path to directory connectivity files are located (ie '/path_to_connectivity_directory/connectivity_directory')
   
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
    
    cd(main_directory)
    [common_mode_num, average_freqs] = most_common_mode_bands(frames_threshold);  
    % gets most common number of mode bands across all subjects and the average frequencies of these mode bands
    
    cd(connectivity_loc); % cd to directory where connectivity files are located
    S = dir('sub*'); % gets all subjects with connectivity files in directory
    
    
	freqs_all = cell(length(S),4); % cell array where a cell is present for each run of each subject
    % the contents of the cell array are the average (median) frequencies of each
    % mode band for that run for that particular subject
    
    
    medians = NaN(length(S),4,length(average_freqs)); % vector that gets median connectivity for each mode band for each run of each subject
    extremes = NaN(length(S),4,length(average_freqs));% vector that gets extreme outlier connectivity for each mode band for each run of each subject
    moderates = NaN(length(S),4,length(average_freqs));% vector that gets moderate outlier connectivity for each mode band for each run of each subject
    quart3s = NaN(length(S),4,length(average_freqs));% vector that gets third quartile connectivity for each mode band for each run of each subject
    
    parfor i = 1:length(S) % goes through each subject
        filename = S(i).name;    
        RS_connectivity = load(filename); % loads subject specific connectivity
        if isfield(RS_connectivity,'RS_connectivity') == 1
            RS_connectivity = RS_connectivity.RS_connectivity;
        end   
        
        for r = 1:4 % loops through runs
            if length(RS_connectivity) >= r
                if ~isempty(RS_connectivity{r})
                    if RS_connectivity{r}.percentCensored < frames_threshold % if run meets criteria for percent censored
                        if isfield(RS_connectivity{r},'modes') % if modes is a field
                    	    freqs = RS_connectivity{r}.modes.freqs; % extracts mode frequencies for individual subject run
                            
                            % extracts connectivity matrix dependent on
                            % measure
                            if measure == 1
                                cmat = RS_connectivity{r}.modes.rmat;
                            elseif measure == 2
                                cmat = RS_connectivity{r}.modes.cxymat;
                            elseif measure == 3
                                cmat = RS_connectivity{r}.modes.mimat;
                            end
                            
                            medians_conn = NaN(1,length(average_freqs)); % vectors for all mode bands to store median connectivity, and other summary statistics
                            quart3s_conn = NaN(1,length(average_freqs)); 
                            moderates_conn = NaN(1,length(average_freqs));
                            extremes_conn = NaN(1,length(average_freqs));
                            
                            % goes through each mode band for particular
                            % subject run
                            for band = 1:length(cmat)
                                A = ones(size(cmat{band})); 
                                vec = cmat{band}(logical(triu(A,1))); % extracts upper triangular (not including diagonal) connectivity values
                                med_conn = nanmedian(vec); % gets summary statistics for specific mode band of this run
                                medians_conn(band) = med_conn; 
                                quarts = prctile(vec,[25 75]);
                                quart3s_conn(band) = quarts(2);
                                moderates_conn(band) = 1.5*(quarts(2)-quarts(1))+med_conn;
                    	        extremes_conn(band) = 3*(quarts(2)-quarts(1))+med_conn; 
                            end
                            
                            freqs_i = NaN(1,length(average_freqs)); % empty vector size of max number of mode bands across all subjects
                            freqs_i(1:length(nanmedian(freqs,1))) = nanmedian(freqs,1); % populates vector with average frequency for available mode bands
                            if sum(~isnan(freqs_i)) < common_mode_num % if less mode bands than number of most common number of mode bands
                                % need to check to see that mode bands for
                                % run of this particular subject are in the
                                % right location, or whether they need to
                                % be shifted to better match in frequency
                                
                                num_bands = sum(~isnan(freqs_i)) % number of mode bands for subject run
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
                                    temp_diff = nansum(abs(temp_freqs_i - average_freqs).*(1./average_freqs));%.*normalizing_width);
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
                                medians_conn(shift_band+1:shift_band+num_bands) = medians_conn(1:num_bands);
                                extremes_conn(shift_band+1:shift_band+num_bands) = extremes_conn(1:num_bands);
                                quart3s_conn(shift_band+1:shift_band+num_bands) = quart3s_conn(1:num_bands);
                                moderates_conn(shift_band+1:shift_band+num_bands) = moderates_conn(1:num_bands);
                                
                                if shift_band ~= 0 % if values shifted
                                    medians_conn(1:shift_band) = NaN;
                                    extremes_conn(1:shift_band) = NaN;
                                    quart3s_conn(1:shift_band) = NaN;
                                    moderates_conn(1:shift_band) = NaN;
                                    freqs_i(1:shift_band) = NaN;
                                end
                            end
                            
                            % stores run summary statistics in larger
                            % sample arrays
                    	    freqs_all(i,r) = {freqs_i};
                            medians(i,r,:) = medians_conn;
                            extremes(i,r,:) = extremes_conn;
                            quart3s(i,r,:) = quart3s_conn;
                            moderates(i,r,:) = moderates_conn;
                        end
                    end
                
                end
            end
        end
    end
    
    % bootstraps median of each summary statistic for each mode band across
    % all subjects and runs
    medf = @(x) nanmedian(x,'all');
    for i = 1:size(medians,3)
        medci(1:2,i) = bootci(2000,{medf,medians(:,:,i)},'type','bca');
        thirdci(1:2,i) = bootci(2000,{medf,quart3s(:,:,i)},'type','bca');
        modci(1:2,i) = bootci(2000,{medf,moderates(:,:,i)},'type','bca');
        extci(1:2,i) = bootci(2000,{medf,extremes(:,:,i)},'type','bca');
    end
    
    % takes upper confidence interval of these values
    % each of these variables has length of number of mode bands
    % so a statistical threshold for each mode band for these four
    % different summary statistics
    medianThresh = medci(2,:);
    thirdThresh = thirdci(2,:);
    modThresh = modci(2,:);
    extThresh = extci(2,:);
    
    % cds to main directory and saves mode thresholds
    cd(main_directory);
    c = date;
    save(strcat('connectivity_thresholds_mode_bands',identifier,'_',measStr(1:end-1),'_',c,'.mat'),'medianThresh','thirdThresh','modThresh','extThresh','freqs_all','common_mode_num','average_freqs');
end 
