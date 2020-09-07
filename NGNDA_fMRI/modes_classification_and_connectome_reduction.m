function modes_classification_and_connectome_reduction(cortical_atlas, lower_filter_bound, upper_filter_bound, frames_threshold)

% Purpose: To denoise fMRI data by removing signal components (modes) that are outside of the expected frequency range
%          and those with amplitudes outside a statistically estimated range.
%          Based on their characteristic (dominant) frequencies, modes are classified in distinct  
%          mode bands that are consistent across brains, to enable comparison between brains
% Also calls function calculate_mode_thresh.m to calculate statistical thresholds for frequency and amplitude and save reduced_connectomes to save denoised (mode-reduced) data
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: September 6, 2020  
% INPUTS:
%         lower_filter_bound = lower frequency for bandpass filter
%         upper_filter_bound = upper frequency for bandpass filter
%         frames_threshold = threshold for % frames censored above which a data run is excluded


% OUTPUTS: Denoised (mode-reduced) connectome structure for each subject. Structure
%          contains original RS_net structure from RS_connectomes folder
%          and mode-reduced data (resynthesized data is a linear superposition of components (modes) with non-artifactual contributions to the original signal; noise-related components that do not meet criteria are excluded
																			                   %          Mode classification: some mode bands exist only for a small number of fMRI signals in a dataset
%          These are excluded. Only modes that are present in more than 40% of the parcellated fMRI signals (nodes)
%          are kept in the dataset. 
%          Modes are reclassified into these bands for a particular node if they have frequencies
%          within the extreme outliers limits of that mode band.


% Function calculate_mode_characteristics.m includes the estimation of thresholds based on bootstrapping a summary statistic (median)
% A sufficiently large cohort with mode characteristic parameters should be processed so that the estimated thresholds reflect 
% the large-scale dataset (e.g., ABCD data)

    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    cd(main_directory); % cd into main directory
    
    % calculates frequency thresholds based off of bandpass filter design
    % calculates amplitude threshold based on bootstrapping the median mode amplitudes
    % any modes outside the frequency or amplitude thresholds are excluded
    [lower_lim_amp, upper_lim_amp, lower_lim_freq, upper_lim_freq] = calculate_mode_thresh(lower_filter_bound, upper_filter_bound, frames_threshold);
    
    modeinfo_loc = ''; % enter path to directory containing mode info files (ie '/path_to_mode_info_directory/mode_info_directory')
    % paths to mode info
    connectome_loc = '';  % enter path to directory containing connectome files (ie '/path_to_connectome_directory/connectome_directory')
    connectome_loc_new = '';  % enter path to directory where you want to save reduced connectome files (ie '/path_to_reduced_connectome_directory/reduced_connectome_directory')
    
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
    
    cd(modeinfo_loc) % cd into mode info location
    S = dir('sub*'); % gets all subjects with mode info
    
    % Track subjects that have misclassified modes to save reduced connectome
    problematic_subs = cell(length(S),1);

%True parallelization is along the cohort dimension since the run-specific calculations are completed very quickly
    parfor s = 1:length(S)
        cd(modeinfo_loc); 
        RS_modeInfo = load(S(s).name); % loads mode info for subject
        cd(connectome_loc); % cd into directory with connectomes
        k_ind = regexp(S(s).name, '_', 'once')-1;
        subID = S(s).name(1:k_ind); % extracts subject id from mode info file name for subject
        filename = dir(strcat(subID,'*')); % gets filename for connectome for subject
        RS_net = load(filename(1).name); % loads connectome for subject
        
        if isfield(RS_modeInfo,'RS_modeInfo') == 1
            RS_modeInfo = RS_modeInfo.RS_modeInfo;
        end
        
        if isfield(RS_net,'RS_net') == 1
            RS_net = RS_net.RS_net;
        end
        
        isNormalized = 0;
        for j = 1:length(RS_net)
            if ~isempty(RS_net{j})
                if RS_net{j}.percentCensored < frames_threshold % if runs are not empty, are normalized, and meet threshold for percent censored
                    if isfield(RS_net{j},'rawData_normalized')
                        isNormalized = 1;
                        if ~isempty(RS_net{j}.rawData_normalized)
                        % finds modes that meet amplitude threshold and are within
                        % proper frequency range
                        %% frequency thresholding
                            % keeps modes that meet frequency threshold
                            % criteria
                            numNodes = size(RS_modeInfo{j}.freq,1);
                            vals2incFreq = (RS_modeInfo{j}.freq > lower_lim_freq) & (RS_modeInfo{j}.freq < upper_lim_freq);
							freq = RS_modeInfo{j}.freq.*vals2incFreq;
							freq(~vals2incFreq) = NaN;
                
							amp = RS_modeInfo{j}.amp.*vals2incFreq;
							amp(~vals2incFreq) = NaN; 

							newFreq = freq;
							newAmp = amp;
							modenums = 1:size(RS_modeInfo{j}.freq,2);
                
%% Monotonicity checks: Given the unsupervised EMD decomposition approach, spurious modes can be estimated
% Given the the decomposition corresponds to a dyadic filter, mode frequencies need to be monotonically decreasing
% Modes with characteristic frequencies that are not monotonically decreasing are removed

							for node = 1:size(freq,1)
								nodefreqs = freq(node,:);
								nonNaNmodenums = modenums(~isnan(nodefreqs));
								nonNaNnodefreqs = nodefreqs(~isnan(nodefreqs));
                    
								diffs_between_modes = diff(nonNaNnodefreqs);
								modenums_ofDiffs = nonNaNmodenums(2:end);
								% mode frequencies with less than 0.001 Hz difference,
								% considered identical.
								troublemodes_higher = modenums_ofDiffs(diffs_between_modes > 0.001);
								troublemodes_same = modenums_ofDiffs(abs(diffs_between_modes) <= 0.001);
                    
                            % removes higher frequency mode
								if ~isempty(troublemodes_higher)
									newFreq(node,troublemodes_higher) = NaN;
									newAmp(node,troublemodes_higher) = NaN;
								end
                    
                            % removes mode with lower amplitude if two are the same
								if ~isempty(troublemodes_same)
									troublemodes_min_1 = nonNaNmodenums(find(ismember(nonNaNmodenums,troublemodes_same))-1);
									for m = 1:length(troublemodes_same)
										if amp(node,troublemodes_same(m)) > amp(node,troublemodes_min_1(m))
											newFreq(node,troublemodes_min_1(m)) = NaN;
											newAmp(node,troublemodes_min_1(m)) = NaN;
										else
											newFreq(node,troublemodes_same(m)) = NaN;
											newAmp(node,troublemodes_same(m)) = NaN;
										end
									end
								end
							end
                
							freq = newFreq; amp = newAmp;
                
							%% reconstruct raw signal from good modes
                            
                            % reconstructs raw signal from modes meeting
                            % frequency criteria and after monotonicity
                            % correction
                           
							modesKept = [];
							RMSE = [];
							rawData_new = NaN(size(RS_net{j}.rawData_normalized));
							decomposedNodes_new = {};
                    
							for node = 1:size(freq,1)
								nodefreqs = freq(node,:);
								nonNaNmodenums = modenums(~isnan(nodefreqs));
                                % new raw signal for node is sum of all mode time
                                % series meeting requirements
								if ~isempty(nonNaNmodenums)
									rawData_new(node,:) = nansum(RS_net{j}.decomposedNodes{node}(nonNaNmodenums,:),1);
                                end
                                % mode numbers kept are recorded, and
                                % individual mode time series used in
                                % reconstructing the raw signal are stored
								modesKept(node,1:length(nonNaNmodenums)) = nonNaNmodenums;
								decomposedNodes_new(end+1) = {RS_net{j}.decomposedNodes{node}(nonNaNmodenums,:)};
                                % computes root mean squared error between
                                % raw signal and reduced raw signal
								RMSE(end+1) = sqrt(nanmean((RS_net{j}.rawData_normalized(node,:) - rawData_new(node,:)).^2));
							end
                    
							modesKept(modesKept == 0) = NaN;
                
							RS_net{j}.rawData_reduced_modes = rawData_new;
							RS_net{j}.modesKept = modesKept;
							RS_net{j}.decomposedNodes_reduced = decomposedNodes_new;
							RS_net{j}.RMSE = RMSE;
                
 %% Threshold modes based on amplitudes      
                            % keeps modes that meet amplitude threshold
                            % criteria
							vals2incAmp = (log(amp) > lower_lim_amp) & (log(amp) < upper_lim_amp);
							freq = freq.*vals2incAmp;
							freq(~vals2incAmp) = NaN;
                
							amp = amp.*vals2incAmp;
							amp(~vals2incAmp) = NaN; 
                
							newFreq = freq;
							newAmp = amp;
                
                        %% mode characterization                 
							fracData = sum(~isnan(newFreq),1)/numNodes;
							indNaN = find(fracData <= 0.4);
                            % any mode bands with less than 40% 
                            % nodes having a mode in that mode band are
                            % excluded from further analysis
                    % mode bands that are kept
                            % keeps all mode bands meeting above criteria
							bands2keep = setdiff(modenums,indNaN);
							newFreq(:,indNaN) = []; newAmp(:,indNaN) = [];
                            
							freqMedians = nanmedian(newFreq,1);
							freq_first_third_qrt = prctile(newFreq,[25 75],1);
							freqIQRs = freq_first_third_qrt(2,:) - freq_first_third_qrt(1,:);
                            
                            % calculates upper and lower moderate outliers
                            % frequencies
                            % on each mode band to aid in mode
                            % reclassification
							UBound = freqMedians + freqIQRs*1.5;
							LBound = freqMedians - freqIQRs*1.5;

							nodes_list = 1:numNodes;
							modeCharacterized_Numbers = repmat(bands2keep,numNodes,1);
							modeCharacterized_Numbers(isnan(newFreq)) = NaN;
                            
                            % goes through each mode band. If a node is
                            % outside the moderate outlier frequencies for
                            % that mode, or if it is nan, that mode is
                            % considered poorly characterized and must be
                            % addressed
                            
                            % if nan, looks if there is a mode in the mode
                            % band above and below in frequency of its
                            % current mode band location. If either of
                            % these modes can be recharacterized to current
                            % mode band of the nan and fall within the
                            % frequency range, in addition to not being
                            % well characterized within its own frequency
                            % range, that mode is recharacterized to fill
                            % the space of the nan, and the mode band it
                            % came from is then assigned nan
                            % if a mode is present but falls outside the
                            % frequency range for that mode band, checks if
                            % a mode above or below is a better fit, if so
                            % either of those modes are assigned that
                            % mode band as long as they are not well
                            % characterized in their own mode band
                            % Also, if mode is mischaracterized in current
                            % mode band but fits better in the mode band
                            % above or below, and no mode is present there
                            % or if present, the mode is outside the
                            % frequency bounds for that mode band, the
                            % current mode is assigned that value.
                            % 
                            % we proceed from highest frequency mode band
                            % to lowest frequency mode band
                          
							for mode_band = 1:length(bands2keep)
								band_freqs = newFreq(:,mode_band);
								badNodes = nodes_list((band_freqs >= UBound(mode_band) | band_freqs <= LBound(mode_band)) | isnan(band_freqs));
                    
								for node = badNodes
									nodefreqs = freq(node,:);
									nonNaNmodenums = modenums(~isnan(nodefreqs));
									current_mode = bands2keep(mode_band);
									prev_mode = nanmax(nonNaNmodenums(nonNaNmodenums < current_mode));
									next_mode = nanmin(nonNaNmodenums(nonNaNmodenums > current_mode));
                        %% action if nan
									if length(bands2keep) == 1
										if ~isempty(prev_mode) && (freq(node,prev_mode) < UBound(mode_band) && freq(node,prev_mode) > LBound(mode_band))
											modeCharacterized_Numbers(node,mode_band) = prev_mode;
											newFreq(node,mode_band) = freq(node,prev_mode);
											newAmp(node,mode_band) = amp(node,prev_mode);
										elseif ~isempty(next_mode) && (freq(node,next_mode) < UBound(mode_band) && freq(node,next_mode) > LBound(mode_band))
											modeCharacterized_Numbers(node,mode_band) = next_mode;
											newFreq(node,mode_band) = freq(node,next_mode);
											newAmp(node,mode_band) = amp(node,next_mode);
										end

									elseif length(bands2keep) > 1
                            
										if isnan(newFreq(node,mode_band))
											if mode_band == 1
												if ~isempty(prev_mode) && (freq(node,prev_mode) < UBound(mode_band) && freq(node,prev_mode) > LBound(mode_band))
													modeCharacterized_Numbers(node,mode_band) = prev_mode;
													newFreq(node,mode_band) = freq(node,prev_mode);
													newAmp(node,mode_band) = amp(node,prev_mode);
												end
											elseif mode_band == length(bands2keep)
												if ~isempty(next_mode) && (freq(node,next_mode) < UBound(mode_band) && freq(node,next_mode) > LBound(mode_band))
													modeCharacterized_Numbers(node,mode_band) = next_mode;
													newFreq(node,mode_band) = freq(node,next_mode);
													newAmp(node,mode_band) = amp(node,next_mode);
												end
											end
              %% action if not nan              
										else
                            %% action if first mode band
											if mode_band == 1
												if isnan(newFreq(node,mode_band+1)) && (newFreq(node,mode_band) < UBound(mode_band+1) && newFreq(node,mode_band) > LBound(mode_band+1))
													modeCharacterized_Numbers(node,mode_band+1) = current_mode;
													newFreq(node,mode_band+1) = newFreq(node,mode_band);
													newAmp(node,mode_band+1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;
													if ~isempty(prev_mode) && (freq(node,prev_mode) < UBound(mode_band) && freq(node,prev_mode) > LBound(mode_band))
														newAmp(node,mode_band) = amp(node,prev_mode);
														newFreq(node,mode_band) = freq(node,prev_mode);
														modeCharacterized_Numbers(node,mode_band) = prev_mode;
													end
                                            
												elseif (newFreq(node,mode_band) < UBound(mode_band+1) && newFreq(node,mode_band) > LBound(mode_band+1)) && ~(newFreq(node,mode_band+1) < UBound(mode_band+1) && newFreq(node,mode_band+1) > LBound(mode_band+1))
													newFreq(node,mode_band+1:end) = newFreq(node,mode_band:end-1);
													newAmp(node,mode_band+1:end) = newAmp(node,mode_band:end-1);
													modeCharacterized_Numbers(node,mode_band+1:end) = bands2keep(mode_band:end-1);
													modeCharacterized_Numbers(node,mode_band) = NaN;
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													if ~isempty(prev_mode) && (freq(node,prev_mode) < UBound(mode_band) && freq(node,prev_mode) > LBound(mode_band))
														newAmp(node,mode_band) = amp(node,prev_mode);
														newFreq(node,mode_band) = freq(node,prev_mode);
														modeCharacterized_Numbers(node,mode_band) = prev_mode;
													end
                                                
												elseif ~isempty(prev_mode) && (freq(node,prev_mode) < UBound(mode_band) && freq(node,prev_mode) > LBound(mode_band))
													newAmp(node,mode_band) = amp(node,prev_mode);
													newFreq(node,mode_band) = freq(node,prev_mode);
													modeCharacterized_Numbers(node,mode_band) = prev_mode;
												end
                                
                            %% action if not first or last mode band
											elseif (mode_band ~= 1 && mode_band ~= length(bands2keep))
												if isnan(newFreq(node,mode_band-1)) && (newFreq(node,mode_band) < UBound(mode_band-1) && newFreq(node,mode_band) > LBound(mode_band-1))
													modeCharacterized_Numbers(node,mode_band-1) = current_mode;
													newFreq(node,mode_band-1) = newFreq(node,mode_band);
													newAmp(node,mode_band-1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;
												elseif isnan(newFreq(node,mode_band+1)) && (newFreq(node,mode_band) < UBound(mode_band+1) && newFreq(node,mode_band) > LBound(mode_band+1))
													modeCharacterized_Numbers(node,mode_band+1) = current_mode;
													newFreq(node,mode_band+1) = newFreq(node,mode_band);
													newAmp(node,mode_band+1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;
                                                
												elseif (newFreq(node,mode_band) < UBound(mode_band-1) && newFreq(node,mode_band) > LBound(mode_band-1)) &&  ~(newFreq(node,mode_band-1) < UBound(mode_band-1) && newFreq(node,mode_band-1) > LBound(mode_band-1))
													modeCharacterized_Numbers(node,mode_band-1) = current_mode;
													newFreq(node,mode_band-1) = newFreq(node,mode_band);
													newAmp(node,mode_band-1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;

												elseif (newFreq(node,mode_band) < UBound(mode_band+1) && newFreq(node,mode_band) > LBound(mode_band+1)) && ~(newFreq(node,mode_band+1) < UBound(mode_band+1) && newFreq(node,mode_band+1) > LBound(mode_band+1))
													newFreq(node,mode_band+1:end) = newFreq(node,mode_band:end-1);
													newAmp(node,mode_band+1:end) = newAmp(node,mode_band:end-1);
													modeCharacterized_Numbers(node,mode_band+1:end) = bands2keep(mode_band:end-1);
													modeCharacterized_Numbers(node,mode_band) = NaN;
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN; 
												end
                                %% action if last mode band
											elseif mode_band == length(bands2keep)
												if isnan(newFreq(node,mode_band-1)) && (newFreq(node,mode_band) < UBound(mode_band-1) && newFreq(node,mode_band) > LBound(mode_band-1))
													modeCharacterized_Numbers(node,mode_band-1) = current_mode;
													newFreq(node,mode_band-1) = newFreq(node,mode_band);
													newAmp(node,mode_band-1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;
													if ~isempty(next_mode) && (freq(node,next_mode) < UBound(mode_band) && freq(node,next_mode) > LBound(mode_band))
														modeCharacterized_Numbers(node,mode_band) = next_mode;
														newAmp(node,mode_band) = amp(node,next_mode);
														newFreq(node,mode_band) = freq(node,next_mode);
													end
												elseif (newFreq(node,mode_band) < UBound(mode_band-1) && newFreq(node,mode_band) > LBound(mode_band-1)) &&  ~(newFreq(node,mode_band-1) < UBound(mode_band-1) && newFreq(node,mode_band-1) > LBound(mode_band-1))
													modeCharacterized_Numbers(node,mode_band-1) = current_mode;
													newFreq(node,mode_band-1) = newFreq(node,mode_band);
													newAmp(node,mode_band-1) = newAmp(node,mode_band);
													newFreq(node,mode_band) = NaN;
													newAmp(node,mode_band) = NaN;
													modeCharacterized_Numbers(node,mode_band) = NaN;
													if ~isempty(next_mode) && (freq(node,next_mode) < UBound(mode_band) && freq(node,next_mode) > LBound(mode_band))
														modeCharacterized_Numbers(node,mode_band) = next_mode;
														newAmp(node,mode_band) = amp(node,next_mode);
														newFreq(node,mode_band) = freq(node,next_mode);
													end
												elseif ~isempty(next_mode) && (freq(node,next_mode) < UBound(mode_band) && freq(node,next_mode) > LBound(mode_band))
													modeCharacterized_Numbers(node,mode_band) = next_mode;
													newAmp(node,mode_band) = amp(node,next_mode);
													newFreq(node,mode_band) = freq(node,next_mode);
												end
											end
										end
  
									end	
								end
                            end
                            % stores all characterized mode amplitudes and
                            % frequencies for the mode bands that were kept
							RS_net{j}.newFreq = newFreq;
							RS_net{j}.newAmp = newAmp;
                            
                            % the numbers of the modes in the original
                            % decomposed mode structure are stored in their
                            % characterized mode bands
							RS_net{j}.modeCharacterized_Numbers = modeCharacterized_Numbers;
							decomposedNodes_characterized = cell(numNodes,1);
                            
                            % extracts time series for all kept characterized modes and
                            % stores them in the reduced connectome
                            % structure
							for node = 1:numNodes
								modes = modeCharacterized_Numbers(node,:);
								nodeData = NaN(length(bands2keep),size(RS_net{j}.rawData,2));
								for m = 1:length(bands2keep)
									if ~isnan(modes(m))
										nodeData(m,:) = RS_net{j}.decomposedNodes{node}(modes(m),:);
									end
								end
								decomposedNodes_characterized(node) = {nodeData};
							end
							RS_net{j}.decomposedNodes_characterized = decomposedNodes_characterized;
                
						end
                    end
                else % if data doesnt meet percent censored criteria stores empty arrays for that run
                    RS_net{j}.decomposedNodes = {};
                    RS_net{j}.rawData = [];
                    RS_net{j}.rawData_normalized = [];                
                end
            end
        end
        
        cd(main_directory);
        
        str1 = sprintf(strcat('%s_RS_net_new-baselineYear1Arm1',identifier,'.mat'),subID);
        if isNormalized
            try
                save_reduced_connectomes(str1,RS_net,connectome_loc_new); % saves reduced connectome to desired location
		    catch % if run contains nodes with misclassified modes, save subject id in cell
                problematic_subs(s) = {subID};
            end
        end
        
    end
    
    cd(main_directory)
    c = date;
    if sum(~cellfun(@isempty,problematic_subs))
        save(strcat('problematic_subs_mode_classification',c,'.mat'),'problematic_subs'); % if trying to save or load one of the reduced connectomes or normal connectomes is problematic or fails, keep track of these by saving subject id to main directory
    end
    
    
end
