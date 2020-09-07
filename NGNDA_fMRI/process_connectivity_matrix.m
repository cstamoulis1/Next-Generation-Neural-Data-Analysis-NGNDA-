function process_connectivity_matrix(subID, measure, cortical_atlas, frames_threshold)
% Purpose: To compute connectivity matrices for denoised (mode-reduced and resynthesized) raw data and
%          individual mode band data, using multiple approaches (signal cross-correlation, cross-coherence and Mutual Information (MI)
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: September 6, 2020
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
    %         frames_threshold = upper limit of % of frames censored at which a
    %                            run of data is no longer satisfactory

    % OUTPUTS: connectivity matrix (N nodes by N nodes) for raw data and 
    % individual mode bands for each run based off desired measure
    
    if measure == 1
        measStr = 'xcorr/';
    elseif measure == 2
        measStr = 'coherence/';
    elseif measure == 3
        measStr = 'mi/';
        pathToMIToolbox = ''; % string path to MIToolbox directory needed for MI calculation (ie '/path_to_MIToolbox/MIToolbox'
        addpath(pathToMIToolbox);
    end
    
    
    connectome_loc_new = '';  % enter path to directory where reduced connectome files are located (ie '/path_to_reduced_connectome_directory/reduced_connectome_directory')
    connectivity_loc = '';  % enter path to directory where you want to save connectivity files (ie '/path_to_connectivity_directory/connectivity_directory')
    
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
    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    
    cd(connectome_loc_new); % cd into directory where reduced connectome data is located
 
    S = dir(strcat(subID,'*'));
    filename = S(1).name;   % load reduced connectome data for subject
    
    RS_net = load(filename,'RS_net');
    % adjusts structure if necessary
    if isfield(RS_net,'RS_net') == 1
        RS_net = RS_net.RS_net;
    end
    
% loads information on mode frequency and amplitude;
% adjusts structure if necessary

cd(main_directory);
isNormalized = 0; % to check to make sure data was normalized
for run1 = 1:length(RS_net)
% checks that run is not empty
    if ~isempty(RS_net{run1})
        RS_connectivity{run1}.percentCensored = RS_net{run1}.percentCensored;
        if RS_net{run1}.percentCensored < frames_threshold % check that run meets criteria for percent censored
            if isfield(RS_net{run1},'rawData_normalized') % check to make sure data was normalized
                isNormalized = 1;
				if ~isempty(RS_net{run1}.rawData_normalized)
                    fs = 1/RS_net{run1}.TR; % sampling frequency in Hz
					rawData = RS_net{run1}.rawData_reduced_modes; % extract reduced raw regional data
					numNodes = size(rawData,1); % number of regions
					decomposedNodes = RS_net{run1}.decomposedNodes_characterized; % extract characterized modal data for each node
					modeCharacterized_Numbers = RS_net{run1}.modeCharacterized_Numbers; % extracts to size mode connectivity cell arrays
					newFreq = RS_net{run1}.newFreq; 
					newAmp = RS_net{run1}.newAmp;
					RS_connectivity{run1}.modes.amps = newAmp; % saves frequency and amplitude info into connectivity file
					RS_connectivity{run1}.modes.freqs = newFreq;
				    modes_cmat = cell(size(modeCharacterized_Numbers,2),1); % cell array where connectivity for each mode will be stored
				    modes_lagmat = cell(size(modeCharacterized_Numbers,2),1); % cell array where lags for each mode will be stored
                % runs through each mode
				    for ind = 1:size(modeCharacterized_Numbers,2) % goes through each mode band
						cmat = NaN(numNodes); % initializes connectivity and lag arrays
						lagmat = NaN(numNodes);
						parfor i = 1:numNodes % goes through each node
                            if measure == 3
                                addpath(pathToMIToolbox);
                            end
                            
							if ~isnan(modeCharacterized_Numbers(i,ind)) % if node has a mode in this mode band
								sig1 = decomposedNodes{i}(ind,:); % extracts time series for node in this mode band
								if ~isnan(sig1) % makes sure non nan
									for j = 1:numNodes % loops through nodes
										if ~isnan(modeCharacterized_Numbers(j,ind)) % if node has a mode in this mode band
											sig2 = decomposedNodes{j}(ind,:);
											if ~isnan(sig2)
												if measure == 3 % mi
                                                    cmat(i,j) = mi(sig1',sig2'); % no second output for mutual information, only one mi value
                                                elseif measure == 2 % coherence
                                                    [c, lags] = mscohere(sig1,sig2,[],[],[],fs); % vector of coherence values (c) and vector of frequences (lags) 
                                                    % expressed in terms of fs, where magnitude squared coherence is estimated 
                          % mscohere(sig1,sig2,window,noverlap,nfft,fs)
                          % window = [] = hamming
                          % noverlap = [] = 50 %
                          % nfft = [] = max(256,2p), where p = [log2 N]
                          % and N is the length of the input signal sig1
                                                    % 
                                                    cmax = max(c); % takes peak coherence between node i and j
                                                    cmat(i,j) = cmax; % stores in connectivity matrix
                                                    
                                                    if all(c == 1) % if all coherence is 1, node i = nod j
                                                        lagmat(i,j) = 0; % f could be anything so just set to 0
                                                    else
                                                        lagmat(i,j) = lags(c == cmax); % else take frequency at which coherence is max
                                                    end
                                                    
                                                elseif measure == 1 % if xcorr
                                                    [c, lags] = xcorr(sig1,sig2,'coeff'); % takes xcorr between node i and j
													cmax = max(c); % peak cross correlation value stored in connectivity matrix
													cmat(i,j) = cmax;
                                                    lagmat(i,j) = lags(c == cmax); % lag at which peak xcorr also stored in connectivity matrix
                                                end
											end
										end
									end
								end
							end
						end
						modes_cmat(ind) = {cmat};  % stores connectivity and lag matrices for this mode band in cell array
						modes_lagmat(ind) = {lagmat};
						clear cmat lagmat sig1 sig2
				    end

% saves the mode connectivity matrices to a struct for this run
                    if measure == 3
                        RS_connectivity{run1}.modes.mimat = modes_cmat;
                    elseif measure == 1
                        RS_connectivity{run1}.modes.lagmat = modes_lagmat;
                        RS_connectivity{run1}.modes.rmat = modes_cmat;
                    elseif measure == 2
                        RS_connectivity{run1}.modes.wmat = modes_lagmat;
                        RS_connectivity{run1}.modes.cxymat = modes_cmat;
                    end
					clear modes_lagmat modes_cmat
                
% analyzing the raw data now
					cmat = NaN(numNodes);
					lagmat = NaN(numNodes);
% computes simple cross correlation matrix for the raw data
					parfor i = 1:numNodes
                        if measure == 3
                            addpath(pathToMIToolbox);
                        end
                        
						sig1 = rawData(i,:);
						if ~isnan(sig1)
							for j = 1:numNodes
								sig2 = rawData(j,:);
								if ~isnan(sig2)
                                    if measure == 3
                                        cmat(i,j) = mi(sig1',sig2');
                                    elseif measure == 2
                                        [c, lags] = mscohere(sig1,sig2,[],[],[],fs);
                                        cmax = max(c);
                                        cmat(i,j) = cmax;
                                                    
                                        if all(c == 1)
                                            lagmat(i,j) = 0;
                                        else
                                            lagmat(i,j) = lags(c == cmax);
                                        end
                                                    
                                    elseif measure == 1
                                        [c, lags] = xcorr(sig1,sig2,'coeff');
									    cmax = max(c);
									    cmat(i,j) = cmax;
                                        lagmat(i,j) = lags(c == cmax);
                                    end
                                end
							end
						end
					end
% saves to matrix to the same structure as the modes

					if measure == 3
                        RS_connectivity{run1}.raw.mimat = cmat;
                    elseif measure == 1
                        RS_connectivity{run1}.raw.lagmat = lagmat;
                        RS_connectivity{run1}.raw.rmat = cmat;
                    elseif measure == 2
                        RS_connectivity{run1}.raw.wmat = lagmat;
                        RS_connectivity{run1}.raw.cxymat = cmat;
                    end
                    
					clear cmat lagmat sig1 sig2 c lags;
               end
            end
        end
    end
end

if isNormalized
    str1 = strcat(subID,'_RS_connectivity_baselineYear1',identifier,'_',measStr(1:end-1),'.mat'); % filename for connectivity file, if you want to change please keep subject ID as initial
    % component of filename for rest of code
    save(str1,'RS_connectivity','-v7.3'); % saves connectivity and moves into designated connectivity folder
    movefile(str1,connectivity_loc);
end

end
