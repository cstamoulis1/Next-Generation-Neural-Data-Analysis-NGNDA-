function [lower_lim_amp, upper_lim_amp, lower_lim_freq, upper_lim_freq] = calculate_mode_thresh(lower_filter_bound, upper_filter_bound, frames_threshold)
% purpose: to calculate statistical thresholds for mode amplitude and
%          frequency in order to discard modes which have characteristic
%          frequencies that do not meet these requirements
% Called by modes_classification_connectome_reduction.m 
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: September 5, 2020
% INPUTS: upper_filter_bound = upper frequency used in bandpass filter
%         lower_filter_bound = lower frequency used in bandpass filter
%         frames_threshold = upper limit of % of frames censored at which a
%                            run of data is no longer satisfactory
% OUTPUTS: lower_lim_amp = lower amplitude threshold in log based form
%          upper_lim_amp = upper ampltiude threshold in log based form
%          lower_lim_freq = lower frequency threshold
%          upper_lim_freq = upper frequency threshold

 % keep track of which cortical atlas you are using and where files are
 % located if you are trying multiple parcellation schemes

    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    cd(main_directory); % cd into main directory
    
    if upper_filter_bound <= lower_filter_bound
        error('Please enter appropriate bounds for band pass filter. The upper bound must be greater than the lower bound.');
    end
    
    TR = 0.8; % enter TR, for ABCD TR = 0.8

    fs = 1/TR; % sampling and nyquist frequencies
    fN = fs/2;

% constructs 3rd order elliptical bandpass filter used in connectome
% construction step
    n = 3;
    Rp = 0.5;
    Rs = 20;
    ftype = 'bandpass';
    Wp = [lower_filter_bound upper_filter_bound]/fN;
    [b_band,a_band] = ellip(n,Rp,Rs,Wp,ftype);

    [h, w] = freqz(b_band,a_band,10000);
    

    f = w*fN/pi;
    dB = mag2db(abs(h));
    [~, f0] = crossing(dB,f,-20); % find transition zone on filter bounds (where -20 dB attentuation first occurs)

    lower_lim_freq = lower_filter_bound + abs(lower_filter_bound - f0(1)); % transition zone based lower frequency bound
    upper_lim_freq = upper_filter_bound - abs(upper_filter_bound - f0(2)); % transition zone based upper frequency bound
    clearvars -except lower_lim_freq upper_lim_freq cortical_atlas lower_filter_bound upper_filter_bound frames_threshold
    
    % keep track of which cortical atlas you are using and where files are located if you are trying multiple parcellation schemes
    modeinfo_loc = ''; % enter path to directory containing mode info files (ie '/pathtomodeinfodirectory/modeinfodirectory')
    % paths to mode info
    
    cd(modeinfo_loc); % cds into modeinfo directory
    S = dir('sub*'); % extracts all mode info files

    amp_meds = []; % all median amplitudes across all nodes in each mode band for all runs and subjects
    for s = 1:length(S)
        
        RS_modeInfo = load(S(s).name); % loads mode info file
        if isfield(RS_modeInfo,'RS_modeInfo') == 1
            RS_modeInfo = RS_modeInfo.RS_modeInfo;
        end
        
        for j = 1:length(RS_modeInfo)
            if ~isempty(RS_modeInfo{j})
                if RS_modeInfo{j}.percentCensored < frames_threshold 
                    if isfield(RS_modeInfo{j},'freq')
					    % modes with appropriate frequency range
					    vals2include = (RS_modeInfo{j}.freq > lower_lim_freq) & (RS_modeInfo{j}.freq < upper_lim_freq);
                            
						freq = RS_modeInfo{j}.freq.*vals2include; % characteristic frequencies of modes that met the frequency cutoffs
						freq(~vals2include) = NaN;
						newFreq = freq;
						amp = RS_modeInfo{j}.amp.*vals2include;
						% sets modes outside frequency range to nan
						amp(~vals2include) = NaN; 
						newAmp = amp;
						modenums = 1:size(freq,2);
                            
% Check to ensure mode order with monotonically
% decreasing frequency; these checks are necessary to make sure only non-artifactual modes are included in the threshold estimation
                        
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
                            
                        % takes log of amplitudes since underlying
                        % distribution is log normal and thresholds
                        % cannot be negative
						newAmp = log(newAmp);
							
					    meds = nanmedian(newAmp,1); % calculates median amplitude of each mode band
						% appends extreme outliers to vector
						amp_meds = [amp_meds meds];
					end
                end
            end 
        end
    end
    amp_meds(isnan(amp_meds)) = []; % any nans are removed
    
    % calculates CI on moderate outlier thresholds
    func_ext_up = @(x) nanmedian(x) + 1.5*iqr(x);
    func_ext_lo = @(x) nanmedian(x) - 1.5*iqr(x);
    upci = bootci(2000,{func_ext_up,amp_meds},'type','bca');
    lowci = bootci(2000,{func_ext_lo,amp_meds},'type','bca');
    
    % calculates upper and lower thresholds for amplitude
    % these are log based
    upper_lim_amp = upci(1);
    lower_lim_amp = lowci(2);
end
