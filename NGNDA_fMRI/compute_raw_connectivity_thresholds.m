function compute_raw_connectivity_thresholds(measure, cortical_atlas, frames_threshold)
  % purpose: To compute statistical thresholds for median, third quartile, moderate outlying and extreme outlying connectivity from broadband connectomes
    %          Each statistic is calculated from the available sample of runs and bootstrapping is used to estimate corresponding confidence intervals 
    %          For each statistic, the upper confidence interval is used as the threshold   
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
    
    %         frames_threshold = Threshold of % frames with excessive motion above which a data run is excluded
    
    % OUTPUTS: medianThresh = upper confidence interval of median of all median connectivity
    %                        value for each subject
    %         thirdThresh = upper confidence interval of median of all third quartile of connectivity
    %                        value for each subject
    %         modThresh = upper confidence interval of median of all upper moderate outlier connectivity
    %                        value for each subject
    %         extThresh = upper confidence interval of median of all upper extreme outlier connectivity
    %                        value for each subject
    
    main_directory = ''; % path to directory where main code is located
    addpath(main_directory); % adds main directory to matlab path
    connectivity_loc = '';  % enter path to directory where connectivity files are (ie '/path_to_connectivity_directory/connectivity_directory')
   
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
    
    cd(connectivity_loc); % cds into directory where connectivity matrix is located
    S = dir('sub*'); % connectivity files of all subjects in directory
    
    % variables to store median connectivity, 3rd quartile connectivity,
    % moderate outlier connectivity, extreme outlier connectivity for each
    % run of each subject
    medians = NaN(length(S),4);
    quart3s = NaN(length(S),4);
    moderates = NaN(length(S),4);
    extremes = NaN(length(S),4);
    
    parfor i = 1:length(S)
        filename = S(i).name; % subject specific connectivity filename
        RS_connectivity = load(filename); % load RS_connectivity file
        if isfield(RS_connectivity,'RS_connectivity') == 1
            RS_connectivity = RS_connectivity.RS_connectivity;
        end   
        
        % goes through each run of subject
        for r = 1:4
            if length(RS_connectivity) >= r
                if ~isempty(RS_connectivity{r})
                    if RS_connectivity{r}.percentCensored < frames_threshold % if run is below threshold for percent censored
                        if isfield(RS_connectivity{r},'modes') % makes sure modes field is 
                            
                            % extracts connectivity matrix dependent on
                            % measure
                            if measure == 1
                                cmat = RS_connectivity{r}.raw.rmat;
                            elseif measure == 2
                                cmat = RS_connectivity{r}.raw.cxymat;
                            elseif measure == 3
                                cmat = RS_connectivity{r}.raw.mimat;
                            end
                            
                            A = ones(size(cmat));
                            
                            vec = cmat(logical(triu(A,1))); % extracts values in upper triangular portion (excluding diagonal) of connectivity matrix
                            med1 = nanmedian(vec); % takes median of upper triangular 
                            medians(i,r) = med1;
                            quarts = prctile(vec,[25 75]);
                            quart3s(i,r) = quarts(2); % stores 75th percentile connectivity of upper triangular
                            moderates(i,r) = 1.5*(quarts(2)-quarts(1))+med1; % stores upper moderate outlier of upper triangular
                            extremes(i,r) = 3*(quarts(2)-quarts(1))+med1; % stores upper extreme outlier of upper triangular
                        end
                    end
                end
            end
        end 
    end

    medf = @(x) nanmedian(x,'all'); % bootstraps median of each of these groups of summary statistics
    medci = bootci(2000,{medf,medians},'type','bca');
    thirdci = bootci(2000,{medf,quart3s},'type','bca');
    modci = bootci(2000,{medf,moderates},'type','bca');
    extci = bootci(2000,{medf,extremes},'type','bca');
    
    % takes upper confidence interval of each of these statistics and sets
    % these as connectivity thresholds to use
    medianThresh = medci(2);
    thirdThresh = thirdci(2);
    modThresh = modci(2);
    extThresh = extci(2);
    cd(main_directory); % cds into main directory
    c = date;
    % saves all statistical thresholds for raw data in main directory
    save(strcat('connectivity_thresholds_raw',identifier,'_',measStr(1:end-1),'_',c,'.mat'),'medianThresh','thirdThresh','modThresh','extThresh');
    
end
