function [decomposedNodes, rawData, rawData_normalized, percentCensored, badFrames, FD, TR] = construct_connectome_helper(subID,n,sliceOrder,frames2remove,parcellation,frames_threshold)
% Also calls function glmfitResidOnly.m
% purpose: to process an individual run of fMRI data. Perform preprocessing
%          steps (remove initial frames, coregistration, normalization of
%          functional image, motion regression, filtering, and
%          interpolation). Parcellate preprocessed data into specific
%          parcellation scheme with various regions, and use ensemble
%          empirical mode decomposition to decompose each regions signal
%          into multiple modal signals.
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School         
% Last version: September 5, 2020
% INPUTS: subID = identification code for each subject                                   
%                 (ie, sub-NDARINVXXXXXXXX)
%         n = current run being processed
%         sliceOrder = order fMRI slices were acquired in milliseconds
%         frames2remove = number of initial frames to remove based on
%                         scanner type
%         parcellation = 3d array with integer values identifying region
%         frames_threshold = upper limit of % of frames censored, above which a
%                            data run is excluded from analysis
%
  % Following parcellation, each parcellated brain region (node) is decomposed using 
% an unsupervised decomposition approached (a modified Ensemble Empirical Mode Decomposition) 
% resulting in each node (a broadband time series) having a number of modes (a set of narrowband time series) each with 
% a distinct characteristics frequency. The broadband time series are also saved. 

% OUTPUTS: decomposedNodes = cell array with a cell for each node, within
%                            each cell is an array of modes for that node
%          rawData = fMRI data processed and parcellated to X number of nodes (brains regions)
%          The dimension of the raw data saved is number of nodes by number of time points (in the fMRI time series)           
%          rawData_normalized = normalized raw data (necessary for comparison of data coming  from different scanners)
%          percentCensored = percent of frames that were censored due to motion
%          badFrames = frames which exceeded the motion threshold and were censored
%          FD = motion estimates for each time point in mm
%          TR = period between successive time points
%

%Data preprocessing uses the SPM package

spm12path = ''; % spm12 path
addpath(spm12path); % adds spm12 path to matlab session
currentFolder = pwd; % marks current folder (supposed to be the location where all main code is
addpath(currentFolder) % adds current folder containing other codes to matlab path
frontPathT = ''; % path to T1 data (string) where T1 folder containing all subjects is located
frontPathR = ''; % path to RS data (string) where RS folder containing all subjects is located

endPath = '/ses-baselineYear1Arm1'; % path from subject specific folder to func or anat subfolder
%All the ABCD folder organization has the same nomenclature sub-NDARINVXXXXXXXX/ses-baselineYear1Arm1/func/ or /anat/
%For any other studies these paths need to be modified. 
structural_batch_path = ''; % string path to fmri SPM12 processing batch file (i.e., '/pathtofile/filename.mat')
load(structural_batch_path); % loads batch script for RS processing

pathSubStruct = strcat(frontPathT,'/T1/',subID,endPath,'/anat/'); % full path to subject structural data
pathSubFunc = strcat(frontPathR,'/RS/',subID,endPath,'/func/'); % full path to subject functional data

structuralImage = strcat(subID,'_ses-baselineYear1Arm1_run-01_T1w.nii'); % filename for structural T1 image
runstr = sprintf('0%d',n); % str form of run number
filename = strcat(subID,'_ses-baselineYear1Arm1_task-rest_run-',runstr); % filename for RS nifti and motion files
T = tdfread(strcat(pathSubFunc,filename,'_motion.tsv')); % loads realignment parameters
V = niftiread(strcat(pathSubFunc,filename,'_bold.nii')); % loads nifti file for RS data
% Scrub first 16 frames to allow scanner to stabilize
hdrV = niftiinfo(strcat(pathSubFunc,filename,'_bold.nii')); % load nifti header info
TR = hdrV.PixelDimensions(4); % gets TR 

% dimensions of nifti file
dimX = size(V,1);
dimY = size(V,2);
dimZ = size(V,3);
dimT = size(V,4); % number of frames

%% Initial Frame Removal
V(:,:,:,1:frames2remove) = []; % removes initial frames depending on scanner type and software
%For the ABCD data this code follows the recommendation of the ABCD developers/imaging team; frames2remove is chosen by the user

dimT = dimT - frames2remove;
hdrV.ImageSize = [dimX dimY dimZ dimT]; % sets new image size
% creates new vxt without initial frames 
% file created in current folder
newfile = strcat('noinit',filename,'_bold.nii'); % file name for nifti with initial frames removed
niftiwrite(V,newfile,hdrV); % creates file with initial frames removed

% clears variables to make space
clear V; % clears large variables
clear dimX dimY dimZ;
clear hdrV;
% moves file created in working directory to directory where functional
% data is located
movefile(newfile, pathSubFunc); % moves noinit file to subject specific RS folder

%% Calculate Mean Image - necessary prior to co-registration with subject-specific structural MRI
% creates path to each volume of the data with no initial frames
fullfile = strcat(pathSubFunc,newfile);
images = cell(dimT,1);
% file paths must be appended with ,i where i is the volume number in the
% nifti file
for i = 1:dimT
    istr = sprintf(',%d',i);
    images(i,1) = {strcat(fullfile,istr)};
end
% inputs parameters for realignment. Only computes mean image
% .interp is 4th degree B spline interpolation
matlabbatch{1}.spm.spatial.realign.write.data = images;

%% Coregister mean image to structural
% path to mean functional image
sourceImageDir = strcat(pathSubFunc,'mean',newfile,',1');
% path to structural image

refImageDir = strcat(pathSubStruct,structuralImage,',1'); % structural image is reference image
matlabbatch{2}.spm.spatial.coreg.estimate.ref = {refImageDir};
matlabbatch{2}.spm.spatial.coreg.estimate.source = {sourceImageDir};
matlabbatch{2}.spm.spatial.coreg.estimate.other = images;                                                    

%% Slice Time Correction

matlabbatch{3}.spm.temporal.st.scans = {images};
matlabbatch{3}.spm.temporal.st.so = sliceOrder;


%% Normalization of functional images and structural CSF mask
% new file paths to slice time corrected images which begin with filenames that begin with 'a'
fullfile = strcat(pathSubFunc,'a',newfile);
images = cell(dimT,1);
for i = 1:dimT
    istr = sprintf(',%d',i);
    images(i,1) = {strcat(fullfile,istr)};
end
% path to deformation field
deformImagePath = strcat(pathSubStruct,'y_',structuralImage); % path to deformation field
% path to csf mask
matlabbatch{4}.spm.spatial.normalise.write.subj.def = {deformImagePath};
matlabbatch{4}.spm.spatial.normalise.write.subj.resample = images;

%% run matlabbatch
spm('defaults','fmri');
spm_jobman('initcfg');
spm_jobman('run',matlabbatch);
clear matlabbatch

%% Motion regressors + interpolation + bandpass filtering
%These preprocessing steps follow the recommendations of the ABCD developers/imaging team
delete(strcat(pathSubFunc,newfile));
delete(strcat(pathSubFunc,'a',newfile));
delete(strcat(pathSubFunc,'mean',newfile)); % deletes all files created in the above steps besides the normalized file with wa in front
% load normalized and preprocessed brain
V = niftiread(strcat(pathSubFunc,'wa',newfile)); % loads normalized, slice time corrected, initial frames removed file
delete(strcat(pathSubFunc,'wa',newfile)); % deletes file
% delete(strcat(pathSubFunc,'wa',newfile));
% zero NaNs, convert to double and demean
V(isnan(V))=0; % any NaN is set to 0
V = double(V); % converts all values to double
V = V - mean(V,4); % zero-means the data 
dimX = size(V,1);
dimY = size(V,2);
dimZ = size(V,3);

% for ABCD brain

% calculates movement in x,y,z, pitch, yaw and roll
transx = T.trans_x(2:end) - T.trans_x(1:end-1);
transy = T.trans_y(2:end) - T.trans_y(1:end-1);
transz = T.trans_z(2:end) - T.trans_z(1:end-1);
rotx = T.rot_x(2:end) - T.rot_x(1:end-1);
roty = T.rot_y(2:end) - T.rot_y(1:end-1);
rotz = T.rot_z(2:end) - T.rot_z(1:end-1);

% Computes nyquist and sampling frequency
fs = 1/TR;
fN = fs/2;

% constructs 3rd order elliptical bandstop filter for motion filtering
n = 3;
Rp = 0.5;
Rs = 20;
ftype = 'stop';
Wp = [0.28 0.46]/fN; % frequency bounds for motion filtering .28 Hz to 0.46 Hz
[b_stop,a_stop] = ellip(n,Rp,Rs,Wp,ftype);

% filters motion parameters in order to remove effects of breathing on motion
%Filters data on both directions given the non-zero and non-linear phase of the IIR (elliptical) filter
transx = filtfilt(b_stop,a_stop,transx);
transy = filtfilt(b_stop,a_stop,transy);
transz = filtfilt(b_stop,a_stop,transz);
rotx = filtfilt(b_stop,a_stop,rotx);
roty = filtfilt(b_stop,a_stop,roty);
rotz = filtfilt(b_stop,a_stop,rotz);

% Computes frame displacement (FD) for each volume, rotation is transformed to displacement
% by assuming distance from center of brain to cortex is approximately 50 mm
% FD computed according to Power JD, Barnes KA, Snyder AZ, Schlaggar BL, Petersen SE (2012) 
% Spurious but systematic correlations in functional connectivity 
% MRI networks arise from subject motion. Neuroimage 59:2142–2154.
FD = abs(transx)+abs(transy)+abs(transz)+abs(rotx*pi/180*50)+abs(roty*pi/180*50)+abs(rotz*pi/180*50);

% computes derivatives of 6 motion parameters
dtx = (transx(2:end)-transx(1:end-1));
dty = (transy(2:end)-transy(1:end-1));
dtz = (transz(2:end)-transz(1:end-1));
dtrotx = (rotx(2:end)-rotx(1:end-1));
dtrotz = (rotz(2:end)-rotz(1:end-1));
dtroty = (roty(2:end)-roty(1:end-1));

% computes squares of parameters and derivatives
transxSq = transx.^2;
transySq = transy.^2;
transzSq = transz.^2;
rotxSq = rotx.^2;
rotySq = roty.^2;
rotzSq = rotz.^2;
dtxSq = dtx.^2;
dtySq = dty.^2;
dtzSq = dtz.^2;
dtrotxSq = dtrotx.^2;
dtrotySq = dtroty.^2;
dtrotzSq = dtrotz.^2;

% ignores frames corresponding to those initial frames removed in the volume
FD = FD(frames2remove:end);
transx = transx(frames2remove:end);
transy = transy(frames2remove:end);
transz = transz(frames2remove:end);
rotx = rotx(frames2remove:end);
roty = roty(frames2remove:end);
rotz = rotz(frames2remove:end);
transxSq = transxSq(frames2remove:end);
transySq = transySq(frames2remove:end);
transzSq = transzSq(frames2remove:end);
rotxSq = rotxSq(frames2remove:end);
rotySq = rotySq(frames2remove:end);
rotzSq = rotzSq(frames2remove:end);

dtx = dtx(frames2remove-1:end);
dty = dty(frames2remove-1:end);
dtz = dtz(frames2remove-1:end);
dtrotx = dtrotx(frames2remove-1:end);
dtroty = dtroty(frames2remove-1:end);
dtrotz = dtrotz(frames2remove-1:end);
dtxSq = dtxSq(frames2remove-1:end);
dtySq = dtySq(frames2remove-1:end);
dtzSq = dtzSq(frames2remove-1:end);
dtrotxSq = dtrotxSq(frames2remove-1:end);
dtrotySq = dtrotySq(frames2remove-1:end);
dtrotzSq = dtrotzSq(frames2remove-1:end);

% volumes with FDs > 0.3 mm are marked as bad frames
% 0.3 mm threshold is based off of ABCD developers/imaging team suggestion
ts = 1:dimT;
badFrames = ts(FD > 0.3);
goodFrames = ts(FD <= 0.3);
% calculates percentage of bad frames out of total frames
percentCensored = (length(badFrames)/dimT)*100;

if percentCensored < frames_threshold % if run % censored is less than frames threshold, does rest of computation
    % constructs TimePoints x 24 design matrix
    X = [transx transy transz rotx roty rotz transxSq transySq transzSq rotxSq rotySq rotzSq dtx dty dtz dtrotx dtroty dtrotz dtxSq dtySq dtzSq dtrotxSq dtrotySq dtrotzSq];
    % regresses motion parameters R, derivative of motion parameters R',
    % squares of derivatives and motion parameters R^2 and R'^2
    
    % censors frames that have excessive motion
    X(badFrames,:) = 0;
    X = double(X);

    clear transx transy transz rotx roty rotz transxSq transySq transzSq rotxSq rotySq rotzSq dtx dty dtz dtrotx dtroty dtrotz dtxSq dtySq dtzSq dtrotxSq dtrotySq dtrotzSq;
    % clears unneeded space
    zerosTC = zeros(dimT,1); % helps make sure voxel is non zero so we are not wasting computation on empty voxel
    zerosTR = zeros(1,dimT);

    % parses bad frames in order to perform different interpolation schemes
    % makima and spline interpolation perform worse near edges of available
    % points, so instead use nearest neighbor for interpolation near edges

%Interpolation necessary to maintain continuity of signal once bad frames are removed
    last = badFrames(badFrames > (dimT - 4));
    first = badFrames(badFrames < 4);
    mid = badFrames(badFrames >= 4 & badFrames <= (dimT - 4));

    % constructs 3rd order elliptical bandpass filter to filter voxel time
    % series themselves 

    ftype = 'bandpass';
Wp = [0.01 0.25]/fN; % Bounds are chosen based on common practices by community; user can change these bounds as desired
    [b_band,a_band] = ellip(n,Rp,Rs,Wp,ftype);

    % iterates through each voxel, regresses movement parameters, then
    % interpolates censored frames, then performs bandpass filtering
% True paralellization at the spatial dimension (voxel) level
    for i = 1:dimX
        Vi = squeeze(V(i,:,:,:));
        parfor j = 1:dimY
            Vj = squeeze(Vi(j,:,:));
            for k = 1:dimZ
                a = Vj(k,:);
                if ~isequal(a,zerosTC) && ~isequal(a,zerosTR) % makes sure voxel has nonzero timeseries
            % glm regression
                    badVals = a(badFrames);
                    a(badFrames) = 0; % zeros bad timepoints for regression
                    b = glmfitResidOnly(X,a); %fit glm model and retain only the residuals
                    % bandpass filter
                    b(badFrames) = badVals; % adds back bad frames that are later interpolated
                    Vj(k,:) = filtfilt(b_band,a_band,b); % filters time series with band pass filter
                end
            end
            Vi(j,:,:) = Vj;
        end
        V(i,:,:,:) = Vi;
    end

    clear Vi Vj a b
    clear X;

    for i = 1:dimX
        Vg = squeeze(V(i,:,:,goodFrames));
        Vb = zeros(dimY,dimZ,length(badFrames));
        parfor j = 1:dimY
            Vgs = squeeze(Vg(j,:,:));
            for k = 1:dimZ
                if ~isequal(Vgs(k,:),zerosTC) && ~isequal(Vgs(k,:),zerosTR) % checks to make sure we are not performing calculation on empty voxel
                    % Interpolation. 
		    % first and last four frames are interpolated with nearest neighbor (best approach for signal edges) 
                    % all other frames interpolated with makima algorithm
                    Vb(j,k,:) = [interp1(goodFrames,Vgs(k,:),first,'nearest','extrap') interp1(goodFrames,Vgs(k,:),mid,'makima') interp1(goodFrames,Vgs(k,:),last,'nearest','extrap')];
                end
            end
        end
        V(i,:,:,badFrames) = Vb; % bad frames are set to interpolated values
    end
    clear good Vg Vb goodFrames first mid last Vgs

    % clears variables for space

%% Region Parcellation

    V = single(V); % converts to single for parcellation step
    V(isnan(V))=0; 
    numNodes = max(parcellation,[],'all'); % number of nodes or regions in parcellation

    rawData = zeros(numNodes,dimT); % initialized regional time series matrix
    
% parcellates 4D matrix into X regions (by N time points), with X corresponding to 
% the size of the chosen template or the estimated (data-derived) number of regions

    for t = 1:dimT
        Vt = squeeze(V(:,:,:,t));
        parfor reg = 1:numNodes
            rawData(reg,t) = nanmean(Vt(parcellation(:) == reg),'all'); 
        end
    end
    
    rawData = double(rawData);
    clear Vt;
    clear V;
    clear reg;

    rawData_normalized = rawData/nanmedian(nanmedian(abs(rawData),2)); % normalized raw data by median of median amplitude across each region
    % clears matrices to make space
    % creates Matlab structure variable with raw region time series, percent censored
    % frames, nodes with modes by time, and centroid locations for nodes
    decomposedNodes = cell(numNodes,1);

% performs decomposition of node time series into narrowband components (modes)

    ns = 1:numNodes; % vector of nodes

% Depending on the selected template/atlas or parcellation approach, some region time series may be zero 
% These are excluded from the decomposition process 
    goodnodes = ns(any(rawData_normalized,2)); % on non zero regional time series

    meds = nanmedian(abs(rawData_normalized(goodnodes,:)),2);
    prc = prctile(abs(rawData_normalized(goodnodes,:)),[25 75], 2);
    thresh = meds + 3*(prc(:,2)-prc(:,1)); % extreme outlier threshold for amplitude for each regions
    maxthresh = nanmax(thresh); % max outlier threshold across all regions 
    goods = cell(length(goodnodes),1);
%True parallelization at the mode decomposition level for each node
    parfor i = 1:length(goodnodes) %The decomposition step is computationally expensive so it is parallelized
        reg = goodnodes(i);
        raw = rawData_normalized(reg,:);
        
        % Substitutes high-amplitude artifacts by time series statistic
raw(abs(raw) >= maxthresh) = nanmedian(abs(raw(abs(raw) < maxthresh))); 
% any values in regional time series over max (over time series) of the time series-specific extreme outlier threshold is set 
%to median amplitude of region time series
% this is done in order to ease the computation time of the Ensemble
% Empirical Mode Decomposition of the regional time series. Large artifacts substantially increase the computation time of the 
% decomposition
% Here the complete EEMD approach is used
% Torres, M. E., Colominas, M. A., Schlotthauer, G., and Flandrin, P. (2011, May):
% A complete ensemble empirical mode decomposition with adaptive noise
% (ICASSP), 2011 (pp. 4144-4147). DOI: 10.1109/ICASSP.2011.594726

        goods(i,1) = {ceemdan(raw, 0.2, 200, 2000)};
    end

    decomposedNodes(goodnodes) = goods;
else % we are not keeping data for this run so set these fields to zero
    decomposedNodes = {};
    rawData = [];
    rawData_normalized = [];
end
clearvars -except decomposedNodes rawData rawData_normalized percentCensored badFrames FD TR
end
