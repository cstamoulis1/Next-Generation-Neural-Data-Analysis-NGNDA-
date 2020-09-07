function construct_connectome_main(subID,cortical_atlas,frames_threshold)
% Main function that calls function construct_connectome_helper and performs preprocessing of fMRI resting-state data
% Also calls function parcellation_generator
% Author: Sean Parks, Computational Neuroscience Laboratory, Boston Childrens Hospital/Harvard Medical School
% Last version: September 5, 2020
% purpose: to process resting state fMRI data for a particular subject.
%          reduced data from over 400k voxels to lower number for regions
%          in a parcellation. Performs ensemble empirical mode
%          decomposition to decompose signal for each region into multiple
%          modes for later noise processing.

% INPUTS: subID = subject identification
%         cortical_atlas = 1 Gordon
%         cortical_atlas = 2 is Schaefer 1000
%         cortical_atlas = 3 is Schaefer 300
%         cortical_atlas = 4 is Schaefer 800
%         cortical_atlas = 5 is Schaefer 100
%         frames_threshold = upper limit of % of frames censored at which a
%                            run of data is no longer satisfactory

% OUTPUTS: saves raw connectome file in the associated folder for chosen parcellation template/atlas



path_to_mriinfo = ''; % string path to directory containing mri info from ABCD (subject id, scanner manufacturer, software)
cd(path_to_mriinfo); %

% must have downloaded abcd_mri01 from NDA and removed headers and saved as .csv file

mriinfo_filename = ''; % filename.csv of abcd_mri01 csv file
M = readtable(mriinfo_filename);

% extracts the subject id column, manufacturer column, and software column
subjects = table2cell(M(2:end,1));
manufacturer = table2cell(M(2:end,5));
software = table2cell(M(2:end,8));

% modifies scanner manufacturer name
for s = 1:length(subjects)
subjects(s) = {strcat('sub-',erase(subjects{s},'_'))}; % turns subject id into sub-NDARINVXXXXXXXX format
scannertype = manufacturer{s};
    if strcmp(scannertype,'SIEMENS')
        manufacturer(s) = {'Siemens'};
    elseif strcmp(scannertype,'GE MEDICAL SYSTEMS')
        manufacturer(s) = {'GE'};
    elseif strcmp(scannertype,'Philips Medical Systems')
        manufacturer(s) = {'Philips'};
    end
end

main_directory = ''; % path to main folder containing code
cd(main_directory); % cd into main code folder
addpath(main_directory); % adds main code folder to matlab path

pathtoSPM12 = ''; % path to folder containing SPM12 code
addpath(pathtoSPM12); % add SPM12 to your matlab session

scantype = manufacturer{find(strcmp(subjects,subID))}; % gets scanner manufacturer for specific subject being processed

pathtoSliceOrder = ''; % string path to slice order file (ie '/pathtofile/sliceTime_ABCD.mat')
load(pathtoSliceOrder); % loads slice time/order for each scanner 

% depending on scanner, slice order and number of initial frames to remove are selected
% during preprocessing
if strcmp(scantype,'Philips')
    sliceOrder = stP;
    frames2remove = 8;
elseif strcmp(scantype,'GE')
    sliceOrder = stG;
    soft = software{find(strcmp(subjects,subID))};
    if contains(soft,'DV25')
        frames2remove = 5;
    elseif contains(soft,'DV26')
        frames2remove = 16;
    end
elseif strcmp(scantype,'Siemens')
    sliceOrder = stS;
    frames2remove = 8;
end

frontPathT = ''; % path to T1 data (string) where T1 folder containing all subjects is located
frontPathR = ''; % path to RS data (string) where RS folder containing all subjects is located

endPath = '/ses-baselineYear1Arm1'; % path from subject specific folder to func or anat subfolder
%For the ABCD dataset, all folders are organized in the same way 
pathSubStruct = strcat(frontPathT,'/T1/',subID,endPath,'/anat/'); % full path to subject structural data
pathSubFunc = strcat(frontPathR,'/RS/',subID,endPath,'/func/'); % full path to subject functional data

structural_batch_path = ''; % string path to structural mri SPM12 processing batch file (ie '/pathtofile/filename.mat')
load(structural_batch_path); % loads batch script for T1 processing

structuralImage = strcat(subID,'_ses-baselineYear1Arm1_run-01_T1w.nii'); % filename of T1 structural image
cd(pathSubStruct); % cd into folder of structural image

% if the warp field file is present removes the warp field file
if isfile(strcat(pathSubStruct,'y_',structuralImage))
delete(strcat(pathSubStruct,'y_',structuralImage));
end

% deletes any .mat files present
if isfile('*mat')
delete *.mat
end

% cd back to main folder
cd(main_directory);

refImageDir = strcat(pathSubStruct,structuralImage,',1'); % sets structural T1 image as reference image
% input structural image directory
matlabbatch{1}.spm.spatial.preproc.channel.vols = {refImageDir};


% Runs structural preprocessing through SPM
% this preprocessing consists of segmentation, which generates a warp field
% for the structural image and transformation from subject-specific space to MNI space

spm('defaults','fmri');

spm_jobman('initcfg');

spm_jobman('run',matlabbatch);

clear matlabbatch; 

% deletes unnecessary tissues created during segmentation step

delete(strcat(pathSubStruct,'c1',structuralImage));
delete(strcat(pathSubStruct,'c2',structuralImage));
delete(strcat(pathSubStruct,'c3',structuralImage));
delete(strcat(pathSubStruct,'c4',structuralImage));
delete(strcat(pathSubStruct,'c5',structuralImage));

% deletes all .mat files created during segmentation step
cd(pathSubStruct); % cd into subject specific structural folder
delete *.mat;

% deletes all .mat files created in any prior SPM processing of the same data 
cd(pathSubFunc); % cd into subject specific functional folder
if isfile('*.mat')
delete *.mat
end

% deletes mean functional image, images with no initial frames, 
% normalized images, slice time corrected images 
% This is necessary so that spm does not crash when trying to create a file that
% already exists and cannot be overwritten

if isfile('mean*')
delete mean*
end

if isfile('noinit*')
delete noinit*
end

if isfile('wanoinit*')
delete wanoinit*
end

if isfile('anoinit*')
delete anoinit*
end

save_location = ''; % string path to location where you want connectomes saved (ie '/pathToSaveLocation/

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

fstruct = dir('sub*.nii'); % gets info on all runs of RS data
numRuns = length(fstruct); % get number of runs of RS data
cd(main_directory); % cds back to current directory containing all code

if numRuns > 0 % if there are any runs of data present
%some brains do not have the first run or consecutive runs
for p = 1:length(fstruct)
    run_vec(p) = str2num(fstruct(p).name(57:58)); %these strings provide the run number
end


[parcellation, ROITable] = parcellation_generator(cortical_atlas); % function generates a parcellation for subject in MNI space, and a corresponding table that describes each region
% in the parcellation

% goes through each run of RS data and produces a connectome for that run
    for n = 1:length(run_vec) 
        [decomposedNodes, rawData, rawData_norm, percentCensored, timePointsCensored, FD, TR] = construct_connectome_helper(subID,run_vec(n),sliceOrder,frames2remove,parcellation,frames_threshold);
        RS_net{run_vec(n)}.decomposedNodes = decomposedNodes;
        RS_net{run_vec(n)}.rawData = rawData;
        RS_net{run_vec(n)}.rawData_normalized = rawData_norm;
        RS_net{run_vec(n)}.percentCensored = percentCensored;
        RS_net{run_vec(n)}.timePointsCensored = timePointsCensored;
        RS_net{run_vec(n)}.FD = FD;
        RS_net{run_vec(n)}.scannerType = scantype;
        RS_net{run_vec(n)}.ROIInfo = ROITable;
        RS_net{run_vec(n)}.TR = TR;
    end
end

% saves generated file in designated location
cd(save_location);
str1 = sprintf(strcat('%s_RS_net-baselineYear1Arm1',identifier,'.mat'),subID); % set filename for connectome
save(str1,'RS_net'); % saves connectome in designated location

delete(strcat(pathSubStruct,'y_',structuralImage)); % delete warped structural image
cd(pathSubFunc);
delete *.mat % delete any .mat file generated in RS processing with SPM
cd(main_directory)
end
