function [parcellation, ROITable] = parcellation_generator(cortical_atlas)
% Purpose: to generate a parcellation scheme in MNI space by combining specific cortical
%          atlas and Melbourne subcortical atlas
%
% Inputs:
% Cortical_atlas = 1 gives Gordon cortical parcellation
% Cortical_atlas = 2 gives Schaefer 1000 node parcellation
% Cortical_atlas = 3 gives Schaefer 300 node parcellation
% Cortical_atlas = 4 gives Schaefer 800 node parcellation
% Cortical_atlas = 5 gives Schaefer 100 node parcellation
%
% Outputs:
% parcellation = 3d array with each voxel assigned a specific integer value
%                associated with a specific region in the atlas
%                parcellation in MNI space
% ROITable = table of centroid location and name for each region of a particular
%            atlas

parcPath = ''; % path to all parcellations (ie '/pathtoparcellations')

% loads cortex specific parcellation
if cortical_atlas == 1
    cortical_regions = uint16(niftiread(strcat(parcPath,'/Gordon_Parcels_MNI_222.nii')));
    ROIcort = readtable(strcat(parcPath,'/Gordon_MNI_Labels_2mm.csv'));
elseif cortical_atlas == 2
    cortical_regions = uint16(niftiread(strcat(parcPath,'/Schaefer2018_1000Parcels_17Networks_order_FSLMNI152_2mm.nii')));
    ROIcort = readtable(strcat(parcPath,'/Schaefer17Network1000parc_Labels_2mm.csv'));
elseif cortical_atlas == 3
    cortical_regions = uint16(niftiread(strcat(parcPath,'/Schaefer2018_300Parcels_17Networks_order_FSLMNI152_2mm.nii')));
    ROIcort = readtable(strcat(parcPath,'/Schaefer17Network300parc_Labels_2mm.csv'));
elseif cortical_atlas == 4
    cortical_regions = uint16(niftiread(strcat(parcPath,'/Schaefer2018_800Parcels_17Networks_order_FSLMNI152_2mm.nii')));
    ROIcort = readtable(strcat(parcPath,'/Schaefer17Network800parc_Labels_2mm.csv'));
elseif cortical_atlas == 5
    cortical_regions = uint16(niftiread(strcat(parcPath,'/Schaefer2018_100Parcels_17Networks_order_FSLMNI152_2mm.nii')));
    ROIcort = readtable(strcat(parcPath,'/Schaefer17Network100parc_Labels_2mm.csv'));
end

% melbourne subcortical parcellation
subcortical_regions = uint16(niftiread(strcat(parcPath,'/Melbourne_Tian_Subcortex_S4_3T.nii')));
ROIsub = readtable(strcat(parcPath,'/MelbourneParcellationLabels_S4_3T_2mm.csv'));

% cerebellum parcellation
cerebellum_regions = uint16((niftiread(strcat(parcPath,'/Cerebellum-MNIsegment_2mm.nii'))));
ROIcer = readtable(strcat(parcPath,'/Cerebellar_MNI_Segmentation_Labels_2mm.csv'));


regSize = size(cortical_regions);
dimX = regSize(1);
dimY = regSize(2);
dimZ = regSize(3);

numOverlapVox = 0;
parcellation = uint16(zeros(dimX,dimY,dimZ));
maxcort = max(cortical_regions,[],'all');
maxsub = max(subcortical_regions,[],'all');
% maxcereb = max(cerebellum_regions,[],'all');

% checks if there are voxels assigned to multiple regions (ie both
% subcortical areas and cerebellum) and keeps voxel assigned to nothing.
% if no overlap assigns voxel a region number according to the specific
% parcellation it comes from 
% numbering starts with cortex then subcortical regions and then cerebellum 
for i = 1:dimX
    for j = 1:dimY
        for k = 1:dimZ
            assignreg = [cortical_regions(i,j,k) subcortical_regions(i,j,k) cerebellum_regions(i,j,k)];
            if sum(assignreg == 0) < 2 
                numOverlapVox = numOverlapVox + 1;
            else
                if assignreg(1) ~= 0
                    parcellation(i,j,k) = assignreg(1);
                elseif assignreg(2) ~= 0
                    parcellation(i,j,k) = assignreg(2) + maxcort;
                elseif assignreg(3) ~= 0
                    parcellation(i,j,k) = assignreg(3) + maxsub + maxcort;
                end
            end
        end
    end
end

% constructs table with information on each region (name, RAS coordinates)
ROITable = ROIcort;
ROITable = [ROITable; table2cell(ROIsub); table2cell(ROIcer)];
end
