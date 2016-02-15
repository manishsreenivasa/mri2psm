% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_3_segmentBMF
% Load bone, muscle, fat point clouds, and classify them into limb segments using "cutting planes"
% The classified point clouds are used to compute the segment mass, COM and inertia

clear; 

patientRootDir = '../../';
load([patientRootDir,'model/data_bmf.mat']);

bSave = 1;
bPlot = 1;

display('###################### SEGMENT PROPERTIES #####################');
display('###############################################################');
display(' '); display(' ');

% Segment limbs from pointclouds
[pelvis, upperLeg_r, lowerLeg_r, foot_r, upperLeg_l, lowerLeg_l, foot_l] = fnc_extractBMF_IDs (data_bmf, patientRootDir, bPlot);

% Densities in g/mm^3 from Ganley et al. 2004, Gait & Posture 19, 133-140
rho_Bone = 2.5e-3;
rho_Muscle = 1.08e-3;
rho_Fat = 0.9e-3;
rho_bmf= [rho_Bone rho_Muscle rho_Fat];

nSlices = length(data_bmf);
% Mass, center of mass and inertia
[mass_pelvis,     com_pelvis,     inertia_pelvis    ] = fnc_compute_segment_properties(pelvis,     rho_bmf, 100000, 'linear',nSlices);
[mass_upperLeg_r, com_upperLeg_r, inertia_upperLeg_r] = fnc_compute_segment_properties(upperLeg_r, rho_bmf, 100000, 'linear',nSlices);
[mass_upperLeg_l, com_upperLeg_l, inertia_upperLeg_l] = fnc_compute_segment_properties(upperLeg_l, rho_bmf, 100000, 'linear',nSlices);
[mass_lowerLeg_r, com_lowerLeg_r, inertia_lowerLeg_r] = fnc_compute_segment_properties(lowerLeg_r, rho_bmf, 100000, 'linear',nSlices);
[mass_lowerLeg_l, com_lowerLeg_l, inertia_lowerLeg_l] = fnc_compute_segment_properties(lowerLeg_l, rho_bmf, 100000, 'linear',nSlices);
[mass_foot_r,     com_foot_r,     inertia_foot_r ]    = fnc_compute_segment_properties(foot_r,     rho_bmf, 100000, 'linear',nSlices);
[mass_foot_l,     com_foot_l,     inertia_foot_l ]    = fnc_compute_segment_properties(foot_l,     rho_bmf, 100000, 'linear',nSlices);

if bSave
    save([patientRootDir,'model/data_segment_properties.mat'],'mass_*','com_*','inertia_*');
    display('Saved segment properties data_segment_properties file');
end
display(' ');
