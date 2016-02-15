% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_4_createLUA
% This file loads the optimized subject-specific joint kinematics, and segment mass properties and writes the corresponding PSM and GSM LUA files

clear; 

bSave = 1;
bPlot = 1;

patientRootDir = '../../';
psmLuaFileName = [patientRootDir,'model/data_samplePSM.lua'];
gsmLuaFileName = [patientRootDir,'model/data_sampleGSM.lua'];

if bPlot
    figure ('name','PSM / GSM Models', 'position', [1200 100 800 800]);
    viewAngle = [240 20];
    markerSz  = 5;
end

display('###################### CREATE LUA MODELS ######################');
display('###############################################################');
display(' '); display(' ');

% Load patient characteristics
load([patientRootDir,'model/data_opt_joints.mat']);
load([patientRootDir,'model/data_segment_properties.mat']);
patient_info = tdfread([patientRootDir,'sampleData/patient_characteristics.txt'],',');
nSegments  = 17;

% For the GSM model you have to provide all segment lengths, and for the
% PSM only the upper body lengths are mandatory
segmentLength = [
    0.15;   % 1 Pelvis - From functional
    0.299;  % 2 Thigh R - From functional
    0.305;  % 3 Shank R - From functional
    0;      % 4 Virtual Foot R
    0.173;  % 5 Foot R - From functional
    0.310;  % 6 Thigh L - From functional
    0.302;  % 7 Shank L - From functional
    0;      % 8 Virtual Foot L
    0.173;  % 9 Foot L - From functional
    0.31;   % 10 Trunk - From functional
    0.18;   % 11 Head
    0.225;  % 12 Upperarm R
    0.195;  % 13 Lowerarm R
    0.1;    % 14 Hand R
    0.225;  % 15 Upperarm L
    0.195;  % 16 Lowerarm L
    0.1];   % 17 Hand L

% Create default GSM model
gsm = fnc_getDefaultModel (segmentLength, patient_info, nSegments);
gsm = fnc_getModelMarkers(gsm, patient_info);

% Create PSM starting with default GSM values
psm = fnc_getDefaultModel (segmentLength, patient_info, nSegments);

%% Here onwards the PSM values are updated based on patient-specific values %%
% Segment relative transformations
psm(1).rel_transformation = optRel_pelvis_in_root;
psm(2).rel_transformation = optRel_hipR_in_pelvis;
psm(3).rel_transformation = optRel_kneeR_in_hipR;
psm(4).rel_transformation = optRel_talocruralR_in_kneeR;
psm(5).rel_transformation = optRel_talocalcanealR_in_talocruralR;
psm(6).rel_transformation = optRel_hipL_in_pelvis;
psm(7).rel_transformation = optRel_kneeL_in_hipL;
psm(8).rel_transformation = optRel_talocruralL_in_kneeL;
psm(9).rel_transformation = optRel_talocalcanealL_in_talocruralL;

% Segment lengths
psm(1).segmentLength  =     opt_pelvis_len/1000;
psm(2).segmentLength  = opt_upperLeg_r_len/1000;
psm(3).segmentLength  = opt_lowerLeg_r_len/1000;
psm(4).segmentLength  =                       0;
psm(5).segmentLength  = opt_foot_r_len    /1000;
psm(6).segmentLength  = opt_upperLeg_l_len/1000;
psm(7).segmentLength  = opt_lowerLeg_l_len/1000;
psm(8).segmentLength  =                       0;
psm(9).segmentLength  = opt_foot_l_len    /1000;

% Get PSM markers with update segment lengths
psm = fnc_getModelMarkers(psm, patient_info);

% Segment masses
psm(1).mass  = mass_pelvis    *10^-3;
psm(2).mass  = mass_upperLeg_r*10^-3;
psm(3).mass  = mass_lowerLeg_r*10^-3;
psm(4).mass  = 0              *10^-3;
psm(5).mass  = mass_foot_r    *10^-3;
psm(6).mass  = mass_upperLeg_l*10^-3;
psm(7).mass  = mass_lowerLeg_l*10^-3;
psm(8).mass  = 0                    ;
psm(9).mass  = mass_foot_l    *10^-3;

% Relative joints
joint_pelvis           = [0 0 0];
joint_upperLeg_r       = (opt_jointAxis_femur_r(1:3,4)         - opt_jointAxis_pelvis(1:3,4))      ;
joint_upperLeg_l       = (opt_jointAxis_femur_l(1:3,4)         - opt_jointAxis_pelvis(1:3,4))      ;
joint_lowerLeg_r       = (opt_jointAxis_knee_r(1:3,4)          - opt_jointAxis_femur_r(1:3,4))     ;
joint_lowerLeg_l       = (opt_jointAxis_knee_l(1:3,4)          - opt_jointAxis_femur_l(1:3,4))     ;
joint_virtualFoot_r    = (opt_jointAxis_talocrural_r(1:3,4)    - opt_jointAxis_knee_r(1:3,4))      ;
joint_virtualFoot_l    = (opt_jointAxis_talocrural_l(1:3,4)    - opt_jointAxis_knee_l(1:3,4))      ;
joint_foot_r           = (opt_jointAxis_talocalcaneal_r(1:3,4) - opt_jointAxis_talocrural_r(1:3,4));
joint_foot_l           = (opt_jointAxis_talocalcaneal_l(1:3,4) - opt_jointAxis_talocrural_l(1:3,4));

psm(1).rel_joint  = joint_pelvis;
psm(2).rel_joint  = (opt_jointAxis_pelvis(1:3, 1:3)'      *joint_upperLeg_r   )'*10^-3;
psm(3).rel_joint  = (opt_jointAxis_femur_r(1:3, 1:3)'     *joint_lowerLeg_r   )'*10^-3;
psm(4).rel_joint  = (opt_jointAxis_knee_r(1:3, 1:3)'      *joint_virtualFoot_r)'*10^-3;
psm(5).rel_joint  = (opt_jointAxis_talocrural_r(1:3, 1:3)'*joint_foot_r       )'*10^-3;
psm(6).rel_joint  = (opt_jointAxis_pelvis(1:3, 1:3)'      *joint_upperLeg_l   )'*10^-3;
psm(7).rel_joint  = (opt_jointAxis_femur_l(1:3, 1:3)'     *joint_lowerLeg_l   )'*10^-3;
psm(8).rel_joint  = (opt_jointAxis_knee_l(1:3, 1:3)'      *joint_virtualFoot_l)'*10^-3;
psm(9).rel_joint  = (opt_jointAxis_talocrural_l(1:3, 1:3)'*joint_foot_l       )'*10^-3;
psm(10).rel_joint = [0.0, 0.0, segmentLength(1)/2];
psm(11).rel_joint = [0.0, 0.0, segmentLength(10)];
psm(12).rel_joint = patient_info.shoulder_r';
psm(13).rel_joint = [0.0, 0.0, -segmentLength(12)];
psm(14).rel_joint = [0.0, 0.0, -segmentLength(13)];
psm(15).rel_joint = patient_info.shoulder_l';
psm(16).rel_joint = [0.0, 0.0, -segmentLength(15)];
psm(17).rel_joint = [0.0, 0.0, -segmentLength(16)];

% Relative COM positions
psm(1).com  = (opt_jointAxis_pelvis(1:3, 1:3)'         *(com_pelvis     - opt_jointAxis_pelvis(1:3,4))         )'*10^-3;
psm(2).com  = (opt_jointAxis_femur_r(1:3, 1:3)'        *(com_upperLeg_r - opt_jointAxis_femur_r(1:3,4))        )'*10^-3;
psm(3).com  = (opt_jointAxis_knee_r(1:3, 1:3)'         *(com_lowerLeg_r - opt_jointAxis_knee_r(1:3,4))         )'*10^-3;
psm(4).com  = [0,0,0];
psm(5).com  = (opt_jointAxis_talocalcaneal_r(1:3, 1:3)'*(com_foot_r     - opt_jointAxis_talocalcaneal_r(1:3,4)))'*10^-3;
psm(6).com  = (opt_jointAxis_femur_l(1:3, 1:3)'        *(com_upperLeg_l - opt_jointAxis_femur_l(1:3,4))        )'*10^-3;
psm(7).com  = (opt_jointAxis_knee_l(1:3, 1:3)'         *(com_lowerLeg_l - opt_jointAxis_knee_l(1:3,4))         )'*10^-3;
psm(8).com  = [0,0,0];
psm(9).com  = (opt_jointAxis_talocalcaneal_l(1:3, 1:3)'*(com_foot_l     - opt_jointAxis_talocalcaneal_l(1:3,4)))'*10^-3;

% Segment inertias
psm(1).inertia  = opt_jointAxis_pelvis(1:3, 1:3)'         *inertia_pelvis    *opt_jointAxis_pelvis(1:3, 1:3)         *10^-9;
psm(2).inertia  = opt_jointAxis_femur_r(1:3, 1:3)'        *inertia_upperLeg_r*opt_jointAxis_femur_r(1:3, 1:3)        *10^-9;
psm(3).inertia  = opt_jointAxis_knee_r(1:3, 1:3)'         *inertia_lowerLeg_r*opt_jointAxis_knee_r(1:3, 1:3)         *10^-9;
psm(4).inertia  = zeros(3);
psm(5).inertia  = opt_jointAxis_talocalcaneal_r(1:3, 1:3)'*inertia_foot_r    *opt_jointAxis_talocalcaneal_r(1:3, 1:3)*10^-9;
psm(6).inertia  = opt_jointAxis_femur_l(1:3, 1:3)'        *inertia_upperLeg_l*opt_jointAxis_femur_l(1:3, 1:3)        *10^-9;
psm(7).inertia  = opt_jointAxis_knee_l(1:3, 1:3)'         *inertia_lowerLeg_l*opt_jointAxis_knee_l(1:3, 1:3)         *10^-9;
psm(8).inertia  = zeros(3);
psm(9).inertia  = opt_jointAxis_talocalcaneal_l(1:3, 1:3)'*inertia_foot_l    *opt_jointAxis_talocalcaneal_l(1:3, 1:3)*10^-9;

% Update Visuals
shoulderWidth = abs(patient_info.shoulder_l(2)-patient_info.shoulder_r(2));
psm(1).mesh_dimension  = [0.5*shoulderWidth, 0.8*shoulderWidth, psm(1).segmentLength];      psm(1).mesh_center  = [0, 0, 0];
psm(2).mesh_dimension  = [0.35, 0.35, 1   ]*psm(2).segmentLength;      psm(2).mesh_center  = [0, 0, -1/2]*psm(2).segmentLength;
psm(3).mesh_dimension  = [0.25, 0.25, 1   ]*psm(3).segmentLength;      psm(3).mesh_center  = [0  , 0  , -1/2]*psm(3).segmentLength;
psm(4).mesh_dimension  = [0.01, 0.01, 0.01];                           psm(4).mesh_center  = [0  , 0  , 0   ];
psm(5).mesh_dimension  = [1, 0.4, 0.15]*psm(5).segmentLength;          psm(5).mesh_center  = [0.3  , 0  , 0]*psm(5).segmentLength;
psm(6).mesh_dimension  = [0.35, 0.35, 1   ]*psm(6).segmentLength;      psm(6).mesh_center  = [0  , 0  , -1/2]*psm(6).segmentLength;
psm(7).mesh_dimension  = [0.25, 0.25, 1   ]*psm(7).segmentLength;      psm(7).mesh_center  = [0  , 0  , -1/2]*psm(7).segmentLength;
psm(8).mesh_dimension  = [0.01, 0.01, 0.01];                           psm(8).mesh_center  = [0  , 0  , 0   ];
psm(9).mesh_dimension  = [1, 0.4, 0.15]*psm(9).segmentLength;          psm(9).mesh_center  = [0.3, 0  , 0 ]*psm(9).segmentLength;
psm(10).mesh_dimension = [0.5*shoulderWidth,0.8*shoulderWidth,1*psm(10).segmentLength];    psm(10).mesh_center  = [-1/20, 0,  1/2]*psm(10).segmentLength;
psm(11).mesh_dimension = [0.55, 0.65, 1   ]*psm(11).segmentLength;     psm(11).mesh_center = [0  , 0  ,  1/2]*psm(11).segmentLength;
psm(12).mesh_dimension = [0.3 , 0.3 , 1   ]*psm(12).segmentLength;     psm(12).mesh_center = [0  , 0  , -1/2]*psm(12).segmentLength;
psm(13).mesh_dimension = [0.3 , 0.3 , 1   ]*psm(13).segmentLength;     psm(13).mesh_center = [0  , 0  , -1/2]*psm(13).segmentLength;
psm(14).mesh_dimension = [0.4 , 0.7 , 1   ]*psm(14).segmentLength;     psm(14).mesh_center = [0  , 0  , -1/2]*psm(14).segmentLength;
psm(15).mesh_dimension = [0.3 , 0.3 , 1   ]*psm(15).segmentLength;     psm(15).mesh_center = [0  , 0  , -1/2]*psm(15).segmentLength;
psm(16).mesh_dimension = [0.3 , 0.3 , 1   ]*psm(16).segmentLength;     psm(16).mesh_center = [0  , 0  , -1/2]*psm(16).segmentLength;
psm(17).mesh_dimension = [0.4 , 0.7 , 1   ]*psm(17).segmentLength;     psm(17).mesh_center = [0  , 0  , -1/2]*psm(17).segmentLength;

gsm_total_weight_in_kg = sum([gsm(:).mass])
psm_total_weight_in_kg = sum([psm(:).mass])
display(' ');

if bPlot
    subplot(1,2,1); hold on; title ('Patient-Specific Model')
    psm(1).global_axes = [psm(1).rel_transformation psm(1).rel_joint'*1000; 0 0 0 1];
    [numMarkers,~] = size(psm(1).marker_names);
    for markerNo = 1:numMarkers
        psm(1).global_marker_position(markerNo,:) = psm(1).global_axes*[psm(1).marker_values(markerNo,:)*1000 1]';
        plot3(psm(1).global_marker_position(:,1), psm(1).global_marker_position(:,2), psm(1).global_marker_position(:,3), 'or');
        text(psm(1).global_marker_position(markerNo,1), psm(1).global_marker_position(markerNo,2), psm(1).global_marker_position(markerNo,3), psm(1).marker_names(markerNo,:));
    end
    fnc_plotJointAxis(psm(1).global_axes,'k');
    for limb = 2:10
        if limb == 6 || limb == 10
            parentNo = 1;
        else
            parentNo = limb-1;
        end
        psm(limb).global_axes = psm(parentNo).global_axes*[inv(psm(limb).rel_transformation) psm(limb).rel_joint'*1000; 0 0 0 1];
        [numMarkers,~] = size(psm(limb).marker_names);
        for markerNo = 1:numMarkers
            psm(limb).global_marker_position(markerNo,:) = psm(limb).global_axes*[psm(limb).marker_values(markerNo,:)*1000 1]';
            plot3(psm(limb).global_marker_position(:,1), psm(limb).global_marker_position(:,2), psm(limb).global_marker_position(:,3), 'or');
            text(psm(limb).global_marker_position(markerNo,1), psm(limb).global_marker_position(markerNo,2), psm(limb).global_marker_position(markerNo,3), psm(limb).marker_names(markerNo,:));
        end
        fnc_plotJointAxis(psm(limb).global_axes,'k');
    end
    view([300 20]); axis square; axis equal;
    
    subplot(1,2,2); hold on; title ('Generic-Scaled Model')
    gsm(1).global_axes = [gsm(1).rel_transformation gsm(1).rel_joint'*1000; 0 0 0 1];
    [numMarkers,~] = size(gsm(1).marker_names);
    for markerNo = 1:numMarkers
        gsm(1).global_marker_position(markerNo,:) = gsm(1).global_axes*[gsm(1).marker_values(markerNo,:)*1000 1]';
        plot3(gsm(1).global_marker_position(:,1), gsm(1).global_marker_position(:,2), gsm(1).global_marker_position(:,3), 'or');
        text(gsm(1).global_marker_position(markerNo,1), gsm(1).global_marker_position(markerNo,2), gsm(1).global_marker_position(markerNo,3), gsm(1).marker_names(markerNo,:));
    end
    fnc_plotJointAxis(gsm(1).global_axes,'k');
    for limb = 2:10
        if limb == 6 || limb == 10
            parentNo = 1;
        else
            parentNo = limb-1;
        end
        gsm(limb).global_axes = gsm(parentNo).global_axes*[inv(gsm(limb).rel_transformation) gsm(limb).rel_joint'*1000; 0 0 0 1];
        [numMarkers,~] = size(gsm(limb).marker_names);
        for markerNo = 1:numMarkers
            gsm(limb).global_marker_position(markerNo,:) = gsm(limb).global_axes*[gsm(limb).marker_values(markerNo,:)*1000 1]';
            plot3(gsm(limb).global_marker_position(:,1), gsm(limb).global_marker_position(:,2), gsm(limb).global_marker_position(:,3), 'or');
            text(gsm(limb).global_marker_position(markerNo,1), gsm(limb).global_marker_position(markerNo,2), gsm(limb).global_marker_position(markerNo,3), gsm(limb).marker_names(markerNo,:));
        end
        fnc_plotJointAxis(gsm(limb).global_axes,'k');
    end
    view([120 20]); axis square; axis equal;
end

if bSave
    res = fnc_writeLuaFile (psmLuaFileName, psm);
    if res == 0 
        display('Saved PSM LUA model file');
    else
        display(['Error saving PSM LUA model file, Errorcode = ', res]);
    end
    
    res = fnc_writeLuaFile (gsmLuaFileName, gsm);
    if res == 0
        display('Saved GSM LUA model file');
    else
        display(['Error saving GSM LUA model file, Errorcode = ', res]);
    end
    display(' ');
end