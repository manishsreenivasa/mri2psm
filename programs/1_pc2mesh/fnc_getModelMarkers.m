function model = fnc_getModelMarkers (model, patient_info)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% RKNE_md and RANK_md are virtual markers
% LKNE_md and LANK_md are virtual markers
model(1).marker_names = char('LASI','SACR','RASI');
model(2).marker_names = char('RTHI','RKNE','RKNE_md');
model(3).marker_names = char('RTIB' ,'RANK','RANK_md');

model(4).marker_names = [];
model(5).marker_names = char('RANK','RANK_md','RTOE','RHEE');
model(6).marker_names = char('LTHI','LKNE','LKNE_md');
model(7).marker_names = char('LTIB','LANK','LANK_md');
model(8).marker_names = [];
model(9).marker_names = char('LANK','LANK_md','LTOE','LHEE');
model(10).marker_names = char('C7','LSHO','CLAV','RSHO');

% Pelvis
model(1).marker_values(1,:) = [0.0 , patient_info.interAsisDistance(1)/2  ,  0.0]; % LASI
model(1).marker_values(2,:) = [-patient_info.interAsisDistance(1)*1/2 , 0.0 , 0.02]; % SACR
model(1).marker_values(3,:) = [0.0 , -patient_info.interAsisDistance(1)/2 , 0.0]; % RASI

% Thigh R
model(2).marker_values(1,:) = [0.02 , -patient_info.kneewidth_r(1)*2/3 , -model(2).segmentLength*1/2]; % RTHI
model(2).marker_values(2,:) = [0.0 , -patient_info.kneewidth_r(1)/2 , -model(2).segmentLength]; % RKNE
model(2).marker_values(3,:) = [0.0 ,  patient_info.kneewidth_r(1)/2 , -model(2).segmentLength]; % RKNE_md

% Shank R
model(3).marker_values(1,:) = [0.0 , -(0.05+patient_info.kneewidth_r(1)/2) , -model(3).segmentLength*1/2]; % RTIB 
model(3).marker_values(2,:) = [0.0 , -patient_info.anklewidth_r(1)/2 , -model(3).segmentLength]; % RANK
model(3).marker_values(3,:) = [0.0 ,  patient_info.anklewidth_r(1)/2 , -model(3).segmentLength]; % RANK_md

% Virtual Foot R
model(4).marker_values = [ ];

% Foot R
model(5).marker_values(1,:) = [0.0 , -patient_info.anklewidth_r(1)/2 , 0.0]; % RANK
model(5).marker_values(2,:) = [0.0 , patient_info.anklewidth_r(1)/2 , 0.0]; % RANK_md
model(5).marker_values(3,:) = [0.5*model(5).segmentLength , 0.0 , 0.0]; % RTOE
model(5).marker_values(4,:) = [-0.3*model(5).segmentLength , 0.0 , 0.0]; % RHEE

% Thigh L
model(6).marker_values(1,:) = [0.02 , patient_info.kneewidth_l(1)*2/3 , -model(6).segmentLength*1/2]; % LTHI
model(6).marker_values(2,:) = [0.0 , patient_info.kneewidth_l(1)/2 , -model(6).segmentLength]; % LKNE
model(6).marker_values(3,:) = [0.0 , -patient_info.kneewidth_l(1)/2 , -model(6).segmentLength]; % LKNE_md

% Shank L
model(7).marker_values(1,:) = [0.0 , 0.05+patient_info.kneewidth_l(1)/2 , -model(7).segmentLength*1/2]; % LTIB 
model(7).marker_values(2,:) = [0.0 , patient_info.anklewidth_l(1)/2 , -model(7).segmentLength]; % LANK
model(7).marker_values(3,:) = [0.0 , -patient_info.anklewidth_l(1)/2 , -model(7).segmentLength]; % LANK_md

% Virtual Foot L
model(8).marker_values = [ ];

% Foot L
model(9).marker_values(1,:) = [0.0 , patient_info.anklewidth_l(1)/2 , 0.0]; % LANK
model(9).marker_values(2,:) = [0.0 , -patient_info.anklewidth_l(1)/2 , 0.0]; % LANK_md
model(9).marker_values(3,:) = [0.5*model(9).segmentLength , 0.0 , 0.0]; % LTOE
model(9).marker_values(4,:) = [-0.3*model(9).segmentLength , 0.0 , 0.0]; % LHEE

% Trunk
model(10).marker_values(1,:) = [-patient_info.shoulder_l(2)/3 , 0.0 , 1.15*patient_info.shoulder_l(3)]; % C7
model(10).marker_values(2,:) = patient_info.shoulder_l';
model(10).marker_values(3,:) = [patient_info.shoulder_l(2)/3  ,  0.0  ,  0.9*patient_info.shoulder_l(3)]; % CLAV
model(10).marker_values(4,:) = patient_info.shoulder_r';
