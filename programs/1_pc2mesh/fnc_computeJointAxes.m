function joint_axes = fnc_computeJointAxes (meshes, landmarks)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_computeJointAxes

% Pelvis Center 0
jointCenter_pelvis_z = mean([meshes(1).vertices(landmarks(1).LandmarkIndices([1:2 4:5]),:)]);
jointCenter_pelvis_xy = mean([meshes(1).vertices(landmarks(1).LandmarkIndices([1 4]),:)]);
jointCenter_pelvis = [jointCenter_pelvis_xy(1:2) jointCenter_pelvis_z(3)];

% Pelvic positive Y axis = LPASI - RPASI;
vec_OY = meshes(1).vertices(landmarks(1).LandmarkIndices(4),:)-...
    meshes(1).vertices(landmarks(1).LandmarkIndices(1),:);
vec_OY = vec_OY./norm(vec_OY);
% Pelvic ZY axis = LPASI - LPPSIS;
vec_in_pelvicPlane = meshes(1).vertices(landmarks(1).LandmarkIndices(4),:)-...
    meshes(1).vertices(landmarks(1).LandmarkIndices(5),:);
vec_in_pelvicPlane = vec_in_pelvicPlane./norm(vec_in_pelvicPlane);
% Pelvic positive Z axis
vec_OZ = cross(vec_in_pelvicPlane,vec_OY);
vec_OZ = vec_OZ./norm(vec_OZ);
% Pelvic positive X axis
vec_OX = cross(vec_OY,vec_OZ);
vec_OX = vec_OX./norm(vec_OX);
% Column-major rotation matrix from unit vectors
joint_axes(1,:,:) = [[vec_OX';0] [vec_OY';0] [vec_OZ';0] [jointCenter_pelvis';1]];

jointCenter_r_hip = mean(meshes(2).vertices(landmarks(2).femoralHead,:));
femur_pointCd_r = mean(meshes(2).vertices(landmarks(2).femoralCondyles,:));
% Femur positive Z axis = center of condyles - hip center
vec_OZ = jointCenter_r_hip - femur_pointCd_r;
vec_OZ = vec_OZ./norm(vec_OZ);
% Femur positive Y axis
vec_OY = jointCenter_r_hip-meshes(2).vertices(landmarks(2).femur_pointG,:);
vec_OY = vec_OY./norm(vec_OY);
% Femur positive X axis
vec_OX = cross(vec_OY,vec_OZ);
vec_OX = vec_OX./norm(vec_OX);
vec_OY = cross(vec_OZ,vec_OX); % To orthogonalize the frame
% Column-major rotation matrix from unit vectors
joint_axes(2,:,:) = [[vec_OX';0] [vec_OY';0] [vec_OZ';0] [jointCenter_r_hip';1]];

% Ankle has to come first because of the way the knee axis is constructed
% Ankle Center - Axis as per Isman & Inman 1969, joint center assumed as mid Pt of distal ends of tibia and fibula
% talocrural_r_pointLateral = meshes(3).vertices(landmarks(3).LandmarkIndices(4),:);
talocrural_r_pointLateral = meshes(4).vertices(landmarks(4).LandmarkIndices(4),:);
talocrural_r_pointMedial = meshes(3).vertices(landmarks(3).LandmarkIndices(3),:);
jointCenter_r_talocrural = mean([talocrural_r_pointLateral;talocrural_r_pointMedial]);

% Talocrural Positive Y axis
vec_OY_talocrural = talocrural_r_pointMedial-talocrural_r_pointLateral;
vec_OY_talocrural = vec_OY_talocrural./norm(vec_OY_talocrural);
vec_in_Talocrural_XY_plane = jointCenter_r_talocrural - mean([meshes(5).vertices(landmarks(5).LandmarkIndices(4),:);meshes(5).vertices(landmarks(5).LandmarkIndices(5),:)]);
vec_in_Talocrural_XY_plane = vec_in_Talocrural_XY_plane./norm(vec_in_Talocrural_XY_plane);
% Talocrural Positive Z axis
vec_OZ_talocrural = cross(vec_OY_talocrural,vec_in_Talocrural_XY_plane);
vec_OZ_talocrural = vec_OZ_talocrural./norm(vec_OZ_talocrural);
vec_OZ_talocrural_right = vec_OZ_talocrural;
% Talocrural Positive X axis
vec_OX_talocrural = cross(vec_OY_talocrural,vec_OZ_talocrural);
vec_OX_talocrural = vec_OX_talocrural./norm(vec_OX_talocrural);

% Talocalcaneal Positive X axis
vec_OX_talocalcaneal = mean([meshes(5).vertices(landmarks(5).LandmarkIndices(4),:);meshes(5).vertices(landmarks(5).LandmarkIndices(5),:)])-jointCenter_r_talocrural;
vec_OX_talocalcaneal = vec_OX_talocalcaneal./norm(vec_OX_talocalcaneal);
vec_in_Talocalcaneal_XY_plane = jointCenter_r_talocrural-meshes(5).vertices(landmarks(5).LandmarkIndices(4),:);
vec_in_Talocalcaneal_XY_plane = vec_in_Talocalcaneal_XY_plane./norm(vec_in_Talocalcaneal_XY_plane);
% Talocalcaneal Positive Z axis
vec_OZ_talocalcaneal = cross(vec_in_Talocalcaneal_XY_plane,vec_OX_talocalcaneal);
vec_OZ_talocalcaneal = vec_OZ_talocalcaneal./norm(vec_OZ_talocalcaneal);
% Talocalcaneal Positive Y axis
vec_OY_talocalcaneal = cross(vec_OZ_talocalcaneal,vec_OX_talocalcaneal);
vec_OY_talocalcaneal = vec_OY_talocalcaneal./norm(vec_OY_talocalcaneal);

% Column-major rotation matrix from unit vectors
joint_axes(4,:,:) = [[vec_OX_talocrural';0] [vec_OY_talocrural';0] [vec_OZ_talocrural';0] [jointCenter_r_talocrural';1]];
joint_axes(5,:,:) = [[vec_OX_talocalcaneal';0] [vec_OY_talocalcaneal';0] [vec_OZ_talocalcaneal';0] [jointCenter_r_talocrural';1]];

% Femur Right
jointCenter_r_knee = femur_pointCd_r;
% Tibial positive Z axis
vec_OZ = jointCenter_r_knee-jointCenter_r_talocrural;
vec_OZ = vec_OZ./norm(vec_OZ);
vec_OZ_tibia_right = vec_OZ;
% Tibial positive Y axis - From Knee landmarks
vec_OY = meshes(2).vertices(landmarks(2).femur_pointMc,:) - meshes(2).vertices(landmarks(2).femur_pointLc,:);
vec_OY = vec_OY./norm(vec_OY);
% Tibial positive X axis
vec_OX = cross(vec_OY,vec_OZ);
vec_OX = vec_OX./norm(vec_OX);
vec_OY = cross(vec_OZ,vec_OX); % To orthogonalize the frame
% Column-major rotation matrix from unit vectors
joint_axes(3,:,:) = [[vec_OX';0] [vec_OY';0] [vec_OZ';0] [jointCenter_r_knee';1]];

% Hip Center
jointCenter_l_hip = mean(meshes(6).vertices(landmarks(6).femoralHead,:));
% Femur positive Z axis = center of condyles - hip center
femur_pointCd_l = mean(meshes(6).vertices(landmarks(6).femoralCondyles,:));
vec_OZ = jointCenter_l_hip - femur_pointCd_l;
vec_OZ = vec_OZ./norm(vec_OZ);
% Femur positive Y axis
vec_OY = meshes(6).vertices(landmarks(6).femur_pointG,:)-jointCenter_l_hip;
vec_OY = vec_OY./norm(vec_OY);
% Femur positive X axis
vec_OX = cross(vec_OY,vec_OZ);
vec_OX = vec_OX./norm(vec_OX);
vec_OY = cross(vec_OZ,vec_OX); % To orthogonalize the frame
% Column-major rotation matrix from unit vectors
joint_axes(6,:,:) = [[vec_OX';0] [vec_OY';0] [vec_OZ';0] [jointCenter_l_hip';1]];

% Ankle has to come first because of the way the knee axis is constructed
% Ankle Center - Axis as per Isman & Inman 1969, joint center assumed as mid Pt of distal ends of tibia and fibula
% talocrural_l_pointLateral = meshes(7).vertices(landmarks(7).LandmarkIndices(4),:);
talocrural_l_pointLateral = meshes(8).vertices(landmarks(8).LandmarkIndices(4),:);
talocrural_l_pointMedial = meshes(7).vertices(landmarks(7).LandmarkIndices(3),:);
jointCenter_l_talocrural = mean([talocrural_l_pointLateral;talocrural_l_pointMedial]);

% Talocrural Positive Y axis
vec_OY_talocrural = talocrural_l_pointLateral-talocrural_l_pointMedial;
vec_OY_talocrural = vec_OY_talocrural./norm(vec_OY_talocrural);
vec_in_Talocrural_XY_plane = jointCenter_l_talocrural - mean([meshes(9).vertices(landmarks(9).LandmarkIndices(4),:);meshes(9).vertices(landmarks(9).LandmarkIndices(5),:)]);
vec_in_Talocrural_XY_plane = vec_in_Talocrural_XY_plane./norm(vec_in_Talocrural_XY_plane);
% Talocrural Positive Z axis
vec_OZ_talocrural = cross(vec_OY_talocrural,vec_in_Talocrural_XY_plane);
vec_OZ_talocrural = vec_OZ_talocrural./norm(vec_OZ_talocrural);
vec_OZ_talocrural_left = vec_OZ_talocrural;
% Talocrural Positive X axis
vec_OX_talocrural = cross(vec_OY_talocrural,vec_OZ_talocrural);
vec_OX_talocrural = vec_OX_talocrural./norm(vec_OX_talocrural);

% Talocalcaneal Positive X axis
vec_OX_talocalcaneal = mean([meshes(9).vertices(landmarks(9).LandmarkIndices(4),:);meshes(9).vertices(landmarks(9).LandmarkIndices(5),:)])-jointCenter_l_talocrural;
vec_OX_talocalcaneal = vec_OX_talocalcaneal./norm(vec_OX_talocalcaneal);
vec_in_Talocalcaneal_XY_plane = jointCenter_l_talocrural-meshes(9).vertices(landmarks(9).LandmarkIndices(4),:);
vec_in_Talocalcaneal_XY_plane = vec_in_Talocalcaneal_XY_plane./norm(vec_in_Talocalcaneal_XY_plane);
% Talocalcaneal Positive Z axis
vec_OZ_talocalcaneal = cross(vec_OX_talocalcaneal,vec_in_Talocalcaneal_XY_plane);
vec_OZ_talocalcaneal = vec_OZ_talocalcaneal./norm(vec_OZ_talocalcaneal);
% Talocalcaneal Positive Y axis
vec_OY_talocalcaneal = cross(vec_OZ_talocalcaneal,vec_OX_talocalcaneal);
vec_OY_talocalcaneal = vec_OY_talocalcaneal./norm(vec_OY_talocalcaneal);

% Column-major rotation matrix from unit vectors
joint_axes(8,:,:) = [[vec_OX_talocrural';0] [vec_OY_talocrural';0] [vec_OZ_talocrural';0] [jointCenter_l_talocrural';1]];
joint_axes(9,:,:) = [[vec_OX_talocalcaneal';0] [vec_OY_talocalcaneal';0] [vec_OZ_talocalcaneal';0] [jointCenter_l_talocrural';1]];


% Knee Left
jointCenter_l_knee = femur_pointCd_l;
% Tibial positive Z axis
vec_OZ = jointCenter_l_knee-jointCenter_l_talocrural;
vec_OZ = vec_OZ./norm(vec_OZ);
% Tibial positive Y axis - From Knee landmarks
vec_OY = meshes(6).vertices(landmarks(6).femur_pointLc,:) - meshes(6).vertices(landmarks(6).femur_pointMc,:);
vec_OY = vec_OY./norm(vec_OY);
vec_OZ_tibia_left = vec_OZ;
% Tibial positive X axis
vec_OX = cross(vec_OY,vec_OZ);
vec_OX = vec_OX./norm(vec_OX);
vec_OY = cross(vec_OZ,vec_OX); % To orthogonalize the frame
% Column-major rotation matrix from unit vectors
joint_axes(7,:,:) = [[vec_OX';0] [vec_OY';0] [vec_OZ';0] [jointCenter_l_knee';1]];

right_ankle_plantarflexion_offset = acos(vec_OZ_tibia_right*vec_OZ_talocrural_right')
joint_axes(4,1:3,1:3) = squeeze(joint_axes(4,1:3,1:3))*[cos(right_ankle_plantarflexion_offset) 0 -sin(right_ankle_plantarflexion_offset); 0 1 0; sin(right_ankle_plantarflexion_offset) 0 cos(right_ankle_plantarflexion_offset)];

left_ankle_plantarflexion_offset = acos(vec_OZ_tibia_left*vec_OZ_talocrural_left')
joint_axes(8,1:3,1:3) = squeeze(joint_axes(8,1:3,1:3))*[cos(left_ankle_plantarflexion_offset) 0 -sin(left_ankle_plantarflexion_offset); 0 1 0; sin(left_ankle_plantarflexion_offset) 0 cos(left_ankle_plantarflexion_offset)];


