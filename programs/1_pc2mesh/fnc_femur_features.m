function femur_features = fnc_femur_features(mesh_vertices, femoralHead, femur_pointD, femur_pointG, femur_pointP, femur_pointLt, femur_pointLc, femur_pointMc)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_femur_features

% Point H = center of femoral head (RFmHC)
femur_features.femur_pointH = mean(mesh_vertices(femoralHead,:));
% midPoint of vector GD
femur_features.femur_pointS = (mesh_vertices(femur_pointG,:)+mesh_vertices(femur_pointD,:))/2;

femur_features.vec_SH = femur_features.femur_pointH - femur_features.femur_pointS;
femur_features.vec_SP = mesh_vertices(femur_pointP,:) - femur_features.femur_pointS;
femur_features.femoral_neckShaft_angle = atan2(norm(cross(femur_features.vec_SH,femur_features.vec_SP)),dot(femur_features.vec_SH,femur_features.vec_SP))*180/pi;

femur_features.vec_SLt = mesh_vertices(femur_pointLt,:)-femur_features.femur_pointS;
femur_features.vecLcMc = mesh_vertices(femur_pointMc,:)-mesh_vertices(femur_pointLc,:);
femur_features.vecLcMc_2d = [femur_features.vecLcMc(1) femur_features.vecLcMc(2) 0.0];
femur_features.vec_SLt_2d = [femur_features.vec_SLt(1) femur_features.vec_SLt(2) 0.0];
femur_features.vec_SH_2d = [femur_features.vec_SH(1) femur_features.vec_SH(2) 0.0];

femur_features.femoral_anteversion_angle = atan2(norm(cross(femur_features.vecLcMc_2d,femur_features.vec_SH_2d)),dot(femur_features.vecLcMc_2d,femur_features.vec_SH_2d))*180/pi;
femur_features.femoral_lesserTrochanterTorsion_angle = atan2(norm(cross(femur_features.vecLcMc_2d,femur_features.vec_SLt_2d)),dot(femur_features.vecLcMc_2d,femur_features.vec_SLt_2d))*180/pi;


end
