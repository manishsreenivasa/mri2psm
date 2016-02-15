function res_cost = fnc_compute_morph_cost(x, data_femurMorph)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_compute_morph_cost

femur_features = fnc_femur_features(data_femurMorph.mesh_vertices, data_femurMorph.femoralHead, data_femurMorph.femur_pointD, data_femurMorph.femur_pointG, data_femurMorph.femur_pointP, data_femurMorph.femur_pointLt, data_femurMorph.femur_pointLc, data_femurMorph.femur_pointMc);
d_anteversion_angle_rad = x(1);
mesh_vertices = fnc_femur_morph_anteversion(data_femurMorph.mesh_vertices, data_femurMorph.femoralEpiphysisHead, femur_features, d_anteversion_angle_rad);

femur_features = fnc_femur_features(mesh_vertices, data_femurMorph.femoralHead, data_femurMorph.femur_pointD, data_femurMorph.femur_pointG, data_femurMorph.femur_pointP, data_femurMorph.femur_pointLt, data_femurMorph.femur_pointLc, data_femurMorph.femur_pointMc);
d_neckShaft_angle_rad = x(2);
mesh_vertices = fnc_femur_morph_neckshaft(mesh_vertices, data_femurMorph.femoralEpiphysisHead, femur_features, d_neckShaft_angle_rad);

res_cost = fnc_residual_pt2surf(mesh_vertices, data_femurMorph.mesh_faces, data_femurMorph.pointCloud_points(data_femurMorph.femoralEpiphysisHead_PC,:));

end
