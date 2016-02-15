function mesh_vertices = fnc_femur_morphNeckShaft(mesh_vertices, femoralEpiphysisHead, femur_features, d_neckShaft_angle)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Neck Shaft Morphing
vec_normal_neckShaft = cross(femur_features.vec_SP, femur_features.vec_SH)./norm(cross(femur_features.vec_SP, femur_features.vec_SH));

q1 = cos(d_neckShaft_angle/2);
q2 = vec_normal_neckShaft(1)*sin(d_neckShaft_angle/2);
q3 = vec_normal_neckShaft(2)*sin(d_neckShaft_angle/2);
q4 = vec_normal_neckShaft(3)*sin(d_neckShaft_angle/2);

rot_neckShaft = [2*(q1^2+q2^2)-1 2*(q2*q3-q1*q4) 2*(q2*q4+q1*q3);
                2*(q2*q3+q1*q4) 2*(q1^2+q3^2)-1 2*(q3*q4-q1*q2);
                2*(q2*q4-q1*q3) 2*(q3*q4+q1*q2) 2*(q1^2+q4^2)-1];

for i = 1:length(femoralEpiphysisHead)
    vec_local = mesh_vertices(femoralEpiphysisHead(i),:) - femur_features.femur_pointS;
    mesh_vertices(femoralEpiphysisHead(i),:) = rot_neckShaft*vec_local';
    mesh_vertices(femoralEpiphysisHead(i),:) = mesh_vertices(femoralEpiphysisHead(i),:)+femur_features.femur_pointS;
end

