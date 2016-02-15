function mesh_vertices = fnc_femur_morphAnteversion(mesh_vertices, femoralEpiphysisHead, femur_features, d_anteversion_angle)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Anteversion Morphing
vec_normal_femoralShaft = femur_features.vec_SP./norm(femur_features.vec_SP);

q1 = cos(d_anteversion_angle/2);
q2 = vec_normal_femoralShaft(1)*sin(d_anteversion_angle/2);
q3 = vec_normal_femoralShaft(2)*sin(d_anteversion_angle/2);
q4 = vec_normal_femoralShaft(3)*sin(d_anteversion_angle/2);

rot_femoralShaft = [2*(q1^2+q2^2)-1 2*(q2*q3-q1*q4) 2*(q2*q4+q1*q3);
                2*(q2*q3+q1*q4) 2*(q1^2+q3^2)-1 2*(q3*q4-q1*q2);
                2*(q2*q4-q1*q3) 2*(q3*q4+q1*q2) 2*(q1^2+q4^2)-1];

for i = 1:length(femoralEpiphysisHead)
    vec_local = mesh_vertices(femoralEpiphysisHead(i),:) - femur_features.femur_pointS;
    mesh_vertices(femoralEpiphysisHead(i),:) = rot_femoralShaft*vec_local';
    mesh_vertices(femoralEpiphysisHead(i),:) = mesh_vertices(femoralEpiphysisHead(i),:)+femur_features.femur_pointS;
end

end
