function defaultModel = fnc_getDefaultModel (segmentLength, patient_info, nSegments)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

segmentName   = char('pelvis','thigh_r','shank_r','virtualFoot_r','foot_r','thigh_l','shank_l','virtualFoot_l','foot_l','trunk','head','upperarmR','lowerarmR','handR','upperarmL','lowerarmL','handL');
segmentType   = char('lower_trunk','thigh','shank','virtualFoot','foot','thigh','shank','virtualFoot','foot','upper_trunk','head','upperarm','lowerarm','hand','upperarm','lowerarm','hand');
segmentParent = char('ROOT','pelvis','thigh_r','shank_r','virtualFoot_r','pelvis','thigh_l','shank_l','virtualFoot_l','pelvis','trunk','trunk','upperarmR','lowerarmR','trunk','upperarmL','lowerarmL');

% Regression terms as per Jensen "Body segment mass, radius and radius of gyration proportions of children", J. Biomech 1986. 
% The "torso" ratios from Jensen, have been modified here to account for our two segment torso (pelvis and trunk)
mass_ratio_regressions = [
    -0.0006*0.4  0.4246*0.4;% Pelvis
    0.00364 0.06634;    % Thigh
    0.00122 0.03809;    % Shank
    0.0 0.0;            % Virtual Foot
    0.00015 0.0187;     % Foot
    0.00364 0.06634;    % Thigh
    0.00122 0.03809;    % Shank
    0.0 0.0;            % Virtual Foot
    0.00015 0.0187;     % Foot
    -0.0006*0.6  0.4246*0.6;% Trunk
    -0.0114  0.2376;     % Head
    0.00084 0.022;      % Upperarm
    0.00018 0.01469;    % Lowerarm
    -0.00003 0.00898;    % Hand
    0.00084 0.022;      % Upperarm
    0.00018 0.01469;    % Lowerarm
    -0.00003 0.00898];   % Hand

com_ratio_regressions = [
    0.0018  0.5087;     % Pelvis
    -0.00115 0.4758;     % R Thigh
    -0.003   0.4526;     % R Shank
    0.0 0.0;            % Virtual Foot R
    -0.00186 0.4351;     % R Foot
    -0.00115 0.4758;     % L Thigh
    -0.003   0.4526;     % L Shank
    0.0 0.0;            % Virtual Foot L
    -0.00186 0.4351;     % L Foot
    0.0018  0.5087;     % Trunk
    0.0025  0.4833;     % Head
    -0.00028 0.557;      % R Upperarm
    0.0019  0.561;      % R Lowerarm
    0.00019 0.6005;     % R Hand
    -0.00028 0.557;      % L Upperarm
    0.0019  0.561;      % L Lowerarm
    0.00019 0.6005];    % L Hand

rgyration_ratio_regressions = [
    0.0023 0.7972;     % Pelvis
    -0.00133 0.5536;     % R Thigh
    -0.00224 0.5307;     % R Shank
    0.0 0.0;            % Virtual Foot R
    -0.00203 0.5022;     % R Foot
    -0.00133 0.5536;     % L Thigh
    -0.00224 0.5307;     % L Shank
    0.0 0.0;            % Virtual Foot L
    -0.00203 0.5022;     % L Foot
    0.00120 0.661;     % Trunk
    0.0021  0.5748;     % Head
    -0.00098 0.64;       % R Upperarm
    0.00135 0.633;      % R Lowerarm
    0.00011 0.6461      % R Hand
    -0.00098 0.64;       % L Upperarm
    0.00135 0.633;      % L Lowerarm
    0.00011 0.6461];    % L Hand

Inertia_Mesomorph_xyz = [1.0000    0.5508    0.9864; % Pelvis
                        0.9750    1.0000    0.2021;% R Thigh
                        0.9996    1.0000    0.0865;% R Shank
                        1.0000    1.0000    1.0000;% Virtual Foot R
                        0.8621    1.0000    0.3534; % R Foot
                        0.9750    1.0000    0.2021; % L Thigh
                        0.9996    1.0000    0.0865; % L Shank
                        1.0000    1.0000    1.0000;% Virtual Foot L
                        0.8621    1.0000    0.3534; % L Foot
                        1.0000    0.7508    0.524; % Trunk
                        0.8834    1.0000    0.7771; % Head
                        0.9683    1.0000    0.1177; % R Upperarm
                        1.0000    0.9947    0.0963;% R Lowerarm
                        1.0000    0.8571    0.2000; % R Hand
                        0.9683    1.0000    0.1177; % L Upperarm
                        1.0000    0.9947    0.0963; % L Lowerarm
                        1.0000    0.8571    0.2000];% L Hand
                    
segmentMass      = (mass_ratio_regressions(:,1)     *patient_info.age(1) + mass_ratio_regressions(:,2))*patient_info.weight(1);
segmentCom       = (com_ratio_regressions(:,1)      *patient_info.age(1) + com_ratio_regressions(:,2)).*segmentLength;
rgyration_x      = (rgyration_ratio_regressions(:,1)*patient_info.age(1) + rgyration_ratio_regressions(:,2)).*segmentLength;
segmentInertia_x = segmentMass.*(rgyration_x.^2 - segmentCom.^2);
segmentInertia   = [segmentInertia_x segmentInertia_x segmentInertia_x].*Inertia_Mesomorph_xyz;

% The joints to be consider. first three numbers are translations, the next three rotations
joint_free      = [[zeros(3),eye(3)];[[[0,1,0];[0,0,1];[1,0,0]],zeros(3)]];
joint_spherical = joint_free(4:6,:);
joint_y         = joint_free(4,:);
joint_y_z       = joint_free(4:5,:);
joint_x         = joint_free(6,:);
joint_null      = [];

for i=1:nSegments
    defaultModel(i).name   = strtrim(segmentName(i,:));
    defaultModel(i).type   = strtrim(segmentType(i,:));
    defaultModel(i).parent = strtrim(segmentParent(i,:));
    defaultModel(i).rel_transformation = eye(3);
    defaultModel(i).marker_names       = [];
    defaultModel(i).marker_values      = [];
    defaultModel(i).segmentLength  = segmentLength(i);
    defaultModel(i).mass    = segmentMass(i);
    defaultModel(i).com     = [0, 0,-segmentCom(i)];
    if i == 10 || i == 11 % Trunk and Head
        defaultModel(i).com = [0, 0, 1]*segmentCom(i);
    end
    defaultModel(i).inertia = diag(segmentInertia(i,:));
end

defaultModel(1).rel_joint  = [0.0, 0.0, 0.0];
defaultModel(2).rel_joint  = patient_info.hip_r';
defaultModel(3).rel_joint  = [0.0, 0.0, -segmentLength(2)];
defaultModel(4).rel_joint  = [0.0, 0.0, -segmentLength(3)];
defaultModel(5).rel_joint  = [0.0, 0.0, 0.0];
defaultModel(6).rel_joint  = patient_info.hip_l';
defaultModel(7).rel_joint  = [0.0, 0.0, -segmentLength(6)];
defaultModel(8).rel_joint  = [0.0, 0.0, -segmentLength(7)];
defaultModel(9).rel_joint  = [0.0, 0.0, 0.0];
defaultModel(10).rel_joint = [0.0, 0.0, segmentLength(1)/2];
defaultModel(11).rel_joint = [0.0, 0.0, segmentLength(10)];
defaultModel(12).rel_joint = patient_info.shoulder_r';
defaultModel(13).rel_joint = [0.0, 0.0, -segmentLength(12)];
defaultModel(14).rel_joint = [0.0, 0.0, -segmentLength(13)];
defaultModel(15).rel_joint = patient_info.shoulder_l';
defaultModel(16).rel_joint = [0.0, 0.0, -segmentLength(15)];
defaultModel(17).rel_joint = [0.0, 0.0, -segmentLength(16)];

% Joint Types
defaultModel(1).joint_type  = joint_free;
defaultModel(2).joint_type  = joint_spherical;
defaultModel(3).joint_type  = joint_y_z;
defaultModel(4).joint_type  = joint_y;
defaultModel(5).joint_type  = joint_x;
defaultModel(6).joint_type  = joint_spherical;
defaultModel(7).joint_type  = joint_y_z;
defaultModel(8).joint_type  = joint_y;
defaultModel(9).joint_type  = joint_x;
defaultModel(10).joint_type = joint_spherical;
for i=11:nSegments
    defaultModel(i).joint_type  = joint_null;
end

% Colors
red = [1 0 0]; green = [0.2 0.7 0.3]; blue = [0 0 1]; grey = [0.8 0.8 0.8];
defaultModel(1).mesh_color      = green;     defaultModel(1).mesh_obj     = 'meshes/knubbi_lowertrunk.obj'    ;
defaultModel(2).mesh_color      = blue;      defaultModel(2).mesh_obj     = 'meshes/knubbi_upperleg.obj'    ;
defaultModel(3).mesh_color      = blue;      defaultModel(3).mesh_obj     = 'meshes/knubbi_lowerleg.obj'      ;
defaultModel(4).mesh_color      = [0,0,0];   defaultModel(4).mesh_obj     = 'meshes/knubbi_middletrunk.obj'   ;
defaultModel(5).mesh_color      = blue;      defaultModel(5).mesh_obj     = 'meshes/knubbi_middletrunk.obj'   ;
defaultModel(6).mesh_color      = red;       defaultModel(6).mesh_obj     = 'meshes/knubbi_upperleg.obj'    ;
defaultModel(7).mesh_color      = red;       defaultModel(7).mesh_obj     = 'meshes/knubbi_lowerleg.obj'      ;
defaultModel(8).mesh_color      = [0,0,0];   defaultModel(8).mesh_obj     = 'meshes/knubbi_middletrunk.obj'   ;
defaultModel(9).mesh_color      = red;       defaultModel(9).mesh_obj     = 'meshes/knubbi_middletrunk.obj'   ;
defaultModel(10).mesh_color     = grey;      defaultModel(10).mesh_obj    = 'meshes/knubbi_uppertrunk.obj'    ;
defaultModel(11).mesh_color     = grey;      defaultModel(11).mesh_obj    = 'meshes/knubbi_head.obj'          ;
defaultModel(12).mesh_color     = grey;      defaultModel(12).mesh_obj    = 'meshes/knubbi_upperarm.obj'      ;
defaultModel(13).mesh_color     = grey;      defaultModel(13).mesh_obj    = 'meshes/knubbi_lowerarm.obj'      ;
defaultModel(14).mesh_color     = grey;      defaultModel(14).mesh_obj    = 'meshes/knubbi_head.obj'          ;
defaultModel(15).mesh_color     = grey;      defaultModel(15).mesh_obj    = 'meshes/knubbi_upperarm.obj'      ;
defaultModel(16).mesh_color     = grey;      defaultModel(16).mesh_obj    = 'meshes/knubbi_lowerarm.obj'      ;
defaultModel(17).mesh_color     = grey;      defaultModel(17).mesh_obj    = 'meshes/knubbi_head.obj'          ;

shoulderWidth = abs(patient_info.shoulder_l(2)-patient_info.shoulder_r(2));
defaultModel(1).mesh_dimension  = [0.5*shoulderWidth, 0.8*shoulderWidth, defaultModel(1).segmentLength];      defaultModel(1).mesh_center  = [0, 0, 0];
defaultModel(2).mesh_dimension  = [0.35, 0.35, 1   ]*defaultModel(2).segmentLength;      defaultModel(2).mesh_center  = [0, 0, -1/2]*defaultModel(2).segmentLength;
defaultModel(3).mesh_dimension  = [0.25, 0.25, 1   ]*defaultModel(3).segmentLength;      defaultModel(3).mesh_center  = [0  , 0  , -1/2]*defaultModel(3).segmentLength;
defaultModel(4).mesh_dimension  = [0.01, 0.01, 0.01];                                    defaultModel(4).mesh_center  = [0  , 0  , 0   ];
defaultModel(5).mesh_dimension  = [1, 0.4, 0.15]*defaultModel(5).segmentLength;          defaultModel(5).mesh_center  = [0.3  , 0  , 0]*defaultModel(5).segmentLength;
defaultModel(6).mesh_dimension  = [0.35, 0.35, 1   ]*defaultModel(6).segmentLength;      defaultModel(6).mesh_center  = [0  , 0  , -1/2]*defaultModel(6).segmentLength;
defaultModel(7).mesh_dimension  = [0.25, 0.25, 1   ]*defaultModel(7).segmentLength;      defaultModel(7).mesh_center  = [0  , 0  , -1/2]*defaultModel(7).segmentLength;
defaultModel(8).mesh_dimension  = [0.01, 0.01, 0.01];                                    defaultModel(8).mesh_center  = [0  , 0  , 0   ];
defaultModel(9).mesh_dimension  = [1, 0.4, 0.15]*defaultModel(9).segmentLength;          defaultModel(9).mesh_center  = [0.3, 0  , 0 ]*defaultModel(9).segmentLength;
defaultModel(10).mesh_dimension = [0.5*shoulderWidth,0.8*shoulderWidth,1*defaultModel(10).segmentLength];    defaultModel(10).mesh_center  = [-1/20, 0,  1/2]*defaultModel(10).segmentLength;
defaultModel(11).mesh_dimension = [0.55, 0.65, 1   ]*defaultModel(11).segmentLength;     defaultModel(11).mesh_center = [0  , 0  ,  1/2]*defaultModel(11).segmentLength;
defaultModel(12).mesh_dimension = [0.3 , 0.3 , 1   ]*defaultModel(12).segmentLength;     defaultModel(12).mesh_center = [0  , 0  , -1/2]*defaultModel(12).segmentLength;
defaultModel(13).mesh_dimension = [0.3 , 0.3 , 1   ]*defaultModel(13).segmentLength;     defaultModel(13).mesh_center = [0  , 0  , -1/2]*defaultModel(13).segmentLength;
defaultModel(14).mesh_dimension = [0.4 , 0.7 , 1   ]*defaultModel(14).segmentLength;     defaultModel(14).mesh_center = [0  , 0  , -1/2]*defaultModel(14).segmentLength;
defaultModel(15).mesh_dimension = [0.3 , 0.3 , 1   ]*defaultModel(15).segmentLength;     defaultModel(15).mesh_center = [0  , 0  , -1/2]*defaultModel(15).segmentLength;
defaultModel(16).mesh_dimension = [0.3 , 0.3 , 1   ]*defaultModel(16).segmentLength;     defaultModel(16).mesh_center = [0  , 0  , -1/2]*defaultModel(16).segmentLength;
defaultModel(17).mesh_dimension = [0.4 , 0.7 , 1   ]*defaultModel(17).segmentLength;     defaultModel(17).mesh_center = [0  , 0  , -1/2]*defaultModel(17).segmentLength;