% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_2_computeJointAxes
% This code loads the optimized meshes and computes the joint centers 
% and joint axes
% Bone meshes are checked for impingement and text and visual warnings
% given if impingement is detected

clear

bSave = 1;
bPlot = 1;
bPlot_impingement = 1;
bPlot_pointCloud = 0;

viewAngle = [240 20];
PC_plotSampling = 5;

patientRootDir = '../../';
load([patientRootDir,'sampleData/meshes/data_default_female_landmarks.mat']);
load([patientRootDir,'model/data_opt_meshes.mat']);

if bPlot
   figure ('name','Joint Properties', 'position', [1200 400 600 800]); 
end

display('####################### JOINT PROPERTIES ######################');
display('###############################################################');
display(' '); display(' ');

for i = 1:length(opt_meshes)
    opt_meshes(i).vertices = opt_meshes(i).mesh_vertices_opt(:,1:3);
    opt_meshes(i).faces = opt_meshes(i).mesh_faces_opt;
end

%% First check if there is any bony impingement
lower_limb_names = {'pelvis','femur_r','tibia_r','fibula_r','foot_r','femur_l','tibia_l','fibula_l','foot_l'};
collision_pairs = [1 2; % Pelvis to Femur R
    2 3; % Femur R to Tibia R
    3 4; % Tibia R to Fibula R
    3 5; % Tibia R to Foot R
    4 5; % Fibula R to Foot R
    1 6; % Pelvis to Femur L
    6 7; % Femur L to Tibia L
    7 8; % Tibia L to Fibula L
    7 9; % Tibia L to Foot L
    8 9]; % Fibula R to Foot R
[collision_results] = fnc_checkBoneImpingement(opt_meshes, collision_pairs, lower_limb_names);

%% Compute all joint axes
opt_joint_axes = fnc_computeJointAxes (opt_meshes, data_default_female_landmarks);
opt_jointAxis_pelvis = squeeze(opt_joint_axes(1,:,:));
opt_jointAxis_femur_r = squeeze(opt_joint_axes(2,:,:));
opt_jointAxis_knee_r = squeeze(opt_joint_axes(3,:,:));
opt_jointAxis_talocrural_r = squeeze(opt_joint_axes(4,:,:));
opt_jointAxis_talocalcaneal_r = squeeze(opt_joint_axes(5,:,:));
opt_jointAxis_femur_l = squeeze(opt_joint_axes(6,:,:));
opt_jointAxis_knee_l = squeeze(opt_joint_axes(7,:,:));
opt_jointAxis_talocrural_l = squeeze(opt_joint_axes(8,:,:));
opt_jointAxis_talocalcaneal_l = squeeze(opt_joint_axes(9,:,:));

%% Compute optimized limb lengths
opt_pelvis_len = abs(max(opt_meshes(1).vertices(:,3)) - min(opt_meshes(1).vertices(:,3))); % Top of pelvis to bottom
opt_upperLeg_r_len = sqrt(sum([opt_jointAxis_femur_r(1:3,4) - opt_jointAxis_knee_r(1:3,4)].^2));
opt_upperLeg_l_len = sqrt(sum([opt_jointAxis_femur_l(1:3,4) - opt_jointAxis_knee_l(1:3,4)].^2));
opt_lowerLeg_r_len = sqrt(sum([opt_jointAxis_talocrural_r(1:3,4) - opt_jointAxis_knee_r(1:3,4)].^2));
opt_lowerLeg_l_len = sqrt(sum([opt_jointAxis_talocrural_l(1:3,4) - opt_jointAxis_knee_l(1:3,4)].^2));
opt_foot_r_len = sqrt(sum(([opt_meshes(5).vertices(data_default_female_landmarks(5).footAnteriorPt,:)-opt_meshes(5).vertices(data_default_female_landmarks(5).footPosteriorPt,:)]).^2));
opt_foot_l_len = sqrt(sum(([opt_meshes(9).vertices(data_default_female_landmarks(9).footAnteriorPt,:)-opt_meshes(9).vertices(data_default_female_landmarks(9).footPosteriorPt,:)]).^2));

%% Relative transformation matrices
optRel_pelvis_in_root                   = (opt_jointAxis_pelvis(1:3, 1:3))';
optRel_hipR_in_pelvis                   = (opt_jointAxis_femur_r(1:3, 1:3))'*opt_jointAxis_pelvis(1:3, 1:3);
optRel_kneeR_in_hipR                    = (opt_jointAxis_knee_r(1:3, 1:3))'*opt_jointAxis_femur_r(1:3, 1:3);
optRel_talocruralR_in_kneeR             = (opt_jointAxis_talocrural_r(1:3, 1:3))'*opt_jointAxis_knee_r(1:3, 1:3);
optRel_talocalcanealR_in_talocruralR    = (opt_jointAxis_talocalcaneal_r(1:3, 1:3))'*opt_jointAxis_talocrural_r(1:3, 1:3);
optRel_hipL_in_pelvis                   = (opt_jointAxis_femur_l(1:3, 1:3))'*opt_jointAxis_pelvis(1:3, 1:3);
optRel_kneeL_in_hipL                    = (opt_jointAxis_knee_l(1:3, 1:3))'*opt_jointAxis_femur_l(1:3, 1:3);
optRel_talocruralL_in_kneeL             = (opt_jointAxis_talocrural_l(1:3, 1:3))'*opt_jointAxis_knee_l(1:3, 1:3);
optRel_talocalcanealL_in_talocruralL    = (opt_jointAxis_talocalcaneal_l(1:3, 1:3))'*opt_jointAxis_talocrural_l(1:3, 1:3);

disp('Limb Lengths - Optimized');
disp(['Pelvis = ',num2str(opt_pelvis_len)]);
disp(['Upper Leg (R/L) = ',num2str(opt_upperLeg_r_len),' / ', num2str(opt_upperLeg_l_len)]);
disp(['Lower Leg (R/L) = ',num2str(opt_lowerLeg_r_len),' / ', num2str(opt_lowerLeg_l_len)]);
disp(['Foot (R/L) = ',num2str(opt_foot_r_len),' / ', num2str(opt_foot_l_len)]);

if bSave
    save([patientRootDir,'model/data_opt_joints.mat'],'opt_*','optRel_*');
    display('Saved joint axes in data_opt_joints file');
end
display(' ');

if bPlot
    hold on;
    
    for i=1:length(opt_meshes)
        
        patch ('Faces', opt_meshes(i).faces, 'Vertices', opt_meshes(i).vertices, 'EdgeColor', 'k', 'FaceColor', 'none');
        if i == 5 || i == 9
            fnc_plotJointAxis(squeeze(opt_joint_axes(i,:,:)),'b',5);
        else
            fnc_plotJointAxis(squeeze(opt_joint_axes(i,:,:)),'r',5);
        end
        if bPlot_impingement
            patch('faces',collision_results(i).intSurface.faces,'vertices',collision_results(i).intSurface.vertices,'edgecolor','g','facecolor','none','linewidth',5);
        end
        if bPlot_pointCloud
            plot3(opt_meshes(i).pointCloud_points(1:PC_plotSampling:end,1), opt_meshes(i).pointCloud_points(1:PC_plotSampling:end,2), opt_meshes(i).pointCloud_points(1:PC_plotSampling:end,3),'ob')
        end
    end
    view(viewAngle);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on; axis tight; axis off;
end
