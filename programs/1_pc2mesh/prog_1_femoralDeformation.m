% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_1_femoralDeformation
% This code optimizes the femur mesh (only epihysis head) to match pointcloud data
% Two angles are optimized, anteversion angle, and, neck-shaft angle
% Optimization is done using matlab's optimization toolbox (fmincon -> active-set method)

clear; 

bSave = 1;
bPlot = 1;
PC_plotSampling = 5;
viewAngle = [240 20];

patientRootDir = '../../';
load([patientRootDir,'sampleData/meshes/data_default_female_meshes.mat']);
load([patientRootDir,'sampleData/meshes/data_default_female_landmarks.mat']);
load([patientRootDir,'model/data_opt_meshes.mat']);

x0 = [0.0 0.0];
lb = [-20.0*pi/180 -20.0*pi/180];
ub = [20.0*pi/180 20.0*pi/180];
maxFunEvals = 1000;
maxIter = 200;
options = optimoptions('fmincon', 'Algorithm', 'active-set', 'Display', 'iter', 'MaxFunEvals', maxFunEvals, 'MaxIter', maxIter);

if bPlot
    figure ('name','Femoral Deformation', 'position', [100 500 1000 500]);
end

display('##################### FEMORAL DEFORMATION #####################');
display('###############################################################');
display(' '); display(' ');

for femur_n = 1:2
    if femur_n == 1
        limbNo = 2;
    else
        limbNo = 6;
    end
    
    mesh_vertices = opt_meshes(limbNo).mesh_vertices_opt;
    mesh_faces = opt_meshes(limbNo).mesh_faces_opt;
    pointCloud_points = opt_meshes(limbNo).pointCloud_points;
    lowestFemoralEpiphysisZ =  min(opt_meshes(limbNo).mesh_vertices_opt(data_default_female_landmarks(limbNo).femoralEpiphysisHead,3));
    femoralEpiphysisHead_PC = find(opt_meshes(limbNo).pointCloud_points(:,3) > lowestFemoralEpiphysisZ);
    
    %% Compute initial femur features    
    femoralHead = data_default_female_landmarks(limbNo).femoralHead;
    femur_pointD = data_default_female_landmarks(limbNo).femur_pointD;
    femur_pointG = data_default_female_landmarks(limbNo).femur_pointG;
    femur_pointP = data_default_female_landmarks(limbNo).femur_pointP;
    femur_pointLt = data_default_female_landmarks(limbNo).femur_pointLt; 
    femur_pointLc = data_default_female_landmarks(limbNo).femur_pointLc;
    femur_pointMc = data_default_female_landmarks(limbNo).femur_pointMc;
    femoralEpiphysisHead = data_default_female_landmarks(limbNo).femoralEpiphysisHead;
    femur_features = fnc_femur_features(mesh_vertices, femoralHead, femur_pointD, femur_pointG, femur_pointP, femur_pointLt, femur_pointLc, femur_pointMc);
    
    if bPlot
        %% Plot initial femur meshes
        subplot(1,2,1);
        hold on; axis equal; view(viewAngle); box on; grid on; axis off;
        plot3(pointCloud_points(femoralEpiphysisHead_PC,1), pointCloud_points(femoralEpiphysisHead_PC,2), pointCloud_points(femoralEpiphysisHead_PC,3),'ok')
        patch ('Faces', mesh_faces, 'Vertices', mesh_vertices, 'EdgeColor', 'b', 'FaceColor', 'none')
        axis tight
        
        subplot(1,2,2); hold on;
        axis equal; view([0 90]); box off; grid off; axis off;
        plot3(femur_features.femur_pointH(1), femur_features.femur_pointH(2), femur_features.femur_pointH(3), 'or', 'markersize', 10, 'markerfacecolor', 'r');
        plot3(femur_features.femur_pointS(1), femur_features.femur_pointS(2), femur_features.femur_pointS(3), 'ok', 'markersize', 10, 'markerfacecolor', 'k');
        plot3(mesh_vertices(femur_pointLt,1), mesh_vertices(femur_pointLt,2), mesh_vertices(femur_pointLt,3), 'or', 'markersize', 10, 'markerfacecolor', 'r');
        plot3(mesh_vertices(femur_pointMc,1), mesh_vertices(femur_pointMc,2), mesh_vertices(femur_pointMc,3), 'og', 'markersize', 10, 'markerfacecolor', 'g');
        plot3(mesh_vertices(femur_pointLc,1), mesh_vertices(femur_pointLc,2), mesh_vertices(femur_pointLc,3), 'og', 'markersize', 10, 'markerfacecolor', 'g');
        fnc_plot3LineFromTo(femur_features.femur_pointS, femur_features.femur_pointH,'r', 1);
        fnc_plot3LineFromTo(femur_features.femur_pointS, mesh_vertices(femur_pointLt,1:3),'r', 1);
        axis tight
        drawnow;
    end
    original_mesh_vertices = mesh_vertices;
    original_mesh_faces = mesh_faces;
    data_femurMorph(femur_n).mesh_vertices = mesh_vertices;
    data_femurMorph(femur_n).mesh_faces = mesh_faces;
    data_femurMorph(femur_n).pointCloud_points = pointCloud_points;
    data_femurMorph(femur_n).femoralEpiphysisHead_PC = femoralEpiphysisHead_PC;
    data_femurMorph(femur_n).femoralHead = femoralHead;
    data_femurMorph(femur_n).femur_pointD = femur_pointD;
    data_femurMorph(femur_n).femur_pointG = femur_pointG;
    data_femurMorph(femur_n).femur_pointP = femur_pointP;
    data_femurMorph(femur_n).femur_pointLt = femur_pointLt;
    data_femurMorph(femur_n).femur_pointMc = femur_pointMc;
    data_femurMorph(femur_n).femur_pointLc = femur_pointLc;
    data_femurMorph(femur_n).femoralEpiphysisHead = femoralEpiphysisHead;
    
    %% Optimize femur angles
    objFnc_femur = @(x) fnc_compute_femur_cost(x,data_femurMorph(femur_n));
    
    [x] = fmincon (objFnc_femur, x0, [], [], [], [], lb, ub, [], options);
    
    %% Compute optimized femur features
    femur_features = fnc_femur_features(mesh_vertices, femoralHead, femur_pointD, femur_pointG, femur_pointP, femur_pointLt, femur_pointLc, femur_pointMc);
    d_anteversion_angle_rad = x(1);
    mesh_vertices = fnc_femur_morph_anteversion(mesh_vertices, femoralEpiphysisHead, femur_features, d_anteversion_angle_rad);
    
    femur_features = fnc_femur_features(mesh_vertices, femoralHead, femur_pointD, femur_pointG, femur_pointP, femur_pointLt, femur_pointLc, femur_pointMc);
    d_neckShaft_angle_rad = x(2);
    mesh_vertices = fnc_femur_morph_neckshaft(mesh_vertices, femoralEpiphysisHead, femur_features, d_neckShaft_angle_rad);
    
    femur_features.femoral_anteversion_angle = femur_features.femoral_anteversion_angle + d_anteversion_angle_rad*180/pi;
    femur_features.femoral_neckShaft_angle = femur_features.femoral_neckShaft_angle + d_neckShaft_angle_rad*180/pi;
    display(' ');
    
    opt_meshes(limbNo).femur_features = femur_features;
    opt_meshes(limbNo).mesh_vertices_opt = mesh_vertices;
    opt_meshes(limbNo).mesh_faces_opt = mesh_faces;
    
    if bPlot
        %% Plot optimized femur meshes
        subplot(1,2,1);
        hold on; axis equal; view(viewAngle); box off; grid off; axis off;
        patch ('Faces', mesh_faces, 'Vertices', mesh_vertices(:,1:3), 'EdgeColor', 'g', 'FaceColor', 'none')
        patch ('Faces', original_mesh_faces, 'Vertices', original_mesh_vertices(:,1:3), 'EdgeColor', 'r', 'FaceColor', 'none')
        axis tight
        
        subplot(1,2,2); hold on;
        axis equal; view([0 90]);
        plot3(femur_features.femur_pointH(1), femur_features.femur_pointH(2), femur_features.femur_pointH(3), 'or', 'markersize', 10, 'markerfacecolor', 'r');
        plot3(femur_features.femur_pointS(1), femur_features.femur_pointS(2), femur_features.femur_pointS(3), 'ok', 'markersize', 10, 'markerfacecolor', 'k');
        plot3(mesh_vertices(femur_pointLt,1), mesh_vertices(femur_pointLt,2), mesh_vertices(femur_pointLt,3), 'or', 'markersize', 10, 'markerfacecolor', 'r');
        plot3(mesh_vertices(femur_pointMc,1), mesh_vertices(femur_pointMc,2), mesh_vertices(femur_pointMc,3), 'og', 'markersize', 10, 'markerfacecolor', 'g');
        plot3(mesh_vertices(femur_pointLc,1), mesh_vertices(femur_pointLc,2), mesh_vertices(femur_pointLc,3), 'og', 'markersize', 10, 'markerfacecolor', 'g');
        text(femur_features.femur_pointH(1)+5, femur_features.femur_pointH(2), femur_features.femur_pointH(3), 'H');
        text(femur_features.femur_pointS(1)+5, femur_features.femur_pointS(2), femur_features.femur_pointS(3), 'S');
        text(mesh_vertices(femur_pointLt,1)+5, mesh_vertices(femur_pointLt,2), mesh_vertices(femur_pointLt,3), 'Lt');
        text(mesh_vertices(femur_pointMc,1)+5, mesh_vertices(femur_pointMc,2), mesh_vertices(femur_pointMc,3), 'Mc');
        text(mesh_vertices(femur_pointLc,1)+5, mesh_vertices(femur_pointLc,2), mesh_vertices(femur_pointLc,3), 'Lc');
        fnc_plot3LineFromTo(femur_features.femur_pointS, femur_features.femur_pointH,'r', 2);
        fnc_plot3LineFromTo(femur_features.femur_pointS, mesh_vertices(femur_pointLt,1:3),'r', 2);
        axis tight
    end
    clear data_femurMorph femur_features mesh_vertices mesh_faces original_mesh_faces
    clear original_mesh_vertices d_anteversion_angle_rad d_neckShaft_angle_rad
    clear femoralEpiphysisHead_PC x fval exitflag output Obj_fnc
end

if bSave
    save([patientRootDir,'model/data_opt_meshes.mat'],'opt_meshes');
    display('Saved femoral deformation data in data_opt_meshes file');
end
