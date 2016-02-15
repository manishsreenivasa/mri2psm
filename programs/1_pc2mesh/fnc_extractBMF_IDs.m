function [pelvis, upperLeg_r, lowerLeg_r, foot_r, upperLeg_l, lowerLeg_l, foot_l] = fnc_extractBMF_IDs (data_bmf, patientRootDir, bPlot)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_extract_BMF_IDs

load([patientRootDir,'model/data_opt_meshes.mat']);
load([patientRootDir,'model/data_opt_joints.mat']);
load([patientRootDir,'sampleData/meshes/data_default_female_landmarks.mat']);

if bPlot
    h = figure ('name','Compute Segment Properties', 'position', [100 800 1200 800]);
    viewAngle = [300 20];
    markerSz  = 5;
end

% Extract bone position data from data_bmf structure
bone = fnc_extractBoneFromDataBMF (data_bmf);

% Define Pelvic Planes
pelvis_LandmarkPos = opt_meshes(1).mesh_vertices_opt(data_default_female_landmarks(1).LandmarkIndices,:);
pelvis_RPASIS = pelvis_LandmarkPos(1,:);
pelvis_RPPSIS =  pelvis_LandmarkPos(2,:);
pelvis_LPASIS = pelvis_LandmarkPos(4,:);
pelvis_LPPSIS =  pelvis_LandmarkPos(5,:);
femur_r_pointS = opt_meshes(2).femur_features.femur_pointS;
femur_l_pointS = opt_meshes(6).femur_features.femur_pointS;

vec_LPASIS_LPPSIS = pelvis_LPASIS-pelvis_LPPSIS;
vec_LPASIS_LPPSIS = vec_LPASIS_LPPSIS./norm(vec_LPASIS_LPPSIS);
vec_LPASIS_LpointH = pelvis_LPASIS-femur_r_pointS;
vec_LPASIS_LpointH = vec_LPASIS_LpointH./norm(vec_LPASIS_LpointH);
normal_pelvisLplane = cross(vec_LPASIS_LPPSIS,vec_LPASIS_LpointH);
normal_pelvisLplane = normal_pelvisLplane./norm(normal_pelvisLplane);
d_pelvisLplane = -sum(dot(normal_pelvisLplane,femur_l_pointS));

vec_RPASIS_RPPSIS = pelvis_RPASIS-pelvis_RPPSIS;
vec_RPASIS_RPPSIS = vec_RPASIS_RPPSIS./norm(vec_RPASIS_RPPSIS);
vec_RPASIS_RpointH = pelvis_RPASIS-femur_l_pointS;
vec_RPASIS_RpointH = vec_RPASIS_RpointH./norm(vec_RPASIS_RpointH);
normal_pelvisRplane = cross(vec_RPASIS_RPPSIS,vec_RPASIS_RpointH);
normal_pelvisRplane = normal_pelvisRplane./norm(normal_pelvisRplane);
d_pelvisRplane = -sum(dot(normal_pelvisRplane,femur_r_pointS));

% Define Knee Limits
knee_r_Center = opt_jointAxis_knee_r(1:3,4);
knee_r_Center(3) = (knee_r_Center(3) + opt_meshes(3).mesh_vertices_opt(data_default_female_landmarks(3).LandmarkIndices(1),3) + opt_meshes(3).mesh_vertices_opt(data_default_female_landmarks(3).LandmarkIndices(2),3))/3;
knee_l_Center = opt_jointAxis_knee_l(1:3,4);
knee_l_Center(3) = (knee_l_Center(3) + opt_meshes(7).mesh_vertices_opt(data_default_female_landmarks(7).LandmarkIndices(1),3) + opt_meshes(7).mesh_vertices_opt(data_default_female_landmarks(7).LandmarkIndices(2),3))/3;

% Define Ankle Limits
ankle_l_Center = opt_jointAxis_talocrural_l(1:3,4);
ankle_r_Center = opt_jointAxis_talocrural_r(1:3,4);

% Preallocate IDs
nSlices = length(data_bmf);
boneIDs(nSlices).pelvis = [];       muscleIDs(nSlices).pelvis = [];       fatIDs(nSlices).upperLeg_l = [];
boneIDs(nSlices).upperLeg_l = [];   muscleIDs(nSlices).upperLeg_l = [];   fatIDs(nSlices).upperLeg_l = [];
boneIDs(nSlices).upperLeg_r = [];   muscleIDs(nSlices).upperLeg_r = [];   fatIDs(nSlices).upperLeg_r = [];
boneIDs(nSlices).lowerLeg_l = [];   muscleIDs(nSlices).lowerLeg_l = [];   fatIDs(nSlices).lowerLeg_l = [];
boneIDs(nSlices).lowerLeg_r = [];   muscleIDs(nSlices).lowerLeg_r = [];   fatIDs(nSlices).lowerLeg_r = [];
boneIDs(nSlices).foot_l = [];       muscleIDs(nSlices).foot_l = [];       fatIDs(nSlices).foot_l = [];
boneIDs(nSlices).foot_r = [];       muscleIDs(nSlices).foot_r = [];       fatIDs(nSlices).foot_r = [];

if bPlot
    subplot(1,3,1); hold on; grid on;
    set(gca,'fontsize',15); title('BONE');
    xlabel('X');ylabel('Y');zlabel('Z');
    view(viewAngle);
    
    subplot(1,3,2); hold on; grid on;
    set(gca,'fontsize',15); title('MUSCLE');
    xlabel('X');ylabel('Y');zlabel('Z');
    view(viewAngle);
    
    subplot(1,3,3); hold on; grid on;
    set(gca,'fontsize',15); title('FAT');
    xlabel('X');ylabel('Y');zlabel('Z');
    view(viewAngle);
end

for slice = 1:nSlices
    
    display(['Processing Slice ', int2str(slice)]);
    
    bonePos     = data_bmf(slice).bonePos;
    musclePos   = data_bmf(slice).musclePos;
    fatPos      = data_bmf(slice).fatPos;
    sliceZ = nanmean(musclePos(:,3));
    pixelArea = data_bmf(slice).pixelSize(1)*data_bmf(slice).pixelSize(2);
    
    % Bone IDs
    bone_dist_to_pelvisLplane = 0;
    bone_dist_to_pelvisRplane = 0;
    for pt = 1:length(bonePos)
        bone_dist_to_pelvisLplane(pt) = fnc_dist2plane(bonePos(pt,:), normal_pelvisLplane, d_pelvisLplane);
        bone_dist_to_pelvisRplane(pt) = fnc_dist2plane(bonePos(pt,:), normal_pelvisRplane, d_pelvisRplane);
    end
    if ~isempty(bonePos)
        boneIDs(slice).pelvis = find(bone_dist_to_pelvisLplane > 0 & bone_dist_to_pelvisRplane < 0);
        boneIDs(slice).upperLeg_l = find(bone_dist_to_pelvisLplane < 0 & bonePos(:,3)' > knee_l_Center(3) & bonePos(:,2)' < 0);
        boneIDs(slice).upperLeg_r = find(bone_dist_to_pelvisRplane > 0 & bonePos(:,3)' > knee_r_Center(3) & bonePos(:,2)' > 0);
        boneIDs(slice).lowerLeg_l = find(bonePos(:,3)' > ankle_l_Center(3) & bonePos(:,3)' < knee_l_Center(3) & bonePos(:,2)' < 0);
        boneIDs(slice).lowerLeg_r = find(bonePos(:,3)' > ankle_r_Center(3) & bonePos(:,3)' < knee_r_Center(3) & bonePos(:,2)' > 0);
        boneIDs(slice).foot_l = find(bonePos(:,3)' < ankle_l_Center(3) & bonePos(:,2)' < 0);
        boneIDs(slice).foot_r = find(bonePos(:,3)' < ankle_r_Center(3) & bonePos(:,2)' > 0);
    end
    
    % Muscle IDs
    mus_dist_to_pelvisLplane = 0;
    mus_dist_to_pelvisRplane = 0;
    for pt = 1:length(musclePos)
        mus_dist_to_pelvisLplane(pt) = fnc_dist2plane(musclePos(pt,:), normal_pelvisLplane, d_pelvisLplane);
        mus_dist_to_pelvisRplane(pt) = fnc_dist2plane(musclePos(pt,:), normal_pelvisRplane, d_pelvisRplane);
    end
    muscleIDs(slice).pelvis = find(mus_dist_to_pelvisLplane > 0 & mus_dist_to_pelvisRplane < 0);
    muscleIDs(slice).upperLeg_l = find(mus_dist_to_pelvisLplane < 0 & musclePos(:,3)' > knee_l_Center(3) & musclePos(:,2)' < 0);
    muscleIDs(slice).upperLeg_r = find(mus_dist_to_pelvisRplane > 0 & musclePos(:,3)' > knee_r_Center(3) & musclePos(:,2)' > 0);
    muscleIDs(slice).lowerLeg_l = find(musclePos(:,3)' > ankle_l_Center(3) & musclePos(:,3)' < knee_l_Center(3) & musclePos(:,2)' < 0);
    muscleIDs(slice).lowerLeg_r = find(musclePos(:,3)' > ankle_r_Center(3) & musclePos(:,3)' < knee_r_Center(3) & musclePos(:,2)' > 0);
    muscleIDs(slice).foot_l = find(musclePos(:,3)' < ankle_l_Center(3) & musclePos(:,2)' < 0);
    muscleIDs(slice).foot_r = find(musclePos(:,3)' < ankle_r_Center(3) & musclePos(:,2)' > 0);
    
    %% %%%%%%
    % Fat IDs
    fat_dist_to_pelvisLplane = 0;
    fat_dist_to_pelvisRplane = 0;
    for pt = 1:length(fatPos)
        fat_dist_to_pelvisLplane(pt) = fnc_dist2plane(fatPos(pt,:), normal_pelvisLplane, d_pelvisLplane);
        fat_dist_to_pelvisRplane(pt) = fnc_dist2plane(fatPos(pt,:), normal_pelvisRplane, d_pelvisRplane);
    end
    fatIDs(slice).pelvis = find(fat_dist_to_pelvisLplane > 0 & fat_dist_to_pelvisRplane < 0);
    fatIDs(slice).upperLeg_l = find(fat_dist_to_pelvisLplane < 0 & fatPos(:,3)' > knee_l_Center(3) & fatPos(:,2)' < 0);
    fatIDs(slice).upperLeg_r = find(fat_dist_to_pelvisRplane > 0 & fatPos(:,3)' > knee_r_Center(3) & fatPos(:,2)' > 0);
    fatIDs(slice).lowerLeg_l = find(fatPos(:,3)' > ankle_l_Center(3) & fatPos(:,3)' < knee_l_Center(3) & fatPos(:,2)' < 0);
    fatIDs(slice).lowerLeg_r = find(fatPos(:,3)' > ankle_r_Center(3) & fatPos(:,3)' < knee_r_Center(3) & fatPos(:,2)' > 0);
    fatIDs(slice).foot_l = find(fatPos(:,3)' < ankle_l_Center(3) & fatPos(:,2)' < 0);
    fatIDs(slice).foot_r = find(fatPos(:,3)' < ankle_r_Center(3) & fatPos(:,2)' > 0);
    
    %% Segment into limbs
    pelvis(slice).bonePos        = bonePos(boneIDs(slice).pelvis,:);
    pelvis(slice).fatPos         = fatPos(fatIDs(slice).pelvis,:);
    pelvis(slice).musclePos      = musclePos(muscleIDs(slice).pelvis,:);
    pelvis(slice).boneArea       = length(boneIDs(slice).pelvis)*pixelArea;
    pelvis(slice).muscleArea     = length(muscleIDs(slice).pelvis)*pixelArea;
    pelvis(slice).fatArea        = length(fatIDs(slice).pelvis)*pixelArea;
    pelvis(slice).sliceZ         = sliceZ;
    pelvis(slice).pixelArea      = pixelArea;
    
    upperLeg_r(slice).bonePos    = bonePos(boneIDs(slice).upperLeg_r,:);
    upperLeg_r(slice).fatPos     = fatPos(fatIDs(slice).upperLeg_r,:);
    upperLeg_r(slice).musclePos  = musclePos(muscleIDs(slice).upperLeg_r,:);
    upperLeg_r(slice).boneArea   = length(boneIDs(slice).upperLeg_r)*pixelArea;
    upperLeg_r(slice).muscleArea = length(muscleIDs(slice).upperLeg_r)*pixelArea;
    upperLeg_r(slice).fatArea    = length(fatIDs(slice).upperLeg_r)*pixelArea;
    upperLeg_r(slice).sliceZ     = sliceZ;
    upperLeg_r(slice).pixelArea  = pixelArea;
    
    lowerLeg_r(slice).bonePos    = bonePos(boneIDs(slice).lowerLeg_r,:);
    lowerLeg_r(slice).fatPos     = fatPos(fatIDs(slice).lowerLeg_r,:);
    lowerLeg_r(slice).musclePos  = musclePos(muscleIDs(slice).lowerLeg_r,:);
    lowerLeg_r(slice).boneArea   = length(boneIDs(slice).lowerLeg_r)*pixelArea;
    lowerLeg_r(slice).muscleArea = length(muscleIDs(slice).lowerLeg_r)*pixelArea;
    lowerLeg_r(slice).fatArea    = length(fatIDs(slice).lowerLeg_r)*pixelArea;
    lowerLeg_r(slice).sliceZ     = sliceZ;
    lowerLeg_r(slice).pixelArea  = pixelArea;
    
    foot_r(slice).bonePos        = bonePos(boneIDs(slice).foot_r,:);
    foot_r(slice).fatPos         = fatPos(fatIDs(slice).foot_r,:);
    foot_r(slice).musclePos      = musclePos(muscleIDs(slice).foot_r,:);
    foot_r(slice).boneArea       = length(boneIDs(slice).foot_r)*pixelArea;
    foot_r(slice).muscleArea     = length(muscleIDs(slice).foot_r)*pixelArea;
    foot_r(slice).fatArea        = length(fatIDs(slice).foot_r)*pixelArea;
    foot_r(slice).sliceZ         = sliceZ;
    foot_r(slice).pixelArea      = pixelArea;
    
    upperLeg_l(slice).bonePos    = bonePos(boneIDs(slice).upperLeg_l,:);
    upperLeg_l(slice).fatPos     = fatPos(fatIDs(slice).upperLeg_l,:);
    upperLeg_l(slice).musclePos  = musclePos(muscleIDs(slice).upperLeg_l,:);
    upperLeg_l(slice).boneArea   = length(boneIDs(slice).upperLeg_l)*pixelArea;
    upperLeg_l(slice).muscleArea = length(muscleIDs(slice).upperLeg_l)*pixelArea;
    upperLeg_l(slice).fatArea    = length(fatIDs(slice).upperLeg_l)*pixelArea;
    upperLeg_l(slice).sliceZ     = sliceZ;
    upperLeg_l(slice).pixelArea  = pixelArea;
    
    lowerLeg_l(slice).bonePos    = bonePos(boneIDs(slice).lowerLeg_l,:);
    lowerLeg_l(slice).fatPos     = fatPos(fatIDs(slice).lowerLeg_l,:);
    lowerLeg_l(slice).musclePos  = musclePos(muscleIDs(slice).lowerLeg_l,:);
    lowerLeg_l(slice).boneArea   = length(boneIDs(slice).lowerLeg_l)*pixelArea;
    lowerLeg_l(slice).muscleArea = length(muscleIDs(slice).lowerLeg_l)*pixelArea;
    lowerLeg_l(slice).fatArea    = length(fatIDs(slice).lowerLeg_l)*pixelArea;
    lowerLeg_l(slice).sliceZ     = sliceZ;
    lowerLeg_l(slice).pixelArea  = pixelArea;
    
    foot_l(slice).bonePos        = bonePos(boneIDs(slice).foot_l,:);
    foot_l(slice).fatPos         = fatPos(fatIDs(slice).foot_l,:);
    foot_l(slice).musclePos      = musclePos(muscleIDs(slice).foot_l,:);
    foot_l(slice).boneArea       = length(boneIDs(slice).foot_l)*pixelArea;
    foot_l(slice).muscleArea     = length(muscleIDs(slice).foot_l)*pixelArea;
    foot_l(slice).fatArea        = length(fatIDs(slice).foot_l)*pixelArea;
    foot_l(slice).sliceZ         = sliceZ;
    foot_l(slice).pixelArea      = pixelArea;
    
    if bPlot
        subplot(1,3,1);
        if ~isempty(bonePos)
            plot3(bonePos(boneIDs(slice).pelvis,1), bonePos(boneIDs(slice).pelvis ,2), repmat(sliceZ,length(boneIDs(slice).pelvis),1), 'ok','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).upperLeg_l,1), bonePos( boneIDs(slice).upperLeg_l,2), repmat(sliceZ,length(boneIDs(slice).upperLeg_l),1), 'om','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).upperLeg_r,1), bonePos( boneIDs(slice).upperLeg_r,2), repmat(sliceZ,length(boneIDs(slice).upperLeg_r),1), 'or','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).lowerLeg_l,1), bonePos(boneIDs(slice).lowerLeg_l,2), repmat(sliceZ,length(boneIDs(slice).lowerLeg_l),1), 'og','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).lowerLeg_r,1), bonePos(boneIDs(slice).lowerLeg_r,2), repmat(sliceZ,length(boneIDs(slice).lowerLeg_r),1), 'ob','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).foot_l,1), bonePos(boneIDs(slice).foot_l,2), repmat(sliceZ,length(boneIDs(slice).foot_l),1), 'oc','markersize', markerSz);
            plot3(bonePos(boneIDs(slice).foot_r,1), bonePos(boneIDs(slice).foot_r,2), repmat(sliceZ,length(boneIDs(slice).foot_r),1), 'oy','markersize', markerSz);
        end
        axis equal
        
        subplot(1,3,2);
        plot3(musclePos(muscleIDs(slice).pelvis,1), musclePos(muscleIDs(slice).pelvis ,2), repmat(sliceZ,length(muscleIDs(slice).pelvis),1), 'ok','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).upperLeg_l,1), musclePos(muscleIDs(slice).upperLeg_l,2), repmat(sliceZ,length(muscleIDs(slice).upperLeg_l),1), 'om','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).upperLeg_r,1), musclePos(muscleIDs(slice).upperLeg_r,2), repmat(sliceZ,length(muscleIDs(slice).upperLeg_r),1), 'or','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).lowerLeg_l,1), musclePos(muscleIDs(slice).lowerLeg_l,2), repmat(sliceZ,length(muscleIDs(slice).lowerLeg_l),1), 'og','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).lowerLeg_r,1), musclePos(muscleIDs(slice).lowerLeg_r,2), repmat(sliceZ,length(muscleIDs(slice).lowerLeg_r),1), 'ob','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).foot_l,1), musclePos(muscleIDs(slice).foot_l,2), repmat(sliceZ,length(muscleIDs(slice).foot_l),1), 'oc','markersize', markerSz);
        plot3(musclePos(muscleIDs(slice).foot_r,1), musclePos(muscleIDs(slice).foot_r,2), repmat(sliceZ,length(muscleIDs(slice).foot_r),1), 'oy','markersize', markerSz);
        axis equal
        
        subplot(1,3,3);
        plot3(fatPos(fatIDs(slice).pelvis,1), fatPos(fatIDs(slice).pelvis ,2), repmat(sliceZ,length(fatIDs(slice).pelvis),1), 'ok','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).upperLeg_l,1), fatPos(fatIDs(slice).upperLeg_l,2), repmat(sliceZ,length(fatIDs(slice).upperLeg_l),1), 'om','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).upperLeg_r,1), fatPos(fatIDs(slice).upperLeg_r,2), repmat(sliceZ,length(fatIDs(slice).upperLeg_r),1), 'or','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).lowerLeg_l,1), fatPos(fatIDs(slice).lowerLeg_l,2), repmat(sliceZ,length(fatIDs(slice).lowerLeg_l),1), 'og','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).lowerLeg_r,1), fatPos(fatIDs(slice).lowerLeg_r,2), repmat(sliceZ,length(fatIDs(slice).lowerLeg_r),1), 'ob','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).foot_l,1), fatPos(fatIDs(slice).foot_l,2), repmat(sliceZ,length(fatIDs(slice).foot_l),1), 'oc','markersize', markerSz);
        plot3(fatPos(fatIDs(slice).foot_r,1), fatPos(fatIDs(slice).foot_r,2), repmat(sliceZ,length(fatIDs(slice).foot_r),1), 'oy','markersize', markerSz);
        axis equal
        
        drawnow;
    end
    clear bone_dist_to_* mus_dist_to_* fat_dist_to_* bonePos musclePos fatPos
end
