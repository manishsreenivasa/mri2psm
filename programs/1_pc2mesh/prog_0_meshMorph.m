% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_0_meshMorph
% This code morphs the lower body meshes to match pointcloud data
% Morphing is done using alpha shape matching and optimal affine transformation using the "metch" toolbox and matlab's optimization toolbox
% You can get the metch tooolbox here:http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?metch

bSave = 1;
bPlot = 1;
patientRootDir = '../../';
viewAngle = [320 10];

% Choose which bones/limbs to process
limbIdx = [1 9];
% Only the GaussNewton method has been implemented for now (based on the metch library)
method = 'gaussnewton';
maxIter = 20;
verbose = 0;

%% LOAD DATA %%
load([patientRootDir,'sampleData/meshes/data_default_female_meshes.mat']);
load([patientRootDir,'sampleData/meshes/data_default_female_landmarks.mat']);
load([patientRootDir,'model/data_bmf.mat']);
lower_limb_names = {'pelvis','femur_r','tibia_r','fibula_r','foot_r','femur_l','tibia_l','fibula_l','foot_l'};

% Extract bone position data from data_bmf structure
bone_pos = fnc_extractBoneFromDataBMF(data_bmf);

if bPlot
     figure ('name','Mesh Morphing', 'position', [0 0 1000 500]);
end

display('######################## MESH MORPHING ########################');
display('###############################################################');
display(' '); display(' ');

for limb = limbIdx(1):limbIdx(2)
    
    display(['Morphing ', char(lower_limb_names(limb)), ' using ', method, ' method']); display(' ');
    
    mesh_vertices       = data_default_female_meshes(limb).vertices;
    mesh_faces          = data_default_female_meshes(limb).faces;
    pointCloud_points   = unique(bone_pos(limb).bonePos, 'rows');
    
    [mesh_vertices_alpha, pointCloud_shp] = fnc_matchAlphaShapes (mesh_vertices, pointCloud_points);
    
    switch lower(method)
        case {'gaussnewton'}
            [mesh_vertices_opt, mesh_faces_opt] = fnc_opt_gaussNewton (mesh_vertices_alpha, mesh_faces, pointCloud_points, maxIter, verbose);
    end
    
    opt_meshes(limb).mesh_vertices_orig  = mesh_vertices;
    opt_meshes(limb).mesh_faces_orig     = mesh_faces;
    opt_meshes(limb).mesh_vertices_opt   = mesh_vertices_opt;
    opt_meshes(limb).mesh_faces_opt      = mesh_faces_opt;
    opt_meshes(limb).mesh_vertices_alpha = mesh_vertices_alpha;
    opt_meshes(limb).mesh_faces_alpha    = mesh_faces;
    opt_meshes(limb).pointCloud_points   = pointCloud_points;
    opt_meshes(limb).pointCloud_shp      = pointCloud_shp;
    
    if bPlot
        PC_plotSampling = 1;
        xyOffset = [100 300 0];
        
        subplot(1,2,1); hold on; axis equal; grid on; view(viewAngle);
        plot3(pointCloud_points(1:PC_plotSampling:end,1)+xyOffset(1), pointCloud_points(1:PC_plotSampling:end,2)+xyOffset(2), pointCloud_points(1:PC_plotSampling:end,3)+xyOffset(3),'ok')
        patch ('Faces', mesh_faces, 'Vertices', mesh_vertices, 'EdgeColor', 'b', 'FaceColor', 'none')
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis tight;
        
        subplot(1,2,2); hold on; axis equal; view(viewAngle); grid on;
        pc_used = pointCloud_points;
        plot3(pc_used(1:PC_plotSampling:end,1), pc_used(1:PC_plotSampling:end,2), pc_used(1:PC_plotSampling:end,3),'ok');
%         patch ('Faces', mesh_faces, 'Vertices', mesh_vertices_alpha, 'EdgeColor', 'b', 'FaceColor', 'none');
        patch ('Faces', mesh_faces_opt, 'Vertices', mesh_vertices_opt, 'EdgeColor', 'r', 'FaceColor', 'none');
        xlabel('X'); ylabel('Y'); zlabel('Z');
        axis tight;
        
        drawnow;
    end
    
    clear pointCloud_* mesh_* final_Affine_T
end

if bSave
    save([patientRootDir,'model/data_opt_meshes.mat'],'opt_meshes');
    display('Saved optimized meshes in data_opt_meshes file');
end
