function [collision_results,collision_pairs] = fnc_checkBoneImpingement(opt_meshes, collision_pairs, lower_limb_names)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_checkBoneImpingement (opt_meshes, collision_pairs, mesh_names)
% Check if mesh pairs in collosion_pairs impinge (collide) with each other and give text warnings
% Collision results can be plotted using info in output variable "collision_results"

% Preallocate memory
collision_results(length(collision_pairs)).intMatrix = [];
collision_results(length(collision_pairs)).intSurface = [];
collision_results(length(collision_pairs)).intIndices = [];

for collPairNo = 1:length(collision_pairs)
    
    mesh1.vertices = opt_meshes(collision_pairs(collPairNo,1)).vertices;
    mesh1.faces = opt_meshes(collision_pairs(collPairNo,1)).faces;
    
    mesh2.vertices = opt_meshes(collision_pairs(collPairNo,2)).vertices;
    mesh2.faces = opt_meshes(collision_pairs(collPairNo,2)).faces;
    
    [intMatrix,intSurface] = SurfaceIntersection(mesh1,mesh2);
    collision_results(collPairNo).intMatrix = intMatrix;
    collision_results(collPairNo).intSurface = intSurface;
    
    if ~isempty(intSurface.vertices)
        [rows,cols] = find(intMatrix==1);
        collision_results(collPairNo).intIndices = [rows cols];
        display(['Impingement detected between ', lower_limb_names{collision_pairs(collPairNo,1)}, ' and ', lower_limb_names{collision_pairs(collPairNo,2)}]);
    end
    clear intMatrix intSurface mesh1_colliding_vertices mesh2_colliding_vertices mesh1 mesh2
end
