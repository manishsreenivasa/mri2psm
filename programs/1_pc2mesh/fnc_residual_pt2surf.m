function res = fnc_residual_pt2surf (mesh_vertices, mesh_faces, pc_vertices)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

nv = nodesurfnorm (mesh_vertices, mesh_faces);
[dist0, ~] = dist2surf (mesh_vertices, nv, pc_vertices);
res = sum (dist0);
