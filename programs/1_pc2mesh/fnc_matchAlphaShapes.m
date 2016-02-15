function [mesh_vertices_alpha, pointCloud_shp] = fnc_matchAlphaShapes (mesh_vertices, pointCloud_points, alphaRadius)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

if nargin < 3, alphaRadius = 50; end

[pointCloud_com, ~, pointCloud_shp] = fnc_computeAlphaProperties(pointCloud_points, alphaRadius);
[mesh_com, ~, ~]                    = fnc_computeAlphaProperties(mesh_vertices, alphaRadius);

pointCloud_cent_pos = (pointCloud_points-ones(size(pointCloud_points))*diag(pointCloud_com));
mesh_cent_pos       = (mesh_vertices-ones(size(mesh_vertices))*diag(mesh_com));

[~, pointCloud_In, ~] = fnc_computeAlphaProperties(pointCloud_cent_pos, alphaRadius);
[~, mesh_In, ~]       = fnc_computeAlphaProperties(mesh_cent_pos, alphaRadius);

pointCloud_In_norm = pointCloud_In/norm(pointCloud_In);
mesh_In_norm       = mesh_In/norm(mesh_In);

x_curr = [0 0 0];
options = optimoptions('fminunc','Algorithm','quasi-newton','Display','none');
Obj_fnc = @(x_curr) fnc_opt_In (x_curr, pointCloud_In_norm, mesh_In_norm);
[x_curr] = fminunc (Obj_fnc, x_curr, options);
rx = [1 0 0; 0 cos(x_curr(1)) -sin(x_curr(1)); 0 sin(x_curr(1)) cos(x_curr(1))];
ry = [cos(x_curr(2)) 0 sin(x_curr(2)); 0 1 0; -sin(x_curr(2)) 0 cos(x_curr(2))];
rz = [cos(x_curr(3)) -sin(x_curr(3)) 0; sin(x_curr(3)) cos(x_curr(3)) 0; 0 0 1];

mesh_components     = mesh_cent_pos*rx*ry*rz;
pointCloud_min      = min(pointCloud_points);
mesh_min            = min(mesh_components);
pointCloud_dist     = max(pointCloud_points) - pointCloud_min;
mesh_dist           = max(mesh_components) - mesh_min;
mesh_vertices_alpha = repmat(pointCloud_min,length(mesh_components),1) + (mesh_components-repmat(mesh_min,length(mesh_components),1)).*repmat(pointCloud_dist./mesh_dist,length(mesh_components),1);


%-----------------------------------------------------------------------------------------------------
function res = fnc_opt_In (x_curr, pointCloud_In_norm, mesh_In_norm)

rx = [1 0 0; 0 cos(x_curr(1)) -sin(x_curr(1)); 0 sin(x_curr(1)) cos(x_curr(1))];
ry = [cos(x_curr(2)) 0 sin(x_curr(2)); 0 1 0; -sin(x_curr(2)) 0 cos(x_curr(2))];
rz = [cos(x_curr(3)) -sin(x_curr(3)) 0; sin(x_curr(3)) cos(x_curr(3)) 0; 0 0 1];

mesh_mod = (rx*ry*rz)'*mesh_In_norm*rx*ry*rz;
res = norm(pointCloud_In_norm-mesh_mod);
