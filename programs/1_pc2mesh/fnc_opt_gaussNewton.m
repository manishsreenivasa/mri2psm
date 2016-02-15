function [mesh_vertices_opt, mesh_faces_opt] = fnc_opt_gaussNewton (mesh_vertices, mesh_faces, pc_vertices, maxiter, verbose, C_ub, C_lb)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Uses a modified "regpt2surf" function, originally developed as part of the "metch" toolbox
% 

if nargin < 7, C_lb = [0.8 -0.4 -0.4 -0.4 0.8 -0.4 -0.4 -0.4 0.8 -50.0 -50.0 -50.0]'; end
if nargin < 6, C_ub = [1.2  0.4  0.4  0.4 1.2  0.4  0.4  0.4 1.2  50.0  50.0  150.0]'; end
if nargin < 5, verbose = 0; end
if nargin < 4, maxiter = 50; end

if verbose
    options = optimoptions('fmincon','Algorithm','active-set','Display','iter','MaxFunEvals', 1000, 'MaxIter', maxiter);
else
    options = optimoptions('fmincon','Algorithm','active-set','Display','none','MaxFunEvals', 1000, 'MaxIter', maxiter);
end
mesh_vertices_opt = mesh_vertices;
mesh_faces_opt = mesh_faces;

% Homotopy Iterations
C0 = [1 0 0 0 1 0 0 0 1 0 0 0]'; cmask = ones(12,1); pmask = -1*ones(size(pc_vertices,1),1);
[C, trivial_cost, ~] = fnc_opt_regpt2surf (mesh_vertices_opt, mesh_faces_opt, pc_vertices, 1, 1, zeros(size(pc_vertices,1),1), pmask, C0, C_ub, C_lb, cmask, verbose);
for lambda = 0.1:0.1:1.0
    display(['Homotopy iteration with lambda = ', num2str(lambda)]);
    [C, err, flag_constraint_saturation] = fnc_opt_regpt2surf (mesh_vertices_opt, mesh_faces_opt, pc_vertices, maxiter, lambda, trivial_cost, pmask, C, C_ub, C_lb,  cmask, verbose);
    display(['Fval = ', num2str(sum(err)),', X Current = ', mat2str(C,2)]);
    if flag_constraint_saturation
        display('!!! Constraint Saturation !!!');
        break;
    end
end
   
A = reshape(C(1:9),3,3);
b = C(end-2:end);
inv_AffineT = [[inv(A); 0 0 0] [-inv(A)*b;1]];
for i = 1:length(mesh_vertices_opt)
    mesh_vertices_tmp(i,:) = inv_AffineT*[mesh_vertices_opt(i,:) 1]';
end
mesh_vertices_opt = mesh_vertices_tmp(:,1:3);
mesh_faces_opt = mesh_faces;
