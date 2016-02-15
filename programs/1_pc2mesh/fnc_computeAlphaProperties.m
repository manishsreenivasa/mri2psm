function [com, inertia, shp] = fnc_computeAlphaProperties (points, rad, fig)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_computeAlphaProperties
% Compute alpha shape, and COM and inertia of alpha shape

if nargin < 2 || isempty(rad), rad = inf; end
if nargin < 3, fig = 0; end

[~, shp] = alphavol(points, rad);
shp.vertices = points;
shp.faces = shp.bnd;

[com, inertia] = get_com_inertia(shp.tri, points);

if fig
   subplot(1,2,2); hold on;
   trisurf(shp.bnd, points(:,1), points(:,2), points(:,3), 'FaceColor','none','FaceAlpha',1);
   axis equal;
   view([210 20]);   
end

%--------------------------------------------------------------------------
function [com, inertia] = get_com_inertia(T, X)

% This function has been partially derived based on Dobrovolskis 1996, ICARUS 124,
% 698-704, Article No. 0243

% Empty case
if isempty(T)
    vol =[];
    com = [];    
    inertia = [];
    return
end

% Local coordinates
A(1).v = X(T(:,1),:);
A(2).v = X(T(:,2),:);
A(3).v = X(T(:,3),:);
A(4).v = X(T(:,4),:);

BxC = cross(A(2).v-A(1).v,A(3).v-A(1).v,2);
vol = abs(dot(BxC,A(4).v-A(1).v,2)/6);
com = (A(1).v+A(2).v+A(3).v+A(4).v)/4;
O = A(1).v;
D = A(2).v - A(1).v;
E = A(3).v - A(1).v;
F = A(4).v - A(1).v;
dim= size(A(4).v,1);
for n=1:3
    for m=1:3            
P(n).P(m).inertia = (2*D(:,n).*D(:,m) + 2*E(:,n).*E(:,m) + 2*F(:,n).*F(:,m) +...
D(:,n).*E(:,m) + D(:,m).*E(:,n) +...
D(:,n).*F(:,m) + D(:,m).*F(:,n) +...
E(:,n).*F(:,m) + E(:,m).*F(:,n)).*vol/20; 
    end
end
XX = P(2).P(2).inertia +P(3).P(3).inertia;
YY = P(3).P(3).inertia +P(1).P(1).inertia;
ZZ = P(1).P(1).inertia +P(2).P(2).inertia;
XY = P(1).P(2).inertia;
XZ = P(1).P(3).inertia;
YZ = P(2).P(3).inertia;
XX = XX + (2*com(:,2).*O(:,2) + 2*com(:,3).*O(:,3) + O(:,2).*O(:,2) + O(:,3).*O(:,3)).*vol;
YY = YY + (2*com(:,3).*O(:,3) + 2*com(:,1).*O(:,1) + O(:,3).*O(:,3) + O(:,1).*O(:,1)).*vol;
ZZ = ZZ + (2*com(:,1).*O(:,1) + 2*com(:,2).*O(:,2) + O(:,1).*O(:,1) + O(:,2).*O(:,2)).*vol;
XY = XY + (com(:,1).*O(:,2) + com(:,2).*O(:,1) + O(:,1).*O(:,2)).*vol;
XZ = XZ + (com(:,1).*O(:,3) + com(:,3).*O(:,1) + O(:,1).*O(:,3)).*vol;
YZ = YZ + (com(:,2).*O(:,3) + com(:,3).*O(:,2) + O(:,2).*O(:,3)).*vol;

com = ((A(1).v+A(2).v+A(3).v+A(4).v)/4).*[vol,vol,vol];
com          = sum(com,1)/sum(vol);
inertia(1,1) = sum(XX,1);
inertia(2,2) = sum(YY,1);
inertia(3,3) = sum(ZZ,1);
inertia(1,2) = - sum(XY,1);
inertia(1,3) = - sum(XZ,1);
inertia(2,3) = - sum(YZ,1);
inertia(2,1) = inertia(1,2);
inertia(3,1) = inertia(1,3);
inertia(3,2) = inertia(2,3);
