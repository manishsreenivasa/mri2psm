function dist = fnc_dist2plane (point, planeABC, planeD)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_dist2plane

planeABC = planeABC./norm(planeABC);
dist = dot(planeABC, point) + planeD;
