function [matOutput]= fnc_prog1_removeHands(matInput, connComp)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
removeMat = matInput + connComp;
removeMat (removeMat == 1) = 0;
removeMat (removeMat == 2) = 1;
matOutput = removeMat;
end
