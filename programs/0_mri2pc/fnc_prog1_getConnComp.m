function [ connComp ] = fnc_prog1_getConnComp (mat_bmf)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

bw_connectivity = bwconncomp(mat_bmf);
connComp.NumObjects = 0;
connComp.ImageSize = bw_connectivity.ImageSize;
connComp.ObjTogether = zeros(bw_connectivity.ImageSize);
for i=1:bw_connectivity.NumObjects
    if size(bw_connectivity.PixelIdxList{i},1) > 1000
        Obj= zeros(bw_connectivity.ImageSize);
        [X,Y] = ind2sub(bw_connectivity.ImageSize, bw_connectivity.PixelIdxList{i});
        for j=1:size(X,1)
            Obj(X(j),Y(j)) = 1;
        end
        connComp.ObjectsSeparate{connComp.NumObjects+1} = Obj;
        connComp.ObjTogether = Obj + connComp.ObjTogether;
        connComp.NumObjects = connComp.NumObjects + 1;
    end
end
end
