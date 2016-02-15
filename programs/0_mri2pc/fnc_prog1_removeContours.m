function [ matOutput ] = fnc_prog1_removeContours( matInput, connComp )
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

sizeMatInput = size(matInput);
% Temporary Hack, skip this process for T2 images
if size(sizeMatInput) == [256 256]
    MatOutputInside = matInput;
    MatToAdd = zeros(size(matInput));
else
    for j=1:connComp.NumObjects
        remover{j} = connComp.ObjectsSeparate{j};
        for i=1:10
            matBound = bwboundaries (remover{j});
            remover{j} = remover{j} - fnc_common_cell2bin(sizeMatInput(1), sizeMatInput(2), matBound);
        end
    end
    
    % Get MatOutputOutside
    sumRemover = zeros(sizeMatInput);
    for j=1:connComp.NumObjects
        sumRemover = sumRemover + remover{j};
    end
    sumRemover = sumRemover + matInput;
    sumRemover (sumRemover == 1) = 0;
    sumRemover (sumRemover == 2) = 1;
    matOutput=sumRemover;
end
end
