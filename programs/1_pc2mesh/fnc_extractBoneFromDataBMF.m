function bone = fnc_extractBoneFromDataBMF (data_bmf)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_extractBoneFromDataBMF

bone = struct('idx', {[0 0],[0 0],[0 0],[0 0],[0 0],[0 0],[0 0],[0 0],[0 0]}, 'pos', {[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]});

for scanNr = 1:length(data_bmf)
    for boneNr = 1:9
        try
            idxBoundaryInBone = find([data_bmf(scanNr).IndexFromBoundToBone(:,boneNr)] == 1);
        catch
            continue;
        end
        
        for boundaryNr = 1:length(idxBoundaryInBone)
            currentIndexes = cell2mat(data_bmf(scanNr).chosenBoundaries(idxBoundaryInBone(boundaryNr)));
            [a,b]=size(bone(boneNr).idx);
            bone(boneNr).idx(a:a+length(currentIndexes)-1,:) = currentIndexes;
            
            pixelCenter = repmat([data_bmf(scanNr).pixelCenter(2) data_bmf(scanNr).pixelCenter(1)],length(currentIndexes),1);
            currentPositions = currentIndexes.*repmat([data_bmf(scanNr).pixelSize]',length(currentIndexes),1) + pixelCenter;
            
            currentPositions(:,3) = repmat(data_bmf(scanNr).pixelCenter(3),length(currentIndexes),1);
            bone(boneNr).bonePos(a:a+length(currentIndexes)-1,:) = currentPositions;
            clear a b currentIndexes currentPositions
        end
        clear idxBoundaryInBone
    end
end
