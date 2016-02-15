function [ mat_output] = fnc_prog0_completeBone ( matBound, matBound_bin, trabecullarBone)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% This function is used to add the cortical bone to the trabecullar bone
% It takes the Outside-Boundary of the selected bone and add the points
% that are part of the TrabecullarBone
% This function makes the sum of the translated matrix(up,down,left &
% right) and take the boundary of this new ones

detectionParameter = 4;
[dimY,dimX] = size(matBound_bin);
mat_orign = matBound_bin;
for i=1:detectionParameter
    % Create translated matrices
    yR = imfill(fnc_common_cell2bin(dimY,dimX,{[matBound(:,1)+1,matBound(:,2)]}));
    yL = imfill(fnc_common_cell2bin(dimY,dimX,{[matBound(:,1)-1,matBound(:,2)]}));
    xR = imfill(fnc_common_cell2bin(dimY,dimX,{[matBound(:,1),matBound(:,2)+1]}));
    xL = imfill(fnc_common_cell2bin(dimY,dimX,{[matBound(:,1),matBound(:,2)-1]}));
    
    matBound_bin = xR+xL+yR+yL+imfill(matBound_bin);
    matBound_bin(matBound_bin~=0) = 1;
    matBound_cell = bwboundaries(matBound_bin,'noholes');
    matBound = matBound_cell{1};
end

corticalBone = matBound_bin + trabecullarBone;
corticalBone (corticalBone ~= 2) = 0;
corticalBone (corticalBone == 2) = 1;
matFull = corticalBone + imfill(mat_orign);
matFull = fnc_common_bin2array(matFull);

[~, shp_output] = alphavol(matFull(:,1:2), 5);
shp_output_points = matFull(shp_output.bnd(:,1), :);
points_interp = fnc_prog0_interpolateSpaces(shp_output_points);

indx_remove=[];
for i=1:length(shp_output_points)
    point1=matFull(shp_output.bnd(i,1),:);
    point2=matFull(shp_output.bnd(i,2),:);
    if (point1(2)-point2(2)) == 0
        indx_remove=[indx_remove(:) ; find(((abs(points_interp(:,2)-point1(2))==1)+(points_interp(:,1)==point1(1)))==2)] ;
        indx_remove=[indx_remove(:) ; find(((abs(points_interp(:,2)-point2(2))==1)+(points_interp(:,1)==point2(1)))==2)] ;
    end
    if (point1(1)-point2(1)) == 0
        indx_remove=[indx_remove(:) ; find(((abs(points_interp(:,1)-point1(1))==1)+(points_interp(:,2)==point1(2)))==2)] ;
        indx_remove=[indx_remove(:) ; find(((abs(points_interp(:,1)-point2(1))==1)+(points_interp(:,2)==point2(2)))==2)] ;
    end
end

points_interp(unique(indx_remove),:) = [];
mat_outputBin = fnc_common_cell2bin(dimY,dimX,{points_interp});
mat_outputBin = imfill(mat_outputBin);

mat_output_tmp = bwboundaries(mat_outputBin);
matMax = [];
for j=1:length(mat_output_tmp)
    matMax(j) = length(mat_output_tmp{j});
end
[~, b] = max(matMax);
mat_output = mat_output_tmp{b};
end