function [ mm_mat ] = fnc_common_pixel2mm( pixelSize, pixelCenter, pixel_mat)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
%  This function takes the cordinates of a pixel Matrix and transforms
%  in mm
pixel_spacing_xy = [pixelSize(1) 0 ; 0 pixelSize(2)];
pos_x = pixelCenter(1);
pos_y = pixelCenter(2);
pos_z = pixelCenter(3);

if isempty(pixel_mat)
    mm_mat=[];
else
    mm_mat = pixel_mat*pixel_spacing_xy;
    mm_mat (:,1) = mm_mat(:,1) + pos_y;
    mm_mat (:,2) = mm_mat(:,2) + pos_x;
    mm_mat (:,3) = pos_z;
end
end
