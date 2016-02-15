function [ output_mat ] = fnc_common_bin2array( input_mat)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% This function takes an mxn bin-Matrix and returns an px2 array with the
% its coordinates

[m,n] = size (input_mat);
coord = 1;
if input_mat == zeros(m,n)
    output_mat = [];
end
for i = 1:n
    for j = 1:m
        if input_mat(j,i) == 1
            output_mat(coord,1:2) = [j, i];
            coord=coord+1;
        end
    end
end
end
