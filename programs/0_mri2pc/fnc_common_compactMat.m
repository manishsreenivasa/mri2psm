function [ output_mat ] = fnc_common_compactMat( input_mat)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% This function takes one 1xn cell and returns 1 array as an compact Matrix
if isempty(input_mat)
    output_mat = [];
else
    output_mat = [];
    for i=1:1:length(input_mat)
        output_mat(size(output_mat,1) + 1:size(output_mat,1) + size(input_mat{i},1),:) = input_mat{i};
    end
end
end

