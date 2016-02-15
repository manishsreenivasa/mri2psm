function [ mat_bin ] = fnc_common_cell2bin(dim1, dim2, mat_cell)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% This function takes a mx1 cell Matrix and transforms in a dim1xdim2
% binary Matix
mat_bin= zeros (dim1,dim2);
mat_tmp = cell2mat(mat_cell);
if ~isempty(mat_tmp)
    mat_tmp(find(mat_tmp(:,1)>dim1),:)=[[]];
    mat_tmp(find(mat_tmp(:,2)>dim2),:)=[[]];
    dim_tmp = size(mat_tmp,1);
    for l=1:1:dim_tmp
        mat_bin( mat_tmp(l,1), mat_tmp(l,2)) = 1;
    end
end
end
