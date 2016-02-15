function [ mat_output ] = fnc_prog0_interpolateSpaces( mat_input )
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

mat_output=[];
if size(mat_input,1)>9    
    intp_x  = interp(mat_input(:,1),30);
    intp_y  = interp(mat_input(:,2),30);
    intp_xF = floor(intp_x);
    intp_yF = floor(intp_y);
    intp_xC = ceil(intp_x);
    intp_yC = ceil(intp_y);
else
    intp_xF = floor(mat_input(:,1));
    intp_yF = floor(mat_input(:,2));
    intp_xC = ceil(mat_input(:,1));
    intp_yC = ceil(mat_input(:,2));
end

mat_output(:,1) = intp_xF;
mat_output(:,2) = intp_yF;
if isempty(mat_output)
    mat_output = [0,0];
else
    mat_output(size(intp_xF,1)+1:size(intp_xF,1)+size(intp_xC,1),:) = [intp_xC intp_yC];
    mat_output = unique(mat_output,'rows','stable');
    mat_output (find(mat_output(:,1)<=0),:) = [[]];
    mat_output (find(mat_output(:,2)<=0),:) = [[]];
end

end

