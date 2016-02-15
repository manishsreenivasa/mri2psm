function [ matNames ] = fnc_prog1_assignNamesToIndex( mat_bmf, mat_selectedBound )
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

for i=1:size(mat_bmf,1)
    
    if isempty(find(mat_selectedBound(i,:)))
        currentIndex = find(mat_bmf(i,:));
    else
        currentIndex = find(mat_selectedBound(i,:));
    end
    
    if currentIndex == 1
        matNames{i}='Pelvis';
    elseif currentIndex==2
        matNames{i}='FemurR';
    elseif currentIndex==3
        matNames{i}='TibiaR';
    elseif currentIndex==4
        matNames{i}='FibulaR';
    elseif currentIndex==5
        matNames{i}='FootR';
    elseif currentIndex==6
        matNames{i}='FemurL';
    elseif currentIndex==7
        matNames{i}='TibiaL';
    elseif currentIndex==8
        matNames{i}='FibulaL';
    elseif currentIndex==9
        matNames{i}='FootL';
    else
        matNames{i} =[];
    end
end
end
