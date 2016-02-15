function fnc_prog1_boneSelection (object, eventdata)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Call Workspace-variables
Button_None    = evalin('base','Button_None');
Button_Pelvis  = evalin('base','Button_Pelvis');
Button_FemurR  = evalin('base','Button_FemurR');
Button_TibiaR  = evalin('base','Button_TibiaR');
Button_FibulaR = evalin('base','Button_FibulaR');
Button_FootR   = evalin('base','Button_FootR');
Button_FemurL  = evalin('base','Button_FemurL');
Button_TibiaL  = evalin('base','Button_TibiaL');
Button_FibulaL = evalin('base','Button_FibulaL');
Button_FootL   = evalin('base','Button_FootL');
boundOfBone    = evalin('base','boundOfBone');
nBoundaries    = evalin('base','nBoundaries');

% Store the selection
if object==Button_Pelvis
    boundOfBone(nBoundaries,1)=1;
elseif object==Button_FemurR
    boundOfBone(nBoundaries,2)=1;
elseif object==Button_TibiaR
    boundOfBone(nBoundaries,3)=1;
elseif object==Button_FibulaR
    boundOfBone(nBoundaries,4)=1;
elseif object==Button_FootR
    boundOfBone(nBoundaries,5)=1;
elseif object==Button_FemurL
    boundOfBone(nBoundaries,6)=1;
elseif object==Button_TibiaL
    boundOfBone(nBoundaries,7)=1;
elseif object==Button_FibulaL
    boundOfBone(nBoundaries,8)=1;
elseif object==Button_FootL
    boundOfBone(nBoundaries,9)=1;
elseif object==Button_None
    boundOfBone(nBoundaries,10)=1;
end

assignin('base','boundOfBone', boundOfBone);
uiresume(gcbf)
clear;
end
