function fnc_prog0_mouseMove (object, eventdata)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

currentMousePointer = evalin('base','currentMousePointer');
pos = get (gca, 'CurrentPoint');

% See the cursor as *
plot((pos(1,1)),(pos(1,2)),'*g');

if isempty(currentMousePointer)
    currentMousePointer = (pos(1,1:2));
end

currentMousePointer(size(currentMousePointer,1)+1,:) = (pos(1,1:2));
assignin('base', 'currentMousePointer', currentMousePointer);
end
