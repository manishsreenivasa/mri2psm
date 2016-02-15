% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
nScans = evalin('base','nScans');
if begin_scanNr == nScans
    begin_scanNr = 1;
else
    begin_scanNr = begin_scanNr + 1;
end
save(['./begin_scanNr.mat'], 'begin_scanNr')
uiresume(gcbf)
