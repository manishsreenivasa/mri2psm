function [l_step,r_step] = fnc_getSteppingEvents(c3d_fileName)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

h1 = btkReadAcquisition(c3d_fileName);
events = btkGetEvents(h1);
samplingFreq = btkGetPointFrequency(h1);
try
    l_step = ceil(events.Left_Foot_Strike*samplingFreq)-btkGetFirstFrame(h1);
    if length(l_step)<2
        l_step = [NaN NaN];
    end
    try
        btkGetPoint(h1,'LGroundReactionForce');
%         display('Valid Left Step');
    catch
        l_step = [NaN NaN];
    end
catch
    l_step = [NaN NaN];
end

try
    r_step = ceil(events.Right_Foot_Strike*samplingFreq)-btkGetFirstFrame(h1);
    if length(r_step)<2
        r_step = [NaN NaN];
    end
    try
        btkGetPoint(h1,'RGroundReactionForce');
%         display('Valid Right Step');
    catch
        r_step = [NaN NaN];
    end
catch
    r_step = [NaN NaN];
end
