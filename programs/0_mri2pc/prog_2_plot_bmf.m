% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_2_bmf.m
% Plot the segmented bone, muscle and fat pointclouds
%

clf; clear; clc;

patientRootDir = '../../';
load([patientRootDir,'/model/data_bmf']);

bPlot_bmf_2D = 1;
bPlot_bmf_3D = 0;


if bPlot_bmf_2D
    for scanNr = 1:length(data_bmf)
        pause(0.1); clf;

        bone   = data_bmf(scanNr).bonePos;
        muscle = data_bmf(scanNr).musclePos;
        fat    = data_bmf(scanNr).fatPos;
        
        subplot(2,2,1); hold on; box on; grid on; title('Combined BMF');
        if (~isempty(bone))
            plot (bone(:,1), bone(:,2), 'ok', 'markersize', 1);
        end
        plot (muscle(:,1), muscle(:,2), 'or', 'markersize', 1);
        plot (fat(:,1), fat(:,2), 'oy', 'markersize', 1);
        axis([-200 200 -200 200]);
        axis square;
        
        subplot(2,2,2); hold on; box on; grid on; title('Bone');
        if (~isempty(bone))
            plot (bone(:,1), bone(:,2), 'ok', 'markersize', 1);
        end
        axis([-200 200 -200 200]);
        axis square;
        
        subplot(2,2,3); hold on; box on; grid on; title('Muscle');
        plot (muscle(:,1), muscle(:,2), 'or', 'markersize', 1);
        axis([-200 200 -200 200]);
        axis square;
        
        subplot(2,2,4); hold on; box on; grid on; title('Fat');
        plot (fat(:,1), fat(:,2), 'oy', 'markersize', 1);
        axis([-200 200 -200 200]);
        axis square;
        
        drawnow;
    end
end

if bPlot_bmf_3D
    viewAngle = [300 10];
    for scanNr = 1:length(data_bmf)
        bone    = data_bmf(scanNr).bonePos;
        muscle  = data_bmf(scanNr).musclePos;
        fat     = data_bmf(scanNr).fatPos;
        
        subplot(1,3,1); hold on; view(viewAngle); grid on;
        if (~isempty(bone))
            plot3(bone(:,1), bone(:,2), bone(:,3), 'ok', 'markersize', 1);
        end
        axis([-200 200 -200 200 0 800]);
        
        subplot(1,3,2); hold on; view(viewAngle); grid on;
        plot3(muscle(:,1), muscle(:,2), muscle(:,3), 'or', 'markersize', 1);
        axis ([-200 200 -200 200 0 800]);
        
        subplot(1,3,3); hold on; view(viewAngle); grid on;
        plot3 (fat(:,1), fat(:,2), fat(:,3), 'oy', 'markersize', 1);
        axis ([-200 200 -200 200 0 800]);
        
        drawnow;
        clear bone muscle fat
    end
end
