function fnc_prog0_nextScan (varargin)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Get some varibales from the workspace
currScanNr = evalin('base','currScanNr');
nScans     = evalin('base','nScans');
T1_scans   = evalin('base','T1_scans');
T2_scans   = evalin('base','T2_scans');

if currScanNr == nScans
    currScanNr = 1; % Recycle to first scan
else
    currScanNr = currScanNr + 1;
end

assignin ('base', 'currScanNr', currScanNr);
evalin ('base', 'boundaries = [];');

% Call variables from workspace
data_bmf             = evalin('base','data_bmf');
dataDir_segmentation = evalin('base','dataDir_segmentation');
dataDir_image        = evalin('base','dataDir_image');
dataset              = evalin('base','dataset');
fileNameEnding       = evalin('base','fileNameEnding');

if isempty(data_bmf)
    nmbBound=0;
else
    try
        nmbBound = length(data_bmf(currScanNr).chosenBoundaries);
    catch
        nmbBound = 0;
    end
end
assignin('base','nmbBound',nmbBound);

% Load Data
data_segm  = h5read ([dataDir_segmentation, int2str(currScanNr), fileNameEnding], dataset);
data_image = imread ([dataDir_image, int2str(currScanNr), '.tiff']);

% Assign data values
muscle_segm     = double(squeeze(data_segm(1,:,:)));
background_segm = double(squeeze(data_segm(1,:,:)));
boneFat_segm    = double(squeeze(data_segm(1,:,:)));
hardbone_segm   = double(squeeze(data_segm(1,:,:)));
muscle_segm (muscle_segm ~= 1) = 0;
muscle_segm (muscle_segm == 1) = 1;
background_segm (background_segm ~= 2) = 0;
background_segm (background_segm == 2) = 1;
boneFat_segm (boneFat_segm ~= 3) = 0;
boneFat_segm (boneFat_segm == 3) = 1;
hardbone_segm (hardbone_segm ~= 4) = 0;
hardbone_segm (hardbone_segm == 4) = 1;

% Save variables on the workspace
assignin ('base', 'hardbone_segm', hardbone_segm)
assignin ('base', 'boneFat_segm', boneFat_segm)
assignin ('base', 'muscle_segm', muscle_segm)

% Define distribution of T1 & T2 scans & Make scans pictures brighter
if find(currScanNr == T1_scans)
    data_image(:,:,:) = data_image(:,:,:)*4;
else
    data_image(:,:,:) = data_image(:,:,:)*15;
end

% Plot Data
clf;
figure1 = subplot (1, 10, [8 10]);
imshow (data_image);

figure2 = subplot (1, 10, [2 7]);
axis equal; axis tight; hold on;
spy(boneFat_segm, 'y');
spy(muscle_segm, 'r');
spy(hardbone_segm, 'k');
title(['Scan Number: ', num2str(currScanNr)])

% Plot Existing Boundaries
if ~isempty(data_bmf)
    try 
        for i=1:length(data_bmf(currScanNr).chosenBoundaries)
            plot (data_bmf(currScanNr).chosenBoundaries{i}(:,2), data_bmf(currScanNr).chosenBoundaries{i}(:,1), 'ob', 'markers', 5, 'MarkerFaceColor', 'b')
        end
    catch
        display ('Unable to plot boundaries');
    end
end

assignin('base', 'figure1', figure1)
assignin('base', 'figure2', figure2)

% Define Buttons uicontrol
DrawSelection= uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[40 500 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Draw/Select',...
    'Callback',@fnc_prog0_drawSelection);
SaveDrawSelection= uicontrol('style','text',...
    'unit','pix',...
    'position',[40 450 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Save');
Cancel= uicontrol('style','text',...
    'unit','pix',...
    'position',[40 400 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Cancel');
NextScan= uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[40 350 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Next scan',...
    'Callback',@fnc_prog0_nextScan);
PreviousScan= uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[40 300 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Previous scan',...
    'Callback',@fnc_prog0_previousScan);

if isempty(data_bmf)
    DeleteSelection= uicontrol('style','text',...
        'unit','pix','position',[40 250 280 40],...
        'fontsize',12,...
        'fontweight','bold',...
        'string','Delete Selection');
else
    try
        if isempty(data_bmf(currScanNr).chosenBoundaries)
            DeleteSelection= uicontrol('style','text',...
                'unit','pix','position',[40 250 280 40],...
                'fontsize',12,...
                'fontweight','bold',...
                'string','Delete Selection');
        else
            DeleteSelection= uicontrol('style','pushbutton',...
                'unit','pix','position',[40 250 280 40],...
                'fontsize',12,...
                'fontweight','bold',...
                'string','Delete Selection', ...
                'Callback',@fnc_prog0_deleteSelection);
        end
    catch
        DeleteSelection= uicontrol('style','text',...
            'unit','pix','position',[40 250 280 40],...
            'fontsize',12,...
            'fontweight','bold',...
            'string','Delete Selection');
    end
end

end

