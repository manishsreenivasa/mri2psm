% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_0_chooseBoneBoundaries
% This code reads Ilastik output and asks user to manually choose/draw bone boundaries
% from automatic BoneFat, HardBone, and Muscle segmentation
% Cortical bone segmentation located around choosen bone boundaries are added to the final bone boundaries
%

clear; clf; clc; 
set (gcf,'toolbar','figure')

% Flag to save modifications to memory
bSave = 0;

% Patient Information
patientRootDir  = '../../';
dataDir_image   = [patientRootDir,'sampleData/images/'];
dataset         = '/segmentation'; 
fileNameEnding  = '_Segmentation.h5';

% Load data_bmf if it exists
dataExist_bmf = [patientRootDir,'model/data_bmf.mat'];
if exist(dataExist_bmf) == 2
    load(dataExist_bmf)
    disp('Loaded existing data_bmf')
else
    data_bmf = [];
    disp('Load data_bmf: No data_bmf file found, creating a new one')
end

% Read segmentations and probability densities from Ilastik
dataDir_segmentation = [patientRootDir,'sampleData/segmentation/'];

% Number of available scans
currScanNr = 5;
nScans = length(dir([patientRootDir 'sampleData/images/']))-2;

% Define distribution of T1 & T2 scans
T1_scans = [2:18 22:46 50:66];
T2_scans = [1 19:21 47:49];

% Load Data and assign them values
data_image      = imread([dataDir_image, int2str(currScanNr),'.tiff']);
data_segm       = h5read ([dataDir_segmentation,int2str(currScanNr),fileNameEnding], dataset);
muscle_segm     = double(squeeze(data_segm(1,:,:)));
background_segm = double(squeeze(data_segm(1,:,:)));
boneFat_segm    = double(squeeze(data_segm(1,:,:)));
hardbone_segm   = double(squeeze(data_segm(1,:,:)));
muscle_segm (muscle_segm ~= 1)         = 0;
muscle_segm (muscle_segm == 1)         = 1;
background_segm (background_segm ~= 2) = 0;
background_segm (background_segm == 2) = 1;
boneFat_segm (boneFat_segm ~= 3)       = 0;
boneFat_segm (boneFat_segm == 3)       = 1;
hardbone_segm (hardbone_segm ~= 4)     = 0;
hardbone_segm (hardbone_segm == 4)     = 1;

% Define variables
currentMousePointer = []; % Pointer used while Drawing
boundaries          = []; % Selected Boundaries 
deleteSelection     = 0;  % Variable that defines if selection should be erased

% Make scans pictures brighter
if find(currScanNr == T1_scans)
    data_image(:,:,:) = data_image(:,:,:)*4;
else
    data_image(:,:,:) = data_image(:,:,:)*15;
end

% Plot Data
figure1 = subplot (1, 10, [8 10]); hold on;
imshow (data_image);

figure2 = subplot (1, 10, [2 7]);
axis equal; axis tight; hold on; 
spy(boneFat_segm, 'y');
spy(muscle_segm, 'r');
spy(hardbone_segm, 'k');
if ~isempty(data_bmf)
    try
        for i=1:length(data_bmf(currScanNr).chosenBoundaries)
            plot (data_bmf(currScanNr).chosenBoundaries{i}(:,2), data_bmf(currScanNr).chosenBoundaries{i}(:,1), 'ob', 'markers', 5, 'MarkerFaceColor', 'b');
        end
    catch
        display ('Unable to plot boundaries');
    end
end
title(['Scan Number:  ' num2str(currScanNr)])

% Get number of existing boundaries
if isempty(data_bmf) % If data_bmf is empty
    nmbBound = 0;
else
    try
        nmbBound = length(data_bmf(currScanNr).chosenBoundaries);
    catch
        nmbBound = 0;
    end
end

% Define Buttons for uicontrol
DrawSelection= uicontrol('style','pushbutton',...
    'unit','pix','position',[40 500 280 40],...
    'fontsize',12,'fontweight','bold',...
    'string','Draw/Select','Callback',@fnc_prog0_drawSelection);
SaveDrawSelection= uicontrol('style','text',...
    'unit','pix','position',[40 450 280 40],...
    'fontsize',12,'fontweight','bold',...
    'string','Save');
Cancel= uicontrol('style','text',...
    'unit','pix','position',[40 400 280 40],...
    'fontsize',12,'fontweight','bold',...
    'string','Cancel');
NextScan= uicontrol('style','pushbutton',...
    'unit','pix','position',[40 350 280 40],...
    'fontsize',12,'fontweight','bold',...
    'string','Next scan',...
    'Callback',@fnc_prog0_nextScan);
PreviousScan= uicontrol('style','pushbutton',...
    'unit','pix','position',[40 300 280 40],...
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
