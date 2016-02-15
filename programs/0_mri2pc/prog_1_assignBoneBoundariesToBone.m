% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: prog_1_assignBoneBoundariesToBone
% This code reads choosen bone boundaries and asks user to assign them to bones
%

clear; clf; clc;
set(gcf, 'toolbar', 'figure');

bSave = 0;

% Patient Information
patientRootDir       = '../../';
dataDir_dicom        = [patientRootDir,'sampleData/dicom/'];
dataDir_segmentation = [patientRootDir,'sampleData/segmentation/'];
dataDir_image        = [patientRootDir,'sampleData/images/'];
dataset              = '/segmentation';
fileNameEnding       = '_Segmentation.h5';

% load data_bmf
dataExist_bmf = [patientRootDir,'model/data_bmf.mat'];
if exist(dataExist_bmf) == 2
    load(dataExist_bmf)
    disp('Loaded existing data_bmf')
else
    data_bmf = [];
    error ('No existing data_bmf found. Quitting.')
end

% define where the program begins and the last scan
begin_scanNr = 10;
nScans = length(dir(dataDir_image))-2; % because it starts in 0 and /. ,  /.. are directories

% Define distribution of T1 & T2 scans
T1_scans = [2:18 22:46 50:66];
T2_scans = [1 19:21 47:49];

while begin_scanNr~=nScans+1
    currScanNr = begin_scanNr;
    
    clf;
    
    % Define uicontrols-buttons
    Button_None = uicontrol('style','pushbutton',...
        'unit','pix','position',[400 440 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','None','Callback',@fnc_prog1_boneSelection);
    Button_Pelvis = uicontrol('style','pushbutton',...
        'unit','pix','position',[400 400 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','Pelvis','Callback',@fnc_prog1_boneSelection);
    Button_FemurR = uicontrol('style','pushbutton',...
        'unit','pix','position',[550 360 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FemurR','Callback',@fnc_prog1_boneSelection);
    Button_TibiaR = uicontrol('style','pushbutton',...
        'unit','pix','position',[550 280 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','TibiaR','Callback',@fnc_prog1_boneSelection);
    Button_FibulaR = uicontrol('style','pushbutton',...
        'unit','pix','position',[550 200 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FibulaR','Callback',@fnc_prog1_boneSelection);
    Button_FootR = uicontrol('style','pushbutton',...
        'unit','pix','position',[550 120 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FootR','Callback',@fnc_prog1_boneSelection);
    Button_FemurL = uicontrol('style','pushbutton',...
        'unit','pix','position',[250 360 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FemurL','Callback',@fnc_prog1_boneSelection);
    Button_TibiaL = uicontrol('style','pushbutton',...
        'unit','pix','position',[250 280 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','TibiaL','Callback',@fnc_prog1_boneSelection);
    Button_FibulaL = uicontrol('style','pushbutton',...
        'unit','pix','position',[250 200 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FibulaL','Callback',@fnc_prog1_boneSelection);
    Button_FootL = uicontrol('style','pushbutton',...
        'unit','pix','position',[250 120 150 30],...
        'fontsize',12,'fontweight','bold',...
        'string','FootL','Callback',@fnc_prog1_boneSelection);
    Button_NextScan = uicontrol('style','pushbutton',...
        'unit','pix','position',[1300 350 280 40],...
        'fontsize',12,'fontweight','bold',...
        'string','Next scan',...
        'Callback','fnc_prog1_nextScan');
    Button_PreviousScan = uicontrol('style','pushbutton',...
        'unit','pix','position',[1300 300 280 40],...
        'fontsize',12,'fontweight','bold',...
        'string','Previous scan',...
        'Callback','fnc_prog1_previousScan');
    Button_Exit = uicontrol('style','pushbutton',...
        'unit','pix','position',[1300 250 280 40],...
        'fontsize',12,'fontweight','bold',...
        'string','Exit',...
        'Callback','fnc_prog1_exit');
    
    % Read Data
    data_segm  = h5read ([dataDir_segmentation, int2str(currScanNr), fileNameEnding], dataset);
    data_image = imread([dataDir_image, int2str(currScanNr), '.tiff']);
    data_dicom = dicominfo ([dataDir_dicom num2str(currScanNr)]);
    
    % Create Matrix of the current selected boundaries
    try
        boundaries = data_bmf(currScanNr).chosenBoundaries;
    catch
        boundaries = [];
        data_bmf(currScanNr).chosenBoundaries = [];
    end
    
    % Create Matrix with Indexes for the Bones that are going to be chosen
    boundOfBone = zeros(length(boundaries),10);
    if isempty (data_bmf(currScanNr).IndexFromBoundToBone)
        data_bmf(currScanNr).IndexFromBoundToBone = boundOfBone(:,1:9);
    end
    
    % Make scans pictures brighter
    if find(currScanNr == T1_scans)
        data_image(:,:,:) = data_image(:,:,:)*3;
    else
        data_image(:,:,:) = data_image(:,:,:)*15;
    end
    
    % Make segmented matrices for BoneFat, Muscle, Background and HardBone
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
    
    % Find Connected Components
    connComp = fnc_prog1_getConnComp (hardbone_segm + boneFat_segm + muscle_segm);
    
    % Get the Dicom Information
    dim_height = data_dicom.Height;
    dim_width = data_dicom.Width;
    
    % Insert pixelSize, pixelCenter and pixelDim
    data_bmf(currScanNr).pixelSize   = data_dicom.PixelSpacing;
    data_bmf(currScanNr).pixelCenter = data_dicom.ImagePositionPatient;
    data_bmf(currScanNr).pixelDim    = [dim_height dim_width];
    
    % Insert bonePx and bonePos
    bone_segm = imfill(fnc_common_cell2bin(dim_height, dim_width, {fnc_common_compactMat(data_bmf(currScanNr).chosenBoundaries)}));
    boneFilled = fnc_common_bin2array(bone_segm);
    data_bmf(currScanNr).bonePx = boneFilled;
    data_bmf(currScanNr).bonePos = fnc_common_pixel2mm (data_bmf(currScanNr).pixelSize, data_bmf(currScanNr).pixelCenter, boneFilled);
    
    % Remove Hands
    muscle_segm   = fnc_prog1_removeHands (muscle_segm, connComp.ObjTogether);
    hardbone_segm = fnc_prog1_removeHands (hardbone_segm, connComp.ObjTogether);
    boneFat_segm  = fnc_prog1_removeHands (boneFat_segm, connComp.ObjTogether);
    
    % Insert musclePx and musclePos
    muscle_segm = muscle_segm + hardbone_segm;
    muscle_segm (muscle_segm == 2) = 1;
    muscle_segm = muscle_segm - bone_segm;
    muscle_segm (muscle_segm == -1)=0;
    
    % Remove Countours of muscle
    muscle_segm = fnc_prog1_removeContours (muscle_segm, connComp);
    data_bmf(currScanNr).musclePx = fnc_common_bin2array (muscle_segm);
    data_bmf(currScanNr).musclePos = fnc_common_pixel2mm (data_bmf(currScanNr).pixelSize, data_bmf(currScanNr).pixelCenter, data_bmf(currScanNr).musclePx);
    
    % Insert fatPx and fatPos
    fat_segm = boneFat_segm - bone_segm;
    fat_segm (fat_segm == -1) = 0;
    data_bmf(currScanNr).fatPx = fnc_common_bin2array (fat_segm);
    data_bmf(currScanNr).fatPos = fnc_common_pixel2mm (data_bmf(currScanNr).pixelSize, data_bmf(currScanNr).pixelCenter, data_bmf(currScanNr).fatPx);
    
    % Plot figures
    subplot(2,2,1); title(['Scan No.',int2str(currScanNr)],'FontWeight','bold'); hold on;
    if find(currScanNr == T1_scans)
        axisLims = [0 512 0 432];
    else
        axisLims = [0 256 0 256];
    end
    imshow(data_image);
    set(gca,'YDir','reverse')
    axis equal; axis(axisLims);
    
    if isempty(boundaries)
        subplot(2,2,2); hold on;
        title('Segmentation')
        spy(fat_segm, 'y');
        spy(muscle_segm, 'r');
        delete ([Button_None,Button_Pelvis,Button_FemurR,Button_FemurL,...
            Button_TibiaR,Button_TibiaL,Button_FibulaR,Button_FibulaL,...
            Button_FootR,Button_FootL])
        uiwait(gcf) % Wait till one button is pressed
    else
        for nBoundaries = 1:length(boundaries)
            
            boneNamesIndex = fnc_prog1_assignNamesToIndex(data_bmf(currScanNr).IndexFromBoundToBone, boundOfBone);
            
            if length(boneNamesIndex) < length(boundaries)
                boneNamesIndex(length(boneNamesIndex)+1:length(boundaries)) = {''};
            end
            currentBoundary = boundaries{nBoundaries};
            currentBoundary_Bin = imfill(fnc_common_cell2bin(dim_height, dim_width, {currentBoundary}));
            
            % Plot Segmentation
            subplot(2,2,2); hold on;
            title('Segmentation');
            spy(fat_segm, 'y');
            spy(muscle_segm, 'r');
            for i = 1:length(boundaries)
                spy(fnc_common_cell2bin(dim_height, dim_width, boundaries(i)), 'k')
            end
            spy(currentBoundary_Bin,'b')
            axis equal; axis(axisLims);
            
            % Plot selected Bones
            subplot(2, 2, [3 4]); hold on;
            title('L <--------      Selected bones      --------> R');
            spy(currentBoundary_Bin, 'b')
            text(50, 90, '-- Current Boundary', 'Color', 'b');
            text(50, 70, '-- Identified Boundaries', 'Color', 'r');
            text(50, 50, '-- Non-identified Boundaries', 'Color', 'k');
            
            for i = 1:length(boundaries)
                if ~isempty(boneNamesIndex{i})
                    spy(fnc_common_cell2bin (dim_height, dim_width, boundaries(i)), 'r')
                    if ~isempty (boneNamesIndex{i})
                        xy = min(boundaries{i});
                        text(xy(2) - 5, xy(1) - 5, cell2mat(boneNamesIndex(i)), 'Color', 'r');
                    end
                else
                    spy(fnc_common_cell2bin (dim_height, dim_width, boundaries(i)), 'k')
                end
            end
            axis equal; axis(axisLims);
            
            clear currentBoundary
            uiwait(gcf) % Wait till one button is pressed
            delete(subplot(2,2,2)); delete(subplot(2,2,3));
            
            % if NextScan or PreviousScan was pressed
            if begin_scanNr > currScanNr || begin_scanNr < currScanNr
                break;
            end
            
        end
    end % End of all boundaries in slice
    
    if ~(begin_scanNr > currScanNr || begin_scanNr < currScanNr)
        begin_scanNr=begin_scanNr+1;
        save(['./begin_scanNr.mat'], 'begin_scanNr')
        data_bmf(currScanNr).IndexFromBoundToBone = boundOfBone(:,1:9);
        display('Saving changes to local memory');
    end
    if bSave
        save(dataExist_bmf, 'data_bmf')
        display('Saving changes to data_bmf file');
    end
end
clearvars -except data_bmf
