function fnc_prog0_drawSelection (varargin)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Load Variables
tmp_hardbone_segm = evalin('base', 'hardbone_segm');
tmp_boneFat_segm  = evalin('base', 'boneFat_segm');
tmp_muscle_segm   = evalin('base', 'muscle_segm');
boundaries        = evalin('base', 'boundaries');
currScanNr        = evalin('base', 'currScanNr');
data_bmf          = evalin('base', 'data_bmf');
nmbBound          = evalin('base', 'nmbBound');
figure2           = evalin('base', 'figure2');

click = 0;
figure2 = subplot (1, 10, [2 7]);
text(50, 10, 'Click LEFT mouse button to start drawing. Click again to stop','color','r');
text(50, 30, 'Click RIGHT mouse button anywhere *inside* closed shape to choose it','color','r');
[x_start, y_start, click] = ginput(1);
set (gcf, 'WindowButtonMotionFcn', @fnc_prog0_mouseMove); % Set mouseMove to get positions of the mouse

if click == 1   % Left Mouse Button
    currentSelection=[];
    currentSelection_bin = [];
    mouse_wait = 1;
    
    while mouse_wait ~= 0
        mouse_wait = waitforbuttonpress;
        currentMousePointer = evalin('base','currentMousePointer');
        
        if mouse_wait == 0; % When Click the mouse for the second-time
            set (gcf, 'WindowButtonMotionFcn', []); % Stop Mouse Move
            currentSelection(:,1) = currentMousePointer(:,2);
            currentSelection(:,2) = currentMousePointer(:,1);
            currentSelection = fnc_prog0_interpolateSpaces(currentSelection); %interpolate the spaces left on the selection
            currentSelection_bin = fnc_common_cell2bin(size(tmp_hardbone_segm,1), size(tmp_hardbone_segm,2), {currentSelection});
            
            % Update segmentation based on selection
            tmp_hardbone_segm = tmp_hardbone_segm + currentSelection_bin;
            tmp_hardbone_segm (tmp_hardbone_segm == 2) = 1;
            
            % Make sure both muscle and fat pixels are not overlapping with newly drawn HardBone selection
            tmp_boneFat_segm = tmp_boneFat_segm - currentSelection_bin;
            tmp_boneFat_segm(tmp_boneFat_segm == -1) = 0;
            tmp_muscle_segm = tmp_muscle_segm - currentSelection_bin;
            tmp_muscle_segm(tmp_muscle_segm == -1) = 0;
            
            % Refresh Plot
            cla(figure2)
            figure2 = subplot (1, 10, [2 7]); hold on;
            spy(tmp_boneFat_segm, 'y');
            spy(tmp_muscle_segm, 'r');
            spy(tmp_hardbone_segm, 'k');
            spy(currentSelection_bin,'g');
            
            % Plot Existing Boundaries
            if ~isempty(data_bmf)
                try
                    for i=1:length(data_bmf(currScanNr).chosenBoundaries)
                        plot(data_bmf(currScanNr).chosenBoundaries{i}(:,2), data_bmf(currScanNr).chosenBoundaries{i}(:,1),'ob','markers',5,'MarkerFaceColor','b');
                    end
                catch
                    display ('Unable to plot boundaries');
                end
                
            end
            title(['Scan Number:  ' num2str(currScanNr)])
        end
    end
    
end

if click == 3 % (right mouse click)
    set (gcf, 'WindowButtonMotionFcn', []); % Stop Mouse Move
    boundaries_bone = bwboundaries(tmp_boneFat_segm);
    boundaries_bone_bin = fnc_common_cell2bin(size(tmp_boneFat_segm,1), size(tmp_boneFat_segm,2), boundaries_bone);
    r1_x = round(x_start); r1_y = round(y_start);
    
    % Find first x-coordinate on boundary
    try
        while boundaries_bone_bin(r1_y, r1_x) == 0
            r1_x = r1_x + 1;
        end
    catch
        r1_x = 0; r1_y = 0;
        
    end
    
    if r1_x ~= 0
        % Scroll through all boundaries and search for one that contains current X,Y ordinates
        r1_x_search = []; r1_y_search = 1;
        while isempty (r1_x_search)
            r1_x_search = find (boundaries_bone{r1_y_search,:}(:,1) == r1_y & boundaries_bone{r1_y_search,:}(:,2) == r1_x);
            r1_y_search = r1_y_search + 1;
        end
        selectedBoneBoundary = boundaries_bone{r1_y_search-1};        
        
        % Add hard bone to current selection
        selectedBoneBoundary_bin = fnc_common_cell2bin (size(tmp_boneFat_segm,1), size(tmp_boneFat_segm,2), {selectedBoneBoundary});
        selectedBoneBoundary = fnc_prog0_completeBone (selectedBoneBoundary, selectedBoneBoundary_bin, tmp_hardbone_segm);
        
        % Highlight current selection
        figure2 = subplot (1,10,[2 7]);
        plot(selectedBoneBoundary(:,2), selectedBoneBoundary(:,1), 'og');
        boundaries = selectedBoneBoundary;
        nmbBound = nmbBound + 1;
    end
end


% Save varibles on the workspace
assignin('base','tmp_hardbone_segm',tmp_hardbone_segm)
assignin('base','tmp_boneFat_segm',tmp_boneFat_segm)
assignin('base','tmp_muscle_segm',tmp_muscle_segm)
assignin('base','boundaries',boundaries)
assignin('base','nmbBound',nmbBound)
assignin('base','figure2',figure2)

% Define Buttons uicontrol
DrawSelection= uicontrol('style','text',...
    'unit','pix',...
    'position',[40 500 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Draw/Select');
SaveDrawSelection= uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[40 450 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Save Selection',...
    'Callback',@fnc_prog0_saveDrawSelection);
Cancel= uicontrol('style','pushbutton',...
    'unit','pix',...
    'position',[40 400 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','Cancel',...
    'Callback',@fnc_prog0_cancel);
NextScan= uicontrol('style','text',...
    'unit','pix',...
    'position',[40 350 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','NextScan');
PreviousScan= uicontrol('style','text',...
    'unit','pix',...
    'position',[40 300 280 40],...
    'fontsize',12,...
    'fontweight','bold',...
    'string','PreviousScan');

% Clear temporary variables
assignin('base','currentMousePointer',[]);
clear boundaries_bone boundaries_bone_bin r1_x r1_y r1_x_search r1_y_search selectedBoneBoundary selectedBoneBoundary_bin
end
