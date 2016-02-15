function fnc_prog0_cancel (object, eventdata)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

% Call variables from workspace
hardbone_segm   = evalin('base', 'hardbone_segm');
boneFat_segm    = evalin('base', 'boneFat_segm');
muscle_segm     = evalin('base', 'muscle_segm');
currScanNr      = evalin('base', 'currScanNr');
data_bmf        = evalin('base', 'data_bmf');
deleteSelection = evalin('base', 'deleteSelection');
figure2         = evalin('base', 'figure2');
evalin('base', 'boundaries = [];');

% Reset values
if isempty(data_bmf)
    nmbBound = 0;
else
    try
        nmbBound = length(data_bmf(currScanNr).chosenBoundaries);
    catch
        nmbBound = 0;
    end
end
assignin('base', 'nmbBound', nmbBound);
assignin('base', 'deleteSelection', 0);

% Refresh Plot
cla(figure2);
figure2 = subplot (1, 10, [2 7]);
axis equal; axis tight; hold on;
spy(boneFat_segm, 'y');
spy(muscle_segm, 'r');
spy(hardbone_segm, 'k');
% Plot Existing Boundaries
if ~isempty(data_bmf)
    try 
        for i=1:length(data_bmf(currScanNr).chosenBoundaries)
            plot(data_bmf(currScanNr).chosenBoundaries{i}(:,2), data_bmf(currScanNr).chosenBoundaries{i}(:,1),'ob', 'markers',5,'MarkerFaceColor','b')
        end
    catch
        display ('Unable to plot boundaries');
    end
end
title(['Scan Number:  ' num2str(currScanNr)])

% Save variable in the worksapce
assignin('base', 'figure2', figure2);

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
    'string','Previous scan', ...
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
