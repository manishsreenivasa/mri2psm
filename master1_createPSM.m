clear; clc; close all;

% Master program to load segmented bone, muscle and fat data and create lua
% files. Note: The code runs a lot faster if you disable plotting (set
% "bPlot" flags to 0 in each of the programs called below).

bSaveConsoleOutput = 1;
consoleOutputFile = './matlabConsoleOutput.txt';

if bSaveConsoleOutput
    diary(consoleOutputFile)
    diary ON
end

tic
% Launch program for mesh morphing
run('./programs/1_pc2mesh/prog_0_meshMorph'); drawnow; clear;

% Launch program for deforming femur
run('./programs/1_pc2mesh/prog_1_femoralDeformation'); drawnow; clear;

% Launch program for computing joint properties
run('./programs/1_pc2mesh/prog_2_computeJointAxes'); drawnow; clear; pause(2.0);

% Launch program for computing segment properties
run('./programs/1_pc2mesh/prog_3_segmentBMF'); drawnow; clear;

% Launch program for writing PSM and GSM LUA files
run('./programs/1_pc2mesh/prog_4_createLUA.m'); drawnow; clear;

toc
diary OFF