clear; clc; close all;

% Master program to plot >> some << of the more interesting results from the
% sample data set. You will need to install BTK with Matlab wrappers
% installed in order to run this code

% Launch program for plotting sagittal plane joint angles and torques
run('./programs/2_dynamicsComputations/matlab_plot_functions/plot_sagittalIKID.m');

% Note: if you have "puppeteer" installed, you can also view the animation using the following shell command
% puppeteer model/data_samplePSM.lua sampleData/c3d/sampleData.c3d sampleResults/ik_psm_sampleData.csv