%% Clean up
clear; close all; clc;
restoredefaultpath;

%% My directories
[filepath, scriptName, ~] = fileparts(mfilename('fullpath'));
fprintf('[%s] Adding paths...\n', scriptName);
addpath(genpath('GaussianExtrapolation'));
addpath(genpath('ManifoldLearn'));
addpath(genpath('data'));
fprintf('[%s] Ready.\n', scriptName);
clear scriptName

%% Belkin's
start_ManifoldLearn;

%% GSPBox
addpath('gspbox/');
gsp_start;

%% Wrap up
cd(filepath);
clear;