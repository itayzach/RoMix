%% Clean up
clear; close all; clc;
restoredefaultpath;

%% My directories
[filepath, scriptName, ~] = fileparts(mfilename('fullpath'));
fprintf('[%s] Adding paths...\n', scriptName);
addpath(genpath(filepath));
fprintf('[%s] Ready.\n', scriptName);
clear scriptName

%% Belkin's
% addpath(genpath(fullfile('..','ManifoldLearn')));
% start_ManifoldLearn;

%% GSPBox
% addpath(fullfile('..','gspbox'));
% gsp_start;

%% Wrap up
cd(filepath);
clear;