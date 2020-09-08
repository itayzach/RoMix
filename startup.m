clear; close all; clc;
[filepath, scriptName, ~] = fileparts(mfilename('fullpath'));
fprintf('[%s] Adding paths...\n', scriptName);
addpath(genpath('GaussianExtrapolation'));
addpath(genpath('ManifoldLearn'));
addpath(genpath('data'));
fprintf('[%s] Ready.\n', scriptName);
clear scriptName

start_ManifoldLearn;

cd('../Downloaded/gspbox-git/');
gsp_start;
cd(filepath);
clear;