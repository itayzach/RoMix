%% Restart
clc; clear; close all; 

%% Run illustrative examples
PlotGaussianKernelEigenfunsExample(); %Illustrate the first eigenfunctions of the 1-D Guassian kernel

%% Run toy examples
b_saveFigures = true;
b_saveResults = true;
b_clearLastRun = false;
selectedEigsPresets = {'Get1DUniformPreset', 'Get2DUniformPreset', 'GetSwissRollPreset', 'GetTwoMoonsPreset'};
for presetInd = 1:numel(selectedEigsPresets)
   funcName = selectedEigsPresets{presetInd};
   Main(funcName, b_saveFigures, b_saveResults, b_clearLastRun);
end

%% Run real world examples
b_saveFigures = true;
b_saveResults = true;
b_clearLastRun = false;
selectedGraphSigPresets = {'GetBulgariBeaconsPreset', 'GetMnistVaePreset', 'GetMnistPreset'};
for presetInd = 1:numel(selectedGraphSigPresets)
   funcName = selectedGraphSigPresets{presetInd};
   Main(funcName, b_saveFigures, b_saveResults, b_clearLastRun);
end
