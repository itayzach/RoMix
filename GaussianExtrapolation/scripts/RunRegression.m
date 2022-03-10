%% Restart
clc; clear; close all; 
set(0,'DefaultFigureWindowStyle','normal')

%% Run selected presets
b_clearLastRun = false;
selectedPresets = {'GetBrazilWeatherPreset', 'Get1DUniformPreset', 'GetSwissRollPreset', 'GetMnistPreset'};
for presetInd = 1:numel(selectedPresets)
    funcName = selectedPresets{presetInd};
    Main(funcName, b_clearLastRun);
end

%% Run all other presets
b_clearLastRun = true;
presetsDir = 'presets';
sDir = dir(fullfile('presets','*.m'));
for presetInd = 1:numel(sDir)
    [~, funcName] = fileparts(fullfile(presetsDir,sDir(presetInd).name));
    if ismember(funcName, selectedPresets)
        continue;
    end
    Main(funcName, b_clearLastRun);
    %keyboard;
    %close all;
end
%% Run additional scripts
PlotGaussianKernelEigenfunsExample(); %Illustrate the first eigenfunctions of the 1-D Guassian kernel
TwoMoonsClassifier;
DigitsClassifier;

%% Run MNIST
b_clearLastRun = true;
selectedPresets = {'GetMnistPreset'};
for presetInd = 1:numel(selectedPresets)
    funcName = selectedPresets{presetInd};
    Main(funcName, b_clearLastRun);
end