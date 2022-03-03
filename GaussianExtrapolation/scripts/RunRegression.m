%% Restart
clc; clear; close all; 
set(0,'DefaultFigureWindowStyle','normal')

%% Run selected presets
selectedPresets = {'GetBrazilWeatherPreset', 'Get1DUniformPreset', 'GetSwissRollPreset', 'GetMnistPreset'};
for presetInd = 1:numel(selectedPresets)
    funcName = selectedPresets{presetInd};
    Main(funcName);
end

%% Run all other presets
presetsDir = 'presets';
sDir = dir(fullfile('presets','*.m'));
for presetInd = 1:numel(sDir)
    [~, funcName] = fileparts(fullfile(presetsDir,sDir(presetInd).name));
    if ismember(funcName, selectedPresets)
        continue;
    end
    Main(funcName);
    %keyboard;
    %close all;
end
%% Run additional scripts
PlotGaussianKernelEigenfunsExample(); %Illustrate the first eigenfunctions of the 1-D Guassian kernel
TwoMoonsClassifier;
DigitsClassifier;

%% Run MNIST
selectedPresets = {'GetMnistPreset'};
for presetInd = 1:numel(selectedPresets)
    funcName = selectedPresets{presetInd};
    Main(funcName);
end