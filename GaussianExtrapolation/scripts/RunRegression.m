%% Restart
clc; clear; close all; 

%% Run top examples
b_clearLastRun = false;
selectedEigsPresets = {'Get1DUniformPreset', 'Get2DUniformPreset', 'GetSwissRollPreset'};
for presetInd = 1:numel(selectedEigsPresets)
    funcName = selectedEigsPresets{presetInd};
    Main(funcName, b_clearLastRun);
end

%% Run real world examples
b_clearLastRun = false;
selectedGraphSigPresets = {'GetBrazilWeatherPreset', 'GetMnistPreset'};
for presetInd = 1:numel(selectedGraphSigPresets)
    funcName = selectedGraphSigPresets{presetInd};
    Main(funcName, b_clearLastRun);
end

%% Run all other presets
b_clearLastRun = true;
presetsDir = 'presets';
sDir = dir(fullfile('presets','*.m'));
for presetInd = 1:numel(sDir)
    [~, funcName] = fileparts(fullfile(presetsDir,sDir(presetInd).name));
    if ismember(funcName, [selectedEigsPresets, selectedGraphSigPresets])
        continue;
    end
    Main(funcName, b_clearLastRun);
end
%% Run additional scripts
PlotGaussianKernelEigenfunsExample(); %Illustrate the first eigenfunctions of the 1-D Guassian kernel
TwoMoonsClassifier;
DigitsClassifier;
