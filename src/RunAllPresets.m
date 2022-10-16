%% Restart
clc; clear; close all; 

%% Run illustrative examples
PlotSimpleSin();
PlotGaussianKernelEigenfunsExample(); %Illustrate the first eigenfunctions of the 1-D Guassian kernel

%% Swiss roll illustrations
b_vRandn = [true, false];
b_vGmmLatent = [true,  false];
for b_randn = b_vRandn
    for b_gmmLatent = b_vGmmLatent
        b_saveFigures = ( b_randn &&  b_gmmLatent) || ... Gaussian, latent
                        (~b_randn && ~b_gmmLatent);   ... Uniform, data
        PlotGmmSwissRollExample(b_randn, b_gmmLatent, b_saveFigures);
    end
end

%% Run toy examples
b_saveFigures = true;
b_saveResults = false;
b_clearLastRun = false;
selectedEigsPresets = {'Get1DUniformPreset', 'Get2DUniformPreset', 'GetSwissRollPreset'};
for presetInd = 1:numel(selectedEigsPresets)
   funcName = selectedEigsPresets{presetInd};
   Main(funcName, b_saveFigures, b_saveResults, b_clearLastRun);
end

%% Run real world examples
b_saveFigures = true;
b_saveResults = true;
b_clearLastRun = false;
selectedGraphSigPresets = {'GetBulgariBeaconsPreset', 'GetMnistVaePreset'};
for presetInd = 1:numel(selectedGraphSigPresets)
   funcName = selectedGraphSigPresets{presetInd};
   Main(funcName, b_saveFigures, b_saveResults, b_clearLastRun);
end

%% Run all other presets
b_saveFigures = false;
b_saveResults = false;
b_clearLastRun = true;
presetsDir = 'presets';
sDir = dir(fullfile('presets','*.m'));
for presetInd = 1:numel(sDir)
    [~, funcName] = fileparts(fullfile(presetsDir,sDir(presetInd).name));
    if ~ismember(funcName, [selectedEigsPresets, selectedGraphSigPresets])
        Main(funcName, b_saveFigures, b_saveResults, b_clearLastRun);
    end
end

%% Run additional scripts
TwoMoonsClassifier;
DigitsClassifier;
