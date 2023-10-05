%% Restart
clc; clear; close all; 

%% Paper presets
RunAllPaperPresets;

%% Simple sin
PlotSimpleSin();

%% Run all additional presets
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