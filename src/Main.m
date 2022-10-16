function sWorkspace = Main(presetName, b_saveFigures, b_saveResults, b_clearLastRun)
%% Restart
if ~exist('b_clearLastRun', 'var') || b_clearLastRun
    clc; close all; 
end
if ~exist('b_saveFigures', 'var')
    b_saveFigures = false;
end
if ~exist('b_saveResults', 'var')
    b_saveResults = false;
end
rng('default'); 
%% GMM / Spectral Clustering
clusterMethod = 'GMM'; % 'GMM' / 'SC'
%% Get perset
sPreset = GetPreset(presetName);
%% Verify preset
VerifyPresetParams(sPreset, clusterMethod);
%% Get plot params
sPlotParams = GetPlotParams(sPreset, b_saveFigures);
%% Run loop
if sPreset.b_interpEigenvecs
    b_interpEigenvecs = true;
    sEigResults = RunInterpGraphSignal(sPreset, sPlotParams, b_interpEigenvecs);
    if b_saveResults
        SaveResultsToFile(sPreset, sPlotParams, sEigResults, b_interpEigenvecs);
    end
end
if sPreset.b_runGraphSignals
    b_interpEigenvecs = false;
    sSigResults = RunInterpGraphSignal(sPreset, sPlotParams, b_interpEigenvecs);
    if b_saveResults
        SaveResultsToFile(sPreset, sPlotParams, sSigResults, b_interpEigenvecs);
    end
end
%% Save workspace
sWorkspace = SaveWorkspaceToStruct();
ReadWorkspaceStructToBase(sWorkspace);
end