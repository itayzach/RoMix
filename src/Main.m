function sWorkspace = Main(presetName, b_saveFigures, b_clearLastRun)
%% Restart
if ~exist('b_clearLastRun', 'var') || b_clearLastRun
    clc; close all; 
end
if ~exist('b_saveFigures', 'var')
    b_saveFigures = false;
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
    sEigResults = RunGraphSignalInterp(sPreset, sPlotParams, b_interpEigenvecs);
end
if sPreset.b_runGraphSignals
    b_interpEigenvecs = false;
    sSigResults = RunGraphSignalInterp(sPreset, sPlotParams, b_interpEigenvecs);
end
%% Save workspace
sWorkspace = SaveWorkspaceToStruct();
ReadWorkspaceStructToBase(sWorkspace);
end