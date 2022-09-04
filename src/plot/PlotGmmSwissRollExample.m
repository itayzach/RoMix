function PlotGmmSwissRollExample()
close all; clc; rng('default')

%% Generate Swiss Roll
b_interpEigenvecs = true;
b_saveFigures = true;
sPreset = GetSwissRollPreset(true);
sPreset.n = 1000;
sPreset.nLabeled = 1000;
sPreset.N = 5000;
sPlotParams = GetPlotParams(sPreset, b_saveFigures);
sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs);

%% Run GMM
sDistParams = EstimateDistributionParameters(sDataset.sData.x, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);

%% Plot
nGmmPoints = 2000;
pltTitle = []; %['Dataset with n = ', num2str(sPreset.n), ' points'];
plt2Title = []; %['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(sPreset.gmmNumComponents)];
windowStyle = 'normal';
b_transform = true;
PlotDataset(sPlotParams, sPreset, sDataset.sData.x, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle);
PlotDataset(sPlotParams, sPreset, sDataset.sData.S, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform);
end