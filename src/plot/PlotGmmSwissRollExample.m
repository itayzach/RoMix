function PlotGmmSwissRollExample()
close all; clc; rng('default')

%% Generate Swiss Roll
b_interpEigenvecs = true;
b_saveFigures = true;
b_randn = true;
nLabeled = 1000;
n        = 1000;
N        = 5000;
sPreset = GetSwissRollPreset(b_randn);
sPreset = UpdatePreset(sPreset,b_interpEigenvecs,nLabeled);
sPreset.n = n;
sPreset.N = N;
sPlotParams = GetPlotParams(sPreset, b_saveFigures);
sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs);

%% Run GMM
sDistParams = EstimateDistributionParameters(sDataset.sData.x, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);

%% Plot
nGmmPoints = 2000;
pltTitle = []; %['Dataset with n = ', num2str(sPreset.n), ' points'];
plt2Title = []; %['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(sPreset.gmmNumComponents)];
windowStyle = 'normal';
cXAxisLabels = {'$\theta$', '$t$'};
PlotDataset(sPlotParams, sPreset, sDataset.sData.x, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, false, cXAxisLabels); % latent space
PlotDataset(sPlotParams, sPreset, sDataset.sData.S, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, true, cXAxisLabels); % data space 
PlotGmmResultWithDataset(sPlotParams, sDataset.sData.x, sDistParams, windowStyle, true); % b_plotCovMean = true
PlotGmmResultWithDataset(sPlotParams, sDataset.sData.x, sDistParams, windowStyle, false); % b_plotCovMean = false
end