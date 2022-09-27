function PlotGmmSwissRollExample(b_randn, b_gmmLatent, b_saveFigures)
%
% b_randn = false; b_gmmLatent = false; b_saveFigures = false;
% PlotGmmSwissRollExample(b_randn, b_gmmLatent, b_saveFigures); 
% 
close all; clc; rng('default')

%% Generate Swiss Roll
b_interpEigenvecs = true;

nLabeled = 1000;
n        = 1000;
N        = 5000;
sPreset = GetSwissRollPreset(b_randn, b_gmmLatent);
sPreset = UpdatePreset(sPreset,b_interpEigenvecs,nLabeled);
sPreset.n = n;
sPreset.N = N;
sPlotParams = GetPlotParams(sPreset, b_saveFigures);
sDataset = GenerateDataset(sPlotParams, sPreset, b_interpEigenvecs);
x = sDataset.sData.x;
S = sDataset.sData.S;
%% Run GMM
sDistParams = EstimateDistributionParameters(x, sPreset.gmmNumComponents, sPreset.gmmRegVal, sPreset.gmmMaxIter);

%% Plot
nGmmPoints = 2000;
pltTitle = []; %['Dataset with n = ', num2str(sPreset.n), ' points'];
plt2Title = []; %['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(sPreset.gmmNumComponents)];
windowStyle = 'normal';
cXAxisLabels = {'$\theta$', '$t$'};

% GMM space
b_transform = false; PlotDataset(sPlotParams, sPreset, x, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform, cXAxisLabels);

% other space
if b_gmmLatent   
    b_transform = true; PlotDataset(sPlotParams, sPreset, S, sDataset.sData.y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle, b_transform, cXAxisLabels);
else
    b_transform = false; PlotDataset(sPlotParams, sPreset, S, sDataset.sData.y, pltTitle, [], nGmmPoints, plt2Title, windowStyle, b_transform, cXAxisLabels);
end
if ~b_randn && ~b_gmmLatent
    cXAxisLabels = {'$x_1$','$x_2$','$x_3$'};
else
    cXAxisLabels = [];
end
b_plotCovMean = true; PlotGmmResultWithDataset(sPlotParams, x, sDistParams, b_plotCovMean, cXAxisLabels);
b_plotCovMean = false; PlotGmmResultWithDataset(sPlotParams, x, sDistParams, b_plotCovMean, cXAxisLabels);
end