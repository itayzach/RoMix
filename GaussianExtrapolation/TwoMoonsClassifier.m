%% Restart run
close all; clear; clc;
rng('default');
set(0,'DefaultFigureWindowStyle','normal')

%% Set params
M                = 20;
omega            = 0.3;
gmmNumComponents = 8;
gmmRegVal        = 1e-3;
gmmMaxIter       = 2000;

%% Load dataset
actualDataDist = 'TwoMoons';
nTrain = 200;
nTest = 400;
sDatasetParams.b_loadTwoMoonsMatFile = false;
interpMethod = 'NewPoints';

sDataset = GenerateDataset(actualDataDist, [], [], nTrain, nTest, interpMethod, sDatasetParams);
sSimParams = GetSimParams(M);
sSimParams.sDataset = sDataset;
PlotTwoMoons(sSimParams, sDataset)

%%
sDistParams = EstimateDistributionParameters(sDataset.sData.x, gmmNumComponents, gmmRegVal, gmmMaxIter);
nGmmPoints = 1000;

PlotGMM('GMM', sDistParams.GMModel, nGmmPoints);
nGmmPoints = nTrain;
pltTitle = ['Dataset with n = ', num2str(nTrain), ' points'];
plt2Title = ['Generated ' num2str(nGmmPoints), ' points from GMM with nEstComp = ' num2str(gmmNumComponents)];

windowStyle = 'normal';
PlotDataset(sSimParams, sDataset.sData.x, sDataset.sData.y, pltTitle, sDistParams.GMModel, nGmmPoints, plt2Title, windowStyle);

sKernelParams = GetKernelParams(sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
                = CalcAnalyticEigenvalues(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.dim, gmmNumComponents);
%% LapRLS
gamma_A_laprls = 0.03125;
gamma_I_laprls = 1;
kernel_sigma_laprls = 1/(6*sqrt(2)); 

options_laprls_orig = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma_laprls, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
orig_laprlsc_classifier = experiment_moon(sDataset.sData.x,sDataset.sData.y,sDataset.sData.xt,sDataset.sData.yt,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls_orig);
fprintf('---------------------------------------------------\n')
plot_classifier(orig_laprlsc_classifier, sSimParams, sDataset.sData.xt, zeros(size(sDataset.sData.yt)))


%% EigRLS params
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0.1;  

sClassifier = BuildClassifier(sSimParams, sDataset, sKernelParams, 'eigrls', gamma_A_eigrls, gamma_I_eigrls);
PlotCoefficients(sSimParams, sClassifier.c, sClassifier.vLambdaAnalytic);
PlotClassifier(sSimParams, sDataset, sClassifier)
