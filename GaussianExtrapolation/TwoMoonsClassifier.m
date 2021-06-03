%% Restart run
close all; clear; clc;
rng('default');

%% Set params
M = 8;
sSimParams = GetSimParams(M);

%% Load dataset
actualDataDist = 'TwoMoons';
nComponents = 2;
dataDim = 2;
nTrain = 200;
nTest = 400;
sDatasetParams.b_loadTwoMoonsMatFile = false;
interpMethod = 'NewPoints';

sDataset = GenerateDataset(actualDataDist, dataDim, nComponents, nTrain, nTest, interpMethod, sDatasetParams);
PlotTwoMoons(sSimParams, sDataset)

%%
sDistParams = EstimateDistributionParameters(sDataset, nComponents, 0.1);
sKernelParams = GetKernelParams(sDataset, sDistParams);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
                = CalcAnalyticEigenvalues(sSimParams.CalcEigenFuncsM, sKernelParams, sDataset.dim, nComponents);
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
plot_classifier(orig_laprlsc_classifier, sSimParams.outputFolder, sDataset.sData.xt, zeros(size(sDataset.sData.yt)))


%% EigRLS params
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0.1;  

sClassifier = BuildClassifier(sSimParams, sDataset, sKernelParams, 'eigrls', gamma_A_eigrls, gamma_I_eigrls);
PlotClassifier(sSimParams, sDataset, sClassifier)

% PlotCoefficients(sSimParams, c, vLambdaAnalytic);
