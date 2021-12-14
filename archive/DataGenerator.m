%% Restart
clc;
clear;
rng('default');
close all;
set(0,'DefaultFigureWindowStyle','docked')

%% Plot flags
b_plotCdf = true;

%% Method parameters
b_saturateT   = true;
pCdfDegree    = 5;
invpCdfDegree = 5;

%% Dataset parameters
verticesPDF         = 'Grid'; % 'Gaussian' / 'Uniform' / 'Grid' / 'TwoMoons' / 'TwoSpirals' / 'SwissRoll'
sDatasetParams.xMin = -1;
sDatasetParams.xMax = 1;
dim                 = 1;
nComponents         = 1;
n                   = 10000;
N                   = 0;
%% Generate data
sDataset = GenerateDataset(verticesPDF, dim, nComponents, n, N, 'NewPoints', sDatasetParams);
x = sDataset.sData.x;

%% Learn pCdf
muTilde = 0;
sigmaTilde = 1;
nEvalPoints = n;min(700, round(n/10));
% [xTrainGrid, estMarginalCdf_xTrain, pCdf, invpCdf] = ...
%     PolyfitEstCdf(x, nEvalPoints, pCdfDegree, invpCdfDegree, b_plotCdf);
% PlotPolyCdfDemonstration1(sDataset.xMin, sDataset.xMax, pCdf, xTrainGrid, estMarginalCdf_xTrain, muTilde, sigmaTilde)

%% Learn pInvCdf
kdeCdf = ksdensity(x, x, 'Function', 'cdf');
invpCdf = polyfit(kdeCdf, x, invpCdfDegree);

%% T inverse
nTilde = 10000;
dimTilde = 1;
nBins = nTilde/100;
zTilde = sigmaTilde*randn(nTilde,dimTilde) + muTilde;
z = Tinv(invpCdf, muTilde, sigmaTilde, zTilde);

%% Reference probability distribution, P
ba = sDatasetParams.xMax-sDatasetParams.xMin; % (b - a)
p = (1/ba)*(ba/nBins)*ones(nBins,dim);

%% KLDiv of Tinv vs analytic
h_Tinv = histcounts(z,nBins)';
q_Tinv = h_Tinv/nTilde;
TinvKLDiv = nansum( p .* log2( p./q_Tinv ) );

%% KLDiv of histcounts vs analytic
h_x = histcounts(x,nBins)';
q_hist = h_x/nTilde;
hist_KLDiv = nansum( p .* log2( p./q_hist ) );

%% Plot and calc KL divergence
figure('Name', 'Compare histograms');
subplot(3,1,1)
histfit(zTilde)
title('histogram($\tilde{Z}$)', ...
    'interpreter', 'latex', 'FontSize', 16)
subplot(3,1,2)
histogram(z,nBins);
title(['histogram(Z), $Z=T^{-1}(\tilde{Z})$, KLDiv = ', num2str(TinvKLDiv)], ...
    'interpreter', 'latex', 'FontSize', 16)
subplot(3,1,3)
histogram(x, nBins);
title(['histogram(X), KLDiv = ', num2str(hist_KLDiv)], ...
    'interpreter', 'latex', 'FontSize', 16)


