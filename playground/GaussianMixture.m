clc; close all; clear;
rng('default'); % For reproducibility

nComponents = 2;
nPoints = 5000;
nEigsSingleComponent = 3;
nEigs = nComponents*nEigsSingleComponent;
sSimParams.PlotSpectM = nEigs;
sSimParams.CalcEigenFuncsM = nEigs;
sSimParams.outputFolder = 'figs';

mu = [3; 7];                % Means
sigma = [0.4 0.5]; cat(3,[0.5],[0.3]); % Covariances
% prop = [0.5 0.5];           % Mixing proportions
% gm = gmdistribution(mu,sigma,prop);
% xPoints = random(gm, nPoints);
vSel = rand(nPoints, 1) < 0.5;
% xPoints = (1-vSel).*(sigma(1)*randn(nPoints, 1) + mu(1)) + ...
%           vSel.*(sigma(2)*randn(nPoints, 1) + mu(2));
xPoints = [(sigma(1)*randn(floor(nPoints/2), 1) + mu(1)); ...
           (sigma(2)*randn(ceil(nPoints/2), 1) + mu(2)) ];
sDataset.sData.x = xPoints;
sDataset.dim = 1;
sDataset.actualDataDist = 'GaussianMixture';
      
sKernelParams.kernelType = 'gaussian';
sKernelParams.constsType = 2;
sKernelParams.omega = 0.3;
sKernelParams.sDistParams.u{1} = 1;
sKernelParams.sDistParams.u{2} = 1;

beta{1} = 2*sigma(1).^2/sKernelParams.omega^2;
beta{2} = 2*sigma(2).^2/sKernelParams.omega^2;

%% Analytic eigenfunctions
tPhiAnalytic = zeros(nComponents, nPoints, nEigsSingleComponent);
mLambda = zeros(nComponents, nEigsSingleComponent);
for c = 1:nComponents
    for m = 0:nEigsSingleComponent-1
        sKernelParams.beta{c} = beta{c};
        sKernelParams.sDistParams.mu_1D{c} = mu(c);
        sKernelParams.sDistParams.sigma{c} = sigma(c);
        tPhiAnalytic(c,:,m+1) = sqrt(nComponents/nPoints)*phi(sKernelParams, c, m, xPoints);
        mLambda(c,m+1) = lambda(sKernelParams,c,m);
    end
end

[ vLambdaAnalytic, vIdx ] = sort(mLambda(:), 'descend');
[row, col] = ind2sub([nComponents nEigs], vIdx);
mPhiAnalytic = zeros(nPoints, nEigs);
for i = 1:nEigs
    mPhiAnalytic(:,i) = tPhiAnalytic(row(i),:,col(i));
end

%% Numeric eigenvectors
dist = pdist2(xPoints, xPoints);
A = exp(-dist.^2/(2*sKernelParams.omega^2));

[mPhiNumeric, mLambdaNumeric] = eigs(A, nEigs);
[vLambdaNumeric, idx] = sort(diag(mLambdaNumeric), 'descend');
vLambdaNumeric = (nComponents/nPoints)*vLambdaNumeric;
mPhiNumeric = mPhiNumeric(:,idx);

PlotSpectrum(sSimParams, [], [], vLambdaAnalytic, vLambdaNumeric, []);
%% Flip sign
mPhiNumeric = FlipSign(mPhiAnalytic, mPhiNumeric);

%% RMSE
tPhiAnalyticRMSE(1,:,:) = mPhiAnalytic;
tPhiNumericRMSE(1,:,:) = mPhiNumeric;
vRMSEAnaVsNum = CalcRMSE(tPhiAnalyticRMSE, tPhiNumericRMSE, 'Analytic');
PlotRMSE(sSimParams, sDataset, [], vRMSEAnaVsNum, [])
%% Plot
PlotEigenfunctionsEigenvectors(sSimParams, sDataset, [], 0, nEigs-1, mPhiAnalytic, mPhiNumeric, 'Analytic');
