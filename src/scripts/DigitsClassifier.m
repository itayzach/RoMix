clc; clear; close all;
rng('default')
% sPreset = GetUspsPreset(); b_transpose = false; sPlotParams.actualDataDist = 'USPS';
sPreset = GetMnistPreset(); b_transpose = true; sPlotParams.actualDataDist = 'MNIST';
sDataset = GenerateDataset(sPreset.verticesPDF, sPreset.dim, sPreset.nGenDataCompnts, sPreset.n, sPreset.N, sPreset.dataGenTechnique, sPreset.sDatasetParams);
X = sDataset.sData.x;
sPlotParams.outputFolder = 'figs';
sPlotParams.dim = size(X,2);

%% Plot samples
vSamples = round(linspace(1,sPreset.n,50));
figTitle = ['Given $n = ', num2str(length(vSamples)), '$ points'];
figName = ['MNIST_digits_n_', num2str(length(vSamples))];
PlotDigits(sPlotParams, X(vSamples,:), [], b_transpose, figTitle, figName)
%% GMM on data
% 20  --> AIC = 2548438.46
% 100 --> AIC = 4390387.71
% 150 --> AIC = 80972218.64
% 200 --> AIC = 10878840.01
% 300 --> AIC = 173361967.27
gmmDataComp = 300;
sDistParams = EstimateDistributionParameters(X,gmmDataComp, sPreset.gmmRegVal, sPreset.gmmMaxIter);
nGmmPoints = 50;
xRand = sDistParams.GMModel.random(nGmmPoints);
figTitle = ['Generated $n = ', num2str(nGmmPoints), '$ points with k = ', num2str(sDistParams.GMModel.NumComponents), ' components'];
figName = ['MNIST_GMM_generated_digits_n_', num2str(nGmmPoints), '_k_', num2str(sDistParams.GMModel.NumComponents)];
PlotDigits(sPlotParams, xRand, [], b_transpose, figTitle, figName)


%% Our method
sKernelParams = CalcKernelParams(sDistParams, sPreset.omegaTilde);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
     = CalcAnalyticEigenvalues(sPreset.MTilde, sKernelParams);
[ Phi, lambdaPhi ] = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, X);
mSigHatPhi = RoMix(Phi, sPreset.gamma1, sPreset.gamma2, lambdaPhi, [], sDataset.sData.ymasked, sPreset.b_maskDataFitTerm);
mSigRecPhi = Phi*mSigHatPhi;
vSigCnvrtRecPhi = ConvertSignalByDataset(sPreset.verticesPDF, mSigRecPhi);
vSigCnvrt       = ConvertSignalByDataset(sPreset.verticesPDF, sDataset.sData.y);

trainAccRecPhi = 100*sum(vSigCnvrtRecPhi == vSigCnvrt) / sPreset.n;
fprintf('\nOurs train accuracy    = %.2f%%\n', trainAccRecPhi);
figTitle = ['Prediction on given $n = ', num2str(length(1:201:sPreset.n)), '$ points'];
PlotDigits([], X(1:201:sPreset.n,:), vSigCnvrtRecPhi(1:201:sPreset.n)-b_transpose, b_transpose, figTitle)
%% Generate new examples from GMM and predict digit
nNew = 100;
xRand = sDistParams.GMModel.random(nNew);
PhiNew = CalcAnalyticEigenfunctions(sPreset.MTilde, sKernelParams, xRand);
mSigNew = PhiNew*mSigHatPhi;
mSigCnvrtNew = ConvertSignalByDataset(sPreset.verticesPDF, mSigNew);
figTitle = [num2str(nGmmPoints), ' generated points'];
PlotDigits([], xRand, mSigCnvrtNew-b_transpose, b_transpose,figTitle)

%% Test data
% vPredictionsInt = mSigCnvrtInt;
% mSigCnvrt       = ConvertSignalByDataset(sPreset.verticesPDF, mSig);
% testAccInt = 100*sum(vPredictionsInt == mSigCnvrtRef) / N;
% fprintf('Ours test accuracy     = %.2f%%\n', testAccInt);


%% GMM after PCA
% [eigenvectors,scores,latent,~,explained,mu] = pca(X);
% K = find(cumsum(explained)>95,1);
% z = scores(:,1:K);
% xRec = z*eigenvectors(:,1:K).' + mu;
% ind = randi(size(X,1),20,1);
% PlotDigits(xRec(ind,:),b_transpose)
% 
% gmmLatentComp = 5:20;
% sDistParamsPCA = EstimateDistributionParameters(z,gmmLatentComp, sPreset.gmmRegVal, sPreset.gmmMaxIter);
% nNew = 20;
% zRand = sDistParamsPCA.GMModel.random(nNew);
% xRand = zRand*eigenvectors(:,1:K).' + mu;
% PlotDigits(xRand,b_transpose)





