function [Phi, lambdaAnalytic] = SimpleCalcAnalyticEigenfunctions(X, omega, nEigs)
GMMRegVal = 0.1;
nComponents = 1;
dim = size(X,2);
sDataset.estDataDist = 'Gaussian';
sDataset.sData.x = X;
sDistParams = EstimateDistributionParameters(sDataset, nComponents, GMMRegVal);
sKernelParams = GetKernelParams(sDataset, sDistParams, omega);
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, dim, nComponents);
[Phi, lambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, X, true);
end