function [Phi, lambdaAnalytic] = SimpleCalcAnalyticEigenfunctions(X, omegaTilde, sigmaTilde, muTilde, nEigs)
dim = size(X,2);
nComponents = 1;

%% Distribution parameters
sDistParams.cov{1} = sigmaTilde;
sDistParams.sigma{1} = diag(sigmaTilde);
sDistParams.sigma_eigv{1} = sigmaTilde;
[sDistParams.u{1}, sDistParams.sigma_eigv{1}] = eig(sDistParams.cov{1});
sDistParams.mu{1} = muTilde;
sDistParams.mu_1D{1} = sDistParams.mu{1}*sDistParams.u{1};
sDistParams.estNumComponents = nComponents;

%% Kernel parameters
sKernelParams.kernelType = 'gaussian';
sKernelParams.beta{1} = 2*sDistParams.sigma{1}.^2./omegaTilde^2;
sKernelParams.constsType = 2;
sKernelParams.omega = omegaTilde;
sKernelParams.sDistParams = sDistParams;

%% Calculate eigenvalues (and matching indices)
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams);
%% Calculate eigenfunctions
[Phi, lambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, X, true);
end