function [Phi, lambdaAnalytic] = SimpleCalcAnalyticEigenfunctions(X, omega, sigma, mu, nEigs)
%% Kernel parameters
sKernelParams.kernelType = 'gaussian';
sKernelParams.beta{1} = 2*sigma^2/omega^2;
sKernelParams.constsType = 2;
sKernelParams.omega = omega;

%% Distribution parameters
sKernelParams.sDistParams.mu{1} = mu;
sKernelParams.sDistParams.mu_1D{1} = mu;
sKernelParams.sDistParams.cov{1} = sigma;
sKernelParams.sDistParams.sigma{1} = sigma;
sKernelParams.sDistParams.sigma_eigv{1} = sigma;
sKernelParams.sDistParams.u{1} = 1;
sKernelParams.sDistParams.estNumComponents = 1;

nComponents = 1;
dim = size(X,2);

%% Calculate eigenvalues (and matching indices)
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, dim, nComponents);
%% Calculate eigenfunctions
[Phi, lambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, X, true);
end