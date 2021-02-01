function [Phi, lambdaAnalytic] = SimpleCalcAnalyticEigenfunctions(X, omega, sigma, mu, nEigs)
dim = size(X,2);
nComponents = 1;

%% Distribution parameters
sDistParams.cov{1} = sigma;
sDistParams.sigma{1} = diag(sigma);
sDistParams.sigma_eigv{1} = sigma;
[sDistParams.u{1}, sDistParams.sigma_eigv{1}] = eig(sDistParams.cov{1});
sDistParams.mu{1} = mu;
sDistParams.mu_1D{1} = sDistParams.mu{1}*sDistParams.u{1};
sDistParams.estNumComponents = nComponents;

%% Kernel parameters
sKernelParams.kernelType = 'gaussian';
sKernelParams.beta{1} = 2*sDistParams.sigma{1}.^2./omega^2;
sKernelParams.constsType = 2;
sKernelParams.omega = omega;
sKernelParams.sDistParams = sDistParams;

%% Calculate eigenvalues (and matching indices)
[sKernelParams.vLambdaAnalytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, dim, nComponents);
%% Calculate eigenfunctions
[Phi, lambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, X, true);
end