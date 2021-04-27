function vPr = SimpleEstPorbablityArea(X, sigma, mu)
%% Distribution parameters
sDistParams.cov{1} = sigma;
sDistParams.sigma{1} = diag(sigma);
sDistParams.sigma_eigv{1} = sigma;
[sDistParams.u{1}, sDistParams.sigma_eigv{1}] = eig(sDistParams.cov{1});
sDistParams.mu{1} = mu;
sDistParams.mu_1D{1} = sDistParams.mu{1}*sDistParams.u{1};
sDistParams.estNumComponents = 1;
sDistParams.componentProportion(1) = 1;
sDistParams.estDataDist = 'Gaussian';
sDistParams.dim = size(X,2);
vPr = EstimateProbabilityArea(sDistParams, X);
end