function vPhi_m = phiD(sKernelParams, c, I, x)

u     = sKernelParams.sDistParams.u{c};
mu    = sKernelParams.sDistParams.mu_1D{c};
sigma = sKernelParams.sDistParams.sigma{c};
beta  = sKernelParams.beta{c};

xu = x*u;
normFactor = (1+2*beta).^(1/8)./sqrt(2.^I.*factorial(I));
% hArg = hermiteArg(beta,xu,mu,sigma);
% hermiteArg = (1/4 + beta./2).^(1/4).*(xu-mu)./sigma;
vHm = hermiteD(I, beta,xu,mu,sigma);
mPhi_m = normFactor .* exp( -((xu-mu).^2./(2*sigma.^2)) .* ((sqrt(1+2*beta)-1)./2) ) .* vHm;
assert(~any(isnan(mPhi_m(:))), 'Phi contains NaN');
assert(isreal(mPhi_m), 'Phi is complex');
% vPhi_m = prod(sort(mPhi_m,2),2); % sort to avoid NaNs...
vPhi_m = prod(mPhi_m,2); 
assert(~any(isnan(vPhi_m)), 'Phi contains NaN');

end