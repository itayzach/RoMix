function [vPhi_m] = phi_d(sKernelParams, c, m, x, d)
if ~isvector(x)
    error('x has to be a vector');
end

if strcmp(sKernelParams.kernelType, 'gaussian')
    mu = sKernelParams.sDistParams.mu_1D{c}(d);
    sigma = sKernelParams.sDistParams.sigma{c}(d);
    beta = sKernelParams.beta{c}(d);
    assert(sigma > 0)
    normFactor = (1+2*beta)^(1/8)/sqrt(2^m*factorial(m));
    if m <= 50
        % factorial(171) = Inf, but hermite() is more efficient than hermiteH... 
        % note that for high m it might be numerically unstable
        vHm = hermite(m, (1/4 + beta/2)^(1/4)*(x-mu)/sigma);
     else
        if m == 50
            warning('m >= 50. using MATLAB''s hermiteH instead of factorial hermite. calculation will take time...')
        end
        vHm = hermiteH(m, (1/4 + beta/2)^(1/4)*(x-mu)/sigma);
    end
    vPhi_m = normFactor * exp( -((x-mu).^2/(2*sigma^2)) * ((sqrt(1+2*beta)-1)/2) ) .* vHm;
    assert(~any(isnan(vPhi_m)), 'Phi %d contains NaN', m);
    assert(isreal(vPhi_m), 'Phi %d is complex', m);

elseif strcmp(sKernelParams.kernelType, 'sinc')
%     vPhi_m = real(sqrt(pi./(2*sParams.a*x)) .* besselj(m+0.5, sParams.a*x));
%     vPhi_m(x < 0) = -vPhi_m(x < 0);
%     if m == 0
%         vPhi_m(x == 0) = 1;
%     else
%         vPhi_m(x == 0) = 0;
%     end
    vPhi_m = besselj(m, sKernelParams.a*x);
else
    error('unknown kernelType');
end
end