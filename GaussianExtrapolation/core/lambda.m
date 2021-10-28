function [lambda_m] = lambda(sKernelParams, c, m)

dim = length(m);

if strcmp(sKernelParams.kernelType, 'gaussian')
    vLambda_m = zeros(dim, 1);
    for d = 1:dim
        beta = sKernelParams.beta{c}(d);
        vLambda_m(d) = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m(d);
    end
    lambda_m = prod(vLambda_m);
elseif strcmp(sKernelParams.kernelType, 'sinc')
%     lambda_m = ((m+1)^2*pi^2)/(4*sParams.a^2);
    lambda_m = sKernelParams.a/pi * (2*m + 1);
else
    error('unknown kernelType');
end
end