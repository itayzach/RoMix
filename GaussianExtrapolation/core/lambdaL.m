function [lambdaL_m] = lambdaL(sParams, m)

assert(length(m) == sParams.dim, 'index dimensions mismatch');

if strcmp(sParams.kernelType, 'gaussian')
    vLambdaL_m = zeros(sParams.dim, 1);
    for d = 1:sParams.dim
        if sParams.constsType == 1
            ell = sParams.ell;
            error('Unknown constsType');
        elseif sParams.constsType == 2
            omega = sParams.omega;
            t = sParams.t;
            beta = sParams.beta(d);
            
            lambdaK_m = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m(d);
            
%             vLambdaL_m(d) = (sParams.dim*log(2*pi*omega^2) - 2*log(lambdaK_m))/omega^2;
            vLambdaL_m(d) = (1/t)*log(1/lambdaK_m);
        elseif sParams.constsType == 3
            assert('TODO...')
        else
            error('Unknown constsType');
        end
    end
    lambdaL_m = prod(vLambdaL_m);
elseif strcmp(sParams.kernelType, 'sinc')
    error('sinc is not supported');
else
    error('unknown kernelType');
end
end