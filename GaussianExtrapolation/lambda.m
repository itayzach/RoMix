function [lambda_m] = lambda(sKernelParams, m)

dim = length(m);

if strcmp(sKernelParams.kernelType, 'exp')
    vLambda_m = zeros(dim, 1);
    for d = 1:dim
        if sKernelParams.constsType == 1
            a = sKernelParams.a(d);
            b = sKernelParams.b;

            % Calculate parameters
            c = sqrt(a^2 + 2*a*b);
            A = a + b + c;
            B = b/A;

            % m-th eigenvalue
            vLambda_m(d) = sqrt(2*a/A) * B^m(d);
        elseif sKernelParams.constsType == 2
            beta = sKernelParams.beta(d);
            vLambda_m(d) = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m(d);
        elseif sKernelParams.constsType == 3
            eps = sKernelParams.eps;
            alpha = sKernelParams.alpha(d);

            vLambda_m(d) = (alpha*eps^(2*m(d))) / (( (alpha^2/2)*(1+sqrt(1+(2*eps/alpha)^2))+eps^2 )^(m(d)+1/2));
        else
            error('Unknown constsType');
        end
    end
    lambda_m = prod(vLambda_m);
elseif strcmp(sKernelParams.kernelType, 'sinc')
%     lambda_m = ((m+1)^2*pi^2)/(4*sParams.a^2);
    lambda_m = sKernelParams.a/pi * (2*m + 1);
else
    error('unknown kernelType');
end
end