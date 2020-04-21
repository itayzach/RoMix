function [lambda_m] = lambda(sParams, m)

assert(length(m) == sParams.dim, 'index dimensions mismatch');

if strcmp(sParams.kernelType, 'exp')
    vLambda_m = zeros(sParams.dim, 1);
    for d = 1:sParams.dim
        if sParams.constsType == 1
            a = sParams.a(d);
            b = sParams.b;

            % Calculate parameters
            c = sqrt(a^2 + 2*a*b);
            A = a + b + c;
            B = b/A;

            % m-th eigenvalue
            vLambda_m(d) = sqrt(2*a/A) * B^m(d);
        elseif sParams.constsType == 2
            beta = sParams.beta(d);
            vLambda_m(d) = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m(d);
        elseif sParams.constsType == 3
            eps = sParams.eps;
            alpha = sParams.alpha(d);

            vLambda_m(d) = (alpha*eps^(2*m(d))) / (( (alpha^2/2)*(1+sqrt(1+(2*eps/alpha)^2))+eps^2 )^(m(d)+1/2));
        else
            error('Unknown constsType');
        end
    end
    lambda_m = prod(vLambda_m);
elseif strcmp(sParams.kernelType, 'sinc')
%     lambda_m = ((m+1)^2*pi^2)/(4*sParams.a^2);
    lambda_m = sParams.a/pi * (2*m + 1);
else
    error('unknown kernelType');
end
end