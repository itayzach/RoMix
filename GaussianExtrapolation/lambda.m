function [lambda_m] = lambda(sParams, m)

if sParams.constsType == 1
    a = sParams.a;
    b = sParams.b;

    % Calculate parameters
    c = sqrt(a^2 + 2*a*b);
    A = a + b + c;
    B = b/A;

    % m-th eigenvalue
    lambda_m = sqrt(2*a/A) * B^m;
elseif sParams.constsType == 2
    beta = sParams.beta;
    lambda_m = sqrt(2./(1+beta+sqrt(1+2*beta))) .* (beta./(1+beta+sqrt(1+2*beta))).^m;
else
    error('Unknown constsType');
end

end