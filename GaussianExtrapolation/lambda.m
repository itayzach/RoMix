function [lambda_m] = lambda(sParams, m)

if sParams.constsType == 1
    a = sParams.a;
    b = sParams.b;

    % Calculate parameters
    c = sqrt(a.^2 + 2*a.*b);
    A = a + b + c;
    B = b./A;

    % m-th eigenvalue
    lambda_m = sqrt(2*a./A) .* B.^m;
elseif sParams.constsType == 2
    beta = sParams.beta;
    lambda_m = sqrt(2./(1+beta+sqrt(1+2*beta))) .* (beta./(1+beta+sqrt(1+2*beta))).^m;
elseif sParams.constsType == 3
    eps = sParams.eps;
    alpha = sParams.alpha;
    
    lambda_m = (alpha*eps^(2*m)) ./ (( (alpha.^2/2).*(1+sqrt(1+(2*eps./alpha).^2))+eps^2 ).^(m+1/2));
else
    error('Unknown constsType');
end

end