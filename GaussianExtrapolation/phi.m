%% phi (Squared Exponentional)
function [vPhi_m] = phi(sParams, m, x, d)
if ~isvector(x)
    error('x has to be a vector');
end

if sParams.constsType == 1
    a = sParams.a(d);
    b = sParams.b;

    % Calculate parameters
    c = sqrt(a^2 + 2*a*b);

    % m-th eigenfunction
    if m <= 170
        % factorial(171) = Inf, but hermite() is more efficient than hermiteH...
        vHm = hermite(m, sqrt(2*c)*x);
    else
        vHm = hermiteH(m, sqrt(2*c)*x);
    end

    normFactor = (1/sqrt(2^m*factorial(m)*sqrt(a/c)));
    vPhi_m = normFactor * exp( -(c-a)*x.^2 ) .* vHm;

elseif sParams.constsType == 2
    if exist('d', 'var')
        mu = sParams.mu(d);
        sigma = sParams.sigma(d);
        beta = sParams.beta(d);
    else
        mu = sParams.mu;
        sigma = sParams.sigma;
        beta = sParams.beta;
    end

    normFactor = (1+2*beta)^(1/8)/sqrt(2^m*factorial(m));
    vHm = hermite(m, (1/4 + beta/2)^(1/4)*(x-mu)/sigma);
    vPhi_m = normFactor * exp( -((x-mu).^2/(2*sigma^2)) * ((sqrt(1+2*beta)-1)/2) ) .* vHm;
elseif sParams.constsType == 3
    if exist('d', 'var')
        mu = sParams.mu(d);
        sigma = sParams.sigma(d);     
        alpha = sParams.alpha(d); % Should be here, or alpha(d)?
    else
        mu = sParams.mu;
        sigma = sParams.sigma;
    end    
    eps = sParams.eps;
    

    normFactor = (1+(2*eps/alpha)^2)^(1/8) / sqrt(2^m*factorial(m));
    vHm = hermite(m, (1+(2*eps/alpha)^2).^(1/4).*alpha.*x );
    vPhi_m = normFactor * exp( -(sqrt(1+(2*eps/alpha)^2)-1)*alpha^2*x.^2/2 ) .* vHm;
else
    error('Unknown constsType');
end
end