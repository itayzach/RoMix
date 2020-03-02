%% phi (Squared Exponentional)
function [vPhi_m] = phi(sParams, m, x)
if ~isvector(x)
    error('x has to be a vector');
end

% a = sParams.a;
% b = sParams.b;
% 
% % Calculate parameters
% c = sqrt(a^2 + 2*a*b);
% 
% % m-th eigenfunction
% if m <= 170
%     % factorial(171) = Inf, but hermite() is more efficient than hermiteH...
%     vHm = hermite(m, sqrt(2*c)*x);
% else
%     vHm = hermiteH(m, sqrt(2*c)*x);
% end
% 
% normFactor = (1/sqrt(2^m*factorial(m)*sqrt(a/c)));
% vPhi_m = normFactor * exp( -(c-a)*x.^2 ) .* vHm;

beta = sParams.beta;
mu = sParams.mu;
sigma = sParams.sigma;

normFactor = (1+2*beta)^(1/8)/sqrt(2^m*factorial(m));
vHm = hermite(m, (1/4 + beta/2)^(1/4)*(x-mu)/sigma);
vPhi_m = normFactor * exp( -((x-mu).^2/(2*sigma^2)) * ((sqrt(1+2*beta)-1)/2) ) .* vHm;

end