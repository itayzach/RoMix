%% phi (Squared Exponentional)
function [vPhi_m] = phi(a, b, m, x)

if ~isvector(x)
    error('x has to be a vector');
end

% Calculate parameters
c = sqrt(a^2 + 2*a*b);

% m-th eigenfunction
if m <= 170
    % factorial(171) = Inf, but hermite() is more efficient than hermiteH...
    vHm = hermite(m, sqrt(2*c)*x);
else
    vHm = hermiteH(m, sqrt(2*c)*x);
end

vPhi_m = (1/sqrt(2^m*factorial(m)*sqrt(a/c))) * exp( -(c-a)*x.^2 ) .* vHm;

end