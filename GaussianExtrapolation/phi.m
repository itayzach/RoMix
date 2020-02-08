%% phi (Squared Exponentional)
function [vPhi_m, lambda_m] = phi(a, b, m, x)

if ~isvector(x)
    error('x has to be a vector');
end

% Calculate parameters
c = sqrt(a^2 + 2*a*b);
A = a + b + c;
B = b/A;

% m-th eigenvalue
lambda_m = sqrt(2*a/A) * B^m;

% m-th eigenfunction
% vHm = hermiteH(m, sqrt(2*c)*x);
vHm = hermite(m, sqrt(2*c)*x);
vPhi_m = (1/sqrt(2^m*factorial(m)*sqrt(a/c))) * exp( -(c-a)*x.^2 ) .* vHm;

end