function [vPhi_m, lambda_m] = SqExpEig(a, b, m, x)

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
if m <= 170
    % factorial(171) = Inf, but hermite() is more efficient than hermiteH...
    vHm = hermite(m, sqrt(2*c)*x);
else
    vHm = hermiteH(m, sqrt(2*c)*x);
end
vPhi_m = (1/sqrt(2^m*factorial(m)*sqrt(a/c))) * exp( -(c-a)*x.^2 ) .* vHm;

end

function a=hermite(n,x)
%generate the nth hermite polynomial of x
m = 0:floor(n/2);
[q1,q2] = size(m);
s = ndims(x);
[g1, g2] = size(x);
X = repmat(x,[ones(1,s), q2]);
m = permute(repmat(m,[g1,1,g2]),[1,3,2]);
a = factorial(n)*sum((-1).^m./(factorial(m).*factorial(n-2*m)).*(2*X).^(n-2*m),3);
end