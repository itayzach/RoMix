function [ res ] = hermite (n, x)
if n == 0
    res = ones(size(x,1),1);
    return
end
%% Get Hermite coefficients
h = zeros(1,n+1);
n_fact = factorial(n);
m = 0:floor(n/2);
h(2*m+1) = n_fact .* (-1).^m ./ (factorial(m) .* factorial(n-2.*m)) .* 2.^(n-2.*m);

%% Calculate polynomial
P = length(h)-1;
p = P:-1:0;
res = sum(h.*(x.^p),2);
% isalmostequal(res,polyval(h, x),1e-13)

end
