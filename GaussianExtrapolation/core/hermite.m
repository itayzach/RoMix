function [ res ] = hermite (n, x)
h = zeros(1,n+1);
n_fact = factorial(n);
m = 0:floor(n/2);
h(2*m+1) = n_fact .* (-1).^m ./ (factorial(m) .* factorial(n-2.*m)) .* 2.^(n-2.*m);

% res = polyval(h, x);

%   [ x1^P, x1^(P-1), ..., x1^0 ][ h(P)   ]
%   [ x2^P, x2^(P-1), ..., x2^0 ][ h(P-1) ]
%   [ ... ,   ...   , ..., ...  ][  ...   ]
%   [ xN^P, x2^(P-1), ..., x2^0 ][ h(0)   ]
P = length(h)-1;
Xp = zeros(length(x), P+1);
for p=P:-1:0
    Xp(:,P-p+1)=x(:).^p;
end
res = Xp*h(:);
% isalmostequal(res,polyval(h, x),1e-13)
end
