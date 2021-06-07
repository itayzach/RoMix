function x = Tinv(invpCdf, mu, sigma, xTilde)
assert(~any(isnan(xTilde(:))),'xTilde contain NaNs...');
dim = size(xTilde,2);
x = zeros(size(xTilde));
for d = 1:dim
    xTildeCdf = cdf('Normal', xTilde(:,d), mu(d), sigma(d,d));
    x(:,d) = polyval(invpCdf(d,:), xTildeCdf);
end
assert(~any(isnan(x(:))),'x contain NaNs...');
end