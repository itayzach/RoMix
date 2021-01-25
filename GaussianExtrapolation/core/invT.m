function x = invT(invpCdf, mu, sigma, xTilde)
assert(~any(isnan(xTilde)),'xTilde contain NaNs...');
xTildeCdf = cdf('Normal', xTilde, mu, sigma);
x = polyval(invpCdf, xTildeCdf);
assert(~any(isnan(x)),'x contain NaNs...');
end