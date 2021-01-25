function xTilde = T(pCdf, b_saturate, mu, sigma, x)
polyCdf = polyval(pCdf, x);
if (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1]...');
    if b_saturate
        polyCdf(polyCdf > 0.99) = 0.99; % saturate
        polyCdf(polyCdf < 0.01) = 0.01; % saturate
    end
end
xTilde = icdf('Normal',polyCdf,mu,sigma);
assert(~any(isnan(xTilde)),['xTilde contain NaNs since because of x = ', num2str(x(isnan(xTilde))')]);
end