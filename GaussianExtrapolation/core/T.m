function xTilde = T(pCdf, b_saturate, mu, sigma, x)
dim = size(x,2);
xTilde = zeros(size(x));
for d = 1:dim
    polyCdf = polyval(pCdf(d,:), x(:,d));
    if (any(polyCdf > 1) || any(polyCdf < 0))
        warning('CDF must be in [0,1]...');
        if b_saturate
            polyCdf(polyCdf > 0.99) = 0.99; % saturate
            polyCdf(polyCdf < 0.01) = 0.01; % saturate
        end
    end
    xTilde(:,d) = icdf('Normal',polyCdf,mu(d),sigma(d,d));
end
assert(~any(isnan(xTilde(:))),['xTilde contain NaNs since because of x = ', num2str(x(isnan(xTilde(:)))')]);
end