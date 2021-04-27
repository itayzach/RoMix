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
    nanInd = find(isnan(xTilde(:,d)));
    assert(~any(isnan(xTilde(:,d))), ...
        ['xTilde contain NaNs because of x(', num2str(nanInd), ',', num2str(d), ') = ', num2str(x(nanInd,d)), ...
        newline '    ---> polyCdf(', num2str(nanInd), ',', num2str(d), ') = ', num2str(polyCdf(nanInd))]);
    
    infInd = find(isinf(xTilde(:,d)));
    assert(~any(isinf(xTilde(:,d))), ...
        ['xTilde contain infs because of x(', num2str(infInd), ',', num2str(d), ') = ', num2str(x(infInd,d)), ...
        newline '    ---> polyCdf(', num2str(infInd), ',', num2str(d), ') = ', num2str(polyCdf(infInd))]);
end
assert(~any(isnan(xTilde(:))),['xTilde contain NaNs because of x = ', num2str(x(isnan(xTilde(:)))')]);
assert(~any(isinf(xTilde(:))),['xTilde contain Infs because of x = ', num2str(x(isinf(xTilde(:)))')]);
end