function [xTilde, polyvalCdf] = T(pCdf, b_saturate, muTilde, sigmaTilde, x, b_debugUseNormalCDF, xTrainGrid, estMarginalCdf_xTrain)
if ~exist('b_debugUseNormalCDF', 'var')
    b_debugUseNormalCDF = false;
end
dim = size(x,2);
xTilde = zeros(size(x));
for d = 1:dim
    polyvalCdf = polyval(pCdf(d,:), x(:,d));
    if b_debugUseNormalCDF
        normalCdf = cdf('Normal',x(:,d),muTilde, sigmaTilde);
        figure; 
        plot(x(:,d), normalCdf, '.');
        hold on;
        plot(xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), '.');
        plot(x(:,d), polyvalCdf, '.');
        legend('cdf(normal)', 'ecdf', 'ours', 'Location', 'southeast');
        set(gca,'FontSize', 14);
        polyvalCdf = normalCdf; % override our polyval with normal cdf
    end
    if (any(polyvalCdf > 1) || any(polyvalCdf < 0))
        warning('CDF must be in [0,1]...');
        if b_saturate
            polyvalCdf(polyvalCdf > 0.99) = 0.99; % saturate
            polyvalCdf(polyvalCdf < 0.01) = 0.01; % saturate
        end
    end
    xTilde(:,d) = icdf('Normal',polyvalCdf,muTilde(d),sigmaTilde(d,d));
    CheckForNaNandInf(xTilde, x, polyvalCdf, d)
end
end

function CheckForNaNandInf(xTilde, x, polyvalCdf, d)
nanInd = find(isnan(xTilde(:,d)));
assert(~any(isnan(xTilde(:,d))), ...
    ['xTilde contain NaNs because of x(', num2str(nanInd), ',', num2str(d), ') = ', num2str(x(nanInd,d)), ...
    newline '    ---> polyCdf(', num2str(nanInd), ',', num2str(d), ') = ', num2str(polyvalCdf(nanInd))]);

infInd = find(isinf(xTilde(:,d)));
assert(~any(isinf(xTilde(:,d))), ...
    ['xTilde contain infs because of x(', num2str(infInd), ',', num2str(d), ') = ', num2str(x(infInd,d)), ...
    newline '    ---> polyCdf(', num2str(infInd), ',', num2str(d), ') = ', num2str(polyvalCdf(infInd))]);
end