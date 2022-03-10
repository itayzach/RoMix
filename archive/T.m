function [xTilde, polyvalCdfMatrix] = T(pCdf, b_saturate, muTilde, sigmaTilde, x, xTrain, b_kde, kde_bw)
if ~exist('b_kde', 'var')
    b_kde = false;
end
dim = size(x,2);

if b_kde
    xTilde = zeros(size(x));
    polyvalCdfMatrix = zeros(size(x));
%     kdeCdf = mvksdensity(xTrain, x, 'Function', 'cdf','bandwidth', 1e-6);
%     kdePdf = mvksdensity(xTrain, x, 'Function', 'pdf','bandwidth', 0.7);
    for d = 1:dim
        if exist('kde_bw', 'var') && ~isempty(kde_bw) && kde_bw > 0
            kdeCdf = ksdensity(xTrain(:,d), x(:,d), 'Function', 'cdf', 'Bandwidth', kde_bw);
        else
            kdeCdf = ksdensity(xTrain(:,d), x(:,d), 'Function', 'cdf');
        end
        xTilde(:,d) = icdf('Normal',kdeCdf,muTilde(d),sigmaTilde(d,d));
        polyvalCdfMatrix(:,d) = kdeCdf;
    end
else    
    xTilde = zeros(size(x));
    polyvalCdfMatrix = zeros(size(x));
    for d = 1:dim
        polyvalCdf = polyval(pCdf(d,:), x(:,d));
        
        if any(polyvalCdf > 1)
            warning([ num2str(sum(polyvalCdf > 1)), ' elements are > 1: ', newline, ...
                '       --->   ', num2str(polyvalCdf(polyvalCdf > 1)')]);
            if b_saturate
                warning('Saturating...')
                polyvalCdf(polyvalCdf > 1) = 1-1e-10; % saturate
                %       polyvalCdf(polyvalCdf > 0.99) = 0.99; % saturate
            end
        end
        if (any(polyvalCdf < 0))
            warning([ num2str(sum(polyvalCdf < 0)), ' elements are < 0: ', newline, ...
                '       --->   ', num2str(polyvalCdf(polyvalCdf < 0)')]);
            if b_saturate
                warning('Saturating...')
                polyvalCdf(polyvalCdf < 0) = 1e-10; % saturate
                %       polyvalCdf(polyvalCdf < 0.01) = 0.01; % saturate
            end
        end
        xTilde(:,d) = icdf('Normal',polyvalCdf,muTilde(d),sigmaTilde(d,d));
        CheckForNaNandInf(xTilde(:,d), x(:,d), polyvalCdf, d)
        polyvalCdfMatrix(:,d) = polyvalCdf;
    end
end
end

function CheckForNaNandInf(xTilde, x, polyvalCdf, d)
nanInd = find(isnan(xTilde));
assert(~any(isnan(xTilde)), ...
    ['xTilde contain NaNs because of x(', num2str(nanInd'), ',', num2str(d), ') = ', num2str(x(nanInd)'), ...
    newline '    ---> polyCdf(', num2str(nanInd'), ',', num2str(d), ') = ', num2str(polyvalCdf(nanInd)')]);

infInd = find(isinf(xTilde));
assert(~any(isinf(xTilde)), ...
    ['xTilde contain infs because of x(', num2str(infInd), ',', num2str(d), ') = ', num2str(x(infInd)), ...
    newline '    ---> polyCdf(', num2str(infInd), ',', num2str(d), ') = ', num2str(polyvalCdf(infInd))]);
end