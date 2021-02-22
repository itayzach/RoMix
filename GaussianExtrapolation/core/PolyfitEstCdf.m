function [xGrid, estMarginalCdf, pCdf, invpCdf] = PolyfitEstCdf(x, bins, pCdfDegree, invpCdfDegree)
[n, dim] = size(x);

%% Estimate d-Dim PDF using histcn (TODO: check if runtime can be improved)
[pdf_xGrid, xGridCells, ~, ~] = histcn(x, bins*ones(1,dim)-1);
pdf_xGrid = (1/n)*pdf_xGrid;

%% Initialize
xGrid = zeros(bins(1),dim);
estMarginalCdf = zeros(bins(1),dim);
pCdf = zeros(dim,pCdfDegree+1);
invpCdf = zeros(dim,invpCdfDegree+1);
%% Estimate the marginal CDFs for each dimension
for d = 1:dim
    xGrid(:,d) = (xGridCells{d})';
    pdf_1d = pdf_xGrid;
    for k = 1:dim
        if k == d % the marginal pdf comes from summing over all dimensions other than the k'th
            continue;
        end
        pdf_1d = sum(pdf_1d,k);
    end
    estMarginalCdf(:,d) = cumsum(pdf_1d);  % CDF_x1
    
    % to have analytic expression for the cdf:
    pCdf(d,:) = polyfit(xGrid(:,d), estMarginalCdf(:,d), pCdfDegree); 
    invpCdf(d,:) = polyfit(estMarginalCdf(:,d), xGrid(:,d), invpCdfDegree);
end
end