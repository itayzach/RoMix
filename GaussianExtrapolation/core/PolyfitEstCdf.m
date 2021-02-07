function [xGrid, estMarginalCdf, pCdf, invpCdf] = PolyfitEstCdf(x, bins, pCdfDegree, invpCdfDegree)
[hist_xGrid, xGridCells, mid, loc] = histcn(x, bins-1);

[n, dim] = size(x);
pdf_xGrid = (1/n)*hist_xGrid;

xGrid = zeros(n,dim);
estMarginalCdf = zeros(n,dim);
pCdf = zeros(dim,pCdfDegree+1);
invpCdf = zeros(dim,invpCdfDegree+1);
for d = 1:dim
    xGrid(:,d) = (xGridCells{d})';
    pdf_1d = pdf_xGrid;
    for k = 1:dim
        if k == d
            continue;
        end
        pdf_1d = sum(pdf_1d,k);
    end
    estMarginalCdf(:,d) = cumsum(pdf_1d);  % CDF_x1

    pCdf(d,:) = polyfit(xGrid(:,d), estMarginalCdf(:,d), pCdfDegree); % to have analytic expression for the cdf
    invpCdf(d,:) = polyfit(estMarginalCdf(:,d), xGrid(:,d), invpCdfDegree); % to have analytic expression for the cdf
end
end