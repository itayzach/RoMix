function [xGrid, estMarginalCdf, pCdf, invpCdf] = PolyfitEstCdf(x, nEvalPoints, pCdfDegree, invpCdfDegree, b_kde, b_plotCdf)
[n, dim] = size(x);
if ~exist('b_plotCdf', 'var')
    b_plotCdf = false;
end


%% Initialize
xGrid = zeros(nEvalPoints,dim);
estMarginalCdf = zeros(nEvalPoints,dim);
pCdf = zeros(dim,pCdfDegree+1);
invpCdf = zeros(dim,invpCdfDegree+1);
%% Estimate the pdf from x
if ~exist('b_kde', 'var')
    b_kde = false;
end
if b_kde
    for d = 1:dim
        xGrid(:,d) = linspace(min(x(:,d)),max(x(:,d)),nEvalPoints);
    end
    if dim == 1
        xMeshGrid = xGrid;
    elseif dim == 2
        [X, Y] = meshgrid(xGrid(:,1), xGrid(:,2));
        xMeshGrid(:,1) = X(:);
        xMeshGrid(:,2) = Y(:);
    elseif dim == 3
        [X, Y, Z] = meshgrid(xGrid(:,1), xGrid(:,2), xGrid(:,3));
        xMeshGrid(:,1) = X(:);
        xMeshGrid(:,2) = Y(:);
        xMeshGrid(:,3) = Z(:);
    end
    pdf_xGrid = mvksdensity(x, xMeshGrid, 'Function', 'pdf');
else
    [pdf_xGrid, xGridCells, ~, ~] = histcn(x, nEvalPoints*ones(1,dim)-1);
    pdf_xGrid = (1/n)*pdf_xGrid;
    for d = 1:dim
        xGrid(:,d) = (xGridCells{d})';
    end
end
%% Estimate the marginal CDFs for each dimension   
for d = 1:dim
    pdf_1d = pdf_xGrid;
    for k = 1:dim
        if k == d % the marginal pdf comes from summing over all dimensions other than the k'th
            continue;
        end
        pdf_1d = sum(pdf_1d,k);
    end
    if b_kde
        dx = (max(x(:,d))-min(x(:,d)))/nEvalPoints;
    else
        dx = 1;
    end
    estMarginalCdf(:,d) = cumsum(pdf_1d*dx);  % CDF_x1
    
    % to have analytic expression for the cdf:
    pCdf(d,:) = polyfit(xGrid(:,d), estMarginalCdf(:,d), pCdfDegree); 
    invpCdf(d,:) = polyfit(estMarginalCdf(:,d), xGrid(:,d), invpCdfDegree);
    
    if b_plotCdf
        figure; 
        plot(xGrid(:,d), estMarginalCdf(:,d), 'LineWidth', 2)
        
        
    end
end
end
