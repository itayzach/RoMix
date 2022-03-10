function [xGrid, estMarginalCdf, pCdf, invpCdf] = PolyfitEstCdf(x, nEvalPoints, pCdfDegree, invpCdfDegree, b_plotCdf)
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
[pdf_xGrid, xGridCells, ~, ~] = histcn(x, nEvalPoints*ones(1,dim)-1);
pdf_xGrid = (1/n)*pdf_xGrid;
for d = 1:dim
    xGrid(:,d) = (xGridCells{d})';
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
    dx = 1; % (max(x(:,d))-min(x(:,d)))/nEvalPoints;
    estMarginalCdf(:,d) = cumsum(pdf_1d*dx);  % CDF_x1
    
    % to have analytic expression for the cdf:
    pCdf(d,:) = polyfit(xGrid(:,d), estMarginalCdf(:,d), pCdfDegree); 
    invpCdf(d,:) = polyfit(estMarginalCdf(:,d), xGrid(:,d), invpCdfDegree);
    
    if b_plotCdf
        figure('Name', 'Est. marginal CDF'); 
        plot(xGrid(:,d), estMarginalCdf(:,d), '.', 'LineWidth', 2)
        title(['Estimated marginal CDF (d = ', num2str(d), ')'], 'interpreter', 'latex', 'fontsize', 14)
        set(gca, 'FontSize', 14)
    end
end
end
