function PlotWeightsMatrix(sPlotParams, W, xTrain, origGraphAdjacency, verticesPDF, omega, k)
figure('Name', 'Weight matrix');
imagesc(W);
colorbar();
title('Weight matrix', 'FontSize', 14, 'interpreter', 'latex')
set(gca,'FontSize', 14);

n = length(W);
title('Weights matrix', 'interpreter', 'latex', 'FontSize', 14)
if strcmp(origGraphAdjacency, 'GaussianKernel')
    figTitle = [ '${\bf W}(:,i)$ (Gaussian kernel) for first $i=1:4$ nodes with $\omega = ' num2str(omega) ...
        '\quad n = ' num2str(n) '$'];
elseif strcmp(origGraphAdjacency, 'NearestNeighbor')
    figTitle = [ '${\bf W}(:,i)$ (k-NN) for first $i=1:4$ nodes with $k = ' num2str(k) ...
        '\quad n = ' num2str(n) '$'];
end
dim = size(xTrain,2);
if dim <= 3
    figName = 'W';
    PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], 0, 3, W(:,1:4), [], [], [], ...
        figTitle, figName, 'w')
end
end