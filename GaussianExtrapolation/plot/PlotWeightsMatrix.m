function PlotWeightsMatrix(sPlotParams, W, dist, xTrain, adjacencyType, verticesPDF, omega, k)

figure('Name', 'Distance matrix');
imagesc(dist);
colorbar();
title(['Distance matrix', newline, ...
    'omega = ', num2str(omega, '%.2f'), ', max(dist) = ', num2str(max(dist(:)), '%.2f')], ...
    'FontSize', 14, 'interpreter', 'latex')
set(gca,'FontSize', 14);

figure('Name', 'Weight matrix');
imagesc(W);
colorbar();
title('Weight matrix', 'FontSize', 14, 'interpreter', 'latex')
set(gca,'FontSize', 14);


[n, dim] = size(xTrain);
if dim <= 3
    if strcmp(adjacencyType, 'GaussianKernel')
        figTitle = [ '${\bf W}(:,i)$ (Gaussian kernel) for first $i=1:4$ nodes with $\omega = ' num2str(omega) ...
            '\quad n = ' num2str(n) '$'];
    elseif strcmp(adjacencyType, 'NearestNeighbor')
        figTitle = [ '${\bf W}(:,i)$ (k-NN) for first $i=1:4$ nodes with $k = ' num2str(k) ...
            '\quad n = ' num2str(n) '$'];
    end

    figName = 'W';
    PlotEigenfuncvecScatter(sPlotParams, verticesPDF, xTrain, [], 0, 3, W(:,1:4), [], [], [], ...
        figTitle, figName, 'w')
end
end