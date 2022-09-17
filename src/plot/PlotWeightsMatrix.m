function PlotWeightsMatrix(sPlotParams, W, dist, D, L, xTrain, adjacencyType, omega, k)

if ~isempty(dist)
    figure('Name', 'Distance matrix');
    imagesc(dist);
    colorbar('TickLabelInterpreter', 'latex');
    title(['Distance matrix', newline, ...
        'omega = ', num2str(omega, '%.2f'), ', max(dist) = ', num2str(max(dist(:)), '%.2f')], ...
        'FontSize', 14, 'interpreter', 'latex')
    set(gca,'FontSize', 14);
end

figure('Name', 'Weight matrix');
imagesc(W);
colorbar('TickLabelInterpreter', 'latex');
title('Weight matrix', 'FontSize', 14, 'interpreter', 'latex')
set(gca,'FontSize', 14);

if size(D,1) == size(D,2)
    d = diag(D); % D is a matrix
else
    d = D; % D is a vector
end
figure('Name', 'Node degrees');
plot(d, '.', 'MarkerSize', 10);
title('Node degrees', 'FontSize', 14, 'interpreter', 'latex')
set(gca,'FontSize', 14);

figure('Name', 'Laplacian matrix');
imagesc(L);
colorbar('TickLabelInterpreter', 'latex');
title('$L_n = I - D^{-0.5}WD^{-0.5}$', 'FontSize', 14, 'interpreter', 'latex')
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
    if isempty(sPlotParams)
        actualDataDist = [];
    else
        actualDataDist = sPlotParams.actualDataDist;
    end
    PlotEigenfuncvecScatter(sPlotParams, actualDataDist, xTrain, [], 0, 3, W(:,1:4), [], [], [], ...
        figTitle, figName, 'w');
end
end