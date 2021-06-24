function fig = PlotGMM(pltTitle, GMModel, nPoints)

[x,compIdx] = random(GMModel, nPoints);
dim = size(x,2);
fig = figure('Name', 'GMM');
if dim == 1
    scatter(x, zeros(1,nPoints), 50, ones(1,nPoints), 'filled')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
elseif dim == 2
    scatter(x(:,1),x(:,2),[],compIdx,'filled') % Scatter plot with points of size 10
    colormap(parula(GMModel.NumComponents));
    colorbar();
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
elseif dim == 3
    scatter3(x(:,1), x(:,2), x(:,3),[],compIdx, 'filled');
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$z$', 'interpreter', 'latex', 'FontSize', 16);
    view(10,5);
    set(gca,'FontSize', 14);
end

title(pltTitle, 'interpreter', 'latex', 'FontSize', 16);
end