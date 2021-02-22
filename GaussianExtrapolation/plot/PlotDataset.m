function fig = PlotDataset(sSimParams, x, actualDataDist, plt_title)

[n, dim] = size(x);

fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
if dim == 1
    scatter(x, zeros(1,n), 50, cmap, 'filled')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
elseif dim == 2
    scatter(x(:,1), x(:,2), 'filled');
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'FontSize', 14);
elseif dim == 3
    scatter3(x(:,1), x(:,2), x(:,3), 'filled');
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$z$', 'interpreter', 'latex', 'FontSize', 16);
    view(10,5);
    set(gca,'FontSize', 14);
end
title(strcat(plt_title, " (", actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end

saveas(fig,strcat(sSimParams.outputFolder, filesep, actualDataDist, num2str(dim), 'd', '_histogram'), 'epsc');
end