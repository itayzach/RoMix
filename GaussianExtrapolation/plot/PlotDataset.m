function fig = PlotDataset(sSimParams, x, actualDataDist, pltTitle, GMModel, nGmmPoints, plt2Title, windowStyle)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

[n, dim] = size(x);

fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
%% GMM
if exist('GMModel', 'var')
    subplot(1,2,2)
    [xGmm,compIdx] = random(GMModel, nGmmPoints);
    if dim == 1
        scatter(xGmm, zeros(1,nGmmPoints), 50, ones(1,nGmmPoints), 'filled')
        xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
        set(gca,'YTick',[],'FontSize', 14);
    elseif dim == 2
        scatter(xGmm(:,1),xGmm(:,2),[],compIdx,'filled') % Scatter plot with points of size 10
        colormap(parula(GMModel.NumComponents));
        colorbar();
        xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
        set(gca,'FontSize', 14);
    elseif dim == 3
        scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[],compIdx, 'filled');
        xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$y$', 'interpreter', 'latex', 'FontSize', 16);
        zlabel('$z$', 'interpreter', 'latex', 'FontSize', 16);
%         view(10,5);
        view(30,70);
        set(gca,'FontSize', 14);
    end
    title(plt2Title, 'Interpreter', 'latex', 'FontSize', 14)
    subplot(1,2,1)
end
%% Dataset
if dim == 1
    scatter(x, zeros(1,n), 50, ones(1,n), 'filled')
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
%     view(10,5);
    view(30,70);
    set(gca,'FontSize', 14);
end
title(strcat(pltTitle, " (", actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)

%% Size
if strcmp(windowStyle, 'normal')
    x0     = 400;
    y0     = 400;
    height = 400;
    width  = 1200;
    set(gcf,'Position', [x0 y0 width height])
end
set(0,'DefaultFigureWindowStyle',prevWindowStyle)
%% Save
if isfield(sSimParams, 'outputFolder')
    if ~exist(sSimParams.outputFolder, 'dir')
        mkdir(sSimParams.outputFolder)
    end
    
    saveas(fig,strcat(sSimParams.outputFolder, filesep, actualDataDist, num2str(dim), 'd', '_histogram'), 'epsc');
end
end