function fig = PlotDataset(sSimParams, x, y, pltTitle, GMModel, nGmmPoints, plt2Title, windowStyle)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

[n, dim] = size(x);
fig = figure('Name', sprintf('%d-D %s', dim, sSimParams.sDataset.actualDataDist));
%% GMM
if exist('GMModel', 'var')
    ax(1) = subplot(1,2,2);
    [xGmm,compIdx] = random(GMModel, nGmmPoints);
    xMax = max([max(xGmm); max(x)]);
    xMin = min([min(xGmm); min(x)]);
    if dim == 1
        scatter(xGmm, zeros(1,nGmmPoints), 50, compIdx, 'filled')
        xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
        set(gca,'YTick',[],'FontSize', 14);
        xlim([xMin(1), xMax(1)])
    elseif dim == 2
        scatter3(xGmm(:,1),xGmm(:,2),compIdx,[],compIdx,'filled');
        colormap(ax(1), jet(GMModel.NumComponents));
        colorbar();
        xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
        view(2);
        set(gca,'FontSize', 14);
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
    elseif dim == 3
        scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[],compIdx, 'filled');
        xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
        ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
        zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
%         view(10,5);
        view(30,70);
        set(gca,'FontSize', 14);
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
        zlim([xMin(3), xMax(3)])
    end
    title(plt2Title, 'Interpreter', 'latex', 'FontSize', 14)
    ax(2) = subplot(1,2,1);
end
%% Dataset
if dim == 1
    scatter(x, zeros(1,n), 50, ones(1,n), 'filled')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
    xlim([xMin(1), xMax(1)])
elseif dim == 2
    if ~isempty(y) && isequal(y, floor(y))
        scatter3(x(:,1), x(:,2), y, [], y, 'filled');
        if exist('ax', 'var')
            colormap(ax(2), jet(length(unique(y)))); 
        else
            colormap(jet(length(unique(y))));
        end
        colorbar;
    else
        scatter(x(:,1), x(:,2), 'filled');
    end
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    view(2);
    set(gca,'FontSize', 14);
    if exist('xMin', 'var') && exist('xMax', 'var')
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
    end
elseif dim == 3
    scatter3(x(:,1), x(:,2), x(:,3), 'filled');
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
%     view(10,5);
    view(30,70);
    set(gca,'FontSize', 14);
    if exist('xMin', 'var') && exist('xMax', 'var')
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
        zlim([xMin(3), xMax(3)])
    end
end
title(strcat(pltTitle, " (", sSimParams.sDataset.actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)

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
    
    saveas(fig,strcat(sSimParams.outputFolder, filesep, sSimParams.sDataset.actualDataDist, num2str(dim), 'd', '_dataset'), 'epsc');
end
end