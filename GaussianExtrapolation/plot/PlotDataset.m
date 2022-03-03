function fig = PlotDataset(sPlotParams, x, y, pltTitle, sDistParams, nGmmPoints, plt2Title, windowStyle)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

[n, dim] = size(x);
if exist('sPlotParams', 'var') && ~isempty(sPlotParams)
    actualDataDist = sPlotParams.actualDataDist;
else
    actualDataDist = '';
end
fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
compIdx(:,1) = ones(n,1);
%% GMM
if exist('sDistParams', 'var')
    b_plotDistModel = true;
    b_spectclust = isfield(sDistParams,'SCcompIdx');
    tiledlayout(1,2 + b_spectclust);
    vAx(2) = nexttile(2);
    [xGmm,compIdx(:,2)] = random(sDistParams.GMModel, nGmmPoints);
    xMax = max([max(xGmm); max(x)]);
    xMin = min([min(xGmm); min(x)]);
    for i=1:1+b_spectclust
        if i == 2
            xGmm = x;
            compIdx(:,3) = sDistParams.SCcompIdx;
            vAx(3) = nexttile(3);
        end
        if dim == 1
            scatter(xGmm, zeros(1,nGmmPoints), 50, compIdx(:,i+1), 'filled')
            xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
            set(gca,'YTick',[],'FontSize', 14);
            xlim([xMin(1), xMax(1)])
        elseif dim == 2
            scatter3(xGmm(:,1),xGmm(:,2),compIdx(:,i+1),[],compIdx(:,i+1),'filled');
            colormap(jet(sDistParams.GMModel.NumComponents));
            colorbar();
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
            view(2);
    %         hold on
    %         gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
    %         fcontour(gmPDF,[xMin(1), xMax(1), xMin(2), xMax(2)],'LevelList',[1e-5:1e-4:1e-2],'LineColor','black');
            set(gca,'FontSize', 14);
            xlim([xMin(1), xMax(1)])
            ylim([xMin(2), xMax(2)])
        elseif dim == 3
            scatter3(xGmm(:,1), xGmm(:,2), xGmm(:,3),[],compIdx(:,i+1), 'filled');
%             hold on;
%             % eigenvectors of covariance matrix
%             for c = 1:sDistParams.GMModel.NumComponents 
%                 quiver3(sDistParams.mu{c}(1),sDistParams.mu{c}(2),sDistParams.mu{c}(3),...
%                     sDistParams.u{c}(1,1),sDistParams.u{c}(2,1),sDistParams.u{c}(3,1),sDistParams.sigma{c}(1),'k','LineWidth',5);
%                 quiver3(sDistParams.mu{c}(1),sDistParams.mu{c}(2),sDistParams.mu{c}(3),...
%                     sDistParams.u{c}(1,2),sDistParams.u{c}(2,2),sDistParams.u{c}(3,2),sDistParams.sigma{c}(2),'k','LineWidth',5);
%                 quiver3(sDistParams.mu{c}(1),sDistParams.mu{c}(2),sDistParams.mu{c}(3),...
%                     sDistParams.u{c}(1,3),sDistParams.u{c}(2,3),sDistParams.u{c}(3,3),sDistParams.sigma{c}(3),'k','LineWidth',5);
%             end
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
            zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
    %         view(10,5);
            view(30,70);
            colormap(jet(sDistParams.GMModel.NumComponents));
            colorbar();
            set(gca,'FontSize', 14);
            xlim([xMin(1), xMax(1)])
            ylim([xMin(2), xMax(2)])
            zlim([xMin(3), xMax(3)])
        end
        title(plt2Title, 'Interpreter', 'latex', 'FontSize', 14)
    end
    vAx(1) = nexttile(1);
    UpdateCursorDataTip(fig, vAx, compIdx);
    
end
%% Dataset
if dim == 1
    scatter(x, zeros(1,n), 50, ones(1,n), 'filled')
    xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
    set(gca,'YTick',[],'FontSize', 14);
    xlim([xMin(1), xMax(1)])
elseif dim == 2
%     if ~isempty(y) && isequal(y, floor(y)) && sum(y(:) > 0)
%         scatter3(x(:,1), x(:,2), y, [], y, 'filled');
%         if exist('ax', 'var')
%             colormap(ax(2), jet(length(unique(y)))); 
%         else
%             colormap(jet(length(unique(y))));
%         end
%         colorbar;
%     else
        scatter(x(:,1), x(:,2), 'filled');
%     end
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
    hold on;
%     scatter3(sDistParams.GMModel.mu(:,1), sDistParams.GMModel.mu(:,2), sDistParams.GMModel.mu(:,3), 100, 'r', 'filled');
%     for c = 1:sDistParams.GMModel.NumComponents 
%         quiver3(sDistParams.u{c}(:,1), sDistParams.u{c}(:,2), sDistParams.u{c}(:,3),zeros(3,1),zeros(3,1),zeros(3,1));
%         hold on;
%     end
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
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    if exist('b_plotDistModel', 'var') && b_plotDistModel
        figName = 'dataset_vs_dist';
    else
        figName = 'dataset';
    end
    set(fig,'renderer','Painters')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end