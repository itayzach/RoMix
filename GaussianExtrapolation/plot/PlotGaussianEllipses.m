function PlotGaussianEllipses(sPlotParams, sDistParams)
dim = numel(sDistParams.mu{1});
if dim == 2
    fig = figure;
    gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
    xMax = 1.5*max(cell2mat(sDistParams.sigma')) + max(cell2mat(sDistParams.mu'));
    xMin = 1.5*min(cell2mat(sDistParams.sigma')) + min(cell2mat(sDistParams.mu'));
    fcontour(gmPDF, [xMin(1), xMax(1), xMin(2), xMax(2)]);
elseif dim == 3
    fig = figure;
    cmap = colormap(jet(sDistParams.GMModel.NumComponents));
    for c = 1:sDistParams.GMModel.NumComponents
        nPoints = 20;
        mu = sDistParams.mu{c};
        sigma = sDistParams.sigma{c};
        u = sDistParams.u{c};
        cov = sDistParams.cov{c};

        [elv, elf] = ellipsedata3(cov, mu, nPoints, 1.5);
        patch('Faces',elf,'Vertices',elv,'FaceColor',cmap(c,:),'EdgeColor','none','FaceAlpha',.5);
        hold on;
        quiver3(repmat(mu(1),3,1), repmat(mu(2),3,1), repmat(mu(3),3,1), ...
            (sigma.*u(1,:))', (sigma.*u(2,:))', (sigma.*u(3,:))','k','LineWidth',5);
    end
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
    colorbar;
    grid on;
    view(30,75);
    set(gca,'FontSize', 14);
end

%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'Ellipses';
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end