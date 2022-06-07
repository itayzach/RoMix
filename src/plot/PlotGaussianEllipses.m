function PlotGaussianEllipses(sPlotParams, sDistParams)
dim = numel(sDistParams.mu{1});
if dim == 2
    %fig = figure;
    %gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
    %xMax = 1.5*max(cell2mat(sDistParams.sigma')) + max(cell2mat(sDistParams.mu'));
    %xMin = 1.5*min(cell2mat(sDistParams.sigma')) + min(cell2mat(sDistParams.mu'));
    %fcontour(gmPDF, [xMin(1), xMax(1), xMin(2), xMax(2)]);
elseif dim == 3
    fig = figure;
    cmap = colormap(jet(sDistParams.GMModel.NumComponents));
    alphaVal = 0.3;
    for c = 1:sDistParams.GMModel.NumComponents
        mu = sDistParams.mu{c};
        sigma = sDistParams.sigma{c};
        u = sDistParams.u{c};
        cov = sDistParams.cov{c};
        
        nPoints = 20;
        [elv, elf] = ellipsedata3(cov, mu, nPoints, 1.5);
        patch('Faces',elf,'Vertices',elv,'FaceColor',cmap(c,:),'EdgeColor','none','FaceAlpha',alphaVal);
        hold on;
        quiver3(repmat(mu(1),3,1), repmat(mu(2),3,1), repmat(mu(3),3,1), ...
            (sigma.*u(1,:))', (sigma.*u(2,:))', (sigma.*u(3,:))','k','LineWidth',2);
    end
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
    grid on;
    view(30,75);
    set(gca,'FontSize', 14);
    set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5. 
    set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5. 
    set(fig,'renderer','Painters')
    h = colorbar('TickLabelInterpreter', 'latex');
    set(get(h,'label'),'string','$k$','interpreter','latex','Rotation',0,'FontSize', 16);
    h.Label.Position(1) = 0.5;
    h.Label.Position(2) = 0;
    h.Label.Position(3) = 0;
    h.Label.FontSize = 14;
    h.TickLabels = 1:sDistParams.GMModel.NumComponents;
    h.Ticks = (0.5:sDistParams.GMModel.NumComponents)/sDistParams.GMModel.NumComponents;
    h.Face.Texture.ColorType = 'truecoloralpha';
    pause(0.3); % I have no idea why we need this pause between the previous and next commands >:(
    h.Face.Texture.CData(4,:) = (alphaVal+0.2) * h.Face.Texture.CData(4,:);
    
    
end

%% Save
if (dim == 3) && ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'Ellipses';
    SaveFigure(sPlotParams, fig, figName, {'pdf', 'png'});
end
end