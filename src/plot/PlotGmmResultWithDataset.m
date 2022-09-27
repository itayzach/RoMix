function fig = PlotGmmResultWithDataset(sPlotParams, x, sDistParams, b_plotCovMean, cXAxisLabels)
[n, dim] = size(x);
if dim <= 3
    if exist('sPlotParams', 'var') && ~isempty(sPlotParams)
        actualDataDist = sPlotParams.actualDataDist;
    else
        actualDataDist = '';
    end
    
    fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
    if dim == 2
        % Data
        scatter(x(:,1), x(:,2), 'filled');
        
        % GMM
        gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
        xMax = max(cell2mat(sDistParams.mu')) + 3*max(cell2mat(sDistParams.sigma'));
        xMin = min(cell2mat(sDistParams.mu')) - 3*min(cell2mat(sDistParams.sigma'));
        mMu = cell2mat(sDistParams.mu');
        if b_plotCovMean
            hold on;
            fcontour(gmPDF, [xMin(1), xMax(1), xMin(2), xMax(2)]);
            scatter(mMu(:,1),mMu(:,2),[],'filled','red')
        end
        xlim([xMin(1), xMax(1)])
        ylim([xMin(2), xMax(2)])
    else
        if ~b_plotCovMean
            % Data
            scatter3(x(:,1), x(:,2), x(:,3), 'filled');
        else
            % Data
            scatter3(x(:,1), x(:,2), x(:,3));
            % GMM
            hold on;
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
    %         xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 16);
    %         ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 16);
    %         zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 16);
    %         grid on;
    %         view(30,80);
    %         set(gca,'FontSize', 14);
    %         set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5. 
    %         set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5. 
    %         set(fig,'renderer','Painters')
    %         h = colorbar('TickLabelInterpreter', 'latex');
    %         set(get(h,'label'),'string','$k$','interpreter','latex','Rotation',0,'FontSize', 16);
    %         h.Label.Position(1) = 0.5;
    %         h.Label.Position(2) = 0;
    %         h.Label.Position(3) = 0;
    %         h.Label.FontSize = 14;
    %         h.TickLabels = 1:sDistParams.GMModel.NumComponents;
    %         h.Ticks = (0.5:sDistParams.GMModel.NumComponents)/sDistParams.GMModel.NumComponents;
    %         h.Face.Texture.ColorType = 'truecoloralpha';
    %         pause(0.3); % I have no idea why we need this pause between the previous and next commands >:(
    %         h.Face.Texture.CData(4,:) = (alphaVal+0.2) * h.Face.Texture.CData(4,:);
        end
    end
    if exist('cXAxisLabels','var') && ~isempty(cXAxisLabels)
        if numel(cXAxisLabels) >= 1, xlabel(cXAxisLabels{1}, 'interpreter', 'latex', 'FontSize', 16); end
        if numel(cXAxisLabels) >= 2, ylabel(cXAxisLabels{2}, 'interpreter', 'latex', 'FontSize', 16); end
        if numel(cXAxisLabels) == 3, zlabel(cXAxisLabels{3}, 'interpreter', 'latex', 'FontSize', 16); end
    else
        axis off
    end
    set(gca,'FontSize', 14);
    
    %% Size
    if dim == 2
        hRatio = 2*(xMax(2)-xMin(2))/(xMax(1)-xMin(1));
    else
        hRatio = 1;
        %view(30,80);
        view(30,75);
    end
    x0     = 400;
    y0     = 400;
    height = hRatio*400;
    width  = 600;
    set(gcf,'Position', [x0 y0 width height])
    %% Save
    if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder') && exist('fig','var')
        if b_plotCovMean
            figName = 'gmm_cov_mean';
        else
            figName = 'gmm_data';
        end
        set(fig,'renderer','Painters')
        if exist('alphaVal', 'var')
            SaveFigure(sPlotParams, fig, figName, {'pdf', 'png'});
        else
            SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
        end
    end
end
end