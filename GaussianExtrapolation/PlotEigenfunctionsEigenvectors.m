function [] = PlotEigenfunctionsEigenvectors(sSimParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A)

if sSimParams.dim == 1
    fig = figure;
    %% 1-D Analytic
    subplot(2,1,1);
    for m = 0:sSimParams.PlotEigenFuncsM-1
        plot(sSimParams.x, mPhi_K(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlim([sSimParams.xMin sSimParams.xMax]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sSimParams.PlotEigenFuncsM - 1
            vP_x = p(sSimParams, sSimParams.x, 1);
            plot(sSimParams.x, vP_x, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '$p(x)$');
            hold off
            title('Eigenfunctions of (analytic) kernel operator')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
        end
    end
    
    %% 1-D Numeric
    subplot(2,1,2);
    for m = 0:sSimParams.PlotEigenFuncsM-1
        plot(sSimParams.x_rand, mPhi_A(:,m+1), 'o', 'LineWidth', 1, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlim([sSimParams.xMin sSimParams.xMax]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sSimParams.PlotEigenFuncsM - 1
            histogram(sSimParams.x_rand, 'Normalization', 'pdf', 'LineStyle', ':', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', '$\hat{p}(x)$');
            hold off
            title('Eigenvectors of (numeric) Gram matrix A')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
        end
    end
    %% Save
    set(gcf,'Position',[100 100 600 800])
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_1d'], 'epsc');
elseif sSimParams.dim == 2
    %% 2-D Analytic
    x0     = 10;
    y0     = 50;
    width  = 1800;
    height = 900;
%     sgtitle('Kernel (analytic) Eigenfunctions')
    [mX1, mX2] = meshgrid(sSimParams.x(:,1), sSimParams.x(:,2));
    
    nRows = floor(sqrt(sSimParams.PlotEigenFuncsM+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(sSimParams.PlotEigenFuncsM/nRows);
    
    % Plot pdf
    fig = figure;
    vP_x = p(sSimParams, [mX1(:) mX2(:)]);
    mP_x = reshape(vP_x, size(sSimParams.x,1), size(sSimParams.x,1));
    if sSimParams.b_plot_contourf
        contourf(mX1, mX2, mP_x);
    else
        surf(mX1, mX2, mP_x, 'edgecolor', 'none');
    end
    hold on
    maxIdx = find(imregionalmax(mP_x));
    plot3(mX1(maxIdx), mX2(maxIdx), mP_x(maxIdx), 'ro', 'MarkerSize', 5)
    text(mX1(maxIdx), mX2(maxIdx), mP_x(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
    colormap(gca, 'hot')
    colorbar()
    view(2)
    xlim([ min(sSimParams.x(:,1)) max(sSimParams.x(:,1))])
    ylim([ min(sSimParams.x(:,2)) max(sSimParams.x(:,2))])
    title('$p({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
    set(gcf,'Position', [x0 y0 600 400])
    saveas(fig,[sSimParams.outputFolder filesep 'fig_pdf_2d'], 'epsc');
    
    % Plot eigenfunctions
    fig = figure;
    for m = 0:sSimParams.PlotEigenFuncsM-1  
        subplot(nRows, nCols,m+1);
        if sSimParams.b_plot_contourf
            contourf(mX1, mX2, reshape(mPhi_K(:,m+1), size(sSimParams.x,1), size(sSimParams.x,1)));
        else
            surf(mX1, mX2, reshape(mPhi_K(:,m+1), size(sSimParams.x,1), size(sSimParams.x,1)), 'edgecolor', 'none');
        end
        colormap(gca, 'default')
        colorbar()
        view(2)
        if m == 0
            hold on
            mPhi_0 = reshape(mPhi_K(:,1), size(sSimParams.x,1), size(sSimParams.x,1));
%             maxIdx = find(imregionalmax(mPhi_0));
%             plot3(mX1(maxIdx), mX2(maxIdx), mPhi_0(maxIdx), 'ro', 'MarkerSize', 5)
%             text(mX1(maxIdx), mX2(maxIdx), mPhi_0(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
%                 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
        end
        xlim([ min(sSimParams.x(:,1)) max(sSimParams.x(:,1))])
        ylim([ min(sSimParams.x(:,2)) max(sSimParams.x(:,2))])
        title(['$\phi_{' num2str(m) '}({\bf x}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda_K(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
    end
    set(gcf,'Position', [x0 y0 width height])
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], 'epsc');
    
    %% 2-D Numeric
    if ~isempty(mPhi_A)
        
        % Plot pdf
        fig = figure;
        hist3(sSimParams.x_rand,'CdataMode','auto', 'Nbins', [30 30], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        xlim([ min(sSimParams.x(:,1)) max(sSimParams.x(:,1))])
        ylim([ min(sSimParams.x(:,2)) max(sSimParams.x(:,2))])
        title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
        set(gcf,'Position', [x0 y0 600 400])
        saveas(fig,[sSimParams.outputFolder filesep 'fig_histogram_2d'], 'epsc');
%                 hold on
%                 [ ~, maxIdx] = max(mHist3(:));
%                 plot3(sSimParams.x(maxIdx,1), sSimParams.x(maxIdx,2), mP_x(maxIdx), 'ro', 'MarkerSize', 5)
%                 text(sSimParams.x(maxIdx,1), sSimParams.x(maxIdx,2), mP_x(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
%                     'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
        
        
        
        fig = figure;
%         sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', size(sSimParams.x_rand,1)))
        for m = 0:sSimParams.PlotEigenFuncsM-1
            subplot(nRows, nCols,m+1);
            scatter3(sSimParams.x_rand(:,1), sSimParams.x_rand(:,2), mPhi_A(:,m+1), [], mPhi_A(:,m+1), 'filled');
            colormap(gca, 'default')
            colorbar()
            view(2)
%             if m == 0
%                 hold on
%                 [ ~, maxIdx] = max(mPhi_A(:,1));
%                 plot3(mX1(maxIdx), mX2(maxIdx), mPhi_A(maxIdx), 'ro', 'MarkerSize', 5)
%                 text(mX1(maxIdx), mX2(maxIdx), mPhi_A(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
%                     'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
%             end
            xlim([ min(sSimParams.x(:,1)) max(sSimParams.x(:,1))])
            ylim([ min(sSimParams.x(:,2)) max(sSimParams.x(:,2))])
            title(['$\hat{\phi}_{' num2str(m) '}({\bf x}),$ $\hat{\lambda}_{' num2str(m)  '} = ' num2str(vLambda_A(m+1), '%.4f') '$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end

        set(gcf,'Position', [x0 y0 width height])
        saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_2d'], 'epsc');
    end
else
    error('Not supporting plots for more than 2-D')
end

end