function [] = PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A)

if sParams.dim == 1
    fig = figure;
    %% 1-D Analytic
    subplot(2,1,1);
    for m = 0:sParams.PlotEigenFuncsM-1
        plot(sParams.x, mPhi_K(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlim([sParams.xMin sParams.xMax]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sParams.PlotEigenFuncsM - 1
            vP_x = p(sParams, sParams.x, 1);
            plot(sParams.x, vP_x, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '$p(x)$');
            hold off
            title('Eigenfunctions of (analytic) kernel operator')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
        end
    end
    
    %% 1-D Numeric
    subplot(2,1,2);
    for m = 0:sParams.PlotEigenFuncsM-1
        plot(sParams.x_rand, mPhi_A(:,m+1), 'o', 'LineWidth', 1, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlim([sParams.xMin sParams.xMax]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sParams.PlotEigenFuncsM - 1
            histogram(sParams.x_rand, 'Normalization', 'pdf', 'LineStyle', ':', 'DisplayStyle', 'stairs', 'LineWidth', 2, 'DisplayName', '$\hat{p}(x)$');
            hold off
            title('Eigenvectors of (numeric) Gram matrix A')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
        end
    end
    %% Save
    set(gcf,'Position',[100 100 600 800])
    saveas(fig,[sParams.sSim.outputFolder filesep 'fig_eigenvectors_1d.png']);
elseif sParams.dim == 2
    %% 2-D Analytic
    x0     = 10;
    y0     = 250;
    width  = 1000;
    height = 700;
    fig = figure;
%     sgtitle('Kernel (analytic) Eigenfunctions')
    [mX1, mX2] = meshgrid(sParams.x(:,1), sParams.x(:,2));
    for m = 0:sParams.PlotEigenFuncsM-1  
        if m == 0
            vP_x = p(sParams, [mX1(:) mX2(:)]);
            mP_x = reshape(vP_x, size(sParams.x,1), size(sParams.x,1));
            subplot(ceil(sqrt(sParams.PlotEigenFuncsM+1)),ceil(sqrt(sParams.PlotEigenFuncsM+1)),1);
            if sParams.sSim.b_plot_contourf
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
            xlim([ min(sParams.x(:,1)) max(sParams.x(:,1))])
            ylim([ min(sParams.x(:,2)) max(sParams.x(:,2))])
            title('$p(x)$', 'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end
        
        subplot(ceil(sqrt(sParams.PlotEigenFuncsM+1)),ceil(sqrt(sParams.PlotEigenFuncsM+1)),m+2);
        if sParams.sSim.b_plot_contourf
            contourf(mX1, mX2, reshape(mPhi_K(:,m+1), size(sParams.x,1), size(sParams.x,1)));
        else
            surf(mX1, mX2, reshape(mPhi_K(:,m+1), size(sParams.x,1), size(sParams.x,1)), 'edgecolor', 'none');
        end
        colormap(gca, 'default')
        colorbar()
        view(2)
        if m == 0
            hold on
            mPhi_0 = reshape(mPhi_K(:,1), size(sParams.x,1), size(sParams.x,1));
            maxIdx = find(imregionalmax(mPhi_0));
            plot3(mX1(maxIdx), mX2(maxIdx), mPhi_0(maxIdx), 'ro', 'MarkerSize', 5)
            text(mX1(maxIdx), mX2(maxIdx), mPhi_0(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
        end
        xlim([ min(sParams.x(:,1)) max(sParams.x(:,1))])
        ylim([ min(sParams.x(:,2)) max(sParams.x(:,2))])
        title(['$\phi_{' num2str(m) '}(x),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda_K(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
    end
    set(gcf,'Position', [x0 y0 width height])
    saveas(fig,[sParams.sSim.outputFolder filesep 'fig_eigenfunctions_2d.png']);
    
    %% 2-D Numeric
    if ~isempty(mPhi_A)
        fig = figure;
%         sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', size(sParams.x_rand,1)))
        for m = 0:sParams.PlotEigenFuncsM-1
            if m == 0
                subplot(ceil(sqrt(sParams.PlotEigenFuncsM+1)),ceil(sqrt(sParams.PlotEigenFuncsM+1)),1);
                hist3(sParams.x_rand,'CdataMode','auto', 'Nbins', [30 30], 'edgecolor', 'flat');
                colormap(gca, 'hot')
                colorbar()
                view(2)
                xlim([ min(sParams.x(:,1)) max(sParams.x(:,1))])
                ylim([ min(sParams.x(:,2)) max(sParams.x(:,2))])
                title('$\hat{p}(x)$', 'Interpreter', 'latex', 'FontSize', 14)
                set(gca,'FontSize', 14);
%                 hold on
%                 [ ~, maxIdx] = max(mHist3(:));
%                 plot3(sParams.x(maxIdx,1), sParams.x(maxIdx,2), mP_x(maxIdx), 'ro', 'MarkerSize', 5)
%                 text(sParams.x(maxIdx,1), sParams.x(maxIdx,2), mP_x(maxIdx), sprintf('$(%.2f, %.2f)$', mX1(maxIdx), mX2(maxIdx)), ...
%                     'VerticalAlignment', 'bottom', 'Interpreter', 'latex', 'FontSize', 14)
            end
            
            subplot(ceil(sqrt(sParams.PlotEigenFuncsM+1)),ceil(sqrt(sParams.PlotEigenFuncsM+1)),m+2);
            scatter3(sParams.x_rand(:,1), sParams.x_rand(:,2), mPhi_A(:,m+1), [], mPhi_A(:,m+1), 'filled');
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
            xlim([ min(sParams.x(:,1)) max(sParams.x(:,1))])
            ylim([ min(sParams.x(:,2)) max(sParams.x(:,2))])
            title(['$\hat{\phi}_{' num2str(m) '}(x),$ $\hat{\lambda}_{' num2str(m)  '} = ' num2str(vLambda_A(m+1), '%.4f') '$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end

        set(gcf,'Position', [x0+width y0 width height])
        saveas(fig,[sParams.sSim.outputFolder filesep 'fig_eigenvectors_2d.png']);
    end
else
    error('Not supporting plots for more than 2-D')
end

end