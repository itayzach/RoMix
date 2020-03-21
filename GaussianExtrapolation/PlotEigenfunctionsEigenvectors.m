function [] = PlotEigenfunctionsEigenvectors(sParams, sSimParams, mPhi_K, mPhi_A)

if sParams.dim == 1
    fig = figure;
    %% 1-D Analytic
    subplot(2,1,1);
    for m = 0:sParams.PlotEigenFuncsM-1
        plot(sParams.x, mPhi_K(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sParams.PlotEigenFuncsM - 1
            vP_x = p(sParams, sParams.x, 1);
            plot(sParams.x, vP_x, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '$p(x)$');
            hold off
            title('Eigenfunctions of (analytic) kernel operator')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
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
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
        end
    end
    %% Save
    set(gcf,'Position',[100 100 600 800])
    %     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_1d.png']);
elseif sParams.dim == 2
    %% 2-D Analytic
    fig = figure;
    sgtitle('Kernel (analytic) Eigenfunctions')
    [mX1, mX2] = meshgrid(sParams.x(:,1), sParams.x(:,2));
    for m = 0:sParams.PlotEigenFuncsM-1  
        subplot(3,3,m+1);
        surf(mX1, mX2, reshape(mPhi_K(:,m+1), size(sParams.x,1), size(sParams.x,1)), 'edgecolor', 'none');
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
        set(gcf,'Position',[10 250 1900 700])
    %     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_2d.png']);
    %% 2-D Numeric
    fig = figure;
    sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', size(sParams.x_rand,1)))
%     [mX1_rand, mX2_rand] = meshgrid(sParams.x_rand(:,1), sParams.x_rand(:,2));
    for m = 0:sParams.PlotEigenFuncsM-1
%         subplot(2,sParams.PlotEigenFuncsM,sParams.PlotEigenFuncsM+m+1);
        subplot(3,3,m+1);
%         surf(mX1_rand, mX2_rand, mPhi_A(:,m+1), 'edgecolor', 'none')
        scatter3(sParams.x_rand(:,1), sParams.x_rand(:,2), mPhi_A(:,m+1), [], mPhi_A(:,m+1), 'filled');
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
    
    set(gcf,'Position',[10 250 1900 700])
    %     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_2d.png']);
else
    error('Not supporting plots for more than 2-D')
end

end