function [] = PlotEigenfunctionsEigenvectors(sParams, sSimParams, mPhi_K, mPhi_A)

fig = figure;
if sParams.dim == 1
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
        plot(sParams.x_rand, mPhi_A(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
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
    sgtitle('Kernel (analytic) Eigenfunctions')
    x1 = x.';
    x2 = x.';
    [mX1, mX2] = meshgrid(x1, x2);
    for m = 0:sParams.PlotEigenFuncsM-1  
        vPhi_m_x1 = phi(sParams, m, x1, 1);
        vPhi_m_x2 = phi(sParams, m, x2, 2);

        % outter product since phi(x1,x2)=phi(x1)phi(x2)
        mPhi_m_x1x2 = vPhi_m_x1.' * vPhi_m_x2; 

        subplot(2,1,m+1);
        surf(mX1, mX2, mPhi_m_x1x2, 'edgecolor', 'none')
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
    %% 2-D Numeric
    subplot(2,1,2);
    sgtitle(sprintf('Eigenvectors of (numeric) A; n = %d', n))
    for m = 0:sParams.PlotEigenFuncsM-1
        subplot(2,2,m+1);
        surf(mX1, mX2, mPhi_Am, 'edgecolor', 'none')
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
    
    %% Save
    set(gcf,'Position',[100 100 600 800])
    %     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_2d.png']);
else
    error('Not supporting plots for more than 2-D')
end

end