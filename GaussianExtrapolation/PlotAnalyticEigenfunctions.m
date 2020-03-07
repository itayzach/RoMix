function [mPhi_m_x, vLambda] = PlotAnalyticEigenfunctions(sParams, sSimParams)
dx = 0.01;
x = (sParams.xMin:dx:sParams.xMax-dx).';
fig = figure;

mPhi_m_x = zeros(length(x), sParams.PlotEigenFuncsM);
vLambda = zeros(sParams.PlotEigenFuncsM, 1);

if sParams.dim == 1     
    for m = 0:sParams.PlotEigenFuncsM-1  
        mPhi_m_x(:,m+1) = phi(sParams, m, x, 1);
        vLambda(m+1) = lambda(sParams, m);
        % subplot(floor(sParams.PlotEigenFuncsM/2), floor(sParams.PlotEigenFuncsM/2)+1, pltIdx)
        plot(x, mPhi_m_x(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        if m == sParams.PlotEigenFuncsM - 1
            vP_x = p(sParams, x, 1);
            % subplot(floor(sParams.PlotEigenFuncsM/2), floor(sParams.PlotEigenFuncsM/2)+1, pltIdx)
            plot(x, vP_x, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '$p(x)$');
            hold off
            title('Kernel (analytic) Eigenfunctions')
            legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
%             print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
            saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_1d.png']);
        end
    end
elseif sParams.dim == 2
    sgtitle('Kernel (analytic) Eigenfunctions')
    x1 = x.';
    x2 = x.';
    [mX1, mX2] = meshgrid(x1, x2);
    for m = 0:sParams.PlotEigenFuncsM-1  
        vPhi_m_x1 = phi(sParams, m, x1, 1);
        vPhi_m_x2 = phi(sParams, m, x2, 2);

        % outter product since phi(x1,x2)=phi(x1)phi(x2)
        mPhi_m_x1x2 = vPhi_m_x1.' * vPhi_m_x2; 

        subplot(2,sParams.PlotEigenFuncsM/2,m+1);
        surf(mX1, mX2, mPhi_m_x1x2, 'edgecolor', 'none')
    %         view(2)
        colorbar()
        xlabel('$x_1$', 'Interpreter', 'latex')
        ylabel('$x_2$', 'Interpreter', 'latex')
        zlabel(['$\phi_' num2str(m) '(x_1,x_2)$'], 'Interpreter', 'latex')
    end
%     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_2d'], '-depsc')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_2d.png']);
else
    error('cannot plot for dim > 2')
end
end