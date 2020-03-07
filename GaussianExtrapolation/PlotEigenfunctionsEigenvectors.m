function [] = PlotEigenfunctionsEigenvectors(sParams, sSimParams, mPhi_K, mPhi_A)

fig = figure;
subplot(2,1,1);
for m = 0:sParams.PlotEigenFuncsM-1            
% subplot(floor(sParams.PlotEigenFuncsM/2), floor(sParams.PlotEigenFuncsM/2)+1, pltIdx)
    plot(sParams.x, mPhi_K(:,m+1), 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
    hold on
    xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
    if m == sParams.PlotEigenFuncsM - 1
        vP_x = p(sParams, sParams.x, 1);
        % subplot(floor(sParams.PlotEigenFuncsM/2), floor(sParams.PlotEigenFuncsM/2)+1, pltIdx)
        plot(sParams.x, vP_x, 'LineStyle', ':', 'LineWidth', 2, 'DisplayName', '$p(x)$');
        hold off
        title('Eigenfunctions of (analytic) kernel operator')
        legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
%         print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
%         saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenfunctions_1d.png']);
    end
end

% fig = figure;
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
set(gcf,'Position',[100 100 600 800])
%     print(fig, [sSimParams.outputFolder filesep 'fig_eigenfunctions_1d'], '-depsc')
saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvectors_1d.png']);


end