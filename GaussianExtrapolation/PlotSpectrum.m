function [] = PlotSpectrum(sParams, vLambda_K, vLambda_A)

fig = figure;
subplot(2,1,1);
stem(0:sParams.PlotSpectM-1, vLambda_K, 'DisplayName', '$\lambda_m$');
hold on;
stem(0:sParams.PlotSpectM-1, vLambda_A, 'DisplayName', '$\hat{\lambda}_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Analytic $(\lambda_m)$ vs. Numeric $(\hat{\lambda}_m)$ eigenvalues', ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

subplot(2,1,2);
stem(0:sParams.PlotSpectM-1, log(vLambda_K),'LineStyle','none', 'DisplayName', '$\log(\lambda_m)$');
hold on;
stem(0:sParams.PlotSpectM-1, log(vLambda_A),'LineStyle','none', 'DisplayName', '$\log(\hat{\lambda}_m)$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% print(fig, [sParams.sSim.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], '-depsc')
title('Analytic $(\lambda_m)$ vs. Numeric $(\hat{\lambda}_m)$ eigenvalues', ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(gcf,'Position',[100 100 600 500])
saveas(fig,[sParams.sSim.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], 'epsc');
end