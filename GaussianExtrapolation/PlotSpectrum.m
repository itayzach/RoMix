function [] = PlotSpectrum(sParams, vLambda_K, vLambda_A)

fig = figure;
subplot(2,1,1);
stem(0:sParams.PlotSpectM-1, vLambda_K, 'DisplayName', '$\lambda_m(K)$');
hold on;
stem(0:sParams.PlotSpectM-1, vLambda_A, 'DisplayName', '$\lambda_m(A)$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Analytic vs. Numeric eigenvalues');


subplot(2,1,2);
stem(0:sParams.PlotSpectM-1, log(vLambda_K),'LineStyle','none', 'DisplayName', '$\log(\lambda_m(K))$');
hold on;
stem(0:sParams.PlotSpectM-1, log(vLambda_A),'LineStyle','none', 'DisplayName', '$\log(\lambda_m(A))$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% print(fig, [sParams.sSim.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], '-depsc')
title('Analytic vs. Numeric eigenvalues log scale');
set(gcf,'Position',[100 100 600 500])
saveas(fig,[sParams.sSim.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd.png']);
end