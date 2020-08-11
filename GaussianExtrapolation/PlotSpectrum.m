function [] = PlotSpectrum(sSimParams, dataDim, vLambdaAnalytic, vLambdaNumeric, vLambdaNystrom)

fig = figure;
subplot(2,1,1);
stem(0:sSimParams.PlotSpectM-1, vLambdaAnalytic, 'DisplayName', '$\lambda^{\phi}_m$');
hold on;
stem(0:sSimParams.PlotSpectM-1, vLambdaNumeric, 'DisplayName', '$\lambda^{v}_m$');
stem(0:sSimParams.PlotSpectM-1, vLambdaNystrom, 'DisplayName', '$\lambda^{\hat{v}}_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Eigenvalues', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);

subplot(2,1,2);
stem(0:sSimParams.PlotSpectM-1, log(vLambdaAnalytic),'LineStyle','none', 'DisplayName', '$\log(\lambda^{\phi}_m)$');
hold on;
stem(0:sSimParams.PlotSpectM-1, log(vLambdaNumeric),'LineStyle','none', 'DisplayName', '$\log(\lambda^{v}_m)$');
stem(0:sSimParams.PlotSpectM-1, log(vLambdaNystrom),'LineStyle','none', 'DisplayName', '$\log(\lambda^{\hat{v}}_m)$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
% print(fig, [sSimParams.sSim.outputFolder filesep 'fig_eigenvalues_' num2str(sSimParams.dim) 'd'], '-depsc')
title('Eigenvalues log-scale', ...
    'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(gcf,'Position',[100 100 600 500])
saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(dataDim) 'd'], 'epsc');
end