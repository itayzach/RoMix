function vLambda = PlotSpectrum(sParams, sSimParams, vLambda_A)

vLambda = zeros(sParams.PlotSpectM, 1);
for PlotSpectM = 0:sParams.PlotSpectM-1
    vLambda(PlotSpectM+1) = lambda(sParams, PlotSpectM);
end

fig = figure;
subplot(2,1,1);
stem(0:sParams.PlotSpectM-1, vLambda);
xlabel('$PlotSpectM$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Analytic eigenvalues');


subplot(2,1,2);
stem(vLambda_A);
xlabel('$M$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14);
title('Numeric eigenvalues');

% print(fig, [sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], '-depsc')
saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd.png']);
end