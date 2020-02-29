function vLambda = PlotSpectrum(sParams, sSimParams, vLambda_A)

vLambda = zeros(sParams.M, 1);
for m = 0:sParams.M-1
    vLambda(m+1) = lambda(sParams.a, sParams.b, m);
end

fig = figure;
subplot(2,1,1);
stem(0:sParams.M-1, vLambda);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Analytic eigenvalues');


subplot(2,1,2);
stem(vLambda_A);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14);
title('Numeric eigenvalues');

% print(fig, [sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], '-depsc')
saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd.png']);
end