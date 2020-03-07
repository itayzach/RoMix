function vLambda_K = PlotSpectrum(sParams, sSimParams, vLambda_A)

vLambda_K = zeros(sParams.PlotSpectM, 1);
for m = 0:sParams.PlotSpectM-1
    vLambda_K(m+1) = prod(lambda(sParams, m), 2);
end


vLambda_An = (1/length(vLambda_A)) * vLambda_A(1:sParams.PlotSpectM);
nom = abs(vLambda_K - vLambda_An);
denom = abs(vLambda_K);

fig = figure;
subplot(2,1,1);
stem(0:sParams.PlotSpectM-1, vLambda_K, 'DisplayName', '$\lambda_m(K)$');
hold on;
stem(0:sParams.PlotSpectM-1, vLambda_An, 'DisplayName', '$\lambda_m(A)$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
xlabel('$M$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
title('Analytic eigenvalues');


subplot(2,1,2);
% plot(0:sParams.PlotSpectM-1, nom./denom);
% xlabel('$M$', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('$E(\lambda(K), \lambda(A))$', 'Interpreter', 'latex', 'FontSize', 14);
% title('Numeric eigenvalues');
plot(0:sParams.PlotSpectM-1, log(vLambda_K), 'DisplayName', '$\log(\lambda_m(K))$');
hold on;
plot(0:sParams.PlotSpectM-1, log(vLambda_An), 'DisplayName', '$\log(\lambda_m(A))$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'best');
% print(fig, [sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd'], '-depsc')
saveas(fig,[sSimParams.outputFolder filesep 'fig_eigenvalues_' num2str(sParams.dim) 'd.png']);
end