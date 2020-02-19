function [] = PlotSpectrum(sParams, sSimParams)
dx = 0.01;
x = (-5:dx:5-dx).';

vLambda = zeros(sParams.M, 1);
for m = 0:sParams.M-1
    [~, vLambda(m+1)] = phi(sParams.a, sParams.b, m, x);
end

fig1 = figure();
plot(0:sParams.M-1, vLambda);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$\lambda_m$', 'Interpreter', 'latex', 'FontSize', 14)
print(fig1, [sSimParams.outputFolder filesep 'fig2_eigenvalues'], '-depsc')
end