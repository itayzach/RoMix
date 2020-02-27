function vLambda = PlotSpectrum(sParams, vLambda_A)

vLambda = zeros(sParams.M, 1);
for m = 0:sParams.M-1
    vLambda(m+1) = lambda(sParams.a, sParams.b, m);
end

figure;
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
end