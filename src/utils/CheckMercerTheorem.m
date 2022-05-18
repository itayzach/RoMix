function CheckMercerTheorem(sDistParams, Phi, lambdaPhi, W, x)
[n, MTilde] = size(Phi);
gmmNumComponents = sDistParams.estNumComponents;
dim = sDistParams.dim;
%%
LambdaPhi = diag(lambdaPhi);
WPhi = Phi*LambdaPhi*Phi.';
WPhiDiff = abs(W - WPhi);

nMercer = min(250,n);
figure('Name', 'Mercer''s Thm check'); tiledlayout('flow')
nexttile; imagesc(W(1:nMercer,1:nMercer)); colorbar; title ('W', 'Interpreter','latex')
nexttile; imagesc(WPhi(1:nMercer,1:nMercer)); colorbar; title('$W_\Phi = \Phi \Lambda \Phi^T$', 'Interpreter','latex')
nexttile; imagesc(WPhiDiff(1:nMercer,1:nMercer)); colorbar; title(['$|W - W_\Phi|$', newline, '$||W - W_\Phi||$ = ', num2str(norm(W-WPhi,'fro')) '.  max = ' num2str(max(WPhiDiff(:))) ], 'Interpreter','latex')

%%
[V, lambdaV] = eigs(W,MTilde);
lambdaV = diag(lambdaV)/n;
nexttile; 
plot(lambdaPhi,'o'); 
hold on; 
plot(lambdaV,'x'); 
legend('$\lambda_\Phi$', '$\frac{1}{n}\lambda_V$', 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'SouthOutside', 'NumColumns', 2)
title('Numeric vs. analytic eigenvalues')
sgtitle('Mercer''s Thm check (only first 50 points)')



end