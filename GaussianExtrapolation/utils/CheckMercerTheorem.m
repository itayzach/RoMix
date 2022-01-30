function CheckMercerTheorem(Phi, lambdaPhi, gmmNumComponents, W)
[n, MTilde] = size(Phi);
WPhi = Phi*diag(lambdaPhi)*Phi.';
WPhiDiff = abs(W - WPhi);

nMercer = 50;
figure; tiledlayout('flow')
nexttile; imagesc(W(1:nMercer,1:nMercer)); colorbar; title ('W')
nexttile; imagesc(WPhi(1:nMercer,1:nMercer)); colorbar; title('$W_\Phi = \Phi \Lambda \Phi^T$', 'Interpreter','latex')
nexttile; imagesc(WPhiDiff(1:nMercer,1:nMercer)); colorbar; title(['$|W - W_\Phi|$', newline, '$||W - W_\Phi||$ = ', num2str(norm(W-WPhi,'fro')) '.  max = ' num2str(max(WPhiDiff(:))) ], 'Interpreter','latex')
[V2, Lambda2] = eigs(W,MTilde);
lambda2 = diag(Lambda2)/(n/gmmNumComponents);
figure; plot(lambdaPhi,'o'); hold on; plot(lambda2,'x'); legend('$\lambda_\Phi$', '$\lambda_V$','interpreter', 'latex','fontsize',14)  
end