MMercer = M;
MMercer = min(MTilde, MMercer);
mercer = N*(PhiTildeInt(:,1:MMercer)*diag(lambdaAnalyticTilde(1:MMercer))*PhiTildeInt(:,1:MMercer)');
ATilde = SimpleCalcAdjacency(xTildeInt,'GaussianKernel', omega);
isalmostequal(ATilde, mercer, 1e-4, [], false)

coraLatent = load('data\coraLatentVGAE.mat');

figure;
tiledlayout(2,2);
ax(1) = nexttile();
    imagesc(mercer);
    colorbar()
    caxis([0 max(1,max(mercer(:)))])
    title('$\Phi \Lambda \Phi^T$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
ax(2) = nexttile();
    imagesc(W);
    colorbar()
    caxis([0 1])
    title('${\bf W} = \exp\big(\|x-y\|^2/\omega^2\big)$', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
ax(3) = nexttile();
    imagesc(coraLatent.adj_norm);
    colorbar();
    title('$A$ (norm)', 'interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
linkaxes(ax,'xy')

