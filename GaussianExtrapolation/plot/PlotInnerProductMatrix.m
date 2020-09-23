function fig = PlotInnerProductMatrix(sSimParams, sDistParams, sDataset, nysRatio, mPhi, figName)

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

fig = figure('Name', sprintf('%s IP matrix', figName));

if strcmp(figName, 'Numeric')
    mInnerProduct = mPhi.'*mPhi;
    pltTitle = [figName ' - $V^T V$'];
elseif strcmp(figName, 'Analytic')
    n = length(sDataset.sData.x);
    mInnerProduct = n^sDistParams.dim * (mPhi.' * diag(sDistParams.vPr)* mPhi);
    pltTitle = [figName ' - $\int \phi_i(x) \phi_j(x) p(x) dx = \Phi^T$diag(Pr)$\Phi$'];
elseif strcmp(figName, 'Nystrom')
    mInnerProduct = mPhi.'*mPhi;
    pltTitle = [figName ' - $\hat{V}^T \hat{V}$'];    
end
imagesc(mInnerProduct)
colormap(gca, 'hot')
colorbar()
title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
if isempty(nysRatio)
    simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd');
else
    simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
end
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix,  '_inner_product_matrix_', figName), 'epsc');
end