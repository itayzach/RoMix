function [] = PlotEigenfuncvecScatter(sSimParams, sDataset, nysRatio, firstEigenIdx, lastEigIdx, mPhi, vLambda, figName)

assert(sDataset.dim == 2, 'Not supported')

%% Plot params
x0     = 10;
y0     = 50;
width  = 1800;
height = 900;
nRows = floor(sqrt(sSimParams.PlotEigenFuncsM+1));
if nRows > 4
    nRows = 4;
end
nCols = ceil(sSimParams.PlotEigenFuncsM/nRows);
%% Plot
fig = figure;
for m = firstEigenIdx:lastEigIdx
    subplot(nRows, nCols,m-firstEigenIdx+1);
    scatter3(sDataset.mData.x(:,1), sDataset.mData.x(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
    colormap(gca, 'default')
    colorbar()
    caxis([min(mPhi(:,m+1)) max(mPhi(:,m+1))])
    view(2)
    xlim([ min(sDataset.mData.x(:,1)) max(sDataset.mData.x(:,1))])
    ylim([ min(sDataset.mData.x(:,2)) max(sDataset.mData.x(:,2))])
    if strcmp(figName, 'Analytic')
        title(['$\phi^{Ana}_{' num2str(m) '}({\bf x}),$ $\lambda^{Ana}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    elseif strcmp(figName, 'Numeric')
        title(['$\phi^{Num}_{' num2str(m) '}({\bf x}),$ $\lambda^{Num}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    elseif strcmp(figName, 'Nystrom')
        title(['$\phi^{Nys}_{' num2str(m) '}({\bf x}),$ $\lambda^{Nys}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    else
        assert('invalid option')
    end
    set(gca,'FontSize', 14);
end
set(gcf,'Position', [x0 y0 width height])

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
simPrefix = strcat(sDataset.actualDataDist, '_', num2str(nysRatio*100, '%d'), 'prec');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_eigenvectors_2d_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'epsc');


end