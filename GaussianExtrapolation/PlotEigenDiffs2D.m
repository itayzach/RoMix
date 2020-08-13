function [] = PlotEigenDiffs2D(sSimParams, sDataset, nysRatio, firstEigenIdx, lastEigIdx, mPhiToCompare, mPhiNumeric, figName)

if sDataset.dim == 2
    fig = figure('Name', 'Error');
    %% Plot params
    x0     = 10;
    y0     = 50;
    width  = 1800;
    height = 900;
    nEig = lastEigIdx - firstEigenIdx + 1;
    nRows = floor(sqrt(nEig+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(nEig/nRows);
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        scatter3(sDataset.sData.x(:,1), sDataset.sData.x(:,2), mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), [], mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), 'filled');
        colormap(gca, 'default')
        colorbar()
        caxis([min(mPhiToCompare(:,m+1)) max(mPhiToCompare(:,m+1))])
        view(2)
        xlim([ min(sDataset.sData.x(:,1)) max(sDataset.sData.x(:,1))])
        ylim([ min(sDataset.sData.x(:,2)) max(sDataset.sData.x(:,2))])
        if strcmp(figName, 'Analytic')
            title(['$\phi_{' num2str(m) '}({\bf x_i}) - v_{' num2str(m) '}({\bf i})$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
        elseif strcmp(figName, 'Nystrom')
            title(['$\hat{v}_{' num2str(m) '}({\bf_i}) - v_{' num2str(m) '}({\bf i})$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
        else
            assert('invalid option')
        end
        set(gca,'FontSize', 14);
    end
    set(gcf,'Position', [x0 y0 width height])
else
    error('Not supporting')
end

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
if isempty(nysRatio)
    simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd');
else
    simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
end
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_eigen_diffs_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'epsc');

end