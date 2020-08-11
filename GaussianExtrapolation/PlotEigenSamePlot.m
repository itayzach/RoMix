function [] = PlotEigenSamePlot(sSimParams, sDataset, nysRatio, firstEigenIdx, lastEigIdx, mPhiAnalytic, mPhiNumeric, mPhiNystrom)

if sDataset.dim == 1
    fig = figure;
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
        plot(sDataset.sData.x, mPhiAnalytic(:,m+1), '.', 'DisplayName', [ '$\phi_{' num2str(m) '}$' ]);
        hold on
        plot(sDataset.sData.x, mPhiNumeric(:,m+1), '.', 'DisplayName', [ '$v_{' num2str(m) '}$' ]);
        plot(sDataset.sData.x, mPhiNystrom(:,m+1), '.', 'DisplayName', [ '$\hat{v}_{' num2str(m) '}$' ]);
%         title(['Blue: ' '$\phi_{' num2str(m) '}$' '; Red: $v_{' num2str(m) '}$; Yellow:  - $\hat{v}_{' num2str(m) '}$'], 'Interpreter', 'latex', 'FontSize', 14)
        xlim([min(sDataset.sData.x) max(sDataset.sData.x)]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14)
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
simPrefix = strcat(sDataset.actualDataDist, '_', num2str(nysRatio*100, '%d'), 'prec');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_same_plot_eigenvectors_1d_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx)), 'epsc');

end