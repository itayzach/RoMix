function [] = PlotEigenfunctionsEigenvectors(sSimParams, sDataset, nysRatio, firstEigenIdx, lastEigIdx, mPhiToCompare, mPhiNumeric, figName)

if sDataset.dim == 1
    fig = figure('Name', 'Eigenfuncs/vecs');
    %% 1-D To Compare
    subplot(2,2,1);
    for m = firstEigenIdx:lastEigIdx
        plot(sDataset.sData.x, mPhiToCompare(:,m+1), '.', 'DisplayName', [ '$\phi_{' num2str(m) '}(x_i)$' ]);
        hold on
        xlim([min(sDataset.sData.x) max(sDataset.sData.x)]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
    end
    set(gca,'FontSize', 14);
    title(figName, 'Interpreter', 'latex', 'FontSize', 14)
%     legend('Interpreter', 'latex', 'FontSize', 14)
    %% 1-D Numeric
    subplot(2,2,2);
    for m = firstEigenIdx:lastEigIdx
        plot(sDataset.sData.x, mPhiNumeric(:,m+1), '.', 'LineWidth', 1, 'DisplayName', [ '$v_{' num2str(m) '}(i)$' ]);
        hold on
        xlim([min(sDataset.sData.x) max(sDataset.sData.x)]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
    end
    set(gca,'FontSize', 14);
    title('Numeric', 'Interpreter', 'latex', 'FontSize', 14)
    %% Diff
    subplot(2,1,2);
    for m = firstEigenIdx:lastEigIdx
        plot(sDataset.sData.x, mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), '.', 'LineWidth', 1, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
        hold on
        xlim([min(sDataset.sData.x) max(sDataset.sData.x)]);
        xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
    end
    set(gca,'FontSize', 14);
    title('Error', 'Interpreter', 'latex', 'FontSize', 14)
    set(gcf,'Position',[100 100 1000 800])

else
    error('Not supporting')
end

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
simPrefix = strcat(sDataset.actualDataDist, num2str(sDataset.dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_eigenvectors_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'epsc');

end