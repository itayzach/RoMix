function [] = PlotRMSE(sPlotParams, vRmseInt, vRmseNys)

M = length(vRmseInt);

fig = figure('Name', 'RMSE');
plot((0:M-1)', vRmseInt.', '-o', 'LineWidth', 2, 'DisplayName', ...
    'RMSE$(v^{{\bf int}}_m, v^{{\bf ref}}_m)$'), hold on;
plot((0:M-1)', vRmseNys.', '-x', 'LineWidth', 2, 'DisplayName', ...
    'RMSE$(v^{{\bf nys}}_m, v^{{\bf ref}}_m)$')
rmseylim = max([vRmseInt vRmseNys]);
ylim([0 rmseylim]);
xlim([0 M-1]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location',  'SouthOutside', 'NumColumns', 2)
set(gca,'FontSize', 14);
if M < 20
    set(gca, 'XTick', 0:M-1);
end
x0     = 10;
y0     = 50;
height = 350;
width  = 600;
set(gcf,'Position', [x0 y0 width height])
%% Save
if isfield(sPlotParams, 'outputFolder')
    if ~exist(sPlotParams.outputFolder, 'dir')
        mkdir(sPlotParams.outputFolder)
    end
    simPrefix = strcat(sPlotParams.sDataset.actualDataDist, num2str(sPlotParams.sDataset.dim), ...
        'd', '_', sPlotParams.matrixForEigs);

    saveas(fig,strcat(sPlotParams.outputFolder, filesep, simPrefix, '_RMSE_eigs_0_to', ...
        num2str(M)), 'epsc');
end


end