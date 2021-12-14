function [] = PlotAccuracy(sPlotParams, mAcc, cDispName)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

[M, nMethods] = size(mAcc);
cLineStyle = {'-', '--', '-.'};
fig = figure('Name', 'Accuracy');
for methodId = 1:nMethods
    vAcc = mAcc(:,methodId);
    plot((0:M-1)', vAcc, 'Marker', 'o', 'LineStyle', cLineStyle{methodId}, 'LineWidth', 2, 'DisplayName', cDispName{methodId});
    hold on;
end
ylim([50 100]);
xlim([0 M-1]);
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('Accuracy [$\%$]', 'Interpreter', 'latex', 'FontSize', 14)
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

    saveas(fig,strcat(sPlotParams.outputFolder, filesep, simPrefix, '_Acc_eigs_0_to_', ...
        num2str(M-1)), 'epsc');
end

set(0,'DefaultFigureWindowStyle',windowStyle)
end