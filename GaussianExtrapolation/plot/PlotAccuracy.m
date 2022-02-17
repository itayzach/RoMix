function [] = PlotAccuracy(sPlotParams, mAcc, mStd, cDispName, figName, minAcc, xTickNames, pltTitle)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

[M, nMethods] = size(mAcc);
cLineStyle = {'-', '--', '-.'};
fig = figure('Name', 'Accuracy');
for methodId = 1:nMethods
    vAcc = mAcc(:,methodId);
    vStd = mStd(:, methodId);
    errorbar((0:M-1)', vAcc, vStd, 'Marker', 'o', 'LineStyle', cLineStyle{methodId}, 'LineWidth', 2, 'DisplayName', cDispName{methodId});
    hold on;
end

if ~exist('minAcc', 'var') || isempty(minAcc)
    minAcc = min(mAcc(:))-10;
end
ylim([minAcc, 100]);
xlim([0 M-1]);
ylabel('Accuracy [$\%$]', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location',  'SouthOutside', 'NumColumns', 2)
if exist('pltTitle', 'var')
    title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
end
if exist('xTickNames', 'var') && ~isempty(xTickNames)
    set(gca,'xtick',(0:M-1)','xticklabel',xTickNames)
    xtickangle(45)
else
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
end
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
    simPrefix = strcat(sPlotParams.actualDataDist, num2str(sPlotParams.dim), ...
        'd', '_', sPlotParams.matrixForEigs);

    saveas(fig,strcat(sPlotParams.outputFolder, filesep, simPrefix, '_', figName), 'epsc');
end

set(0,'DefaultFigureWindowStyle',windowStyle)
end