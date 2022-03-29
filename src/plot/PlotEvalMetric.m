function [] = PlotEvalMetric(sPlotParams, mMetrics, cDispName, figName, minMaxMetric, xTickNames, pltTitle)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

[M, nMethods] = size(mMetrics);
cLineStyle = {'-', '--', '-.'};
fig = figure('Name', 'EvalMetric');
for methodId = 1:nMethods
    vMetric = mMetrics(:,methodId);
    plot((0:M-1)', vMetric, 'Marker', 'o', 'LineStyle', cLineStyle{methodId}, 'LineWidth', 2, 'DisplayName', cDispName{methodId});
%     errorbar()
    hold on;
end

if ~exist('minMaxMetric', 'var') || isempty(minMaxMetric)
    minMaxMetric(1) = min(mMetrics(:));
    minMaxMetric(2) = max(mMetrics(:));
    
end
ylim(minMaxMetric);
xlim([0 M-1]);
% ylabel('Coherence', 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location',  'SouthOutside', 'NumColumns', 2)
if exist('pltTitle', 'var')
    title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
end
if exist('xTickNames', 'var')
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
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end

set(0,'DefaultFigureWindowStyle',windowStyle)
end