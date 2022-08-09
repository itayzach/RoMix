function [] = PlotAccuracy(sPlotParams, mAcc, mStd, cDispName, figName, minAcc, xTickNames, pltTitle, ylab, xlab, b_zoomIn)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

[M, nMethods] = size(mAcc);
fig = figure('Name', 'Accuracy');
if ~exist('xTickNames', 'var') || isempty(xTickNames)
    xVals = (0:M-1)';
else
    xVals = str2double(xTickNames)';
    if any(isnan(xVals))
        xVals = (0:M-1)';
    end
end
for methodId = 1:nMethods
    vAcc = mAcc(:,methodId);
    vStd = mStd(:, methodId);
    errorbar(xVals, vAcc, vStd, 'Marker', 'o', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', cDispName{methodId});
    hold on;
end

if ~exist('minAcc', 'var') || isempty(minAcc)
    minAcc = min(mAcc(:))-2;
end
if ~exist('ylab', 'var')
    ylab = 'Accuracy [$\%$]';
end
%ylim([max(0,minAcc), 100]);
%xlim([0 M-1]);
xlim([min(xVals), max(xVals)]);
ylabel(ylab, 'Interpreter', 'latex', 'FontSize', 14)
legend('Interpreter', 'latex', 'FontSize', 14, 'Location',  'SouthOutside', 'NumColumns', 2)
if exist('pltTitle', 'var')
    title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
end
if exist('xTickNames', 'var') && ~isempty(xTickNames)
    set(gca,'xtick',xVals,'xticklabel',xTickNames)
    xtickangle(45)
    if exist('xlab','var')
        xlabel(xlab, 'Interpreter', 'latex', 'FontSize', 14)
    end
else
    xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14)
end
set(gca,'FontSize', 14);
%if M < 20
%    set(gca, 'XTick', 0:M-1);
%end
x0     = 10;
y0     = 50;
height = 350;
width  = 600;
set(gcf,'Position', [x0 y0 width height])
if exist('b_zoomIn', 'var') && ~isempty(b_zoomIn) && b_zoomIn
    p_ax=gca;
    srcPos =[45 92 150 93.5];
    targetPos = [50 72 145 80];    
    az = zoomPlot(p_ax, srcPos, targetPos);
end
%% Save
SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});

set(0,'DefaultFigureWindowStyle',windowStyle)
end