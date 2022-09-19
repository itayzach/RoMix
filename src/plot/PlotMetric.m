function [] = PlotMetric(sPlotParams, mVals, mStd, cDispName, figName, vylim, xTickNames, pltTitle, ylab, xlab, b_zoomIn)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

[M, nMethods] = size(mVals);
fig = figure('Name', figName);
if ~exist('xTickNames', 'var') || isempty(xTickNames)
    xVals = (0:M-1)';
else
    xVals = str2double(xTickNames)';
    if any(isnan(xVals))
        xVals = (0:M-1)';
    end
end
for methodId = 1:nMethods
    vAcc = mVals(:,methodId);
    vStd = mStd(:, methodId);
    errorbar(xVals, vAcc, vStd, 'Marker', 'o', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', cDispName{methodId});
    hold on;
end

if ~exist('vylim', 'var') || isempty(vylim)
    vylim = [min(mVals(:)), max(mVals(:))];
end
if ~exist('ylab', 'var') || isempty(ylab)
    ylab = 'Accuracy [$\%$]';
end
ylim(vylim);
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

%% Size
x0     = 10;
y0     = 50;
height = 350;
width  = 600;
set(gcf,'Position', [x0 y0 width height])

%% Zoom in
if exist('b_zoomIn', 'var') && ~isempty(b_zoomIn) && b_zoomIn
    p_ax=gca;

    top3min = sort(min(mVals(2:end,:)));
    top3min = top3min(3);

    top3max = sort(max(mVals(2:end,:)));
    top3max = top3max(end);

    x0src = 45;  y0src = top3min-0.1;
    x1src = 150; y1src = top3max+0.1;
    srcPos =[x0src y0src x1src y1src];
    
    minValsSorted = sort(mVals(:));

    x0tar = 60;  y0tar = min(mVals(:))+1;
    x1tar = 145; y1tar = minValsSorted(2)+2;
    targetPos = [x0tar y0tar x1tar y1tar];

    az = zoomPlot(p_ax, srcPos, targetPos);
end
%% Save
SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});

set(0,'DefaultFigureWindowStyle',windowStyle)
end