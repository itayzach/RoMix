function [] = PlotMetric(sPlotParams, mVals, mStd, cDispName, figName, vylim, xTickNames, pltTitle, ylab, xlab, b_zoomIn, b_ylog)
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
    p(methodId) = errorbar(xVals, vAcc, vStd, 'Marker', 'o', 'LineStyle', '-', 'LineWidth', 2, 'DisplayName', cDispName{methodId});
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
lgd = legend(p, 'Interpreter', 'latex', 'FontSize', 14, 'Location',  'SouthOutside', 'NumColumns', 3);
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
    
    %v2ndRunSorted = sort(mVals(2,:));
    %top3min = v2ndRunSorted(1);
    %top3max = v2ndRunSorted(end-1);

    %x0src = 25;  y0src = top3min-abs(top3min/10);
    %x1src = 150; y1src = top3max+abs(top3max/50);
    %srcPos =[x0src y0src x1src y1src];
    
    [mValsMax, maxInd] = max(mVals,[],2);
    keepind = setdiff(repmat(1:nMethods,M,1)', maxInd','rows');
    mValsNoMax = zeros(size(mVals,1),size(mVals,2)-1);
    for ii=1:M
        mValsNoMax(ii,:) = mVals(ii,keepind(:,ii));
    end
    x0src = 25;  y0src = min(mValsNoMax(:));
    x1src = 150; y1src = max(mValsNoMax(:));
    srcPos =[x0src y0src x1src y1src];


    %minValsSorted = sort(mVals(:));
    %x0tar = 60;  y0tar = min(mVals(:))+1;
    %x1tar = 145; y1tar = minValsSorted(2)+2;
    %vValsSorted = sort(mVals(:));
    range = min(mValsMax)-max(mValsNoMax(:));
    x0tar = 60;  y0tar = max(mValsNoMax(:)) + range*0.1;
    x1tar = 145; y1tar = min(mValsMax) - range*0.1;   
    targetPos = [x0tar y0tar x1tar y1tar];

    az = zoomPlot(p_ax, srcPos, targetPos, lgd);
end
if exist('b_ylog', 'var') && ~isempty(b_ylog) && b_ylog
    set(gca, 'YScale', 'log');
    set(gca, 'YTickLabel', {'10^{0}', '10^{1}'})
    set(gca, 'YTick', [1 10])
end

%% Save
SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});

set(0,'DefaultFigureWindowStyle',windowStyle)
end