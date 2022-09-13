function fig = PlotGmmResultWithDataset(sPlotParams, x, sDistParams, windowStyle, b_plotCovMean)
prevWindowStyle = get(0,'DefaultFigureWindowStyle');
if ~exist('windowStyle', 'var')
    windowStyle = prevWindowStyle;
end
set(0,'DefaultFigureWindowStyle',windowStyle)

[n, dim] = size(x);
assert(dim == 2);
if exist('sPlotParams', 'var') && ~isempty(sPlotParams)
    actualDataDist = sPlotParams.actualDataDist;
else
    actualDataDist = '';
end

fig = figure('Name', sprintf('%d-D %s', dim, actualDataDist));
% Data
scatter(x(:,1), x(:,2), 'filled');

% GMM
gmPDF = @(x,y) arrayfun(@(x0,y0) pdf(sDistParams.GMModel,[x0 y0]),x,y);
xMax = max(cell2mat(sDistParams.mu')) + 3*max(cell2mat(sDistParams.sigma'));
xMin = min(cell2mat(sDistParams.mu')) - 3*min(cell2mat(sDistParams.sigma'));
mMu = cell2mat(sDistParams.mu');
if b_plotCovMean
    hold on;
    fcontour(gmPDF, [xMin(1), xMax(1), xMin(2), xMax(2)]);
    scatter(mMu(:,1),mMu(:,2),[],'filled','red')
end
xlim([xMin(1), xMax(1)])
ylim([xMin(2), xMax(2)])
axis off
set(gca,'FontSize', 14);

%% Size
if strcmp(windowStyle, 'normal')
    hRatio = 2*(xMax(2)-xMin(2))/(xMax(1)-xMin(1));
    x0     = 400;
    y0     = 400;
    height = hRatio*400;
    width  = 600*(1 + (exist('b_plotDistModel', 'var')));
    set(gcf,'Position', [x0 y0 width height])
end
set(0,'DefaultFigureWindowStyle',prevWindowStyle)
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder') && exist('fig','var')
    if b_plotCovMean
        figName = 'gmm_cov_mean';
    else
        figName = 'gmm_data';
    end
    set(fig,'renderer','Painters')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end