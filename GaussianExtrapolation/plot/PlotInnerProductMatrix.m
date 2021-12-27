function fig = PlotInnerProductMatrix(sPlotParams, mPhi, vPr, pltTitle, figName)

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

fig = figure('Name', sprintf('%s IP matrix', figName));

if ~exist('vPr', 'var') || isempty(vPr)
    mInnerProduct = mPhi.'*mPhi;
else
%     mInnerProduct = N^(dim-1) * (mPhi.' * diag(vPr)* mPhi);
    mInnerProduct = mPhi.' * diag(vPr)* mPhi;
end
imagesc(mInnerProduct)
% colormap(gca, 'hot')
colorbar()
title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])
set(0,'DefaultFigureWindowStyle',windowStyle)
%% Save
if isfield(sPlotParams, 'outputFolder') && ~isempty(sPlotParams)
    if ~exist(sPlotParams.outputFolder, 'dir')
        mkdir(sPlotParams.outputFolder)
    end
    simPrefix = strcat(sPlotParams.actualDataDist, num2str(sPlotParams.dim), ...
        'd', '_', sPlotParams.matrixForEigs);

    saveas(fig,strcat(sPlotParams.outputFolder, filesep, simPrefix, '_ip_matrix_', figName), 'epsc');
end
end
