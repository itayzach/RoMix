function fig = PlotInnerProductMatrix(sSimParams, dim, vPr, graphName, nysRatio, mPhi, pltTitle, figName)

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

N = size(mPhi,1);

fig = figure('Name', sprintf('%s IP matrix', figName));

if isempty(vPr)
    mInnerProduct = mPhi.'*mPhi;
else
    mInnerProduct = N^dim * (mPhi.' * diag(vPr)* mPhi);
end
imagesc(mInnerProduct)
% colormap(gca, 'hot')
colorbar()
title(pltTitle, 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])
set(0,'DefaultFigureWindowStyle',windowStyle)
%% Save
if isfield(sSimParams, 'outputFolder')
    if ~exist(sSimParams.outputFolder, 'dir')
        mkdir(sSimParams.outputFolder)
    end
    if isempty(nysRatio)
        simPrefix = strcat(graphName, num2str(dim), 'd');
    else
        simPrefix = strcat(graphName, num2str(dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
    end
    saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix,  '_inner_product_matrix_', figName), 'epsc');
end
end
