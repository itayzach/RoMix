function PlotDigits(sPlotParams, X, y, b_transpose, figTitle, figName)
fig = figure; 
tiledlayout('flow')
[n,d] = size(X);
imsize = sqrt(d);

for i = 1:n
    nexttile;
    image = reshape(X(i,:),imsize,imsize) ;
    if b_transpose
        imshow(image.')
    else
        imshow(image)
    end
    if exist('y','var') && ~isempty(y)
        title(y(i),'Interpreter','latex','FontSize',14)
    end
end
sgtitle(figTitle, 'interpreter', 'latex', 'fontsize', 14)
nRows = floor(sqrt(n));
nCols = ceil(n/nRows);
x0     = 10;
y0     = 50;
height = 100*nRows;
width  = 200*nCols;
set(gcf,'Position', [x0 y0 width height])

if isfield(sPlotParams, 'outputFolder')
    if ~exist(sPlotParams.outputFolder, 'dir')
        mkdir(sPlotParams.outputFolder)
    end
    saveas(fig, strcat(sPlotParams.outputFolder, filesep, figName), 'epsc');
end
end