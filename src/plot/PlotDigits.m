function fig = PlotDigits(sPlotParams, X, y, b_transpose, figTitle, figName)
fig = figure('Name', figName); 
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
if exist('figTitle', 'var')
    sgtitle(figTitle, 'interpreter', 'latex', 'fontsize', 14)
end
nRows = floor(sqrt(n));
nCols = ceil(n/nRows);
x0     = 10;
y0     = 50;
height = 100*nRows;
width  = 200*nCols;
set(gcf,'Position', [x0 y0 width height])

if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end