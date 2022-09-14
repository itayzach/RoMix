function fig = PlotDigits(sPlotParams, X, y, b_transpose, figTitle, figName)
[n,d] = size(X);
imsize = sqrt(d);

fig = figure('Name', figName); 
if exist('y','var') && ~isempty(y)
    %% With title and white spaces
    tiledlayout('flow','TileSpacing','none','Padding','none')
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
else
    %% No titles, and tight (concate all images to one image)
    nRows = 2;%floor(sqrt(n));
    nCols = ceil(n/nRows);
    while nCols > 10
        nRows = nRows + 1;
        nCols = ceil(n/nRows);
    end
    
    images = zeros(imsize,imsize,n);
    for i = 1:n
        image = reshape(X(i,:),imsize,imsize) ;
        if b_transpose
            image = image.';
        end
        images(:,:,i) = image;
    end
    concateImages = imtile(images,[],"GridSize",[nRows, nCols]);
    imshow(concateImages);
    
%     x0     = 10;
%     y0     = 50;
%     height = 100*nRows;
%     width  = 200*nCols;
%     set(gcf,'Position', [x0 y0 width height])
end

%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end