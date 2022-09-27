function localCmap = PlotGraphSignals(sPlotParams, suptitle, figName, cData, cSignals, cSigStr, cCirclesInd, cMarkers, cColors, xylim, cmap, grid, tx, sense)
cData = reshape(cData,[],1);
cSignals = reshape(cSignals,[],1);
cSigStr = reshape(cSigStr,[],1);
nSignals = numel(cData);
if exist('cNumCircles', 'var')
    cCirclesInd = reshape(cCirclesInd,[],1);
end
if exist('cMarkers', 'var')
    cMarkers = reshape(cMarkers,[],1);
else
    cMarkers = cell(nSignals,1); cMarkers(:) = {'.'};
end
if exist('cColors', 'var')
    cColors = reshape(cColors,[],1);
end
dim = size(cData{1}, 2);
assert(dim <= 3, 'Not supported')
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

if dim == 1
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 800;
    height = 400;
    fig = figure('Name', 'Graph signals');
    for m = 1:nSignals
        dispName = cSigStr{m};
        if exist('cColors', 'var') && ~isempty(cColors)
            plot(cData{m}, cSignals{m}, cMarkers{m}, 'Color', cColors{m}, 'DisplayName', dispName);
        else
            plot(cData{m}, cSignals{m}, cMarkers{m}, 'DisplayName', dispName);
        end
        xlim([ min(cData{m}) max(cData{m}) ])
        hold on;
        set(gca,'FontSize', 14);
    end
    if exist('suptitle', 'var') && ~isempty(suptitle)
        title(suptitle,'Interpreter', 'latex', 'FontSize', 14);
    end
    nRows = floor(sqrt(2*nSignals+1));
    if nRows > 2
        nRows = 2;
    end
    nCols = ceil(nSignals/nRows);
    if nCols == 1
        nCols = 2;
    end
    legend('Interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',nCols)
    set(gcf,'Position', [x0 y0 width height])
    localCmap = [];
elseif dim == 2 || dim == 3
    %% Plot params
%     nRows = floor(sqrt(nSignals));
%     if nRows > 4
%         nRows = 4;
%     end
%     nCols = ceil(nSignals/nRows);
    if nSignals <= 2
        nRows = 1; 
        nCols = nSignals;
        x0     = 10;
        y0     = 50;
        height = 400;
        width  = 500*nCols;
    else
        nRows = 2;
        nCols = ceil(nSignals/nRows);
        x0     = 10;
        y0     = 50;
        if dim == 2
            xMax = max(cell2mat(cellfun(@max, cData, 'UniformOutput', false)));
            xMin = max(cell2mat(cellfun(@min, cData, 'UniformOutput', false)));
            hRatio = (xMax(2) - xMin(2))/(xMax(1) - xMin(1));
            hRatio = 2.5*hRatio;
            %hRatio = 1.5*hRatio;
            assert(hRatio > 0)
        else
            hRatio = 1;
        end
        height = 300*hRatio*nRows;
        width  = 400*nCols;
    end
   
    %% Plot
    fig = figure('Name', 'Graph signals');
    tiledlayout(nRows, nCols);
    ax = zeros(nSignals,1);
    if exist('cmap', 'var') && ~isempty(cmap)
        cMapMax = cmap(2);
        cMapMin = cmap(1);
    else
        cMapMax = max(cSignals{1});
        cMapMin = min(cSignals{1});
    end
    for m = 1:nSignals
        mData = cData{m};
        vSignal = cSignals{m};
        sigStr = cSigStr{m};
        if ~exist('cCirclesInd', 'var')
            error('should be here')
            vCirclesInd = 1:length(vSignal);
            vRectInd = [];
        else
            vCirclesInd = cCirclesInd{m};
            vRectInd = setdiff(1:length(vSignal),vCirclesInd);
        end
        ax(m) = nexttile;
        if exist('tx', 'var') % Meaning: strcmp(sPlotParams.actualDataDist, 'BulgariBeacons') == 1
            colormap(jet);
            %imagesc(grid.lon, (grid.lat), (reshape(vSignal, length(grid.lat), length(grid.lon))));
            scatter3(mData(vCirclesInd,1), mData(vCirclesInd,2), vSignal(vCirclesInd), 10, vSignal(vCirclesInd), 'filled','s');
            view(2);
            c=colorbar('TickLabelInterpreter', 'latex');
            c.Label.String='dBm';
            if cMapMin < cMapMax
               caxis([cMapMin cMapMax]);
            end
            xlabel('Lon. [deg]'); 
            ylabel('Lat. [deg]');
            set(gca,'YDir','normal')
            %set(get(gcf,'Children'),'YDir','normal');

            hold on;
            plot(tx.lon, tx.lat, 'db','MarkerFaceColor','b','MarkerSize',2.5);
            plot(sense.lon, sense.lat, 'yo','MarkerSize',3);
            
        else
            if dim == 2
                scatter3(mData(vCirclesInd,1), mData(vCirclesInd,2), vSignal(vCirclesInd), 10, ...
                   vSignal(vCirclesInd), 'filled');
                if vCirclesInd < length(mData)
                    hold on;
                    scatter3(mData(vRectInd,1), mData(vRectInd,2), vSignal(vRectInd), 50, ...
                        vSignal(vRectInd), 'filled', 's');
                end
            else % dim == 3
                scatter3(mData(vCirclesInd,1), mData(vCirclesInd,2), mData(vCirclesInd,3), 10, ...
                    vSignal(vCirclesInd), 'filled');
                if vCirclesInd < length(mData)
                    hold on;
                    scatter3(mData(vRectInd,1), mData(vRectInd,2), mData(vRectInd,3), 10, ...
                        vSignal(vRectInd), 'filled', 's');
                end
            end
            colormap(gca, 'jet');
            colorbar('TickLabelInterpreter', 'latex');
            if cMapMin < cMapMax
               caxis([cMapMin cMapMax]);
            end
            xlim([ min(mData(:,1)) max(mData(:,1))])
            ylim([ min(mData(:,2)) max(mData(:,2))])
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 14)
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 14)
            if dim == 2
                view(2); %view(20,40);
            else % dim == 3
                view(30,75);
                zlim([ min(mData(:,3)) max(mData(:,3))])
                zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', 14)
            end
        end
        dispName = sigStr;
        
        title(dispName, 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
        if exist('xylim', 'var') && ~isempty(xylim)
            xlim([xylim(1), xylim(2)]);
            ylim([xylim(3), xylim(4)]);
        end
    end
    if exist('suptitle', 'var')
        sgtitle(suptitle,'Interpreter', 'latex', 'FontSize', 16);
    end
%     linkaxes(ax)
    set(fig,'renderer','Painters')
    set(gcf,'Position', [x0 y0 width height])
    localCmap = [cMapMin; cMapMax];
end

%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = [figName, '_signals'];
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
set(0,'DefaultFigureWindowStyle',windowStyle)

end
