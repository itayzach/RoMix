function localCmap = PlotGraphSignals(sPlotParams, suptitle, figName, cData, cSignals, cSigStr, cNumCircles, cMarkers, cColors, xylim, cmap)
cData = reshape(cData,[],1);
cSignals = reshape(cSignals,[],1);
cSigStr = reshape(cSigStr,[],1);
if exist('cNumCircles', 'var')
    cNumCircles = reshape(cNumCircles,[],1);
end
if exist('cMarkers', 'var')
    cMarkers = reshape(cMarkers,[],1);
end
if exist('cColors', 'var')
    cColors = reshape(cColors,[],1);
end
dim = size(cData{1}, 2);
nSignals = numel(cData);
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
    nRows = 2;
    nCols = ceil(nSignals/nRows);
    
    x0     = 10;
    y0     = 50;
    height = 300*nRows;
    width  = 400*nCols;
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
        if ~exist('cNumCircles', 'var')
            nCircles = length(vSignal);
        else
            nCircles = cNumCircles{m};
        end
        ax(m) = nexttile;
        if dim == 2
            scatter3(mData(1:nCircles,1), mData(1:nCircles,2), vSignal(1:nCircles), [], ...
                vSignal(1:nCircles), 'filled');
            if nCircles < length(mData)
                hold on;
                scatter3(mData(nCircles+1:end,1), mData(nCircles+1:end,2), vSignal(nCircles+1:end), [], ...
                    vSignal(nCircles+1:end), 'filled', 's');
            end
        else % dim == 3
            scatter3(mData(1:nCircles,1), mData(1:nCircles,2), mData(1:nCircles,3), [], ...
                vSignal(1:nCircles), 'filled');
            if nCircles < length(mData)
                hold on;
                scatter3(mData(nCircles+1:end,1), mData(nCircles+1:end,2), mData(nCircles+1:end,3), [], ...
                    vSignal(nCircles+1:end), 'filled', 's');
            end
        end
        colormap(gca, 'jet');
        colorbar();
        caxis([cMapMin cMapMax]);
        xlim([ min(mData(:,1)) max(mData(:,1))])
        ylim([ min(mData(:,2)) max(mData(:,2))])
        if dim == 2
            view(2); %view(20,40);
        else % dim == 3
            view(30,75);
            zlim([ min(mData(:,3)) max(mData(:,3))])
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