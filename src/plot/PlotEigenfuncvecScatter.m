function localCmap = PlotEigenfuncvecScatter(sPlotParams, actualDataDist, mData, cmap, ...
    firstEigenIdx, lastEigIdx, mPhi, vLambda, lambdaStr, G, suptitle, figureName, phiStr, mPhi2, phi2Str, mData2)
dim = size(mData, 2);
assert(dim <= 3, 'Not supported')
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
if ~exist('mData2', 'var')
    mData2 = mData;
end
if dim == 1
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 800;
    height = 400;
    fig = figure('Name', '1D Scatter');
    mPhiForCmap = mPhi(:,firstEigenIdx+1:lastEigIdx+1);
    if exist('mPhi2', 'var')
        nEigenFuncsToPlot = 2*(lastEigIdx-firstEigenIdx+1);
        mPhi2ForCmap = mPhi2(:,firstEigenIdx+1:lastEigIdx+1);
        localCmap = [min([mPhiForCmap(:); mPhi2ForCmap(:)]), max([mPhiForCmap(:); mPhi2ForCmap(:)])];
    else
        nEigenFuncsToPlot = lastEigIdx-firstEigenIdx+1;
        localCmap = [min(mPhiForCmap(:)), max(mPhiForCmap(:))];
    end
    for m = firstEigenIdx:lastEigIdx
        dispName = ['$' phiStr '_{' num2str(m) '}$'];
        if exist('phi2Str', 'var')
            dispName2 = ['$' phi2Str '_{' num2str(m) '}$'];
        end
        if exist('mPhi2', 'var')
            plot(mData2(:), mPhi2(:,m+1), 'o', 'DisplayName', dispName2);
            hold on;
        end
        plot(mData(:), mPhi(:,m+1), '.', 'DisplayName', dispName);
        xlim([ min(mData) max(mData) ])
        hold on;
        set(gca,'FontSize', 14);
    end
    if exist('suptitle', 'var')
        title(suptitle,'Interpreter', 'latex', 'FontSize', 14);
    end
    
    nRows = floor(sqrt(2*nEigenFuncsToPlot+1));
    if nRows > 2
        nRows = 2;
    end
    nCols = ceil(nEigenFuncsToPlot/nRows);
    if isempty(cmap)
        ylim(localCmap);
    else
        ylim(cmap);
    end
    legend('Interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',nCols)
    set(gcf,'Position', [x0 y0 width height])
elseif dim == 2 || dim == 3
    %% Plot params
    nEigenFuncsToPlot = lastEigIdx-firstEigenIdx+1;
    nRows = floor(sqrt(nEigenFuncsToPlot+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(nEigenFuncsToPlot/nRows);
    
    x0     = 10;
    y0     = 50;
    height = 300*nRows;
    width  = 400*nCols;
    %% Plot
    if exist('mPhi2', 'var')
        nLoops = 2;
    else
        nLoops = 1;
    end
    if isempty(cmap)
        mPhiForCmap = mPhi(:,firstEigenIdx+1:lastEigIdx+1);
        if exist('mPhi2', 'var')
            mPhi2ForCmap = mPhi2(:,firstEigenIdx+1:lastEigIdx+1);
        else
            mPhi2ForCmap = [];
        end
        localCmap = [min([mPhiForCmap(:); mPhi2ForCmap(:)]), max([mPhiForCmap(:); mPhi2ForCmap(:)])];
    end
    for i = 1:nLoops
        if i == 2
            mPhi = mPhi2;
            phiStr = phi2Str;
        end
        fig = figure('Name', [ num2str(dim) 'D Scatter']);
        tiledlayout(nRows, nCols)
        vAx = zeros(nRows*nCols);
        for m = firstEigenIdx:lastEigIdx
            vAx(m+1-firstEigenIdx) = nexttile;
            if dim == 2
                scatter3(mData(:,1), mData(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
            else % dim == 3
                scatter3(mData(:,1), mData(:,2), mData(:,3), [], mPhi(:,m+1), 'filled');
                UpdateCursorDataTip(fig, vAx, mPhi);
            end
            colormap(gca, 'jet')
            colorbar('TickLabelInterpreter', 'latex');
            if ~isempty(cmap)
                localCmap(1) = cmap(1);
                localCmap(2) = cmap(2);
            end
            caxis([localCmap(1) localCmap(2)])
            xlim([ min(mData(:,1)) max(mData(:,1))])
            ylim([ min(mData(:,2)) max(mData(:,2))])
            if dim == 2
                view(2); %view(20,40);
            else % dim == 3
                view(30,75);
                zlim([ min(mData(:,3)) max(mData(:,3))])
            end
            
            dispName = ['$' phiStr '_{' num2str(m) '}$'];
            if exist('vLambda', 'var') && ~isempty(vLambda)
                lambda_m_str = ['$' lambdaStr '_{' num2str(m) '} = ' num2str(vLambda(m+1), '%.4f') '$'];
                dispName = strcat(dispName, ', ', lambda_m_str);
            end
            title(dispName, 'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end
        if exist('suptitle', 'var')
            sgtitle(suptitle,'Interpreter', 'latex', 'FontSize', 16);
        end
        set(gcf,'Position', [x0 y0 width height])
        set(fig,'renderer','Painters')
    end
end

%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = [figureName, '_', sPlotParams.matrixForEigs, '_eigenvectors', '_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx)];
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end

set(0,'DefaultFigureWindowStyle',windowStyle)

end