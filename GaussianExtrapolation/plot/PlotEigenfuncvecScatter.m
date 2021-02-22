function [] = PlotEigenfuncvecScatter(sSimParams, actualDataDist, mData, nysRatio, ...
    firstEigenIdx, lastEigIdx, mPhi, vLambda, c, G, suptitle, figureName, phiStr, lambdaStr)
dim = size(mData, 2);
assert(dim <= 3, 'Not supported')
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')
if dim == 1
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 800;
    height = 400;
    fig = figure('Name', '1D Scatter');
    for m = firstEigenIdx:lastEigIdx
        if exist('lambdaStr', 'var')
            dispName = ['$' phiStr '_{' num2str(m) '},$ $' lambdaStr '_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            dispName = ['$' phiStr '_{' num2str(m) '}$'];
        end
        plot(mData(:), mPhi(:,m+1), '.', 'DisplayName', dispName);
%         scatter(mData(:), V(:,m+1), 'filled', 'DisplayName', dispName);
        xlim([ min(mData) max(mData) ])
        hold on;
        set(gca,'FontSize', 14);
    end
    if exist('suptitle', 'var')
        title(suptitle,'Interpreter', 'latex', 'FontSize', 14);
    end
    legend('Interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',min(ceil(length(firstEigenIdx:lastEigIdx)/2),3))
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
    fig = figure('Name', [ num2str(dim) 'D Scatter']);
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        if sSimParams.b_GSPBoxPlots
            param.show_edges = false;
            gsp_plot_signal(G,mPhi(:,m+1),param); 
        else
            if dim == 2
                scatter3(mData(:,1), mData(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
            else % dim == 3
                scatter3(mData(:,1), mData(:,2), mData(:,3), [], mPhi(:,m+1), 'filled');
            end
            colormap(gca, 'default')
            colorbar()
            cMin = min(mPhi(:,m+1));
            cMax = max(mPhi(:,m+1));
            caxis([cMin cMax])
            xlim([ min(mData(:,1)) max(mData(:,1))])
            ylim([ min(mData(:,2)) max(mData(:,2))])
            if dim == 2
                view(2); %view(20,40);
            else % dim == 3
                view(10,5);
                zlim([ min(mData(:,3)) max(mData(:,3))])
            end
        end
        if exist('lambdaStr', 'var')
            dispName = ['$' phiStr '_{' num2str(m) '},$ $' lambdaStr '_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            dispName = ['$' phiStr '_{' num2str(m) '}$'];
        end
        if exist('c', 'var') && ~isempty(c)
            dispName = strcat(dispName, ', $c = ', num2str(c(m+1), '%.4f'), '$');
        end
        title(dispName, 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
    end
    if exist('suptitle', 'var')
        sgtitle(suptitle,'Interpreter', 'latex', 'FontSize', 16);
    end
    set(gcf,'Position', [x0 y0 width height])
end

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
if isempty(nysRatio)
    simPrefix = strcat(actualDataDist, num2str(dim), 'd');
else
    simPrefix = strcat(actualDataDist, num2str(dim), 'd', '_', num2str(nysRatio*100, '%d'), 'prec');
end
if ~exist('figureName', 'var')
    figureName = evecsName;
end
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, ...
    '_', figureName, '_eigenvectors', '_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx)), 'png');
set(0,'DefaultFigureWindowStyle',windowStyle)

end