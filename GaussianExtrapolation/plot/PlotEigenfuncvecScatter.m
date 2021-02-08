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
        if strcmp(evecsName, 'Analytic')
            dispName = ['$\phi_{' num2str(m) '},$ $\lambda^{\phi}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(evecsName, 'Numeric')
            dispName = ['$v_{' num2str(m) '},$ $\lambda^{v}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(evecsName, 'Nystrom')
            dispName = ['$\hat{v}_{' num2str(m) '},$ $\lambda^{\hat{v}}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            assert('invalid option')
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
elseif dim == 2
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 1200; 1800;
    height = 600; 900;
    nEigenFuncsToPlot = lastEigIdx-firstEigenIdx+1;
    nRows = floor(sqrt(nEigenFuncsToPlot+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(nEigenFuncsToPlot/nRows);
    %% Plot
    fig = figure('Name', '2D Scatter');
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        if sSimParams.b_GSPBoxPlots
            param.show_edges = false;
            gsp_plot_signal(G,mPhi(:,m+1),param); 
        else
            scatter3(mData(:,1), mData(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
            colormap(gca, 'default')
            colorbar()
            caxis([min(mPhi(:,m+1)) max(mPhi(:,m+1))])
            view(2); %view(20,40);
            xlim([ min(mData(:,1)) max(mData(:,1))])
            ylim([ min(mData(:,2)) max(mData(:,2))])
        end
        if strcmp(evecsName, 'Analytic')
            dispName = ['$\phi_{' num2str(m) '}({\bf x_i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
            if exist('c', 'var') && ~isempty(c)
                dispName = strcat(dispName, ', $c = ', num2str(c(m+1), '%.4f'), '$');
            end
        elseif strcmp(evecsName, 'Numeric')
            dispName = ['$v_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(evecsName, 'Nystrom')
            dispName = [ '$\hat{v}_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            assert('invalid option')
        end
        title(dispName, 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
    end
    if exist('suptitle', 'var')
        sgtitle(suptitle,'Interpreter', 'latex', 'FontSize', 16);
    end
    set(gcf,'Position', [x0 y0 width height])
elseif dim == 3
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 1200; 1800;
    height = 600; 900;
    nEigenFuncsToPlot = lastEigIdx-firstEigenIdx+1;
    nRows = floor(sqrt(nEigenFuncsToPlot+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(nEigenFuncsToPlot/nRows);
    %% Plot
    fig = figure('Name', '3D Scatter');
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        if sSimParams.b_GSPBoxPlots
            assert('should be enabled')
            param.show_edges = false;
            gsp_plot_signal(G,mPhi(:,m+1),param); 
        else
            scatter3(mData(:,1), mData(:,2), mData(:,3), [], mPhi(:,m+1), 'filled');
            colormap(gca, 'default')
            colorbar()
            caxis([min(mPhi(:,m+1)) max(mPhi(:,m+1))])
            view(10,5);
            xlim([ min(mData(:,1)) max(mData(:,1))])
            ylim([ min(mData(:,2)) max(mData(:,2))])
            zlim([ min(mData(:,3)) max(mData(:,3))])
        end
        if strcmp(evecsName, 'Analytic')
            dispName = ['$\phi_{' num2str(m) '}({\bf x_i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
            if exist('c', 'var') && ~isempty(c)
                dispName = strcat(dispName, ', $c = ', num2str(c(m+1), '%.4f'), '$');
            end
        elseif strcmp(evecsName, 'Numeric')
            dispName = ['$v_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(evecsName, 'Nystrom')
            dispName = [ '$\hat{v}_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            assert('invalid option')
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