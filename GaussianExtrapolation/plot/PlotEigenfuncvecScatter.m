function [] = PlotEigenfuncvecScatter(sSimParams, actualDataDist, mData, nysRatio, ...
    firstEigenIdx, lastEigIdx, mPhi, vLambda, c, G, suptitle, figureName, phiStr, mPhi2, phi2Str)
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
        dispName = ['$' phiStr '_{' num2str(m) '}$'];
        if exist('phi2Str', 'var')
            dispName2 = ['$' phi2Str '_{' num2str(m) '}$'];
        end
        if exist('mPhi2', 'var')
            plot(mData(:), mPhi2(:,m+1), 'o', 'DisplayName', dispName2);
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
    if exist('mPhi2', 'var')
        nLoops = 2;
    else
        nLoops = 1;
    end
    for i = 1:nLoops
        if i == 2
            mPhi = mPhi2;
            phiStr = phi2Str;
        end
        fig = figure('Name', [ num2str(dim) 'D Scatter']);
        for m = firstEigenIdx:lastEigIdx
            subplot(nRows, nCols,m-firstEigenIdx+1);
            if isfield(sSimParams, 'b_GSPBoxPlots') && sSimParams.b_GSPBoxPlots
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
            dispName = ['$' phiStr '_{' num2str(m) '}$'];
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