function [] = PlotEigenDiffs(sSimParams, actualDist, x, nysRatio, firstEigenIdx, lastEigIdx, mPhiToCompare, mPhiNumeric, figName, phiToCompareLhsStr, phiToCompareRhsStr, b_plotErrVsNodeInd)

windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

dim = size(x,2);

if b_plotErrVsNodeInd
    fig = figure('Name', 'Error');
    %% Plot params
    x0     = 10;
    y0     = 100;
    width  = 1200;
    height = 600;
    nEig = lastEigIdx - firstEigenIdx + 1;
    nRows = floor(sqrt(nEig+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(nEig/nRows);
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        plot(1:length(x), abs(mPhiNumeric(:,m+1) - mPhiToCompare(:,m+1)), '.');
        title(['$|' phiToCompareLhsStr '_{' num2str(m) '} - ' phiToCompareRhsStr '_{' num2str(m) '}|$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
        xlabel('Node index', 'Interpreter', 'latex', 'FontSize', 14)
        %legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14)
        set(gca, 'YScale', 'log')
        set(gca,'FontSize', 14);
    end
    set(gcf,'Position', [x0 y0 width height])
else
    if dim == 1
        fig = figure('Name', 'Error');
        %% Plot params
        x0     = 10;
        y0     = 100;
        width  = 1200;
        height = 600;
        nEig = lastEigIdx - firstEigenIdx + 1;
        nRows = floor(sqrt(nEig+1));
        if nRows > 4
            nRows = 4;
        end
        nCols = ceil(nEig/nRows);
        for m = firstEigenIdx:lastEigIdx
            subplot(nRows, nCols,m-firstEigenIdx+1);
            plot(x, mPhiNumeric(:,m+1) - mPhiToCompare(:,m+1), '.');
            title(['$' phiToCompareLhsStr '_{' num2str(m) '} - ' phiToCompareRhsStr '_{' num2str(m) '}$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
            xlim([min(x) max(x)]);
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
            %legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end
        set(gcf,'Position', [x0 y0 width height])
    elseif dim == 2
        fig = figure('Name', 'Error');
        %% Plot params
        x0     = 10;
        y0     = 100;
        width  = 1200;
        height = 600;
        nEig = lastEigIdx - firstEigenIdx + 1;
        nRows = floor(sqrt(nEig+1));
        if nRows > 4
            nRows = 4;
        end
        nCols = ceil(nEig/nRows);
        for m = firstEigenIdx:lastEigIdx
            subplot(nRows, nCols,m-firstEigenIdx+1);
            scatter3(x(:,1), x(:,2), mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), [], mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), 'filled');
            colormap(gca, 'default')
            colorbar()
            caxis([min(mPhiToCompare(:,m+1)) max(mPhiToCompare(:,m+1))])
            view(2)
            xlim([ min(x(:,1)) max(x(:,1))])
            ylim([ min(x(:,2)) max(x(:,2))])
            title(['$' phiToCompareLhsStr '_{' num2str(m) '} - ' phiToCompareRhsStr '_{' num2str(m) '}$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end
        set(gcf,'Position', [x0 y0 width height])
    elseif dim == 3
        fig = figure('Name', 'Error');
        %% Plot params
        x0     = 10;
        y0     = 100;
        width  = 1200;
        height = 600;
        nEig = lastEigIdx - firstEigenIdx + 1;
        nRows = floor(sqrt(nEig+1));
        if nRows > 4
            nRows = 4;
        end
        nCols = ceil(nEig/nRows);
        for m = firstEigenIdx:lastEigIdx
            subplot(nRows, nCols,m-firstEigenIdx+1);
            scatter3(x(:,1), x(:,2), x(:,3), [], mPhiToCompare(:,m+1) - mPhiNumeric(:,m+1), 'filled');
            colormap(gca, 'default')
            colorbar()
            caxis([min(mPhiToCompare(:,m+1)) max(mPhiToCompare(:,m+1))])
            view(10,5);
            xlim([ min(x(:,1)) max(x(:,1))])
            ylim([ min(x(:,2)) max(x(:,2))])
            zlim([ min(x(:,3)) max(x(:,3))])
            title(['$' phiToCompareLhsStr '_{' num2str(m) '} - ' phiToCompareRhsStr '_{' num2str(m) '}$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
        end
        set(gcf,'Position', [x0 y0 width height])    
    else
        error('Not supporting')
    end
end
%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end
saveas(fig,strcat(sSimParams.outputFolder, filesep, actualDist, '_eigen_diffs_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'png');
set(0,'DefaultFigureWindowStyle',windowStyle)
end