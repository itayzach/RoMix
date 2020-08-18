function [] = PlotEigenfuncvecScatter(sSimParams, actualDataDist, mData, nysRatio, firstEigenIdx, lastEigIdx, mPhi, vLambda, figName)

dim = size(mData, 2);
assert(dim <= 2, 'Not supported')
if dim == 1
    %% Plot params
    x0     = 10;
    y0     = 50;
    width  = 1800;
    height = 900;
    fig = figure('Name', '1D Scatter');
    for m = firstEigenIdx:lastEigIdx
        if strcmp(figName, 'Analytic')
            dispName = ['$\phi^{Ana}_{' num2str(m) '}({\bf x}),$ $\lambda^{Ana}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(figName, 'Numeric')
            dispName = ['$\phi^{Num}_{' num2str(m) '}({\bf x}),$ $\lambda^{Num}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        elseif strcmp(figName, 'Nystrom')
            dispName = ['$\phi^{Nys}_{' num2str(m) '}({\bf x}),$ $\lambda^{Nys}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'];
        else
            assert('invalid option')
        end
        plot(mData(:), mPhi(:,m+1), '.', 'DisplayName', dispName);
%         scatter(mData(:), mPhi(:,m+1), 'filled', 'DisplayName', dispName);
        xlim([ min(mData) max(mData) ])
        hold on;
        set(gca,'FontSize', 14);
    end
    legend('Interpreter', 'latex', 'FontSize', 14)
    set(gcf,'Position', [x0 y0 width height])
elseif dim == 2
    %% Plot params
    x0     = 10;
    y0     = 50;
    width  = 1800;
    height = 900;
    nRows = floor(sqrt(sSimParams.PlotEigenFuncsM+1));
    if nRows > 4
        nRows = 4;
    end
    nCols = ceil(sSimParams.PlotEigenFuncsM/nRows);
    %% Plot
    fig = figure('Name', '2D Scatter');
    for m = firstEigenIdx:lastEigIdx
        subplot(nRows, nCols,m-firstEigenIdx+1);
        scatter3(mData(:,1), mData(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
        colormap(gca, 'default')
        colorbar()
        caxis([min(mPhi(:,m+1)) max(mPhi(:,m+1))])
        view(2)
        xlim([ min(mData(:,1)) max(mData(:,1))])
        ylim([ min(mData(:,2)) max(mData(:,2))])
        if strcmp(figName, 'Analytic')
            title(['$\phi_{' num2str(m) '}({\bf x_i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
        elseif strcmp(figName, 'Numeric')
            title(['$v_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
        elseif strcmp(figName, 'Nystrom')
            title(['$\hat{v}_{' num2str(m) '}({\bf i}),$ $\lambda_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
                'Interpreter', 'latex', 'FontSize', 14)
        else
            assert('invalid option')
        end
        set(gca,'FontSize', 14);
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
saveas(fig,strcat(sSimParams.outputFolder, filesep, simPrefix, '_eigenvectors_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'epsc');


end