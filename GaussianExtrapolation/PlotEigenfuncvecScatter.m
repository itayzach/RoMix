function [] = PlotEigenfuncvecScatter(sParams, firstEigenIdx, lastEigIdx, mPhi, vLambda, figName)

assert(sParams.dim == 2, 'Not supported')

%% Plot params
x0     = 10;
y0     = 50;
width  = 1800;
height = 900;
nRows = floor(sqrt(sParams.PlotEigenFuncsM+1));
if nRows > 4
    nRows = 4;
end
nCols = ceil(sParams.PlotEigenFuncsM/nRows);
%% Plot
fig = figure;
for m = firstEigenIdx:lastEigIdx
    subplot(nRows, nCols,m-firstEigenIdx+1);
    scatter3(sParams.x_rand(:,1), sParams.x_rand(:,2), mPhi(:,m+1), [], mPhi(:,m+1), 'filled');
    colormap(gca, 'default')
    colorbar()
    caxis([min(mPhi(:,m+1)) max(mPhi(:,m+1))])
    view(2)
    xlim([ min(sParams.x_rand(:,1)) max(sParams.x_rand(:,1))])
    ylim([ min(sParams.x_rand(:,2)) max(sParams.x_rand(:,2))])
    if strcmp(figName, 'Analytic')
        title(['$\phi^{A}_{' num2str(m) '}({\bf x}),$ $\lambda^{A}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    elseif strcmp(figName, 'Numeric')
        title(['$\phi^{Num}_{' num2str(m) '}({\bf x}),$ $\lambda^{Num}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    elseif strcmp(figName, 'Nystrom')
        title(['$\phi^{Nys}_{' num2str(m) '}({\bf x}),$ $\lambda^{Nys}_{' num2str(m)  '} = ' num2str(vLambda(m+1), '%.4f') '$'], ...
            'Interpreter', 'latex', 'FontSize', 14)
    else
        assert('invalid option')
    end
    set(gca,'FontSize', 14);
end
set(gcf,'Position', [x0 y0 width height])

%% Save
if sParams.sSim.twomoons_dataset
    simPrefix = strcat('Two_moons_', num2str(sParams.nysRatio*100, '%d'), 'prec');
else
    simPrefix = strcat('Gaussian_', num2str(sParams.nysRatio*100, '%d'), 'prec');
end
if ~exist(sParams.sSim.outputFolder, 'dir')
    mkdir(sParams.sSim.outputFolder)
end

saveas(fig,strcat(sParams.sSim.outputFolder, filesep, simPrefix, '_eigenvectors_2d_m_', num2str(firstEigenIdx), '_to_', num2str(lastEigIdx), '_', figName), 'epsc');


end