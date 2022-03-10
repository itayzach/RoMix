function PlotCovEigs(sPlotParams, sDistParams)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

mCovEigs = cell2mat(sDistParams.sigma');
minEig = min(mCovEigs(:));
vHowManyMoreThanMinDim = sum(mCovEigs - minEig > 0.05*minEig);
vPrincipalDims = vHowManyMoreThanMinDim > 0.05*sDistParams.estNumComponents;
mSigma = mCovEigs(:,vPrincipalDims);

fig = figure('Name', 'Covariance eigenvalues');
[C, dim] = size(mSigma);
scatter3(repmat((1:C)',dim,1), mSigma(:), repelem((1:dim)',C), [], repelem((1:dim)',C),'filled');
xlabel('$c$','Interpreter', 'latex', 'FontSize', 14)
ylabel('$\sigma^{(d)}_c$', 'Interpreter', 'latex', 'FontSize', 14)
% colormap(lines(dim));
colormap(jet(dim));
h = colorbar();
set(get(h,'label'),'string','$d$','interpreter','latex','Rotation',0,'FontSize', 16);
h.Label.Position(1) = 0.5;
h.Label.Position(2) = 1;
h.Label.Position(3) = 0;
view(2)

title('Cov eigs', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(0,'DefaultFigureWindowStyle',windowStyle)
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'CovEigs';
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end