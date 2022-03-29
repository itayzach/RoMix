function PlotCovEigs(sPlotParams, sDistParams)
windowStyle = get(0,'DefaultFigureWindowStyle');
set(0,'DefaultFigureWindowStyle','normal')

mCovEigs = cell2mat(sDistParams.sigma');
[origC, origD] = size(mCovEigs);
minEig = min(mCovEigs(:));
lastPrincipalDims = find(mCovEigs(1,:) > 1.05*minEig,1,'last');
lastPrincipalComp = find(mCovEigs(:,1) > 1.05*minEig,1,'last');

mSigma = mCovEigs(1:lastPrincipalComp,1:lastPrincipalDims);
[C, dim] = size(mSigma);

fig = figure('Name', 'Covariance eigenvalues');
% lines() returns 7 unique colors. Make sure dim <= 7:
if dim <= size(unique(lines(dim),'rows'),1)
    cmap = colormap(lines(dim));
else
    cmap = colormap(jet(dim));
end
for d=1:dim
    scatter((1:C)', mSigma(:,d), [], cmap(d,:), 'filled');
    hold on
end
xlabel('$c$','Interpreter', 'latex', 'FontSize', 14)
ylabel('$\sigma^{(d)}_c$', 'Interpreter', 'latex', 'FontSize', 14)
grid on;
h = colorbar('TickLabelInterpreter', 'latex');
h.TickLabels = 1:dim;
h.Ticks = (0.5:dim)/dim;

set(get(h,'label'),'string','$d$','interpreter','latex','Rotation',0,'FontSize', 16);
h.Label.Position(1) = 0.5;
h.Label.Position(2) = 0;
h.Label.Position(3) = 0;
h.Label.FontSize = 14;

title('Cov eigs $> 1.05\min(\sigma)$', 'Interpreter', 'latex', 'FontSize', 14);
set(gca,'FontSize', 14);
set(fig,'renderer','Painters')
set(0,'DefaultFigureWindowStyle',windowStyle)
%% Save
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'CovEigs';
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end