function fig = PlotHistogram(sSimParams, sDataset)

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

x = [sDataset.sData.x; sDataset.sData.xt];

fig = figure('Name', sprintf('%d-D histogram', sDataset.dim));
if sDataset.dim == 1
    histogram(x, 100);
elseif sDataset.dim == 2
    hist3(x,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
    colormap(gca, 'hot')
    colorbar()
    view(2)
    xlim([ min(x(:,1)) max(x(:,1))])
    ylim([ min(x(:,2)) max(x(:,2))])
end
title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end

saveas(fig,strcat(sSimParams.outputFolder, filesep, sDataset.actualDataDist, num2str(sDataset.dim), 'd', '_histogram'), 'epsc');
end