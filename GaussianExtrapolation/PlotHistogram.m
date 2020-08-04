function [] = PlotHistogram(sParams)

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

fig = figure;
hist3(sParams.x_rand,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
colormap(gca, 'hot')
colorbar()
view(2)
xlim([ min(sParams.x_rand(:,1)) max(sParams.x_rand(:,1))])
ylim([ min(sParams.x_rand(:,2)) max(sParams.x_rand(:,2))])
title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])

%% Save
if sParams.sSim.twomoons_dataset
    simPrefix = strcat('Two_moons');
else
    simPrefix = strcat('Gaussian');
end
if ~exist(sParams.sSim.outputFolder, 'dir')
    mkdir(sParams.sSim.outputFolder)
end

saveas(fig,strcat(sParams.sSim.outputFolder, filesep, simPrefix, '_histogram_2d'), 'epsc');
end