function fig = PlotHistogram(sSimParams, v, actualDataDist, plt_title, b_histfit)

x0     = 10;
y0     = 50;
width  = 600;
height = 400;

if ~exist('b_histfit', 'var')
    b_histfit = false;
end

dim = size(v,2);

fig = figure('Name', sprintf('%d-D histogram', dim));
if dim == 1
    if b_histfit
        histfit(v,100);
    else
        histogram(v, 100);
    end
elseif dim == 2
    hist3(v,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
    colormap(gca, 'hot')
    colorbar()
    view(2)
    xlim([ min(v(:,1)) max(v(:,1))])
    ylim([ min(v(:,2)) max(v(:,2))])
end
title(strcat(plt_title, " (", actualDataDist, ")"), 'Interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);
set(gcf,'Position', [x0 y0 width height])

%% Save
if ~exist(sSimParams.outputFolder, 'dir')
    mkdir(sSimParams.outputFolder)
end

saveas(fig,strcat(sSimParams.outputFolder, filesep, actualDataDist, num2str(dim), 'd', '_histogram'), 'epsc');
end