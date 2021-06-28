function [] = PlotTwoMoons(sSimParams, sDataset)
x0     = 10;
y0     = 250;
width  = 600;
height = 400;

fig = figure;
plot2D(sDataset.sData.xt,zeros(size(sDataset.sData.yt)),5,'k*');
plot2D(sDataset.sData.x,sDataset.sData.y,20,'ks');
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
set(gca,'FontSize', 14);
set(gcf,'Position',[x0 y0 width height])
if isfield(sSimParams, 'outputFolder')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_two_moons'], 'epsc');
end
end