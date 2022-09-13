function PlotSimpleSin()
f0 = 1e3;
fs = 100e3;
nCyc = 3;
t = 0:1/fs:nCyc/f0-1/fs;
%nSamples = nCyc*fs/f0;

y = 1 + sin(2*pi*f0*t);% + 0.05*randn(1,nSamples);

fig = figure;
hold on;
%yline(0,'k-','linewidth',2);
plot(t,y,'k','linewidth',2); 
h = stem(t(6:10:end),y(6:10:end),'filled','LineStyle','none');
h.Color = 'black';
h.MarkerFaceColor = 'white';
set(gca,'XColor', 'none','YColor','none')
x0     = 10;
y0     = 50;
height = 200;
width  = 600;
set(gcf,'Position', [x0 y0 width height])

sPlotParams.outputFolder = 'figs';
SaveFigure(sPlotParams,fig,'SimpleSin',{'epsc', 'png'})
end