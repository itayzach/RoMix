function PlotManifoldExample()
xSamples = -3:2:3;
ySamples = -3:2:3;
%xSamples = -4 + (4+4)*rand(3,1);
%ySamples = -4 + (4+4)*rand(3,1);

[XSamples, YSamples] = meshgrid(xSamples,ySamples);
ZSamples = Manifold(XSamples,YSamples);

x = -4:0.1:4;
y = -4:0.1:4;
[X, Y] = meshgrid(x,y);
Z = Manifold(X,Y);
sPlotParams = GetPlotParams();


fig = figure;
surf(X, Y, Z,"EdgeColor","none")
hold on
scatter3(XSamples(:), YSamples(:), ZSamples(:)+0.05,[],(1:length(XSamples(:))),'filled')
oldcmap = colormap("bone");
colormap( flipud(oldcmap) );
line([-4 -4],[-4  4],[-2 -2],'color','black')
line([-4 -4],[ 4  4],[-2 10],'color','black')
line([-4  4],[ 4  4],[10 10],'color','black')
line([ 4  4],[ 4 -4],[10 10],'color','black')
line([ 4  4],[-4 -4],[10 -2],'color','black')
line([ 4 -4],[-4 -4],[-2 -2],'color','black')
zlim([-2 10])
view(-60,36)
text(-3.5,3.5,0,'$\mathcal{X}$','Interpreter','latex', 'FontSize',16)
text(-2.1,-2.8,4,'$\mathcal{M}$','Interpreter','latex', 'FontSize',16)
text(-2.8,3.6,3.8,'$v_1$','Interpreter','latex', 'FontSize',14)
text(3.1,-3,3.7,'$v_2$','Interpreter','latex', 'FontSize',14)
grid off
set(gca,'XColor', 'none','YColor','none','ZColor','none')
set(gcf, 'PaperPosition', [0 0 7 5]); %Position plot at left hand corner with width 5 and height 5. 
set(gcf, 'PaperSize', [7 5]); %Set the paper to have width 5 and height 5. 
%set(fig,'renderer','Painters')

SaveFigure(sPlotParams, fig, 'Manifold', {'eps','png'})
% print('-painters','-depsc',['figs' filesep 'eps' filesep 'Manifold'])

function z = Manifold(x,y)
 z = cos(x) + cos(y)+5;
end
end