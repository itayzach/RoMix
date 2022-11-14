function [] = PlotClassifier(sPlotParams, sDataset, sClassifier)

fig = figure;
%sgtitle(['EigRLS with $M$ = ' num2str(sClassifier.MTilde) newline ...
%    '$\omega = $ ' num2str(sClassifier.omega) ';   $\gamma_A = ' num2str(sClassifier.gamma_A, '%.4f') '$ $\gamma_I = ' num2str(sClassifier.gamma_I, '%.4f') '$'], ...
%    'Interpreter', 'latex', 'FontSize', 14);
%subplot(2,1,1)
surf(sClassifier.XX1, sClassifier.XX2, sClassifier.mPhi_X_c, 'edgecolor', 'none')
colormap('hot')
view(-25,60)
hold on;
pos=find(sDataset.sData.y==1);
neg=find(sDataset.sData.y==-1);
unlab=find(sDataset.sData.y==0);
scatter3(sDataset.sData.x(unlab,1), sDataset.sData.x(unlab,2), sClassifier.vPhi_xTrain_c(unlab), [],'filled');
scatter3(sDataset.sData.x(pos,1),sDataset.sData.x(pos,2), sClassifier.vPhi_xTrain_c(pos),200, 'rd','MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
scatter3(sDataset.sData.x(neg,1),sDataset.sData.x(neg,2), sClassifier.vPhi_xTrain_c(neg),200, 'go' ,'MarkerFaceColor','g','MarkerEdgeColor','k'); hold on;
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
zlim([-1 1]);
h = colorbar('TickLabelInterpreter', 'latex');
%set(fig,'renderer','Painters')
set(gca,'FontSize', 14);
x0     = 10;
y0     = 100;
width  = 600;
height = 350;
set(gcf,'Position',[x0 y0 width height])
if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'classifier_surf'; 
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end

fig = figure;
contourf(sClassifier.XX1,sClassifier.XX2,sClassifier.mPhi_X_c,[0 0]);shading flat;
colormap('hot')
%h = colorbar('TickLabelInterpreter', 'latex');
%h.TickLabels = [-1, 1];
%h.Ticks = [ -0.5, 0.5];
hold on;
plot2D(sDataset.sData.x,sDataset.sData.y,15,'ks');
%plot2D(sDataset.sData.xt,zeros(size(sDataset.sData.xt,1),1),15,'k*');
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
set(fig,'renderer','Painters')
set(gca,'FontSize', 14);
x0     = 10;
y0     = 100;
width  = 450;
height = 300;
set(gcf,'Position',[x0+width y0 width height])

if ~isempty(sPlotParams) && isfield(sPlotParams, 'outputFolder')
    figName = 'classifier'; 
    SaveFigure(sPlotParams, fig, figName, {'epsc', 'png'});
end
end

function plot2D(X,Y,markersize, marker)
pos=find(Y==1);
neg=find(Y==-1);
unlab=find(Y==0);

scatter(X(unlab,1),X(unlab,2),'filled'); hold on;     
plot(X(pos,1),X(pos,2),'rd','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
plot(X(neg,1),X(neg,2),'bo' ,'MarkerSize',markersize,'MarkerFaceColor','g','MarkerEdgeColor','k'); hold on;
end