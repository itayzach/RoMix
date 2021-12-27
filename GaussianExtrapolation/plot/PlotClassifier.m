function [] = PlotClassifier(sSimParams, sDataset, sClassifier)

x0     = 10;
y0     = 250;
width  = 550;
height = 700;

fig = figure;
sgtitle(['EigRLS with $M$ = ' num2str(sClassifier.MTilde) newline ...
    '$\omega = $ ' num2str(sClassifier.omega) ';   $\gamma_A = ' num2str(sClassifier.gamma_A, '%.4f') '$ $\gamma_I = ' num2str(sClassifier.gamma_I, '%.4f') '$' newline ...
    'Test error = ' num2str(sClassifier.error, '%.1f'), '$\%$' ], ...
    'Interpreter', 'latex', 'FontSize', 14);
subplot(2,1,1)
surf(sClassifier.XX1, sClassifier.XX2, sClassifier.mPhi_X_c, 'edgecolor', 'none')
hold on;
pos=find(sDataset.sData.y==1);
neg=find(sDataset.sData.y==-1);
unlab=find(sDataset.sData.y==0);
scatter3(sDataset.sData.x(unlab,1), sDataset.sData.x(unlab,2), sClassifier.vPhi_xTrain_c(unlab), 'ks');
scatter3(sDataset.sData.x(pos,1),sDataset.sData.x(pos,2), sClassifier.vPhi_xTrain_c(pos),200, 'rd','MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
scatter3(sDataset.sData.x(neg,1),sDataset.sData.x(neg,2), sClassifier.vPhi_xTrain_c(neg),200, 'bo' ,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
colorbar;
set(gca,'FontSize', 14);
subplot(2,1,2)
contourf(sClassifier.XX1,sClassifier.XX2,sClassifier.mPhi_X_c,[0 0]);shading flat;
hold on;
plot2D(sDataset.sData.x,sDataset.sData.y,15,'ks');
plot2D(sDataset.sData.xt,zeros(size(sDataset.sData.xt,1),1),15,'k*');
set(gca,'FontSize', 14);
set(gcf,'Position',[x0+width+width y0 width height])
if isfield(sSimParams, 'outputFolder')
    saveas(fig,[sSimParams.outputFolder filesep 'fig_eigrls_M_' num2str(sClassifier.MTilde)], 'epsc');
end
end

function plot2D(X,Y,markersize, marker)
pos=find(Y==1);
neg=find(Y==-1);
unlab=find(Y==0);

plot(X(unlab,1),X(unlab,2),marker); hold on;     
plot(X(pos,1),X(pos,2),'rd','MarkerSize',markersize,'MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
plot(X(neg,1),X(neg,2),'bo' ,'MarkerSize',markersize,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
end