function [] = PlotTwoMoonsLapRLS(sSimParams, sDataset, sDistanceParams, omega, mAlpha,gamma1Rep, gamma2Rep,test_error)

classifier.alpha = mAlpha;
classifier.gammas = [gamma1Rep, gamma2Rep];
classifier.test_error = 100 - test_error;

x0     = 10;
y0     = 250;
width  = 550;
height = 700;

alpha_laprls = classifier.alpha;
xTrain = sDataset.sData.x;
yTrain = sDataset.sData.y;
xTest = sDataset.sData.xt;
xMax = max(max([xTrain; xTest],[],1));
xMin = min(min([xTrain; xTest],[],1));
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);

X=[XX1(:) XX2(:)];
tic;
%     mK_xTrain_X = calckernel(classifier.Kernel,classifier.KernelParam, xTrain_laprls, X);
dist = CalcDistance(X,xTrain, sDistanceParams);
mK_xTrain_X = exp(-dist.^2/(2*omega^2));
run_time = toc;
fprintf('K time = %.2f ms\n', run_time*1e3);

%     mK_xTrain_xTrain = calckernel(classifier.Kernel,classifier.KernelParam, xTrain, xTrain);
dist = CalcDistance(xTrain, xTrain, sDistanceParams);
mK_xTrain_xTrain = exp(-dist.^2/(2*omega^2));

vKa = mK_xTrain_X*alpha_laprls;
vKa_train = mK_xTrain_xTrain*alpha_laprls;
mKa = reshape(vKa,length(x1),length(x2));

fig = figure;
sgtitle(['LapRLS with $n$ = ' num2str(length(classifier.alpha)) newline ...
    '$\gamma_A = ' num2str(classifier.gammas(1), '%.4f') '$ $\gamma_I = ' num2str(classifier.gammas(2), '%.4f') '$' newline ...
    'Test error = ' num2str(classifier.test_error, '%.1f'), '$\%$' ], ...
    'Interpreter', 'latex', 'FontSize', 14);
subplot(2,1,1)
surf(XX1,XX2,mKa, 'edgecolor', 'none')
hold on;
pos=find(yTrain==1);
neg=find(yTrain==-1);
unlab=find(yTrain==0);
scatter3(xTrain(unlab,1), xTrain(unlab,2), vKa_train(unlab), 'ks');
scatter3(xTrain(pos,1),xTrain(pos,2), vKa_train(pos),200, 'rd','MarkerFaceColor','r','MarkerEdgeColor','k'); hold on;
scatter3(xTrain(neg,1),xTrain(neg,2), vKa_train(neg),200, 'bo' ,'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;

xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
colorbar;
set(gca,'FontSize', 14);
subplot(2,1,2)
contourf(XX1,XX2,mKa,[0 0]);shading flat;
hold on;
plot2D(xTrain,yTrain,15,'ks');
plot2D(xTest,zeros(size(sDataset.sData.xt,1),1),15,'k*');

set(gca,'FontSize', 14);
set(gcf,'Position',[x0+width y0 width height])
if isfield(sSimParams, 'output_folder') && ~isempty(sSimParams.output_folder)
    saveas(fig,[sSimParams.output_folder filesep 'fig_' classifier.Name], 'epsc');
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