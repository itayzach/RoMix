function [] = PlotGaussianSanity(xTrain, xTildeTrain, muTilde, sigmaTilde, b_debugUseNormalCDF, ...
    xTrainGrid, estMarginalCdf_xTrain, polyvalCdfMatrix)
dim = size(xTrain,2);
[~, sortInd] = sort(xTildeTrain);
figure;
tiledlayout(2,1)
ax1 = nexttile;
plot(xTrain(sortInd),xTrain(sortInd), 'o')
hold on;
plot(xTrain(sortInd),xTildeTrain(sortInd),'.');
set(gca,'FontSize', 14);
legend('$x$', '$\tilde{x}$', 'interpreter', 'latex', 'fontsize', 14, 'location', 'southeast')
ax2 = nexttile;
plot(xTrain(sortInd),abs(xTildeTrain(sortInd) - xTrain(sortInd)), '.')
set(gca,'FontSize', 14);
linkaxes([ax1 ax2],'xy')
if ~b_debugUseNormalCDF
    for d=1:dim
        normalCdf = cdf('Normal',xTrain(:,d),muTilde(d),sigmaTilde(d,d));
        figure;
        plot(xTrain(:,d), normalCdf, 'o');
        hold on;
        plot(xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), '.');
        plot(xTrain(:,d), polyvalCdfMatrix(:,d), '.');
        legend('cdf(normal)', 'ecdf', 'ours', 'Location', 'southeast');
        set(gca,'FontSize', 14);
        
        figure;
        plot(xTrain(:,d), abs(polyvalCdfMatrix(:,d)-normalCdf), '.');
    end
end
sgtitle('Gaussian sanity check', 'interpreter', 'latex', 'FontSize',16)
end