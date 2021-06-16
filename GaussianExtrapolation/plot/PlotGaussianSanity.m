function [] = PlotGaussianSanity(xTrain, xTildeTrain, muTilde, sigmaTilde, b_debugUseNormalCDF, ...
    xTrainGrid, estMarginalCdf_xTrain, polyvalCdfMatrix, b_kde)
dim = size(xTrain,2);
[~, sortInd] = sort(xTildeTrain);

nFigRows = 1 + ~b_kde*~b_debugUseNormalCDF*dim;

%%
figure('Name', 'Sainty check');
iFig = 1;
subplot(nFigRows,2,iFig)
plot(xTrain(sortInd),xTrain(sortInd), 'o', 'DisplayName', '$x$')
hold on;
plot(xTrain(sortInd),xTildeTrain(sortInd),'.', 'DisplayName', '$\tilde{x}$');
legend('interpreter', 'latex', 'fontsize', 14, 'location', 'southeast')
title('$x$ vs. $\tilde{x}$','interpreter', 'latex', 'fontsize', 14)
set(gca,'FontSize', 14);

%%
iFig = iFig + 1;
subplot(nFigRows,2,iFig)
plot(xTrain(sortInd),abs(xTildeTrain(sortInd) - xTrain(sortInd)), '.')
title('$x$ vs. $\tilde{x}$ error','interpreter', 'latex', 'fontsize', 14)
set(gca,'FontSize', 14);

%%
if ~b_debugUseNormalCDF && ~b_kde
    for d=1:dim
        normalCdf = cdf('Normal',xTrain(:,d),muTilde(d),sigmaTilde(d,d));
        iFig = iFig + 1;
        subplot(nFigRows,2,iFig)
        plot(xTrain(:,d), normalCdf, 'o');
        hold on;
        plot(xTrainGrid(:,d), estMarginalCdf_xTrain(:,d), '.');
        plot(xTrain(:,d), polyvalCdfMatrix(:,d), '.');
        legend('cdf(normal)', 'ecdf', 'ours', 'Location', 'southeast');
        title(['Our CDF vs. normal CDF (d = ', num2str(d), ')'],'interpreter', 'latex', 'fontsize', 14)
        set(gca,'FontSize', 14);
        
        iFig = iFig + 1;
        subplot(nFigRows,2,iFig)
        plot(xTrain(:,d), abs(polyvalCdfMatrix(:,d)-normalCdf), '.');
        title(['Our CDF vs. normal CDF error (d = ', num2str(d), ')'],'interpreter', 'latex', 'fontsize', 14)
        set(gca,'FontSize', 14);
    end
end

sgtitle('Gaussian sanity check', 'interpreter', 'latex', 'FontSize',16)
end