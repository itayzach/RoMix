% clear; clc; close all;
% rng('default')

n = 5000;
mu = 0;
sigma = 1;
pDegree = 10;
x = sigma*randn(n,1) + mu;

%% KDE
pts = linspace(min(x),max(x),n/10);
[cdfKDE,xGridKDE] = ksdensity(x, pts, 'Function', 'cdf'); 
pCdfKDE = polyfit(xGridKDE,cdfKDE,pDegree);

%% ecdf
[cdfECDF, xGridECDF] = ecdf(x);
pCdfECDF = polyfit(xGridECDF,cdfECDF,pDegree);

%% hist
% h = histogram(x,n/10); 
% cdfHist = cumsum(h.Values/n);
% xGridHist = h.BinEdges(1:end-1);
% close;
[pdfHist, histBins] = histcounts(x,n/10);
cdfHist = cumsum(pdfHist/n);
xGridHist = histBins(1:end-1);
pCdfHist = polyfit(xGridHist,cdfHist,pDegree);

%% normal cdf
xAxis = min(x):1e-5:max(x);
normalCdf = cdf('Normal',xAxis, mu, sigma);

figure;
subplot(2,1,1)
plot(xAxis, polyval(pCdfKDE,xAxis));
hold on
plot(xAxis, polyval(pCdfECDF,xAxis));
plot(xAxis, polyval(pCdfHist,xAxis));
plot(xAxis, normalCdf, '.');
legend('normal cdf', 'pCdfKDE', 'pCdfECDF', 'pCdfHist', 'location', 'northwest');

subplot(2,1,2)
plot(xAxis, abs(normalCdf-polyval(pCdfKDE,xAxis)));
hold on
plot(xAxis, abs(normalCdf-polyval(pCdfECDF,xAxis)));
plot(xAxis, abs(normalCdf-polyval(pCdfHist,xAxis)));
legend('pCdfKDE', 'pCdfECDF', 'pCdfHist', 'location', 'northwest');

fprintf('KDE:  %f\n', norm(normalCdf-polyval(pCdfKDE,xAxis))/norm(normalCdf));
fprintf('ECDF: %f\n', norm(normalCdf-polyval(pCdfECDF,xAxis))/norm(normalCdf));
fprintf('Hist: %f\n', norm(normalCdf-polyval(pCdfHist,xAxis))/norm(normalCdf));

% muTilde = mu;
% sigmaTilde = sigma;
% xTilde = icdf('Normal',polyval(pCdfKDE,x),muTilde,sigmaTilde);
% [~, sortInd] = sort(xTilde);
% % sortInd = sortInd(round(n*0.1):round(n*0.9));
% figure;
% tiledlayout(2,1)
% ax1 = nexttile;
% plot(x(sortInd),x(sortInd), 'o')
% hold on;
% plot(x(sortInd),xTilde(sortInd),'.');
% set(gca,'FontSize', 14);
% legend('$x$', '$\tilde{x}$', 'interpreter', 'latex', 'fontsize', 14, 'location', 'northwest')
% ax2 = nexttile;
% plot(x(sortInd),abs(xTilde(sortInd) - x(sortInd)), '.')
% legend('$|x-\tilde{x}|$', 'interpreter', 'latex', 'fontsize', 14, 'location', 'northwest')
% set(gca,'FontSize', 14);
% % linkaxes([ax1 ax2],'xy')