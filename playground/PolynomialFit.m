%% Restart
clc; clear; close all;
rng('default')
%% Parameters
muTilde = 0;
sigmaTilde = 5;
pCdfDegree = 2;
invpCdfDegree = 2;
%% (Source) uniform data
% --------------------------------------------------------------------------------------------------
% Generate data
% --------------------------------------------------------------------------------------------------
n = 2000;
xMax = 10;
xMin = -10;

xTrain = (xMax - xMin)*rand(n,1) + xMin;

figure('Name', 'Histogram of xTrain');
histogram(xTrain,100);
title('Histogram of $x_{{\bf train}}$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

% --------------------------------------------------------------------------------------------------
% Learn CDF(x) from xTrain by fitting ecdf to a polynomial
% --------------------------------------------------------------------------------------------------
[estCdf_xTrainGrid, xTrainGrid] = ecdf(xTrain);
xTrainGrid = xTrainGrid(1:end-1);
estCdf_xTrainGrid = estCdf_xTrainGrid(1:end-1);

pCdf = polyfit(xTrainGrid, estCdf_xTrainGrid, pCdfDegree); % to have analytic expression for the cdf
invpCdf = polyfit(estCdf_xTrainGrid, xTrainGrid, invpCdfDegree); % to have analytic expression for the cdf

% --------------------------------------------------------------------------------------------------
% Verification - Plot estCdf(xTrain) vs. polyCdf(xTest)
% --------------------------------------------------------------------------------------------------
N = 2000;
xTestGrid = linspace(xMin,xMax,N)';
polyCdf_xTestGrid = polyval(pCdf, xTestGrid);
b_saturate = false;
if b_saturate && (any(polyCdf_xTestGrid > 1) || any(polyCdf_xTestGrid < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf_xTestGrid(polyCdf_xTestGrid > 1) = 1-eps; % saturate
    polyCdf_xTestGrid(polyCdf_xTestGrid < 0) = eps; % saturate
end

figure('Name', 'Verify estCdfPoly');
plot(xTrainGrid, estCdf_xTrainGrid,'.');
hold on;
plot(xTestGrid, polyCdf_xTestGrid, '.');
title('${\bf ecdf}(x_{{\bf train}}$) vs. ${\bf polyCdf}(x_{{\bf test}})$', 'interpreter', 'latex', 'FontSize', 16);
legend('${\bf ecdf}(x_{{\bf train}})$', '${\bf polyCdf}(x_{{\bf test}})$', ...
    'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

%% estCdf -> Gaussian iCDF
xTildeTestGrid = icdf('Normal',polyCdf_xTestGrid,muTilde,sigmaTilde);

figure('Name', 'xTilde = T(xTest) = icdf(polyCdf(xTest))');
plot(polyCdf_xTestGrid,xTildeTestGrid,'.')
title('$\tilde{x} = T(x_{{\bf test}}) = {\bf icdf}({\bf polyCdf}(x_{{\bf test}}))$', 'interpreter', 'latex', 'FontSize', 16);
xlabel('${\bf polyCdf}(x_{{\bf test}})$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

figure('Name', 'Histogram of xTildeTestGrid');
histfit(xTildeTestGrid ,100);
title('Histogram of $\tilde{x}_{{\bf test}}$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

%% Show test points transformation
nTestPoints = 50;
xSmallTest = linspace(xMin+0.5,xMax-0.5,nTestPoints)';%(xMax - xMin)*rand(nTestPoints,1) + xMin;

% --------------------------------------------------------------------------------------------------
% T(x)
% --------------------------------------------------------------------------------------------------
xTestTilde = T(pCdf, b_saturate, muTilde, sigmaTilde, xSmallTest);

% --------------------------------------------------------------------------------------------------
% T^{-1}(x)
% --------------------------------------------------------------------------------------------------
xSmallTestEst = invT(invpCdf, muTilde, sigmaTilde, xTestTilde);

% --------------------------------------------------------------------------------------------------
% Plot
% --------------------------------------------------------------------------------------------------
sz = 25;
cmap = xSmallTest;
figure('Name', 'x <-> xTilde');
subplot(2,1,1)
scatter(xSmallTest, zeros(1,nTestPoints), 100, cmap, 'o')
hold on;
scatter(xSmallTestEst, zeros(1,nTestPoints), 50, cmap, 'filled')
colormap('jet')
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
legend('$x_{{\bf test}}$', '$T^{-1}(T(x_{{\bf test}}))$','interpreter', 'latex', 'FontSize', 14);
title('Original nodes','interpreter', 'latex', 'FontSize', 16);
set(gca,'YTick',[],'FontSize', 14);

subplot(2,1,2);
scatter(xTestTilde, zeros(1,nTestPoints), 50, cmap, 'filled')
colormap('jet')
xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
set(gca,'YTick',[],'FontSize', 14);

%% T(x)
function xTilde = T(pCdf, b_saturate, mu, sigma, x)
polyCdf = polyval(pCdf, x);
if (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1]...');
    if b_saturate
        polyCdf(polyCdf > 1) = 1 - eps; % saturate
        polyCdf(polyCdf < 0) = eps; % saturate
    end
end
xTilde = icdf('Normal',polyCdf,mu,sigma);
assert(~any(isnan(xTilde)),'xTilde contain NaNs...');
end

%% invT(x)
function x = invT(invpCdf, mu, sigma, xTilde)
assert(~any(isnan(xTilde)),'xTilde contain NaNs...');
xTildeCdf = cdf('Normal', xTilde, mu, sigma);
x = polyval(invpCdf, xTildeCdf);
assert(~any(isnan(x)),'x contain NaNs...');
end