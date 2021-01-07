%% Restart
clc; clear; close all;
rng('default')
%% (Source) uniform data
% --------------------------------------------------------------------------------------------------
% Generate data
% --------------------------------------------------------------------------------------------------
n = 2000;
xMax = 10;
xMin = -10;

x = (xMax - xMin)*rand(n,1) + xMin;

figure('Name', 'Histogram of X'); 
histogram(x,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

% --------------------------------------------------------------------------------------------------
% Learn CDF(x) from x
% --------------------------------------------------------------------------------------------------
[estCdf, x_ecdf] = ecdf(x);
x_ecdf = x_ecdf(1:end-1);
estCdf = estCdf(1:end-1);

% --------------------------------------------------------------------------------------------------
% Fit ecdf to a polynomial
% --------------------------------------------------------------------------------------------------
pCdfDegree = 10;
invpCdfDegree = 10;
pCdf = polyfit(x_ecdf, estCdf, pCdfDegree); % to have analytic expression for the cdf
invpCdf = polyfit(estCdf, x_ecdf, invpCdfDegree); % to have analytic expression for the cdf

% --------------------------------------------------------------------------------------------------
% Verification - Plot x estimated CDF vs. analytic x poly CDF
% --------------------------------------------------------------------------------------------------
N = 2000;
axis = linspace(xMin,xMax,N)';
polyCdf = polyval(pCdf, axis);
b_saturate = false;
if b_saturate && (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf(polyCdf > 1) = 1-1e-10; % saturate
    polyCdf(polyCdf < 0) = 1e-10; % saturate
end

figure('Name', 'Verify estCdfPoly');
plot(x_ecdf, estCdf,'.');
hold on;
plot(axis, polyCdf, '.');
title('ecdf(x) vs. polynomial estimation of CDF', 'interpreter', 'latex', 'FontSize', 16);
legend('ecdf(x)', 'polyvalCdf(axis)', ...
    'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

%% estCdf -> Gaussian iCDF
mu = 0;
sigma = 5;
xTilde = icdf('Normal',polyCdf,mu,sigma);

figure('Name', 'xTilde = T(x) = icdf(polyCdf(x))'); 
plot(polyCdf,xTilde,'.')
title('$\tilde{x} = T(x) = {\bf icdf}({\bf polyCdf}(x))$', 'interpreter', 'latex', 'FontSize', 16);
xlabel('${\bf polyCdf}(x)$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

figure('Name', 'xTilde hist'); 
histfit(xTilde ,100);
title('$\tilde{x}$ histogram', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

%% Show test points transformation
nTestPoints = 50;
xTest = linspace(xMin+0.1,xMax-0.1,nTestPoints)';%(xMax - xMin)*rand(nTestPoints,1) + xMin;

polyCdf = polyval(pCdf, xTest);
if b_saturate && (any(polyCdf > 1) || any(polyCdf < 0))
    warning('CDF must be in [0,1], saturating...');
    polyCdf(polyCdf > 1) = 1 - eps; % saturate
    polyCdf(polyCdf < 0) = eps; % saturate
end

xTestTilde = icdf('Normal',polyCdf,mu,sigma);

xTildeCdf = cdf('Normal', xTestTilde, mu, sigma);
x_new_est = polyval(invpCdf, xTildeCdf);

sz = 25;
cmap = xTest;

figure('Name', 'x & xTilde');
subplot(2,1,1)
scatter(xTest, zeros(1,nTestPoints), 100, cmap, 'o')
hold on;
scatter(x_new_est, zeros(1,nTestPoints), 50, cmap, 'filled')
colormap('jet')
xlabel('$x$', 'interpreter', 'latex', 'FontSize', 16);
legend('$x$', '$T^{-1}(T(x))$','interpreter', 'latex', 'FontSize', 14);
title('Original nodes','interpreter', 'latex', 'FontSize', 16);
set(gca,'YTick',[],'FontSize', 14);

subplot(2,1,2);
scatter(xTestTilde, zeros(1,nTestPoints), 50, cmap, 'filled')
colormap('jet')
xlabel('$\tilde{x}$', 'interpreter', 'latex', 'FontSize', 16);
title('Transformed nodes','interpreter', 'latex', 'FontSize', 16);
set(gca,'YTick',[],'FontSize', 14);