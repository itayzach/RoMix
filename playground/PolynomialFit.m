%% Restart
clc; clear; close all;
rng('default')
%% Number of points
n = 1000;

%% (Source) uniform data
% --------------------------------------------------------------------------------------------------
% Generate data
% --------------------------------------------------------------------------------------------------
xMax = 4;
xMin = -3;
x = (xMax - xMin)*rand(n,1) + xMin;
figure('Name', 'Histogram of X'); 
histogram(x,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

% --------------------------------------------------------------------------------------------------
% Estimate x CDF
% --------------------------------------------------------------------------------------------------
[estCdf, x_new] = ecdf(x);
x_new = x_new(1:end-1);
estCdf = estCdf(1:end-1);
pCdf = polyfit(x_new,estCdf,5); % to have analytic expression for the cdf
polyvalCdf = polyval(pCdf, x);

% --------------------------------------------------------------------------------------------------
% Sort x points
% --------------------------------------------------------------------------------------------------
[xSorted, xSortedInd] = sort(x);
polyvalCdfSorted = polyval(pCdf, xSorted);

% --------------------------------------------------------------------------------------------------
% Plot x estimated CDF vs. analytic x poly CDF
% --------------------------------------------------------------------------------------------------
figure('Name', 'Verify estCdfPoly');
plot(x_new, estCdf,'.');
hold on;
plot(x, polyvalCdf, 'o');
plot(xSorted, polyvalCdfSorted, '.');

title('ecdf(x) vs. polynomial estimation of CDF(x)', 'interpreter', 'latex', 'FontSize', 16);
legend('ecdf(x)', 'polyvalCdf(x)','polyvalCdfSorted(x)', ...
    'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

%% (Target) Gaussian data
% --------------------------------------------------------------------------------------------------
% Generate data
% --------------------------------------------------------------------------------------------------
sigma = 5; mu = 0;
xTilde = sigma*randn(n,1) + mu;
figure('Name', 'Histogram of xTilde'); 
histogram(xTilde,100);
title('Histogram of $\tilde{X}$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

% --------------------------------------------------------------------------------------------------
% Sort xTilde points
% --------------------------------------------------------------------------------------------------
[xTildeSorted, xTildeSortedInd] = sort(xTilde);
cdfTildeSorted = cdf('Normal',xTildeSorted,mu,sigma);

% --------------------------------------------------------------------------------------------------
% Plot Gaussian CDF
% --------------------------------------------------------------------------------------------------
figure('Name', 'Gaussian cdf(xTilde)'); 
plot(xTildeSorted,cdfTildeSorted,'.');
title('Gaussian cdf($\tilde{X}$)', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

%% Polynomial fit between CDF(x) and CDF(xTilde)
T = polyfit(polyvalCdfSorted,cdfTildeSorted,7);
Tinv = polyfit(cdfTildeSorted,polyvalCdfSorted,7);
%% Plot fit
dx = 0.001;
zeroOneAxis = (0:dx:1-dx)'; % possible values of a CDF are always [0,1]
gaussCdfZeroOne = polyval(T, zeroOneAxis);

figure('Name', 'Polyval(T,[0,1])');
plot(zeroOneAxis,gaussCdfZeroOne,'.')
title('polyval(T, [0,1]) - should be Gaussian cdf', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

icdfZeroOne = icdf('Normal',gaussCdfZeroOne,mu,sigma);

figure('Name', 'icdf([0,1])'); 
histfit(icdfZeroOne,100);
title('icdf([0,1]) - should be Gaussian pdf', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);