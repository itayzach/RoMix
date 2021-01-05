%% Restart
clc; clear; close all;
rng('default')
%% (Source) uniform data
% --------------------------------------------------------------------------------------------------
% Generate data
% --------------------------------------------------------------------------------------------------
n = 1000;
xMax = 4;
xMin = -3;
x = (xMax - xMin)*rand(n,1) + xMin;

figure('Name', 'Histogram of X'); 
histogram(x,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

% --------------------------------------------------------------------------------------------------
% Generate axis
% --------------------------------------------------------------------------------------------------
N = 2000;
axis = linspace(xMin,xMax,N)';
% --------------------------------------------------------------------------------------------------
% Estimate x CDF
% --------------------------------------------------------------------------------------------------
[estCdf, x_new] = ecdf(x);
x_new = x_new(1:end-1);
estCdf = estCdf(1:end-1);
pCdf = polyfit(x_new,estCdf,5); % to have analytic expression for the cdf

polyCdf = polyval(pCdf, axis);
% --------------------------------------------------------------------------------------------------
% Plot x estimated CDF vs. analytic x poly CDF
% --------------------------------------------------------------------------------------------------
figure('Name', 'Verify estCdfPoly');
plot(x_new, estCdf,'.');
hold on;
plot(axis, polyCdf, '.');

title('ecdf(x) vs. polynomial estimation of CDF', 'interpreter', 'latex', 'FontSize', 16);
legend('ecdf(x)', 'polyvalCdf(axis)', ...
    'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14)
set(gca,'FontSize', 14);

%% (Target) Gaussian data
% --------------------------------------------------------------------------------------------------
% Get Gaussian CDF on a grid
% --------------------------------------------------------------------------------------------------
mu = 0;
sigma = 1;
cdfTilde = cdf('Normal',axis,mu,sigma);

% --------------------------------------------------------------------------------------------------
% Plot Gaussian CDF
% --------------------------------------------------------------------------------------------------
figure('Name', 'Gaussian cdf'); 
plot(axis,cdfTilde,'.');
title('Gaussian cdf', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);

%% Polynomial fit between CDF(x) and CDF(xTilde)
T = polyfit(polyCdf, cdfTilde,7);
Tinv = polyfit(cdfTilde, polyCdf,7);
%% Plot fit
dx = 0.001;
zeroOneAxis = (0:dx:1-dx)'; % possible values of a CDF are always [0,1]
gaussCdfZeroOne = polyval(T, zeroOneAxis);

figure('Name', 'Polyval(T,[0,1])');
plot(polyCdf, cdfTilde,'o')
hold on;
plot(zeroOneAxis,gaussCdfZeroOne,'.')
title('polyval(T, [0,1]) - should be Gaussian cdf', 'interpreter', 'latex', 'FontSize', 16);
legend('polyCdf $\to$ cdfTilde', '$T([0,1])$', ...
    'location', 'southeast', 'interpreter', 'latex', 'FontSize', 14); 
set(gca,'FontSize', 14);

icdfZeroOne = icdf('Normal',gaussCdfZeroOne,mu,sigma);

figure('Name', 'icdf([0,1])'); 
% histfit(icdfZeroOne,100);
plot(gaussCdfZeroOne,icdfZeroOne,'.')
title('icdf([0,1])', 'interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);


