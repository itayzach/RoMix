%% Restart
clc; clear; close all;

%% Number of points
n = 1000;

%% Uniform data
xMax = 1;
xMin = -1;
x = (xMax - xMin)*rand(n,1) + xMin;
figure('Name', 'Histogram of Y'); 
histogram(x,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
set(gca,'FontSize', 14);

[xCdf, x2] = ecdf(x);
figure; plot(x2,xCdf);
x2 = x2(2:end);
xCdf = xCdf(2:end);
p = polyfit(x2,xCdf,5);
%% Gaussian data
sigma = 1; mu = 0;
xTilde = sigma*randn(n,1) + mu;

figure('Name', 'Histogram of X'); 
subplot(2,1,1)
histogram(xTilde,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(2,1,2)
plot(xTilde,'.')
xlabel('$i$','interpreter', 'latex')
ylabel('$x(i)$','interpreter', 'latex')
set(gca,'FontSize', 14);

gaussianCDF = cdf('Normal',xTilde,mu,sigma);
figure; plot(xTilde,gaussianCDF,'.');


%% Polynomial fit
p_yCDF_to_gaussianCDF = polyfit(xCdf,gaussianCDF,7);
pinverse = polyfit(gaussianCDF,xCdf,7);
%% Plot fit
y_interp = polyval(p, -1:0.001:1);
figure; plot(-1:0.001:1, y_interp, '.');

gaussianCDF_interp = polyval(p_yCDF_to_gaussianCDF, y_interp);
gaussianNodes_interp = icdf('Normal',gaussianCDF_interp,mu,sigma);

figure; histogram(gaussianNodes_interp);