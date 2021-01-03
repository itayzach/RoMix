%% Restart
clc; clear; close all;

%% Number of points
n = 1000;

%% Gaussian data
sigma = 1; mu = 0;
x = sigma*randn(n,1) + mu;

figure('Name', 'Histogram of X'); 
subplot(2,1,1)
histogram(x,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(2,1,2)
plot(x,'.')
xlabel('$i$','interpreter', 'latex')
ylabel('$x(i)$','interpreter', 'latex')
set(gca,'FontSize', 14);

%% Uniform data
yMax = 1;
yMin = -1;
y = (yMax - yMin)*rand(n,1) + yMin;
figure('Name', 'Histogram of Y'); 
subplot(2,1,1)
histogram(y,100);
title('Histogram of $Y$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(2,1,2)
plot(y,'.')
xlabel('$i$','interpreter', 'latex')
ylabel('$y(i)$','interpreter', 'latex')
set(gca,'FontSize', 14);

%% Polynomial fit
p = polyfit(x,y,7);

%% Plot fit
x1 = linspace(-3*sigma+mu,3*sigma+mu);
y1 = polyval(p,x1);
figure
plot(x,y,'o')
hold on
plot(x1,y1)
hold off