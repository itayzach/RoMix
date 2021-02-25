%% Restart
clc; clear; close all;

%% Uniform data
n = 1000;

yMax = 1;
yMin = 0;
u1 = linspace(yMin,yMax,n);(yMax - yMin)*rand(n,1) + yMin;
u2 = (yMax - yMin)*rand(n,1) + yMin;

%% Box-Muller
z1 = sqrt(-2*log(u1)).*cos(2*pi*u2);
z2 = sqrt(-2*log(u1)).*sin(2*pi*u2);
z = z1 + z2;
%% Plot
figure(); 
subplot(2,2,1)
histogram(u1,100);
title('Histogram of $u_1$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(2,2,2)
histogram(u2,100);
title('Histogram of $u_2$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

subplot(2,2,3)
histogram(z1,100);
title('Histogram of $z_1$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(2,2,4)
histogram(z2,100);
title('Histogram of $z_2$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
% 
% figure;
% histogram(z,100);
% title('Histogram of $z$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);