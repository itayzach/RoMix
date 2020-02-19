close all; clear; clc;

x = (-0.1:0.01:0.1).'; 
l = 0.01; 

exp1 = exp((-x.^2)/(2*l^2));

K = exp1 * exp1.';

[X1, X2] = meshgrid(x);

fig = figure();
fig_left_loc = -1500;
fig_bottom_loc = 400;
fig_width = 1000;
fig_height = 450;
set(fig,'position',[fig_left_loc,fig_bottom_loc,fig_width,fig_height])

surf(X1, X2, K);
view(2);