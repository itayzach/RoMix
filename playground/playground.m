%% Restart
clc; clear; close all;
rng('default');
qID = 2;
Nx = 101;
M = toeplitz([ 0 1 zeros(1, Nx-3) 1]);

x = 1:Nx;
y = zeros(1, length(x));
outputFolder = 'figs';
x0 = 100;
y0 = 100;
width = 700;
height = 400;

%% (a) Construct M
Nx = size(M, 1);

%% (b) Construct D and L
D = diag(sum(M,1));
L = D - M;

%% (c) Compute eigen-decomposition
[mPsi_M, vLambda_M] = eig(M);
vLambda_M = diag(vLambda_M);

[mPsi_L, vLambda_L] = eig(L);
vLambda_L = diag(vLambda_L);

%% ii. eigenvalues
fig = figure;
stem(vLambda_M, 'filled', 'LineStyle', 'none', 'MarkerSize', 1);
xlim([1 Nx]);
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\lambda_k$', 'Interpreter', 'latex', 'FontSize', 16);
title('Adjacency')
set(gca,'FontSize', 14);
set(gcf, 'Position', [x0, y0, width, height])


fig = figure;
stem(vLambda_L, 'filled', 'LineStyle', 'none', 'MarkerSize', 1);
xlim([1 Nx]);
title('Laplacian')
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$\lambda_k$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
set(gcf, 'Position', [x0, y0+height, width, height])
%% iii. eigenvectors
fig = figure;
for i = 1:5
    plot(mPsi_M(:,Nx+1-i), 'LineWidth', 2, 'DisplayName', ['$\psi_{' num2str(Nx+1-i) '}$']);
    hold on;
end
xticks((1:10:Nx)-1);
xlim([0 Nx+1])
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
title('Adjacency')
set(gca,'FontSize', 14);
set(gcf, 'Position', [x0+width, y0, width, height])

fig = figure;
for i = 1:5
    plot(mPsi_L(:,i), 'LineWidth', 2, 'DisplayName', ['$\psi_' num2str(i) '$']);
    hold on;
end
xticks((1:10:Nx)-1);
xlim([0 Nx+1])
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northwest');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 16);
title('Laplacian')
set(gca,'FontSize', 14);
set(gcf, 'Position', [x0+width, y0+height, width, height])