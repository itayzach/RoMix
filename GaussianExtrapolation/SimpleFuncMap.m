clc;
clear;
close all;
rng('default');
set(0,'DefaultFigureWindowStyle','docked')

%% Uniform data
N = 1000;
maxVal = -1;
minVal = 1;
X = (maxVal - minVal)*rand(N,1) + minVal;
figure('Name', 'Histogram of X'); 
histogram(X,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Generate graph
M = 5;
omega = 0.3;
dist = pdist2(X, X);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:4;

figure('Name', 'Eigenvectors of W_G');
plot(X, V(:,vInd),'.');
title('Eigenvectors of $W_G$', 'interpreter', 'latex', 'FontSize', 16);
legend(strcat('$v_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

%% Transfrom to Gaussian
R = (1/sqrt(N))*randn(N,N);
X_tilde = R*X;

figure('Name', 'Histogram of X_tilde');
histogram(X_tilde,100);
title('Histogram of $\tilde{X}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% Build G tilde
M_tilde = 30;
omega = 0.3;
dist_tilde = pdist2(X_tilde, X_tilde);
W_tilde = exp(-dist_tilde.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of W_tilde
[V_tilde, LambdaNumericTilde] = eigs(W_tilde, M_tilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

figure('Name', '(Numeric) Eigenvectors of W_Gtilde');
plot(X_tilde, V_tilde(:,vInd),'.');
legend(strcat('$\tilde{v}_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('(Numeric) Eigenvectors of $W_{\tilde{G}}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Calculate (analytic) eigenfunctions of W_tilde
[Phi_tilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(X_tilde, omega, M_tilde);

figure('Name', '(Anayltic) Eigenfunctions of W_Gtilde');
plot(X_tilde, Phi_tilde(:,vInd),'.');
legend(strcat('$\tilde{\phi}_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('(Anayltic) Eigenfunctions of $W_{\tilde{G}}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Transform eigenvectors
C = V_tilde'*V;
% alphaTilde = C*alpha;

% f_tilde = V_tilde*alphaTilde;
f_tilde = V_tilde*C(:,vInd);

T = (V_tilde*V_tilde')*(V*V');

figure('Name', 'Transformation matrices');
subplot(1,2,1)
imagesc(T);
colorbar();
title('$T$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
subplot(1,2,2)
imagesc(C);
colorbar();
title('$C$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

v_tilde = V_tilde(:,vInd);
Tv_tilde = T*V(:,vInd);

figure('Name', 'f_tilde');
plot(X_tilde, f_tilde,'o');
hold on
plot(X_tilde, Tv_tilde,'.');
title('$\tilde{f}$ on $G$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
legend([strcat('$\tilde{f}_',string(vInd),'$') strcat('$T v_',string(vInd),'$') ]', 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))

% plot(X_tilde, v_tilde,'.');
%% Transform back
X_rec = R\X_tilde;
% V_rec = R\V_tilde*C;
V_rec = pinv(T)*Tv_tilde;

figure('Name', 'Reconstructed eigenvectors of W_G');
plot(X_rec, V_rec,'.');
legend(strcat('$v_{{\bf rec},',string(vInd),'}$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('Reconstructed eigenvectors of $W_G$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

set(0,'DefaultFigureWindowStyle','normal')