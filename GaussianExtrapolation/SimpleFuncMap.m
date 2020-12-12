clc;
clear;
close all;
rng('default');
set(0,'DefaultFigureWindowStyle','normal')
%% Setup
% 'Gaussian_1D' / 'Uniform_1D'
verticesPDF = 'Gaussian_1D';  

% 'randomMatrix' / 'permutation' / 'eye'
verticesTransform = 'permutation';

%% Random data
N = 500;
if strcmp(verticesPDF, 'Uniform_1D')
    maxVal = 1;
    minVal = -1;
    X = (maxVal - minVal)*rand(N,1) + minVal; %linspace(minVal,maxVal,N)';
elseif  strcmp(verticesPDF, 'Gaussian_1D')
    sigma = 1;
    mu = 0;
    X = sigma*randn(N,1) + mu;
else
    error('invalid verticesPDF');
end

set(0,'DefaultFigureWindowStyle','docked')
figure('Name', 'Histogram of X'); 
histogram(X,100);
title('Histogram of $X$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Generate graph
M = 20;
omega = 0.3;
dist = pdist2(X, X);
W = exp(-dist.^2/(2*omega^2));

%% Calculate (numeric) eigenvectors of G
[V, Lambda] = eigs(W,M);
lambda = diag(Lambda);
vInd = 1:3;

figure('Name', 'Eigenvectors of W_G');
plot(X, V(:,vInd),'.');
title('Eigenvectors of $W_G$', 'interpreter', 'latex', 'FontSize', 16);
legend(strcat('$v_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

%% Transform
if strcmp(verticesTransform, 'randomMatrix')
    R = (1/sqrt(N))*randn(N,N);
%     R = diag(randn(N,1));
elseif strcmp(verticesTransform, 'permutation')
    R = eye(N);
    r = randperm(N);
    R = R(r,:);
elseif strcmp(verticesTransform, 'eye')
    R = eye(N);
else
    error('invalid verticesTransform');
end
    
X_tilde = R*X;
sigma_tilde = 1; %std(X_tilde);
mu_tilde = 0; %mean(X_tilde);

figure('Name', 'Histogram of X_tilde');
histfit(X_tilde,100);
title('Histogram of $\tilde{X}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

figure('Name', 'R');
imagesc(R);
colorbar();
title('${\bf R}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Build G tilde
M_tilde = 20;
omega_tilde = 0.3;
dist_tilde = pdist2(X_tilde, X_tilde);
W_tilde = exp(-dist_tilde.^2/(2*omega_tilde^2));

%% Calculate (analytic) eigenfunctions of W_tilde
[Phi_tilde, lambdaAnalyticTilde] = SimpleCalcAnalyticEigenfunctions(X_tilde, omega_tilde, sigma_tilde, mu_tilde, M_tilde);

figure('Name', '(Analytic) Eigenfunctions of W_Gtilde');
plot(X_tilde, Phi_tilde(:,vInd),'.');
legend(strcat('$\tilde{\phi}_',string(vInd),'$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['(Analytic) Eigenfunctions of $W_{\tilde{G}}$' newline ...
    '$\tilde{\omega}$ = ' num2str(omega_tilde, '%.2f') ...
    '; $\tilde{\sigma}$ = ' num2str(sigma_tilde, '%.2f') ...
    '; $\tilde{\mu}$ = ' num2str(mu_tilde, '%.2f') ], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Calculate (numeric) eigenvectors of W_tilde (for comparison with the analytic eigefunctions)
[V_tilde, LambdaNumericTilde] = eigs(W_tilde, M_tilde);
V_tilde = FlipSign(Phi_tilde, V_tilde);
lambdaNumericTilde = diag(LambdaNumericTilde);

figure('Name', '(Numeric) Eigenvectors of W_Gtilde');
plot(X_tilde, V_tilde(:,vInd),'o');
hold on
plot(X_tilde, Phi_tilde(:,vInd),'.');
legend([strcat('$\tilde{v}_',string(vInd),'$') strcat('$\tilde{\phi}_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('Eigenfunctions and eigenvectors of $W_{\tilde{G}}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Debug
b_useAnalyticEigenfuncs = true;
if ~b_useAnalyticEigenfuncs
    Phi_tilde = V_tilde;
end

%% Transform eigenvectors (Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs)
C = pinv(Phi_tilde)*R*V;
T = (Phi_tilde*pinv(Phi_tilde))*R*(V*V');
pinvC = pinv(C);

figure('Name', 'T');
imagesc(T);
colorbar();
title('${\bf T}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
figure('Name', 'C');
imagesc(C);
colorbar();
title('${\bf C}$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
figure('Name', 'pinv(C)');
imagesc(pinvC);
colorbar();
title('${\bf C}^\dagger$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Phi_tilde in terms of V
Phi_tilde_rec = R*V*pinvC;

% f_tilde = Phi_tilde*C(:,vInd);
% Tv_tilde = T*V(:,vInd);

figure('Name', 'Reconstructed eigenfunctions of W_G_tilde');
plot(X_tilde, Phi_tilde(:,vInd),'o');
hold on
plot(X_tilde, Phi_tilde_rec(:,vInd),'.');
legend([strcat('$\tilde{\phi}_',string(vInd),'$') strcat('$\tilde{\phi}_{{\bf rec},',string(vInd),'}$')] , 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Reconstructed eigenvectors on $\tilde{G}$' newline '${\bf \tilde{\Phi}}_{{\bf rec}} = {\bf R}{\bf V} {\bf C}^\dagger$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% V in terms of Phi_tilde
X_rec = R\X_tilde;
V_rec = R\Phi_tilde*C;

figure('Name', 'Reconstructed eigenvectors of W_G');
plot(X, V(:,vInd),'o');
hold on
plot(X_rec, V_rec(:,vInd),'.');
legend([strcat('$v_',string(vInd),'$') strcat('$v_{{\bf rec},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Reconstructed eigenvectors on $G$' newline '${\bf V}_{{\bf rec}} = {\bf R}^{-1}{\bf \tilde{\Phi}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% V in terms of Phi_tilde (take 2)
X_rec = R\X_tilde;
V_rec = pinv(T)*T*V;

figure('Name', 'Reconstructed eigenvectors of W_G');
plot(X, V(:,vInd),'o');
hold on
plot(X_rec, V_rec(:,vInd),'.');
legend([strcat('$v_',string(vInd),'$') strcat('$v_{{\bf rec},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Reconstructed eigenvectors of $W_G$' newline '${\bf V}_{{\bf rec}} = {\bf T}^\dagger {\bf T} {\bf V}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Interpolate V in terms of Phi_tilde
% 'interp' / 'NewRandomPoints'
interpMethod = 'NewRandomPoints'; 
N_int = 2000;

interpRatio = (N+N_int)/N;
if strcmp(interpMethod, 'NewRandomPoints')
    if strcmp(verticesPDF, 'Uniform_1D')
        maxVal = 1;
        minVal = -1;
        X_tilde_int = [X_tilde; (maxVal - minVal)*rand(N_int,1) + minVal];
    elseif  strcmp(verticesPDF, 'Gaussian_1D')
        X_tilde_int = [X_tilde; sigma_tilde*randn(N_int,1)+mu_tilde];
    else
        error('invalid verticesPDF');
    end
elseif strcmp(interpMethod, 'interp')
    X_tilde_int = interp(X_tilde,4);
else
    error('invalid interpMethod')
end

[Phi_tilde_int, ~] = SimpleCalcAnalyticEigenfunctions(X_tilde_int, omega_tilde, sigma_tilde, mu_tilde, M_tilde);

R_int = [R zeros(N,N_int);
         zeros(N_int,N) eye(N_int)];
X_int = R_int\X_tilde_int;
V_int = (1/sqrt(interpRatio))*R_int\Phi_tilde_int*C;

figure('Name', 'Interpolated eigenvectors of W_G');
plot(X, V(:,vInd),'o');
hold on
plot(X_int, V_int(:,vInd),'.');
legend([strcat('$v_',string(vInd),'$') strcat('$v_{{\bf int},',string(vInd),'}$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title(['Interpolated eigenvectors on $G$' newline '${\bf V}_{{\bf int}} = {\bf R_{{\bf int}}}^{-1} {\bf \tilde{\Phi}}_{{\bf int}} {\bf C}$'], 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
