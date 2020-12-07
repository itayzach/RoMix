clc;
clear;
close all;
rng('default');
set(0,'DefaultFigureWindowStyle','docked')

%% Setup
% 'Gaussian_1D' / 'Uniform_1D'
verticesPDF = 'Uniform_1D';  

% 'randomMatrix' / 'permutation' / 'eye'
verticesTransform = 'randomMatrix';

% 'WithR' / 'WithoutR'
funcTransform = 'WithR';

%% Random data
N = 1000;
if strcmp(verticesPDF, 'Uniform_1D')
    maxVal = -1;
    minVal = 1;
    X = (maxVal - minVal)*rand(N,1) + minVal;
elseif  strcmp(verticesPDF, 'Gaussian_1D')
    sigma = 1;
    mu = 0;
    X = sigma*randn(N,1) + mu;
else
    error('invalid verticesPDF');
end
    
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

%% Transfrom
if strcmp(verticesTransform, 'randomMatrix')
    R = (1/sqrt(N))*randn(N,N);
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

figure('Name', 'Histogram of X_tilde');
histfit(X_tilde,100);
title('Histogram of $\tilde{X}$', 'interpreter', 'latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

figure('Name', 'R');
imagesc(R);
colorbar();
title('$R$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

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

%% Transform eigenvectors (Functional maps: C = L2.evecs'*L2.A'*P'*L1.evecs)
if strcmp(funcTransform, 'WithR')
    C = pinv(Phi_tilde)*R*V;
    T = (Phi_tilde*pinv(Phi_tilde))*R*(V*V');
elseif strcmp(funcTransform, 'WithoutR')
    C = pinv(Phi_tilde)*V;
    T = (Phi_tilde*pinv(Phi_tilde))*(V*V');
else
    error('invalid funcTransform');
end


figure('Name', 'T');
imagesc(T);
colorbar();
title('$T$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
figure('Name', 'C');
imagesc(C);
colorbar();
title('$C$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% f_tilde
f_tilde = Phi_tilde*C(:,vInd);
Tv_tilde = T*V(:,vInd);

figure('Name', 'f_tilde');
plot(X_tilde, f_tilde,'o');
hold on
plot(X_tilde, Tv_tilde,'.');
title('$\tilde{f}$ on $G$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);
legend([strcat('$\tilde{f}_',string(vInd),'$') strcat('$T v_',string(vInd),'$') ]', 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))

%% Transform back
X_rec = R\X_tilde;
% V_rec = R\V_tilde*C;
V_rec = pinv(T)*Tv_tilde;

figure('Name', 'Reconstructed eigenvectors of W_G');
plot(X_rec, V_rec,'.');
legend(strcat('$v_{{\bf rec},',string(vInd),'}$'), 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
title('Reconstructed eigenvectors of $W_G$', 'interpreter', 'latex', 'FontSize', 16); set(gca,'FontSize', 14);

set(0,'DefaultFigureWindowStyle','normal')