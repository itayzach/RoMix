%% Init
clc; clear; close all;
rng('default');
set(0,'DefaultFigureWindowStyle','docked')

%% Data
N = 5000;
sigma = 0.5;
mu = 30;
X_I = mu + sigma*randn(N, 1); %sigma*linspace(-1,1,N)';

n = 1000;
X_R = mu + sigma*randn(n, 1); %sigma*linspace(-1,1,n)';

%% Kernel
epsilon = 0.5*sigma^2;
dist = pdist2(X_I, X_R);
A = exp(-dist.^2/epsilon); % "cross"-kernel (Nxn)

W_I = A*A'; % NxN
W_R = A'*A; % nxn
%% Decompose
M = 4; vInd = 1:M;

[Psi_WI,Lambda_WI] = eigs(W_I,M); % Actual eigenvectors of W_I (NxN), entire image
lambda_WI = diag(Lambda_WI);

%% Extension
[Phi,Lambda_WR] = eigs(W_R,M); % Eigenvectors of W_R (nxn), reference set
lambda_WR = diag(Lambda_WR);

PsiExt = A*Phi*diag(1./sqrt(lambda_WR)); % Haddad's extension
PsiExt = FlipSign(Psi_WI,PsiExt);

figure('Name', 'Actual W_I eig vs. Haddad');
plot(X_I, PsiExt,'o');
hold on;
plot(X_I, Psi_WI,'.');
title(['Actual eigenvectors of $W_I = A A^T \in {\bf R}^{N \times N}$, $\psi_j \in {\bf R}^N$' newline 'vs. Haddad''s extension: $\hat{\psi}_j = \frac{1}{\sqrt{\lambda_j}}{\bf A} \phi_j$'], ...
       'interpreter', 'latex', 'FontSize', 16);
legend([strcat('$\hat{\psi}_',string(vInd),'$') strcat('$\psi_',string(vInd),'$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

%% svd of A
[Psi_A,Sigma,Phi_A] = svds(A,M); % svd of A, similarities of image to ref. set
Psi_A = FlipSign(Psi_WI,Psi_A);

figure('Name', 'Actual W_I eigvecs vs. A left sigvecs');
plot(X_I, Psi_A,'o');
hold on;
plot(X_I, Psi_WI,'.');
title(['Actual eigenvectors of $W_I = A A^T \in {\bf R}^{N \times N}$, $\psi_j \in {\bf R}^N$' newline 'vs. $A \in {\bf R}^{N \times n}$  left singular-vectors, $\psi^{(A)}_j \in {\bf R}^N$'], ...
       'interpreter', 'latex', 'FontSize', 16);
legend([ strcat('$\psi^{(A)}_',string(vInd),'$') strcat('$\psi_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

%% Analytic eigenfunctions
omega = sqrt(epsilon/2);
[PsiAnalytic, ~] = SimpleCalcAnalyticEigenfunctions(X_I, omega, sigma, mu, M);
PsiFlipped = FlipSign(PsiAnalytic, Psi_WI);

figure('Name', 'Actual W_I eig vs. analytic eigenfunctions');
plot(X_I, PsiAnalytic,'o');
hold on;
plot(X_I, PsiFlipped,'.');
title(['Actual eigenvectors of $W_I = A A^T \in {\bf R}^{N \times N}$, $\psi_j \in {\bf R}^N$' newline 'vs. analytic eigenvectors, $\psi^{({\bf Ana})}_j$'], ...
       'interpreter', 'latex', 'FontSize', 16);
legend([ strcat('$\psi^{({\bf Ana})}_',string(vInd),'$') strcat('$\psi_',string(vInd),'$') ], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

%% Normalized RW matrix
d1 = sum(W_R,2);
D1 = diag(d1);
W1 = D1^(-1)*W_R*D1^(-1); % normalized kernel
d2 = sum(W1,2);
D2 = diag(d2);
W2 = D2^(-1)*W1;
[Phi_W2,Lambda_W2] = eigs(W2,M);
Phi_W2 = Phi_W2./Phi_W2(:,1);
lambda_W2 = diag(Lambda_W2);

Atilde = A*D1^(-1)*D2^(-1/2);
d = sum(Atilde*D2^(1/2),2); % same as sum(A*D1^(-1),2);
D = diag(d);
A1 = D^(-1)*Atilde*D2^(1/2);

PsiExtNormed = A1*Phi_W2*diag(1./sqrt(lambda_W2)); % Haddad's extension

figure('Name', 'PsiExtNormed');
plot(X_I, PsiExtNormed,'o');
hold on;
plot(X_R,Phi_W2,'.')
title(['PsiExtNormed'], ...
       'interpreter', 'latex', 'FontSize', 16);
legend([strcat('$\psi_',string(vInd),'$') strcat('$\phi^{W_2}_',string(vInd),'$')], 'interpreter', 'latex', 'Location', 'SouthOutside', 'FontSize', 14,'NumColumns',length(vInd))
set(gca,'FontSize', 14);

% figure;
% imagesc(PsiExtNormed'*D.^2*PsiExtNormed); colorbar;
