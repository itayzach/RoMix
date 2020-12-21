%% Init
clc; clear; close all;
rng('default');
set(0,'DefaultFigureWindowStyle','docked')

%% Data
dataset = 'torus'; % 'torus' / 'trefoil'
N = 1000;
b_normalize = true;
b_plotDataset = true;
X = GetDiffMapsDataset(dataset, N, b_normalize, b_plotDataset);
N = length(X);
%% Kernel
epsilon = 0.05;
dist = pdist2(X, X);
K = exp(-dist.^2/epsilon);

figure('Name', 'K');
imagesc(K);
colorbar;
title(strcat('$K_{i,j} = \exp(-\|x_i-x_j\|^2/\varepsilon), \quad \varepsilon = ', num2str(epsilon),'$'),...
    'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);
%% Normalize
q = sum(K)';
Q = diag(q);

alpha = 1;
Ka = Q^-alpha*K*Q^-alpha;
da = sum(Ka)';
Da = diag(da);

if alpha == 0
    alphaDesc = 'Normalized graph Laplacian';
elseif alpha == 0.5
    alphaDesc = 'Fokker-Planck diffusion';
elseif alpha == 1
    alphaDesc = 'Heat kernel';
end
alphaStr = strcat('($\alpha = ', num2str(alpha),'$, ',alphaDesc,')');

figure('Name', 'Ka');
imagesc(Ka);
colorbar;
title(strcat('$K^{(\alpha)} = Q^{-\alpha}KQ^{-\alpha}, \quad \varepsilon = ', num2str(epsilon),'$'),...
    'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);
%% RW matrix
Pa = Da^-1*Ka;
isalmostequal(diag(sum(Pa,2)),eye(N),1e-14) % make sure W is row stochastic

figure('Name', 'Pa');
imagesc(Pa);
colorbar;
title(strcat('$P_\alpha = D^{-1}_\alpha K_\alpha,\quad$', alphaStr),...
    'interpreter','latex', 'FontSize', 16); 
set(gca,'FontSize', 14);
% x0 = 100; y0 = 100; width = 550; height = 400;
% set(gcf,'Position', [x0 y0 width height])

%% Decompose
M = 6; % number of eigenvectors

[Psi, Lambda, Phi] = eigs(Pa,M);

% Psi2 = Da^(-1/2)*Psi;
Psi = Psi./Psi(:,1);
%% Psi^T D Psi
% figure;
% imagesc(Psi'*D*Psi);
% colorbar;
% title('$\Psi^T D \Psi$','interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Plot eigenvectors of K
% [Phi_K, ~] = eigs(K, M);
% figure;
% for i = 1:M
%     subplot(2,3,i)
%     if size(X,2) == 3
%         scatter3(X(:,1), X(:,2), X(:,3), 50, Phi_K(:,i), 'filled');
%     else
%         scatter(X(:,1), X(:,2), 50, Phi_K(:,i), 'filled');
%     end
%     colorbar();
%     title(strcat('$\phi_', num2str(i), '$'),'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);
% end
% sgtitle('Eigenvectors of K');
% x0 = 400; y0 = 200; width = 1200; height = 600;
% set(gcf,'Position', [x0 y0 width height])

%% Plot eigenvectors of P
% figure;
% for i = 1:M
%     psiMax = max(Psi(:));
%     psiMin = min(Psi(:));
%     subplot(2,3,i)
%     if size(X,2) == 3
%         scatter3(X(:,1), X(:,2), X(:,3), 50, Psi(:,i), 'filled');
%     else
%         scatter(X(:,1), X(:,2), 50, Psi(:,i), 'filled');
%     end
%     caxis([ psiMin psiMax]);
%     colorbar();
%     title(strcat('$\psi_', num2str(i), '$, $\lambda_', num2str(i), ' = ', num2str(Lambda(i,i)), '$'),...
%         'interpreter','latex', 'FontSize', 16); 
%     set(gca,'FontSize', 14);
% end
% sgtitle(strcat('Eigenvectors of $P_\alpha,\quad$', alphaStr),...
%     'interpreter','latex', 'FontSize', 16);
% x0 = 400; y0 = 200; width = 1200; height = 600;
% set(gcf,'Position', [x0 y0 width height])

%% Diffusion maps
t = 1;
diffmap = diag(Lambda)'.^t.*Psi;

%% Plot
figure('Name', 'Diffusion map');
cmap = 1:N;
scatter(diffmap(:,2), diffmap(:,3), 50, cmap, 'filled')
xlabel('$\lambda_2 \psi^{(\alpha)}_2$','interpreter','latex')
ylabel('$\lambda_3 \psi^{(\alpha)}_3$','interpreter','latex')
title(strcat('$t = ', num2str(t),',\quad$', alphaStr),'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Dataset
function X = GetDiffMapsDataset(dataset, N, b_normalize, b_plotDataset)
    if strcmp(dataset, 'torus')
        R = 10; r = 4;
        tmin = 0;
        tmax = 1;
        t = linspace(tmin,tmax,sqrt(N))';
        [x, y] = meshgrid(t,t);
        x = x(:);
        y = y(:);
        X = [(R + r*cos(2*pi*y)).*cos(2*pi*x), (R + r*cos(2*pi*y)).*sin(2*pi*x), r*sin(2*pi*y) ];    
    elseif strcmp(dataset, 'trefoil')
        tmin = 0;
        tmax = 2*pi;
        t = linspace(tmin,tmax,N+1)'; t = t(1:end-1);
        X = [sin(t)+2*sin(2*t), cos(t)-2*cos(2*t), -sin(3*t) ];
    else
        error('invalid dataset')
    end

    if b_normalize
        X = X - mean(X);
        X = X./std(X);
    end

    if b_plotDataset
        figure('Name', 'Dataset');
        cmap = 1:length(X);
        scatter3(X(:,1), X(:,2), X(:,3), 50, cmap, 'filled');
        title(strcat('$X$', ' (', dataset, ')'), 'interpreter', 'latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
end