%% Init
clc; clear; close all;
rng('default');
set(0,'DefaultFigureWindowStyle','docked')

%% Data
dataset = 'trefoil'; % 'torus' / 'trefoil'
N = 1000;
b_normalize = true;
b_plotDataset = true;
X = GetDiffMapsDataset(dataset, N, b_normalize, b_plotDataset);
N = length(X);
%% Kernel
epsilon = 0.01;
dist = pdist2(X, X);
K = exp(-dist.^2/epsilon);

figure('Name', 'K');
imagesc(K);
colorbar;
title(strcat('$K_{i,j} = \exp(-\|x_i-x_j\|^2/\varepsilon), \quad \varepsilon = ', num2str(epsilon),'$'),...
    'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);
%% (2) Density (q_epsilon)
q_eps = sum(K)';
Q_eps = diag(q_eps);

figure('Name', 'q');
cmap = 1:N;
scatter(1:N, q_eps, 50, cmap, 'filled');
title(strcat('$q$'),...
    'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% (3) Normalize by Q and alpha
alpha = 1;
Ka = Q_eps^-alpha*K*Q_eps^-alpha;
da = sum(Ka)';
Da = diag(da);

if alpha == 0
    alphaDesc = 'Normalized Laplacian. max influence of $q$';
elseif alpha == 0.5
    alphaDesc = 'Fokker-Planck diffusion';
elseif alpha == 1
    alphaDesc = 'Heat kernel. min influence of $q$';
end
alphaStr = strcat('($\alpha = ', num2str(alpha),':\quad$ ',alphaDesc,')');

figure('Name', 'Da');
cmap = 1:N;
scatter(1:N, da, 50, cmap, 'filled');
title(['$D_\alpha, \quad$', newline ,alphaStr],...
    'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);

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
title(['$P_\alpha = D^{-1}_\alpha K_\alpha,\quad$', newline, alphaStr],...
    'interpreter','latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

%% Decompose
M = 6; % number of eigenvectors
[Psi, Lambda, Phi] = eigs(Pa,M);
assert(~any(isnan(diag(Lambda))), 'You have at least one NaN');
Psi = Psi./Psi(:,1);
%% Psi^T D Psi
PsiT_Da_Psi = Psi'*Da*Psi;
figure('Name', 'Psi^T D Psi');
imagesc(PsiT_Da_Psi);
colorbar;
title('$\Psi^T D_\alpha \Psi$','interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);

%% Plot eigenvectors of K
if b_plotEigenvectors
    [Phi_K, ~] = eigs(Ka, M);
    set(0,'DefaultFigureWindowStyle','normal')
    figure('Name', 'K''s eigenvecs');
    for i = 1:M
        subplot(2,3,i)
        if size(X,2) == 3
            scatter3(X(:,1), X(:,2), X(:,3), 50, Phi_K(:,i), 'filled');
        else
            scatter(X(:,1), X(:,2), 50, Phi_K(:,i), 'filled');
        end
        colorbar();
        title(strcat('$\phi_', num2str(i), '$'),'interpreter','latex', 'FontSize', 16); set(gca,'FontSize', 14);
    end
    sgtitle('Eigenvectors of K');
    x0 = 400; y0 = 200; width = 1200; height = 600;
    set(gcf,'Position', [x0 y0 width height])
    set(0,'DefaultFigureWindowStyle','docked')
end
%% Plot eigenvectors of P
if b_plotEigenvectors
    set(0,'DefaultFigureWindowStyle','normal')
    figure('Name', 'Pa''s eigenvecs');
    for i = 1:M
        psiMax = max(Psi(:));
        psiMin = min(Psi(:));
        subplot(2,3,i)
        if size(X,2) == 3
            scatter3(X(:,1), X(:,2), X(:,3), 50, Psi(:,i), 'filled');
        else
            scatter(X(:,1), X(:,2), 50, Psi(:,i), 'filled');
        end
        caxis([ psiMin psiMax]);
        colorbar();
        title(strcat('$\psi^{(\alpha)}_', num2str(i), '$, $\lambda_', num2str(i), ' = ', num2str(Lambda(i,i)), '$'),...
            'interpreter','latex', 'FontSize', 16); 
        set(gca,'FontSize', 14);
    end
    sgtitle(['Eigenvectors of $P_\alpha,\quad$', newline, alphaStr],...
        'interpreter','latex', 'FontSize', 16);
    x0 = 400; y0 = 200; width = 1200; height = 600;
    set(gcf,'Position', [x0 y0 width height])
    set(0,'DefaultFigureWindowStyle','docked')
end

%% Diffusion maps
t = 1;
diffmap = diag(Lambda)'.^t.*Psi;

figure('Name', 'Diffusion map');
cmap = 1:N;
scatter(diffmap(:,2), diffmap(:,3), 50, cmap, 'filled')
xlabel('$\lambda^t_2 \psi^{(\alpha)}_2$','interpreter','latex')
ylabel('$\lambda^t_3 \psi^{(\alpha)}_3$','interpreter','latex')
title(['$t = ', num2str(t),',\quad$', newline, alphaStr],...
    'interpreter','latex', 'FontSize', 16); 
set(gca,'FontSize', 14);

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
%         t = sort(2*pi*rand(N,1));
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