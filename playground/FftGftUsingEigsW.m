%% Restart
clc; clear; rng('default');
close all; 
set(0,'DefaultFigureWindowStyle','normal')

%% Time
T = 1 - 1e-3;
sigmaT = 0.1;
fs = 1e3;
Ts = 1/fs;

% t = (0:Ts:T)';
% N = length(t);
N = 2000;
t = sigmaT*randn(N,1);

interpRatio = 1;
n = N/interpRatio;

% randInd = randperm(N);
% rand_t = t(sort(randInd(1:n)));
% rand_t = linspace(0,T,n)';
%% Signal
sigma = 0;
mu = 0;
noise = sigma*randn(n,1) + mu;
f0 = (fs/2)/10;
signal = sin(2*pi*f0*t) + noise;

%% Classic DFT
VDft = (1/sqrt(n))*dftmtx(n);
dftSignal = abs(VDft.'*signal); % <v_k,signal>_{\ell_2} --> p(x)=uniform

figure;
subplot(3,1,1)
plot(t,signal,'.');
title('signal');  set(gca,'FontSize', 14);
xlabel('t [ms]')
set(gca,'FontSize', 14);
subplot(3,1,2);
plot([real(VDft(:,1:3)) imag(VDft(:,1:3))])
title('DFT eigenvectors'); set(gca,'FontSize', 14);
subplot(3,1,3)
plot(dftSignal)
title('|DFT^T s|'); set(gca,'FontSize', 14);
x0 = 100; y0 = 100; width = 600; height = 600;
set(gcf,'Position', [x0 y0 width height])

%% GFT with eigs(W)
M = 20;
% omega = 0.5*Ts;
omega = 0.3*sigmaT;
dist = pdist2(t, t);
W = exp(-dist.^2/(2*omega^2));
% W = SimpleCalcAdjacency(t, 'GaussianKernel', omega);
% W = toeplitz([ 0 1 zeros(1, n-2)]); % path graph
% W = toeplitz([ 1 1 zeros(1, n-2)]); % path graph with self loops
% W = toeplitz([ 0 1 zeros(1, n-3) 1]); % ring graph
% W = toeplitz([ 1 1 zeros(1, n-3) 1]); % ring graph with self loops


[VGsp, LambdaGsp] = eigs(W,M);
gftSignal = (VGsp.'*signal);
lambdaGsp = diag(LambdaGsp);

figure;
subplot(3,1,1)
plot(t,signal,'.');
title('signal');  set(gca,'FontSize', 14);
xlabel('t [ms]')
set(gca,'FontSize', 14);
subplot(3,1,2);
plot(t, VGsp(:,1:6), '.')
title('W eigenvectors'); set(gca,'FontSize', 14);
subplot(3,1,3)
plot(lambdaGsp, gftSignal)
title('|GFT^T s|'); set(gca,'FontSize', 14);
x0 = 700; y0 = 100; width = 600; height = 600;
set(gcf,'Position', [x0 y0 width height])

