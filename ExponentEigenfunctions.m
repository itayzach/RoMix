%% Restart run
clear; close all; clc;

%% Parameters
% kernel
a = 1;
b = 3;
mu = 0;

% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);

% num of eigenfunctions
L = 10;

% simulation
b_verifyEigenfunctions = true;
b_plotEigenFigs = true;

%% Check eigenfunctions
dx = 0.01;
x = -10:dx:10-dx;

if b_plotEigenFigs
    if ~exist('fig', 'var')
        fig = figure;
        tg = uitabgroup; % tabgroup
    end
    thistab = uitab(tg, 'Title', 'Eigenfunctions'); % build tab
    axes('Parent',thistab); % somewhere to plot
end
for k = 0:L-1
    [phi_k_x, ~] = SEeig(a, b, k, x, mu);

    if b_plotEigenFigs
        if k <= 2
            plot(x, phi_k_x, 'DisplayName', [ '$\phi_' num2str(k) '$' ]);
            xlim([-5 5])
            hold on
        end
    end

    if b_verifyEigenfunctions
        y = -2+mu:0.5:2+mu;
        [phi_k_y, lambda_k] = SEeig(a, b, k, y, mu);
        rhs = lambda_k * phi_k_y;
        lhs = zeros(1, length(y));
        for j = 1:length(y)
            lhs(j) = sum(kernel(x,y(j)-mu,l).*phi_k_x.*p(x,sigma)*dx);
        end
        fprintf('k = %d\n', k)
        fprintf('lambda * phi = ')
        fprintf('%.6f  ', rhs);
        fprintf('\n');
        fprintf('<k, phi>     = ')
        fprintf('%.6f  ', lhs);
        fprintf('\n');
    end
end
if b_plotEigenFigs
    hold off
    title('Eigenfunctions')
    legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
end

%% Read mnist data
% mnist = load('data/mnist.mat');
% mXTest = single(mnist.testX);
% vTestLabels = mnist.testY.';
% [nTestPoints, nPixels] = size(mXTest);
% N = nTestPoints;

%% Graph signals
% nDigits = 10;
% mS = zeros(N, nDigits);
% for k = 0:9
%     mS(:, k+1) = (vTestLabels == k);
% end

%% Reconstruct functions
N = 1000;
x = linspace(-1, 1, N)'; % [x_min x_max] should be small since gaussians decay to zero
f = [sin(5*x) exp(-7*x).*sin(10*x) exp(-2*x).*sin(5*x)];
nFuncs = size(f, 2);

if ~exist('fig', 'var')
    fig = figure;
    tg = uitabgroup; % tabgroup
end
for i = 1:nFuncs
    fi = f(:, i);
    mPhi = zeros(N, L);
    for k = 0:L-1 
        [phi_k_x, ~] = SEeig(a, b, k, x, mu);
        mPhi(:, k+1) = phi_k_x;
    end
    r = 0.1*N;
    R = 1:r:N;
    vCR = pinv(mPhi(R, :)) * fi(R);
    thistab = uitab(tg, 'Title', 'Reconstruction');
    axes('Parent',thistab);
    p1 = plot(x, mPhi * vCR, 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = \Phi c $ with ' num2str(N/r) ' points']);
    hold on
    p2 = plot(x, fi, '-.', 'DisplayName', ['$f_' num2str(i) '$']);
    p3 = plot(x(R), fi(R), 'o');
    hold off
    legend([p1 p2], 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
end

%% SEeig (Squared Exponentional)
function [phi_k, lambda_k] = SEeig(a, b, k, x, mu)

% Calculate parameters
c = sqrt(a^2 + 2*a*b);
A = a + b + c;
B = b/A;

% k-th eigenvalue
lambda_k = sqrt(2*a/A) * B^k;

% k-th eigenfunction
Hk = hermiteH(k, sqrt(2*c)*x);
phi_k = exp( -(c-a)*(x - mu).^2 ) .* Hk;
end

function pr = p(y, sigma)
pr = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );
end

function k = kernel(x,y,l) 
k = exp(-(x-y).^2./(2*l^2));
end