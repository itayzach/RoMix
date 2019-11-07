%% Restart run
clear; close all; clc;
rng(); % set seed
%% Parameters
% kernel
a = 1;
b = 3;

% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);

% num of eigenfunctions
M = 10;

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
for m = 0:M-1
    [vPhi_m_x, ~] = SqExpEig(a, b, m, x);
    vP_x = p(x,sigma);
    
    if b_plotEigenFigs
        if m <= 2
            if m == 0
                plot(x, vP_x, 'DisplayName', '$p(x)$');
                xlim([-5 5])
                hold on
            end
            plot(x, vPhi_m_x, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
            xlim([-5 5])
            xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
        end
    end

    if b_verifyEigenfunctions
        y = -2:0.5:2;
        [phi_m_y, lambda_m] = SqExpEig(a, b, m, y);
        rhs = lambda_m * phi_m_y;
        lhs = zeros(1, length(y));
        for j = 1:length(y)
            vKernel_x = kernel(x,y(j),l);
            lhs(j) = sum(vKernel_x.*vPhi_m_x.*vP_x*dx);
        end
        fprintf('m = %d\n', m)
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

%% Extrapolate functions
% N = 1000;
% x = linspace(-1, 1, N)'; % [x_min x_max] should be small since gaussians decay to zero

dx = 0.002;
x = (-1:dx:1-dx)';
N = length(x);

dt = 0.002;
t = -10:dt:10-dt;
vP_t = p(t,sigma);

f = [sin(5*x) exp(-x).*sin(2.5*x) exp(-2*x).*sin(5*x)];
nFuncs = size(f, 2);

if ~exist('fig', 'var')
    fig = figure;
    tg = uitabgroup; % tabgroup
end
for i = 1:nFuncs
    fi = f(:, i);
    mPhi = zeros(N, M);
    for m = 0:M-1 
        [vPhi_m_x, ~] = SqExpEig(a, b, m, x);
        mPhi(:, m+1) = vPhi_m_x;
    end
    r = 0.1*N;
    R = randperm(N,N/r);%    1:r:N;
    vCR = pinv(mPhi(R, :)) * fi(R);
    fi_hat = mPhi * vCR;
    
%     for m = 0:M-1 
%         for j = 1:length(x)
%             vKernel_t = kernel(t,x(j),l);
%             [vPhi_m_t, ~] = SqExpEig(a, b, m, t);
%             phi_m_x(j) = sum(vKernel_t.*vPhi_m_t.*vP_t*dt);
%         end
%     end
    
    thistab = uitab(tg, 'Title', ['Extrapolate f' num2str(i)]);
    axes('Parent',thistab);
    p1 = plot(x, fi_hat, 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = \Phi c $ with ' num2str(N/r) ' points']);
    title(['Extrapolate $f_' num2str(i) '$'], 'Interpreter', 'latex', 'FontSize', 12)
    hold on
    p2 = plot(x, fi, '-.', 'DisplayName', ['$f_' num2str(i) '$']);
    p3 = plot(x(R), fi(R), 'o');
    hold off
    legend([p1 p2], 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
end


%% SqExpEig (Squared Exponentional)
function [vPhi_m, lambda_m] = SqExpEig(a, b, m, x)

% Calculate parameters
c = sqrt(a^2 + 2*a*b);
A = a + b + c;
B = b/A;

% m-th eigenvalue
lambda_m = sqrt(2*a/A) * B^m;

% m-th eigenfunction
vHm = hermiteH(m, sqrt(2*c)*x);
vPhi_m = exp( -(c-a)*x.^2 ) .* vHm;
end

function vPr = p(y, sigma)
vPr = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );
end

function vKernel = kernel(x,y,l) 
vKernel = exp(-(x-y).^2./(2*l^2));
end