%% Restart run
clear; close all; clc;
rng(0); % set seed

if ~exist('figures', 'dir')
    mkdir('figures')
end

%% Parameters
% kernel
a = 1;
b = 3;

% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);

% num of eigenfunctions
M = 7;

% simulation
b_verifyEigenfunctions = false;
b_plotEigenFigs = false;

%% Check eigenfunctions
dx = 0.01;
x = -10:dx:10-dx;

if b_verifyEigenfunctions || b_plotEigenFigs
    for m = 0:M-1
        [vPhi_m_x, ~] = SqExpEig(a, b, m, x);
        vP_x = p(x,sigma);

        if b_plotEigenFigs
            fig1 = figure(1);
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
                vKernel_y_x = kernel(y(j),x,l);
                lhs(j) = sum(vKernel_y_x.*vPhi_m_x.*vP_x*dx);
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
end
if b_plotEigenFigs
    hold off
%     title('Eigenfunctions')
    legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
    print(fig1, ['figures' filesep 'fig1_eigenfunctions'], '-dpdf')
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
dx = 0.01;
x = (-1:dx:1-dx)';
N = length(x);

dt = 0.01;
t = x; %(-20:dt:20-dt)';
vP_t = p(t,sigma);

mF      = [10*sin(5*x) + 10*exp(-x).*sin(2.5*x)   10*exp(-2*x).*sin(5*x) ];
mF_awgn = [1*randn(N,1) 1*randn(N,1)];
cFstr = {'10\sin(5x) + 10e^{-x}\sin(2.5x)' '10e^{-2x}\sin(5x)'};
nFuncs = size(mF, 2);

% x_eval = ([-5:0.1:5])';
% vF_eval = [sin(5*x_eval) exp(-x_eval).*sin(2.5*x_eval) exp(-2*x_eval).*sin(5*x_eval)];

cFigs = cell(1, nFuncs);

for i = 1:nFuncs
    vFi = mF(:, i);
    vFi_awgn = mF_awgn(:, i);
    vGi = vFi + vFi_awgn;
    mPhi = zeros(N, M);
    for m = 0:M-1 
        [vPhi_m_x, lambda_m] = SqExpEig(a, b, m, x);
        mPhi(:, m+1) = vPhi_m_x;
    end
    r = 0.1*N;
    vR = randperm(N,N/r); % 1:r:N;
    vCR = pinv(mPhi(vR, :)) * vGi(vR);
    vFi_hat = mPhi * vCR;
    
%     % ---------------------------------------------------------------------
%     % Kernel try
%     vFi_eval = vF_eval(:, i);
%     vPhi_m_x = zeros(M, 1);
%     vPhi_m_x2 = zeros(M, 1);
%     vFi_hat2 = zeros(length(x_eval), 1);
%     for j = 1:length(x_eval)
%         for m = 0:M-1
%             vKernel_x_t = kernel(x_eval(j),t,l);
%             [vPhi_m_t, lambda_m] = SqExpEig(a, b, m, t);
%             I = conv(vKernel_x_t, vPhi_m_t.*vP_t*dt, 'same');
%             [ ~, idx ] = min( abs( t-x_eval(j) ) );
%             vPhi_m_x2(m+1) = (1/lambda_m)*I(idx);
%             vPhi_m_x(m+1) = (1/lambda_m)*sum(vKernel_x_t.*vPhi_m_t.*vP_t*dt);
%         end
%         vFi_hat2(j) = sum(vCR.*vPhi_m_x2);
%         res = sin(5*x_eval(j));
%     end
%     % ---------------------------------------------------------------------
    cFigs{i} = figure(i+1);
    p1 = plot(x, vFi_hat, 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = \Phi c $ with ' num2str(N/r) ' points']);
    hold on
%     title(['Extrapolate $f_' num2str(i) ' = ' cFstr{i} '$'], 'Interpreter', 'latex', 'FontSize', 12)
%     p2 = plot(x_eval, vFi_hat2, '*');
    p3 = plot(x, vGi, '-.', 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = ' cFstr{i} '$']);
    p4 = plot(x, vFi, 'b--', 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = ' cFstr{i} ' + $AWGN']);
    p5 = plot(x(vR), vGi(vR), 'o', 'LineWidth', 2);
    hold off
    legend([p1 p3 p4], 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
    print(cFigs{i}, ['figures' filesep 'fig' num2str(i+1) '_extrapolate_f' num2str(i) '.pdf'], '-dpdf')
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
% vHm = hermiteH(m, sqrt(2*c)*x);
vHm = hermite(m, sqrt(2*c)*x);
vPhi_m = exp( -(c-a)*x.^2 ) .* vHm;
end

function vPr = p(y, sigma)
vPr = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );
end

function vKernel = kernel(x,y,l) 
vKernel = exp(-(x-y).^2./(2*l^2));
end