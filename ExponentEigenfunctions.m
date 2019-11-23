%% Restart run
clear; close all; clc;
rng(0); % set seed

outputFolder = 'figs';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
else
    % Remove all epsc files
    filePattern = fullfile(outputFolder, '*.eps');
    theFiles = dir(filePattern);
    for k = 1 : length(theFiles)
      baseFileName = theFiles(k).name;
      fullFileName = fullfile(outputFolder, baseFileName);
      delete(fullFileName);
    end
end

%% Parameters
% kernel
a = 3;
b = 1;

% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);

% num of eigenfunctions
M = 20;

% simulation
b_verifyEigOrth = true;
b_verifyMercersTheorm = true;
b_verifyEigenfunctions = true;
b_plotEigenFigs = true;

%% Check eigenfunctions
if b_verifyEigenfunctions || b_plotEigenFigs
    dx = 0.01;
    x = -1e3:dx:1e3-dx;
    for m = 0:M-1
        [vPhi_m_x, ~] = SqExpEig(a, b, m, x);
        vP_x = p(x,sigma);
        
        if b_plotEigenFigs
            nEigenFuncsToPlot = 4;
            fig1 = figure(1);
            if m < nEigenFuncsToPlot
                if m == 0
                    plot(x, vP_x, '-.', 'LineWidth', 2, 'DisplayName', '$p(x)$');
                    xlim([-5 5])
                    hold on
                end
                plot(x, vPhi_m_x, 'LineWidth', 2, 'DisplayName', [ '$\phi_' num2str(m) '(x)$' ]);
                xlim([-5 5])
                xlabel('$x$', 'Interpreter', 'latex', 'FontSize', 14)
            end
            if m == nEigenFuncsToPlot + 1
                hold off
            %     title('Eigenfunctions')
                legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'northeast')
                print(fig1, [outputFolder filesep 'fig1_eigenfunctions'], '-depsc')
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


if b_verifyEigOrth    
    dx = 0.01;
    x = -1e3:dx:1e3-dx;
    I = zeros(M, M);
    
    vP_x = p(x,sigma);
    
    mPhi = zeros(length(x), M);
    vLambda_m = zeros(M, 1);
    for m = 0:M-1 
        [vPhi_m_x, lambda_m] = SqExpEig(a, b, m, x);
        vLambda_m(m+1) = lambda_m;
        mPhi(:, m+1) = vPhi_m_x;
    end
    
    for k = 0:M-1
        for m = 0:M-1
%             [vPhi_m_x, lambda_m] = SqExpEig(a, b, m, x);
%             [vPhi_k_x, ~] = SqExpEig(a, b, k, x);
%             I(m+1,k+1) = sum(vPhi_m_x.*vPhi_k_x.*vP_x*dx);
            I(m+1,k+1) = trapz(x, vPhi_m_x(:, m+1).*vPhi_m_x(:, k+1).*vP_x);
            if abs(I(m+1,k+1)) < 1e-12
                I(m+1,k+1) = 0;
            end
        end
    end
    isalmostequal(I, eye(M));
    fprintf('All %d eigenfunctions are orthonormal\n', M);
end    
    
if b_verifyMercersTheorm    
    y = -100:5:100;
    x = y;
    lhs = zeros(length(y), length(x));
    rhs = zeros(length(y), length(x));
    for i = 1:length(x)
        for j = 1:length(y)
            for m = 0:M-1
                [phi_m_x, lambda_m] = SqExpEig(a, b, m, x(i));
                [phi_m_y, ~] = SqExpEig(a, b, m, y(j));
                rhs(i,j) = rhs(i,j) + lambda_m*phi_m_x*phi_m_y;
            end
            lhs(i,j) = kernel(x(i),y(j),l);
        end
    end
    
    isalmostequal(rhs, lhs)
    fprintf('Mercer''s theorm confirmed\n');
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
x = (-5:dx:5-dx)';
N = length(x);

noiseVar1 = 0.4;
noiseVar2 = 0.4;

A1 = 5;
B1 = 0.001;
C1 = 0.001;
A2 = 10;
B2 = 7;
C2 = 6;
% mF      = [ A*exp(-x).*(sin(2.5*x) + sin(2*pi*x))   A*exp(-2*x).*sin(5*x) ];
[phi1, ~] = SqExpEig(a, b, 1, x);
[phi2, ~] = SqExpEig(a, b, 6, x);
[phi3, ~] = SqExpEig(a, b, 10, x);
mF      = [ A1*phi1 + B1*phi2 + C1*phi3   A2*exp(-0.2*x.^2).*sin(pi*x) + B2*exp(-0.5*x.^2).*sin(2*pi*x) + B2*exp(-0.3*x.^2).*sin(1*pi*x) ];
mF_awgn = [sqrt(noiseVar1)*randn(N,1) sqrt(noiseVar2)*randn(N,1)];
% cFstr   = {'10e^{-x}\big(\sin(2.5x) + \sin(2\pi x)\big)' '10e^{-2x}\sin(5x)'};
nFuncs  = size(mF, 2);

% x_eval = ([-1:0.1:1])';
% vF_eval = [10*exp(-2*x).*sin(5*x)];

cFigs = cell(1, nFuncs);

for i = 1:nFuncs
    vFi = mF(:, i);
    vFi_awgn = mF_awgn(:, i);
    SNR = snr(vFi, vFi_awgn);
    vGi = vFi + vFi_awgn;
    mPhi = zeros(N, M);
    for m = 0:M-1 
        [vPhi_m_x, lambda_m] = SqExpEig(a, b, m, x);
        mPhi(:, m+1) = vPhi_m_x;
    end
    step = 0.02*N;
    r = N/step;
    assert(M <= r, 'You cannot have less points than eigenfunctions!');
    vR = 1:step:N; %randperm(N,N/r); % 1:r:N;
    vC = pinv(mPhi) * vGi;
    vCR = pinv(mPhi(vR, :)) * vGi(vR);
    vCR_th = vCR.*(abs(vCR) > 1e-2);
%     disp([(0:M-1)'  vCR_th])
    vFi_hat = mPhi * vCR;
    accuracy = 100*(1 - norm(vFi_hat - vFi)/norm(vFi));
%     % ---------------------------------------------------------------------
%     % With kernel try
%     vFi_eval = vF_eval(:, i);
%     vPhi_m_x = zeros(M, 1);
%     vFi_hat2 = zeros(length(x_eval), 1);
%     for j = 1:length(x_eval)
%         for m = 0:M-1
%             vKernel_x_t = kernel(x_eval(j),t,l);
%             [vPhi_m_t, lambda_m] = SqExpEig(a, b, m, t);
%             vPhi_m_x(m+1) = trapz(t, vKernel_x_t.*vPhi_m_t.*vP_t);
%         end
%         vFi_hat2(j) = sum(vCR.*vPhi_m_x);
%         res = 10*exp(-2*x_eval(j)).*sin(5*x_eval(j));
%     end
%     % ---------------------------------------------------------------------
    cFigs{i} = figure(i+1);
    p1 = plot(x, vGi, 'Color',       '#4DBEEE', ...
                      'LineWidth',    2, ...
                      'DisplayName', ['$g(x)$']);
%     title({['$f_' num2str(i) '(x)$ with SNR = ' num2str(SNR) ' [dB]']; num2str(accuracy)}, 'Interpreter', 'latex', 'FontSize', 14)
    hold on
    p2 = plot(x, vFi, 'Color',       '#0072BD', ...
                      'LineWidth',    3, ...
                      'DisplayName', ['$f(x)$']);
    p3 = plot(x, vFi_hat, 'Color', '#D95319', ...
                          'LineWidth', 3, ...
                          'LineStyle', '-.', ...
                          'DisplayName', ['$\hat{f}(x)$']);

    p4 = plot(x(vR), vGi(vR), 'ro');
    hold off
    legend([p2 p1 p3], 'Interpreter', 'latex', 'FontSize', 14, 'Location', 'best')
    print(cFigs{i}, [outputFolder filesep 'fig' num2str(i+1) '_extrapolate_f' num2str(i)], '-depsc')
    fprintf('f%d : r = %d; M = %d; SNR = %.2f; Accuracy = %.2f%%\n', i, r, M, SNR, accuracy);
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
vPhi_m = (1/sqrt(2^m*factorial(m)*sqrt(a/c))) * exp( -(c-a)*x.^2 ) .* vHm;

end

function vPr = p(y, sigma)
vPr = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );
end

function vKernel = kernel(x,y,l) 
vKernel = exp(-(x-y).^2./(2*l^2));
end