%% Restart run
clear; close all; clc;
rng(0); % set seed

outputFolder = 'figs';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder)
else
    % Remove all pdf files
    filePattern = fullfile(outputFolder, '*.pdf');
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
            nEigenFuncsToPlot = 3;
            fig1 = figure(1);
            if m <= nEigenFuncsToPlot
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
                legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
                print(fig1, [outputFolder filesep 'fig1_eigenfunctions'], '-dpdf')
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
x = (0:dx:5-dx)';
N = length(x);

dt = 0.01;
t = (-20:dt:20-dt)';
vP_t = p(t,sigma);

mF      = [10*exp(-x).*(sin(2.5*x) + sin(2*pi*x))   10*exp(-2*x).*sin(5*x) ];
mF_awgn = [0.1*randn(N,1) 0.1*randn(N,1)];
cFstr = {'10e^{-x}\big(\sin(2.5x) + \sin(2\pi x)\big)' '10e^{-2x}\sin(5x)'};
nFuncs = size(mF, 2);

% x_eval = ([-1:0.1:1])';
% vF_eval = [10*exp(-2*x).*sin(5*x)];

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
    r = 0.01*N;
    vR = 1:r:N; %randperm(N,N/r); % 1:r:N;
    vC = pinv(mPhi) * vGi;
    vCR = pinv(mPhi(vR, :)) * vGi(vR);
    vFi_hat = mPhi * vCR;
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
    p1 = plot(x, vFi_hat, 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = \Phi c $ with ' num2str(N/r) ' points']);
    hold on
    title(['$f_' num2str(i) ' - \Phi c = ' num2str(norm(vFi_hat-vFi)) '$'], 'Interpreter', 'latex', 'FontSize', 12)
%     p2 = plot(x_eval, vFi_hat2);
    p3 = plot(x, vGi, '-.', 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = ' cFstr{i} '$ + AWGN']);
    p4 = plot(x, vFi, 'b--', 'LineWidth', 2, 'DisplayName', ['$f_' num2str(i) ' = ' cFstr{i} '$']);
    p5 = plot(x(vR), vGi(vR), 'o');
    hold off
    legend([p1 p3 p4], 'Interpreter', 'latex', 'FontSize', 12, 'Location', 'northeast')
    print(cFigs{i}, [outputFolder filesep 'fig' num2str(i+1) '_extrapolate_f' num2str(i) '.pdf'], '-dpdf')
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