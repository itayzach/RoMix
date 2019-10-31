%% Restart run
clear; close all; clc;

%% Parameters
a = 1;
b = 3;
L = 10;

b_verifyEigenfunctions = false;
b_plotFigures = true;

%% Eigenfunctions
syms x y
cPhi = cell(L, 1);
cLambda = cell(L, 1);
for i = 0:L-1
    [phi_k, lambda_k] = SEeig(a, b, i, y);
    cPhi{i+1} = phi_k;
    cLambda{i+1} = lambda_k;
end

%% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);
p(y) = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );

%% kernel(x,y)
kernel(x,y) = exp(-(x-y)^2/(2*l^2));

%% Set values
dx = 0.01;
x = -10:dx:10-dx;

%% Substitue
% phi_k0 = subs(phi_k0);
% phi_k1 = subs(phi_k1);
% phi_k2 = subs(phi_k2);
% p = subs(p);

%% Check eigenfunctions
if b_verifyEigenfunctions
    for i = 0:L-1
        p_k(y) = cPhi{i+1};
        rhs(y) = cLambda{i+1} * p_k(y);
        lhs(y) = sum(kernel(x,y).*p_k(x).*p(x)*dx);

        y0 = -2:0.5:2;
        rhsEval = single(subs(rhs(y0)));
        lhsEval = single(subs(lhs(y0)));
        fprintf('k = %d\n', i)
        disp([rhsEval; lhsEval])
    %     norm(rhsEval - lhsEval)/(norm(rhsEval)*norm(lhsEval))
    end
end
%% Plot
figure;
tg = uitabgroup; % tabgroup
if b_plotFigures
    x = -1:dx:1-dx;
    thistab = uitab(tg, 'Title', 'Eigenfunctions'); % build tab
    axes('Parent',thistab); % somewhere to plot
    hold on
    for i = 0:2
        p_k(y) = cPhi{i+1};
        plot(x, subs(p_k(x)), 'DisplayName', [ '$\phi_' num2str(i) '$' ]);

    end
    hold off
    title('Eigenfunctions')
    legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
    legend show
end
%% Read mnist data
mnist = load('data/mnist.mat');
mXTest = single(mnist.testX);
vTestLabels = mnist.testY.';
[nTestPoints, nPixels] = size(mXTest);
N = nTestPoints;

%% Graph signals
nDigits = 10;
mS = zeros(N, nDigits);
for k = 0:9
    mS(:, k+1) = (vTestLabels == k);
end

%% Stam functions
N = 4000;
x = linspace(-2, 2, N)'; % [x_min x_max] should be small since gaussians decay to zero
dx = abs(x(2) - x(1));
% x = (-1:dx:1-dx)';
% f = [sin(-10*x)  cos(-5*x)  exp(-x)  exp(x)];
f = [cos(5*x) cos(2.5*x)];
nFuncs = size(f, 2);
% f = [exp(-x)  exp(x)];
% f = mS(:, 1);



thistab = uitab(tg, 'Title', 'Reconstruction'); % build tab
axes('Parent',thistab); % somewhere to plot
hold on
for i = 1:nFuncs
    fi = f(:, i);
    mPhi = zeros(N, L);
    for k = 0:L-1
        p_k(y) = cPhi{k+1};
        vPhikEval = single(subs(p_k(x)));
        mPhi(:, k+1) = vPhikEval;
    end
    r = 0.025*N;
    R = 1:r:N; %randperm(r);
    vCr = pinv(mPhi(R, :)) * fi(R);
    vC = pinv(mPhi) * (fi);
    disp([vCr vC])

    plot(x, mPhi * vCr, 'LineWidth', 2, 'DisplayName', ['$f_i \Phic \text{with c from ' num2str(N/r) ' points}$']);
    plot(x, mPhi * vC, 'DisplayName', ['$\Phic \text{with c from ' num2str(N) ' points} $']);
    plot(x, fi, '-.', 'DisplayName', ['$f_' num2str(i) '$']);
end
hold off
legend('Interpreter', 'latex', 'FontSize', 12, 'Location', 'best')
legend show


%% SEeig (Squared Exponentional)
function [phi_k, lambda_k] = SEeig(a, b, k, x)

% Calculate parameters
c = sqrt(a^2 + 2*a*b);
A = a + b + c;
B = b/A;

% k-th eigenvalue
lambda_k = sqrt(2*a/A) * B^k;

% k-th eigenfunction
Hk = hermiteH(k, sqrt(2*c)*x);
phi_k = simplify(exp( -(c-a)*x^2 ) * Hk);
end