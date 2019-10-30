%%
clear; close all; clc;

%% Parameters
a = 1;
b = 3;
L = 10;


%% Eigenfunctions
syms x y
phi = cell(L, 1);
lambda = cell(L, 1);
for i = 0:L-1
    [phi_k, lambda_k] = SEeig(a, b, i, y);
    phi{i+1} = phi_k;
    lambda{i+1} = lambda_k;
end

%% p(x)
sigma = 1/sqrt(4*a);
l = 1/sqrt(2*b);
p(y) = (1/sqrt(2*pi*sigma^2)) * exp( -y.^2/(2*sigma^2) );

%% k(x,y)
k(x,y) = exp(-(x-y)^2/(2*l^2));

%% Set values
dx = 0.01;
x = -10:dx:10-dx;

%% Substitue
% phi_k0 = subs(phi_k0);
% phi_k1 = subs(phi_k1);
% phi_k2 = subs(phi_k2);
% p = subs(p);

%% Check eigenfunctions
for i = 0:L-1
    p_k(y) = phi{i+1};
    rhs(y) = lambda{i+1} * p_k(y);
    lhs(y) = sum(k(x,y).*p_k(x).*p(x)*dx);

    y0 = -2:0.5:2;
    rhsEval = single(subs(rhs(y0)));
    lhsEval = single(subs(lhs(y0)));
    fprintf('k = %d\n', i)
    disp([rhsEval; lhsEval])
%     norm(rhsEval - lhsEval)/(norm(rhsEval)*norm(lhsEval))
end

%% Plot
x = -1:dx:1-dx;
figure;
hold on
for i = 0:3
    p_k(y) = phi{i+1};
    plot(x, subs(p_k(x)), 'DisplayName', [ 'k = ' num2str(i) ]);
    
end
hold off
legend show



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