%%
clc; clear; close all;
rng('default');
%%
mu = 0; sigma = 1; 
omega = 0.1; 

% mixture)
N = 5001;
x = mu + sigma*randn(N,1);

%%
[Phi_a, lambda_a] = CalcAnalytical(x, mu, sigma, omega);

[ x_sorted, sorted_ind ] = sort(x);
[~, inv_sorted_ind] = sort(sorted_ind);
dx_sorted = [0; diff(x_sorted)];
dx = dx_sorted(inv_sorted_ind);

p_of_x = normpdf(x, mu, sigma);
% p_of_x =(1/sqrt(2*pi*sigma^2)) * exp( -(x-mu).^2./(2*sigma^2) );


P_of_x = p_of_x .* dx;

figure; imagesc(Phi_a'*diag(P_of_x)*Phi_a); colormap('jet'); colorbar;

function [Phi_a, lambda_a] = CalcAnalytical(x,mu,sigma,omega)
if nargin<4 || isempty(omega)
    omega = 1;
end
if nargin<3 || isempty(sigma)
    sigma = 1;
end
N = length(x);
beta = 2*sigma^2/omega^2;
hermite_arg = (1/4 + beta/2)^(1/4)*(x-mu)/sigma;
M = 5;
Hm = zeros(N,M);
Phi_a = zeros(N,M);
lambda_a = zeros(M,1);
exp_term = exp( -((x-mu).^2/(2*sigma^2)) * ((sqrt(1+2*beta)-1)/2) );
for m=0:M-1
    Hm(:,m+1) = hermiteH(m, hermite_arg);
    normFactor = (1+2*beta)^(1/8)/sqrt(2^m*factorial(m));
    Phi_a(:,m+1) = normFactor * exp_term .* Hm(:,m+1);
    lambda_a = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m;
end
end