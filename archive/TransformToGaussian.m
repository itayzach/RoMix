clear; clc; close all;

sigma1 = 0.5;
mu1 = 5;
sigma2 = 0.1;
mu2 = 4;
N = 1000;
x = [mu1 + sigma1*randn(floor(N/3),1); mu2 + sigma2*randn(floor(N/3),1); (mu1 - mu2)*rand(ceil(N/3),1) + mu2];

%% Random projection?
% dataDim = 1;
% newDim = 1000;
% R = randn(dataDim, newDim);
% for d = 1:dataDim
%    R(d,:) = R(d,:)/norm(R(d,:)); % Normalization
% end
% y = x*R;
% x_hat = y*R';

%% Just a random matrix multiplication.
% linear combination of the rows by the coeffs (values of x) should return a gaussian vector

M = 1000;
R = randn(N, N);
for d = 1:N
   R(d,:) = R(d,:)/norm(R(d,:)); % Normalization
end
y = R*x;
x_hat = R^(-1)*y;


figure;
subplot(2,2,1); histogram(x,100); title('$x$ histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);
subplot(2,2,2); histogram(x_hat,100); title('$\hat{x} = R^{-1}y$  histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);
subplot(2,2,3); histfit(R(1,:),100); title('$R_{[1,:]}$ histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);
subplot(2,2,4); histfit(y,100); title('$y = R x$ histogram', 'Interpreter', 'latex', 'FontSize', 14); set(gca,'FontSize', 14);

set(gcf,'Position', [600 200 800 600])