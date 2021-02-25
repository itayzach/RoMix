%%
% 17.8 will use a small subset of gaussian samples, transform to chi
% square the distances and make new N(0,1) samples.
%
% 13.8 draft. still does not work well
% 10.8 test for non Gaussian data

clc; clear; close all;
rng('default');

%%
mu = 0; sigma = 1; 
omega = 0.5; 
beta = 2*sigma^2/omega^2;

% Generate N points from some distribution (here, this is a gaussian
% mixture)
N = 1000;
mu1 = 1; sigma1 = 0.5; 
mu2 = 3; sigma2 = 0.3;
% v = [mu1 + sigma1*randn(floor(N/2),1); mu2 + sigma2*randn(ceil(N/2),1) ]; 
v = mu1 + sigma1*randn(N,1);
% v = (10 - 5)*rand(N, 1) + 5;

% first N/2 points are nodes close to the point mu1. last N/2 points are
% nodes close to mu2.
% also note that sigma1=sigma2=omega, since I want the Gaussian kernel to
% reflect correctly the closeness of points on the original manifold

d = pdist(v, 'euclidean');
z = d.^2/(2*omega^2); % normalized squared distances. Should be distributed as chi-square, since here, the nodes are Gaussian points in R
figure(1); clf;
h = cdfplot(z); hold on;
k = size(v,2);
plot(z, cdf('Chisquare', z, k), '.');
legend('z estimated CDF', 'Chi squared CDF');
h.LineStyle = ':';
h.Marker = 's';
title('CDF of squared distances on original manifold');

%%
% In the general case, the squared distances on the original manifold are
% not chi-square distributed.
% Transform squared distances to chi-squared distribution using the inverse
% probability integral transform:

cdf_est = h.YData'; %CDF estimated values
ind = cdf_est > 0.01 & cdf_est < 0.99;
cdf_est = cdf_est(ind);
z_icdf = h.XData';
z_icdf = z_icdf(ind);
z_tilde_icdf = icdf('Chisquare', cdf_est, k); %G^{-1} mapping

% Model G^{-1}(F) which maps original squared distances to new squared
% distances
PolyOrder = 5;
A = ModelPolyMatrix(z_icdf, PolyOrder);
options.TolX = 1;
% c_hat = lsqnonneg(A, z_tilde_icdf(:), options); % to assure monotone solution restricting to non negative coefs
c_hat = pinv(A'*A)*(A'*z_tilde_icdf);
c_hat = max(0,c_hat);

figure(2); clf;
plot(z_icdf, A*c_hat, '.');
xlabel('z');
ylabel('tilde z');
grid;

figure;
histogram(z_icdf);

% mle(z_icdf,'pdf',@(x,v,s)chi2pdf(x/s,v)/s,'start',[1,200]),

%% Use this model to obtain square distances on the new manifold:
B = ModelPolyMatrix(z(:), PolyOrder);
z_tilde=B*c_hat;
figure;
histogram(z_tilde); % z_tilde should be chi squared


d_tilde = real(sqrt(z_tilde));
figure;
histogram(d_tilde);

[Y, lam] = cmdscale(d_tilde',15);
v_tilde = Y(:,1);

% v_tilde should be gaussian
figure;
histogram(v_tilde,100);

figure(3); clf;
subplot(211); plot(v, v_tilde, 'x'); xlabel('v'); ylabel('tilde v'); title('mapping from original nodes positions to new nodes positions');
subplot(212); plot(z, d_tilde, 'x', sort(z), sort(z), 'k.-'); 
xlabel('d'); ylabel('tilde d'); title('mapping from the original distances to new distances');

figure();
% histfit(v_tilde,100)
hist3(Y(:,1:2),'CdataMode','auto', 'Nbins', [50 50], 'EdgeColor', 'flat', 'FaceColor', 'interp');
colormap(gca, 'hot')
colorbar()
%% define mapping from w to tilde_w:
w_vec = exp(-z(:)); % kernel
w_tilde_vec = exp(-z_tilde); % kernel
C = ModelPolyMatrix(w_vec(:), PolyOrder);
q = pinv(C'*C)*(C'*w_tilde_vec); % least squares

figure(4); clf;
hold on; plot(w_vec, w_tilde_vec, '.', 'DisplayName', 'tilde w');
plot(w_vec, C*q, 'o', 'DisplayName', 'model');
xlabel('w');
legend();
% unfortunatelly this is not the mapping in the spectrum domain. We need to
% think how to compute that.

%% show numerical eigenvectors of the full set and analytical eigefunctions:
all_ind = 1:N; %use first cluster
z_all = abs((v(all_ind,:)-v(all_ind,:)').^2)/(2*omega^2);
W=exp(-z_all);
M = 5;
[V, La] = eigs(W,M);
figure(5); clf;
% subplot(311); plot(v(all_ind,1),V,'.'); title('numerical on v'); xlabel('v')

sigma_est = std(v_tilde);
b_normalize = true;
[Phi_a, lambda_a] = CalcAnalytical(v_tilde, mu, sigma_est, omega, b_normalize);
subplot(311), plot(v_tilde, Phi_a, '.'); title('analytical on new manifold'); xlabel('v tilde')


W_tilde=exp(-abs((v_tilde-v_tilde').^2)/(2*omega^2));
[tV, tLa] = eigs(W_tilde,M);
tV = FlipSign(Phi_a, tV);
subplot(312); plot(v_tilde, tV, '.'); title('numerical on the new manifold'); xlabel('v tilde')

subplot(313); plot(v_tilde, Phi_a-tV, '.'); title('error on the new manifold'); xlabel('v tilde')

[ v_tilde_sorted, v_tilde_sort_ind ] = sort(v_tilde);
[~, v_tilde_inverse_sort] = sort(v_tilde_sort_ind);
dx_sorted = [0; diff(v_tilde_sorted)];
dx_vec = dx_sorted(v_tilde_inverse_sort);
p_of_x = (1./sqrt(2*pi*sigma_est.^2)) .* exp( -(v_tilde-mu).^2./(2*sigma_est.^2) );
P_of_x = p_of_x .* dx_vec;

inner_product_matrix = Phi_a'*diag(P_of_x)*Phi_a;
figure;
imagesc(inner_product_matrix); colorbar; colormap('hot');

% % normalization issue...
% [rV, rLa] = eig(R); minusSqrtR=rV*diag(1./sqrt(diag(rLa)));
% % orthogonalize:
% nPhi_a = Phi_a*minusSqrtR; % normalized phi analytic
% R2 = nPhi_a'*nPhi_a;
% 
% figure;
% imagesc(R2); colorbar; colormap('hot');


% numeric on the new manifold:


%% make a noisy graph signal
f=V*randn(M,1); % Assume that the signal lies exactly in the lower numerical eigenspace
noisy_f = f+ 0.1*randn(size(f));
% now we will remove the noise by projecting onto the lower eigenspace:
% first use the original eigenvectors:
noisy_f_transform = V'*noisy_f;
f_cleaned = V*noisy_f_transform;

figure(6); clf;
h=plot(v(all_ind,1),f,'*',v(all_ind), noisy_f, 'x', v(all_ind), f_cleaned, 's', 'MarkerSize', 3);
xlabel('node v_i pos'); 
legend('original graph signal', 'measured noisy version', 'cleaned using full numerical eigenspace');
h(1).MarkerSize=7;

% clean using the new eigenvectors
f_cleaned2 = tV*(tV'*noisy_f);
hold on;
% f_cleaned3 = nPhi_a*(nPhi_a'*noisy_f);
f_cleaned3 = Phi_a*(Phi_a'*noisy_f);
plot(v(1:N), f_cleaned3, '.', 'DisplayName', 'cleaned by projecting onto W-tilde analytical and whitened vectors'); legend;

function A = ModelPolyMatrix(pos, PolyOrder)
A = zeros(length(pos), PolyOrder+1);
for n=0:PolyOrder
    A(:,n+1)=pos.^n;
end
end

function [Phi_a, lambda_a] = CalcAnalytical(x,mu,sigma,omega, b_normalize)
if nargin<5 || isempty(b_normalize)
    b_normalize = false;
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
    lambda_a(m+1) = sqrt(2/(1+beta+sqrt(1+2*beta))) * (beta/(1+beta+sqrt(1+2*beta)))^m;
end

if b_normalize
    n = length(x);
    Phi_a = (1/sqrt(n))*Phi_a;
end
end
