% 2 Moons Figure
% For TSVM, uncomment some lines below after correctly installing
% SVM Light and a matlab interface to it
% See http://www.cis.tugraz.at/igi/aschwaig/software.html
% Also, you can uncomment some lines to make this code search for best
% parameters.
%
% Author: Vikas Sindhwani (vikass@cs.uchicago.edu)

%% Restart run
close all; clear; clc;
rng(0);

%% RLS params
gamma_A_rls = 0.03125;
gamma_I_rls = 0;

%% LapRLS params
gamma_A_laprls = 0.03125;
gamma_I_laprls = 1;
kernel_sigma_laprls = 1/sqrt(2)/4; 

%% EigRLS params
gamma_A_eigrls = 0.03125;
gamma_I_eigrls = 1;  

%% sParams
[sParams, sSimParams] = GetParameters();
assert(sParams.dim == 2, 'dim should be 2.')
kernel_sigma_eigrls = sParams.omega;

if kernel_sigma_laprls ~= kernel_sigma_eigrls
    error('LapRLS and EigRLS kernel width is different');
end
if sSimParams.b_plotEigenFigs
    [ mPhi_K, vLambda_K ] = CalcAnalyticEigenfunctions(sParams);
    mPhi_A = [];
    PlotEigenfunctionsEigenvectors(sParams, sSimParams, mPhi_K, mPhi_A);
end

%% Load dataset
load 2moons.mat;

l=2; % number of labeled examples
scale_factor = 1;
x = scale_factor*x;
xt = scale_factor*xt;

xMax = max(max(x,[],1));
xMin = min(min(x,[],1));

pos=find(y==1);
neg=find(y==-1);
ipos=randperm(length(pos));
ineg=randperm(length(neg));
y1=zeros(length(y),1);
y1(pos(ipos(1:l)))=1;
y1(neg(ineg(1:l)))=-1;

%% RLS
options_rls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_laprls, ...0.35, ...
                   'GraphWeightParam', 1, ...
                   'GraphWeights', 'binary', ...
                   'GraphNormalize', true);

experiment_moon(x,y1,xt,yt,'rlsc',gamma_A_rls,gamma_I_rls, options_rls);
fprintf('---------------------------------------------------\n')

%% LapRLS - original
options_laprls_orig = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma_laprls, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
experiment_moon(x,y1,xt,yt,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls_orig);
fprintf('---------------------------------------------------\n')
%% LapRLS - Laplacian and f from the same kernel
options_laprls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_laprls, ...
                   'GraphWeightParam', kernel_sigma_laprls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
           
laprlsc_classifier = experiment_moon(x,y1,xt,yt,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls);
fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions
options_eigrls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_eigrls, ...
                   'GraphWeightParam', kernel_sigma_eigrls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
               
eigrls_classifier = experiment_moon(x,y1,xt,yt,'eigrls',gamma_A_eigrls,gamma_I_eigrls, options_eigrls);
fprintf('Using %d eigenfunctions\n', sParams.ExtrplM);
fprintf('---------------------------------------------------\n')

%% Plot classifier using retreived alpahas
%--------------------------------------------------------------------------
% Plot alpha difference
%--------------------------------------------------------------------------
figure;
plot(abs(laprlsc_classifier.alpha - eigrls_classifier.alpha))
ylim([-0.05 0.05]);
title('$|\alpha(LapRLS)-\alpha(EigRLS)|$', 'Interpreter', 'latex');

%--------------------------------------------------------------------------
% Plot test error difference
%--------------------------------------------------------------------------
fprintf([newline newline 'Test error diff: ' num2str(abs(laprlsc_classifier.test_error - eigrls_classifier.test_error)) '\n'])

%--------------------------------------------------------------------------
% Generate axis
%--------------------------------------------------------------------------
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);

X=[XX1(:) XX2(:)];

%--------------------------------------------------------------------------
% Plot LapRLS
%--------------------------------------------------------------------------
alpha_laprls = laprlsc_classifier.alpha;
XTrain_laprls = laprlsc_classifier.xtrain;
K = calckernel(options_laprls.Kernel,options_laprls.KernelParam, XTrain_laprls, X);
Ka = K*alpha_laprls;
mKa = reshape(Ka,length(x1),length(x2));

figure;
sgtitle(sprintf('LapRLS - Laplacian and f from the same kernel\n gamma_A = %.4f, gamma_I = %.4f\nTest error = %.1f%%', laprlsc_classifier.gammas(1), laprlsc_classifier.gammas(2), laprlsc_classifier.test_error));
subplot(2,1,1)
    surf(XX1,XX2,mKa, 'edgecolor', 'none')
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
subplot(2,1,2)
    contourf(XX1,XX2,mKa,[0 0]);shading flat;
    hold on;
    plot2D(x,y1,15);
set(gcf,'Position',[10+600 250 600 700])
%--------------------------------------------------------------------------
% Plot EigRLS
%--------------------------------------------------------------------------
alpha_eigrls = eigrls_classifier.alpha;
mPhi_m_x = zeros(length(x), sParams.ExtrplM);
mPhi_m_X = zeros(length(X), sParams.ExtrplM);
mPhi_m_xt = zeros(length(xt),sParams.ExtrplM);
vLambda = zeros(1, sParams.ExtrplM);
for i = 0:sParams.ExtrplM-1 
    m = OneDim2TwoDimIndex(i);
    mPhi_m_X(:, i+1) = phi(sParams,m,X);
    vLambda(i+1) = lambda(sParams,m);
    mPhi_m_x(:, i+1) = phi(sParams,m,x);
    mPhi_m_xt(:, i+1) = phi(sParams,m,xt);
end
mPLP_X_x = mPhi_m_X*diag(vLambda)*mPhi_m_x.';
mPLP_xt_x = mPhi_m_xt*diag(vLambda)*mPhi_m_x.';
vPLPa_X_x = mPLP_X_x * alpha_eigrls;
vPLPa_X_xt = mPLP_xt_x * alpha_eigrls;
mPLPa_X_x = reshape(vPLPa_X_x, length(x1),length(x2));

isalmostequal(mPLP_X_x, K, 1e-10, '', false);
figure;
subplot(2,1,1);
    imagesc(K); colorbar;
    title('kernel');

subplot(2,1,2);
    imagesc(mPLP_X_x); colorbar;
    title('mercer');
            
            
figure;
sgtitle(sprintf('EigRLS\n gamma_A = %.4f, gamma_I = %.4f\nTest error = %.1f%%', eigrls_classifier.gammas(1), eigrls_classifier.gammas(2), eigrls_classifier.test_error));
subplot(2,1,1)
    surf(XX1,XX2,mPLPa_X_x, 'edgecolor', 'none')
    hold on;
    scatter3(xt(:,1), xt(:,2), vPLPa_X_xt, 'filled');
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
subplot(2,1,2)
    contourf(XX1,XX2,mPLPa_X_x,[0 0]);shading flat;
    hold on;
    plot2D(x,y1,15);
set(gcf,'Position',[10+600+600 250 600 700])    