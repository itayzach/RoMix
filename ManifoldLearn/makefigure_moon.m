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
rng(1);

%% RLS params
gamma_A_rls = 0.03125;
gamma_I_rls = 0;

%% LapRLS params
gamma_A_laprls = 0.03125;
gamma_I_laprls = 1;
kernel_sigma_laprls = 1/(6*sqrt(2)); 

%% EigRLS params
gamma_A_eigrls = 0.03125;
gamma_I_eigrls = 1;  

%% sParams
sParams = GetParameters();
assert(sParams.dim == 2, 'dim should be 2.')
kernel_sigma_eigrls = sParams.omega;

if kernel_sigma_laprls ~= kernel_sigma_eigrls
    error('LapRLS and EigRLS kernel width is different');
end
if sParams.sSim.b_plotEigenFigs
    [ mPhi_K, vLambda_K ] = CalcAnalyticEigenfunctions(sParams);
    [ mPhi_A, vLambda_A ] = CalcNumericEigenvectors(sParams);
    PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A);
    PlotSpectrum(sParams, vLambda_K, vLambda_A);
end

%% Load dataset
l=1; % number of labeled examples
scale_factor = 1;
xTrain = scale_factor*sParams.sDataset.x(1:end,:);
xTest = scale_factor*sParams.sDataset.xt;
yTrain = sParams.sDataset.y(1:end);
yTest = sParams.sDataset.yt;

xMax = max(max(xTrain,[],1));
xMin = min(min(xTrain,[],1));

pos=find(yTrain==1);
neg=find(yTrain==-1);
ipos=randperm(length(pos));
ineg=randperm(length(neg));
y1=zeros(length(yTrain),1);
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

experiment_moon(xTrain,y1,xTest,yTest,'rlsc',gamma_A_rls,gamma_I_rls, options_rls);
fprintf('---------------------------------------------------\n')

%% LapRLS - original
options_laprls_orig = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma_laprls, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
experiment_moon(xTrain,y1,xTest,yTest,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls_orig);
fprintf('---------------------------------------------------\n')
%% LapRLS - Laplacian and f from the same kernel
options_laprls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_laprls, ...
                   'GraphWeightParam', kernel_sigma_laprls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
           
laprlsc_classifier = experiment_moon(xTrain,y1,xTest,yTest,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls);
fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions
options_eigrls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_eigrls, ...
                   'GraphWeightParam', kernel_sigma_eigrls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
               
eigrls_classifier = experiment_moon(xTrain,y1,xTest,yTest,'eigrls',gamma_A_eigrls,gamma_I_eigrls, options_eigrls);
fprintf('Using %d eigenfunctions\n', sParams.ExtrplM-sParams.FirstM);
fprintf('---------------------------------------------------\n')

%% Plot classifier using retreived alpahas
%--------------------------------------------------------------------------
% Plot alpha difference
%--------------------------------------------------------------------------
alpha_err_diff = abs(laprlsc_classifier.alpha - eigrls_classifier.alpha);
figure;
plot(alpha_err_diff)
ylim([0 max(alpha_err_diff)]);
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
    plot2D(xTrain,y1,15);
set(gcf,'Position',[10+600 250 600 700])


%--------------------------------------------------------------------------
% Plot EigRLS (Phi*c)
%--------------------------------------------------------------------------
c = eigrls_classifier.c;
mPhi_m_x = zeros(length(xTrain), sParams.ExtrplM-sParams.FirstM);
mPhi_m_X = zeros(length(X), sParams.ExtrplM-sParams.FirstM);
mPhi_m_xt = zeros(length(xTest),sParams.ExtrplM-sParams.FirstM);
vLambda = zeros(1, sParams.ExtrplM-sParams.FirstM);
for i = sParams.FirstM:sParams.ExtrplM-1 
    m = OneDim2TwoDimIndex(i);
    mPhi_m_X(:, i+1) = phi(sParams,m,X);
    vLambda(i+1) = lambda(sParams,m);
    mPhi_m_x(:, i+1) = phi(sParams,m,xTrain);
    mPhi_m_xt(:, i+1) = phi(sParams,m,xTest);
end


vPhi_X_c = mPhi_m_X*c;
vPhi_xt_c = mPhi_m_xt*c;
mPhi_X_c = reshape(vPhi_X_c,length(x1),length(x2));
% mPhi_xt_c = reshape(vPhi_xt_c,length(x1),length(x2));

figure;
sgtitle(sprintf('EigRLS (Phi c) \n gamma_A = %.4f, gamma_I = %.4f\nTest error = %.1f%%', eigrls_classifier.gammas(1), eigrls_classifier.gammas(2), eigrls_classifier.test_error));
subplot(2,1,1)
    surf(XX1,XX2,mPhi_X_c, 'edgecolor', 'none')
    hold on;
    scatter3(xTest(:,1), xTest(:,2), vPhi_xt_c, 'filled');
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
subplot(2,1,2)
    contourf(XX1,XX2,mPhi_X_c,[0 0]);shading flat;
    hold on;
    plot2D(xTrain,y1,15);
set(gcf,'Position',[10+600+600 250 600 700]) 

%--------------------------------------------------------------------------
% Plot c
%--------------------------------------------------------------------------
figure; 
stem(c);
title('$c$', 'Interpreter', 'latex');
c_backup = c;
% c(c < 1e-3) = 0;
% %--------------------------------------------------------------------------
% % Plot EigRLS (Phi*Lambda*Phi')alpha
% %--------------------------------------------------------------------------
% alpha_eigrls = eigrls_classifier.alpha;
% 
% mPLP_X_x = mPhi_m_X*diag(vLambda)*mPhi_m_x.';
% mPLP_xt_x = mPhi_m_xt*diag(vLambda)*mPhi_m_x.';
% vPLPa_X_x = mPLP_X_x * alpha_eigrls;
% vPLPa_xt_x = mPLP_xt_x * alpha_eigrls;
% mPLPa_X_x = reshape(vPLPa_X_x, length(x1),length(x2));
% 
% isalmostequal(mPLP_X_x, K, 1e-10, '', false);
% figure;
% subplot(2,1,1);
%     imagesc(K); colorbar;
%     title('kernel');
% 
% subplot(2,1,2);
%     imagesc(mPLP_X_x); colorbar;
%     title('mercer');
%             
%             
% figure;
% sgtitle(sprintf('EigRLS (Phi*Lambda*Phi'')alpha\n gamma_A = %.4f, gamma_I = %.4f\nTest error = %.1f%%', eigrls_classifier.gammas(1), eigrls_classifier.gammas(2), eigrls_classifier.test_error));
% subplot(2,1,1)
%     surf(XX1,XX2,mPLPa_X_x, 'edgecolor', 'none')
%     hold on;
%     scatter3(xt(:,1), xt(:,2), vPLPa_xt_x, 'filled');
%     xlabel('$x_1$', 'Interpreter', 'latex')
%     ylabel('$x_2$', 'Interpreter', 'latex')
%     zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
%     colorbar;
% subplot(2,1,2)
%     contourf(XX1,XX2,mPLPa_X_x,[0 0]);shading flat;
%     hold on;
%     plot2D(x,y1,15);
% set(gcf,'Position',[10+600+600 250 600 700])  


%--------------------------------------------------------------------------
% Plot f difference
%--------------------------------------------------------------------------
figure; 
contourf(XX1,XX2,(abs(mKa - mPhi_X_c)));
colorbar;
title('$ | K \alpha - \Phi c |$', 'Interpreter', 'latex');