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

%% number of labeled examples
l=1; 

%% RLS params
gamma_A_rls = 0.03125;
gamma_I_rls = 0;

%% LapRLS params
gamma_A_laprls = 0.03125;
gamma_I_laprls = 1;
kernel_sigma_laprls = 1/(6*sqrt(2)); 

%% EigRLS params
gamma_A_eigrls = 0.01;
gamma_I_eigrls = 0.1;  

%% sParams
sParams = GetParameters();
assert(sParams.dim == 2, 'dim should be 2.')
kernel_sigma_eigrls = sParams.omega;

% if kernel_sigma_laprls ~= kernel_sigma_eigrls
%     error('LapRLS and EigRLS kernel width is different');
% end
if sParams.sSim.b_plotEigenFigs
    [ mPhi_K, vLambda_K ] = CalcAnalyticEigenfunctions(sParams);
    [ mPhi_A, vLambda_A ] = CalcNumericEigenvectors(sParams);
    PlotEigenfunctionsEigenvectors(sParams, mPhi_K, mPhi_A, vLambda_K, vLambda_A);
    PlotSpectrum(sParams, vLambda_K, vLambda_A);
end

%% Load dataset
xTrain = sParams.sDataset.x(1:end,:);
xTest = sParams.sDataset.xt;
yTrain = sParams.sDataset.y(1:end);
yTest = sParams.sDataset.yt;

pos=find(yTrain==1);
neg=find(yTrain==-1);
unlab=find(yTrain==0);
ipos=randperm(length(pos));
ineg=randperm(length(neg));
y1=zeros(length(yTrain),1);
y1(pos(ipos(1:l)))=1;
y1(neg(ineg(1:l)))=-1;

%% Plot two-moons
x0     = 10;
y0     = 250;
width  = 600;
height = 400;

fig = figure;
plot2D(xTest,yTest,5,'k*');
plot2D(xTrain,y1,20,'ks');
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
set(gca,'FontSize', 14);
set(gcf,'Position',[x0 y0 width height])
saveas(fig,[sParams.sSim.outputFolder filesep 'fig_two_moons'], 'epsc');

%% RLS
% options_rls=ml_options('gamma_A',0.1, ...
%                    'NN',6, ...
%                    'Kernel', 'rbf', ...
%                    'KernelParam', kernel_sigma_laprls, ...0.35, ...
%                    'GraphWeightParam', 1, ...
%                    'GraphWeights', 'binary', ...
%                    'GraphNormalize', true);
% 
% experiment_moon(xTrain,y1,xTest,yTest,'rlsc',gamma_A_rls,gamma_I_rls, options_rls);
% fprintf('---------------------------------------------------\n')

%% LapRLS - original
options_laprls_orig = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma_laprls, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
orig_laprlsc_classifier = experiment_moon(xTrain,y1,xTest,yTest,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls_orig);
fprintf('---------------------------------------------------\n')
plot_classifier(orig_laprlsc_classifier, sParams.sSim.outputFolder, xTest, zeros(size(yTest)))

% Plot classifier using retreived alpahas
% fprintf([newline newline 'Test error diff: ' num2str(abs(orig_laprlsc_classifier.test_error - eigrls_classifier.test_error)) '\n'])
% plot_classifier(laprlsc_classifier, sParams.sSim.outputFolder)
% 
% %--------------------------------------------------------------------------
% % Plot alpha
% %--------------------------------------------------------------------------
% alpha = orig_laprlsc_classifier.alpha;
% figure; 
% stem(alpha);
% hold on;
% small_alpha_idx = find(abs(alpha) < 1e-5);
% pos=find(y1==1);
% neg=find(y1==-1);
% scatter([pos; neg], alpha([pos; neg]), 'filled');
% % yline(1e-5, 'r-', 'Threshold')
% title('$\alpha$', 'Interpreter', 'latex');

%% LapRLS - Laplacian and f from the same kernel
% options_laprls=ml_options('gamma_A',0.1, ...
%                    'NN',6, ...
%                    'Kernel', 'rbf', ...
%                    'KernelParam', kernel_sigma_laprls, ...
%                    'GraphWeightParam', kernel_sigma_laprls, ...
%                    'GraphWeights', 'my_heat', ...
%                    'GraphNormalize', false);
% 
% laprlsc_classifier = experiment_moon(xTrain,y1,xTest,yTest,'laprlsc',gamma_A_laprls,gamma_I_laprls, options_laprls);
% fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions
options_eigrls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_eigrls, ...
                   'GraphWeightParam', kernel_sigma_eigrls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
               
options_eigrls.sParams = sParams;
eigrls_classifier = experiment_moon(xTrain,y1,xTest,yTest,'eigrls',gamma_A_eigrls,gamma_I_eigrls, options_eigrls);
fprintf('Using %d eigenfunctions\n', sParams.ExtrplM-sParams.FirstM);
fprintf('---------------------------------------------------\n')

plot_classifier(eigrls_classifier, sParams.sSim.outputFolder, xTest, zeros(size(yTest)))

%--------------------------------------------------------------------------
% Plot c
%--------------------------------------------------------------------------
lambda = sParams.vLambda_K(1:sParams.ExtrplM);
fig = figure; 
plot(lambda, 'LineWidth', 2, 'DisplayName', '$\lambda_m$');
xlabel('$m$', 'Interpreter', 'latex', 'FontSize', 14);
hold on;
c_sq = eigrls_classifier.c.^2;
plot(c_sq, 'LineWidth', 2, 'DisplayName', '$c_m^2$');
plot(c_sq.^2 ./ lambda, 'LineWidth', 2, 'DisplayName', '$c_m^2/\lambda_m$');
legend('Interpreter', 'latex', 'FontSize', 14, 'Location', 'southwest')
% small_c_idx = find(abs(c_sq) < 1e-5);
% scatter(small_c_idx, c_sq(small_c_idx), 'filled');
set(gca,'FontSize', 14);
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'linear');
set(gcf,'Position', [x0 y0 600 400])
saveas(fig,[sParams.sSim.outputFolder filesep 'fig_rkhs_norm_decay'], 'epsc');
% %--------------------------------------------------------------------------
% % Plot f difference
% %--------------------------------------------------------------------------
% figure; 
% contourf(XX1,XX2,(abs(mKa - mPhi_X_c)));
% colorbar;
% title('$ | K \alpha - \Phi c |$', 'Interpreter', 'latex');