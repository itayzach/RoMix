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
gamma_I_rls = 0.03125;
gamma_A_rls = 0;

%% LapRLS params
gamma_I_laprls = 0.03125;
gamma_A_laprls = 1;
kernel_sigma_laprls = 1/sqrt(2)/10; 

%% EigRLS params
gamma_I_eigrls = 100;
gamma_A_eigrls = 1000;

a_k_eigrls = [ 0.01 0.01];
kernel_sigma_eigrls = 1/sqrt(2); 
b_k_eigrls = 1/(2*kernel_sigma_eigrls^2);

scale_factor = 1;
M = 150;
%% old
% a_k = [0.0001 0.0001];
% kernel_sigma = 1/sqrt(2)/10; 
% scale_factor = 1;
% 
% 
% b_k = 1/(2*kernel_sigma^2);
% M = 150;        

%% sParams
sParams.dim = 2;
sParams.constsType = 1;
sParams.a = a_k_eigrls;
sParams.ell = kernel_sigma_eigrls; % kernel width
sParams.b = 1/(2*kernel_sigma_eigrls^2);
sParams.sigma = 1./(2*sParams.a);
sParams.mu = 0*ones(1, sParams.dim);
sParams.c = sqrt(sParams.a.^2 + 2*sParams.a.*sParams.b);
sParams.A = sParams.a + sParams.b + sParams.c;
sParams.B = sParams.b./sParams.A;

%% Load dataset
load 2moons.mat;

l=10; % number of labeled examples

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


% can run a search over best parameters below

% earlier version
% subplot(2,3,1); plot2D(x,y1,12);axis([-1.5 2.5 -1.5 2.5]);
% subplot(2,3,2); decision_surface(x,y1,xt,yt,'svm',-5:0); hold on;
%                 plot2D(x,y1,12);
% subplot(2,3,3); decision_surface(x,y1,xt,yt,'rlsc',-5:0);
%                 hold on;  plot2D(x,y1,10);
% subplot(2,3,4); decision_surface(x,y1,xt,yt,'tsvm',-5:0);
%                 hold on;  plot2D(x,y1,10);
% subplot(2,3,5); decision_surface(x,y1,xt,yt,'lapsvm',[-5:0]);
%                 hold on;  plot2D(x,y1,10);
% subplot(2,3,6); decision_surface(x,y1,xt,yt,'laprlsc',[-5:0]);
%                  hold on;  plot2D(x,y1,10);

% new version
%subplot(1,3,1); plot2D(x,y1,12);axis([-1.5 2.5 -1.5 2.5]);

%% Figure setup
% fig = figure;
% fig_left_loc = -1500;
% fig_bottom_loc = 100;
% fig_width = 750;
% fig_height = 750;
% set(fig,'position',[fig_left_loc,fig_bottom_loc,fig_width,fig_height])



%% RLS
options_rls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_laprls, ...0.35, ...
                   'GraphWeightParam', 1, ...
                   'GraphWeights', 'binary', ...
                   'GraphNormalize', true);
% subplot(2,2,1); 
experiment_moon(x,y1,xt,yt,'rlsc',gamma_I_rls,gamma_A_rls, options_rls);
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')
%subplot(2,3,3); decision_surface(x,y1,xt,yt,'rlsc',-5:0);
%               hold on;  plot2D(x,y1,10);

% warning('For TSVM, Install SVMLight and a matlab interface to it. See http://www.cis.tugraz.at/igi/aschwaig/software.html');
% Use below with Antons's code
%subplot(2,3,2); experiment_moon(x,y1,xt,yt,'tsvm',0.125,0);
%               hold on;  plot2D(x,y1,10);


%% EigRLS
% subplot(1,2,2);
% experiment_moon(x,y1,xt,yt,'eigrls',0.03125,1);
% hold on;  
% plot2D(x,y1,15);

%% LapRLS - original
options_laprls_orig = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma_laprls, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
% subplot(2,2,2);
experiment_moon(x,y1,xt,yt,'laprlsc',gamma_I_laprls,gamma_A_laprls, options_laprls_orig);
% hold on;  
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')
%% LapRLS - Laplacian and f from the same kernel
options_laprls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_laprls, ...
                   'GraphWeightParam', kernel_sigma_laprls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
           
% subplot(2,2,3);
laprlsc_classifier = experiment_moon(x,y1,xt,yt,'laprlsc',gamma_I_laprls,gamma_A_laprls, options_laprls);
alpha_laprls = laprlsc_classifier.alpha;
XTrain_laprls = laprlsc_classifier.xtrain;
% hold on;  
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions
options_eigrls=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma_eigrls, ...
                   'GraphWeightParam', kernel_sigma_eigrls, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
               
options_eigrls.a_k = a_k_eigrls;
options_eigrls.b_k = b_k_eigrls;
options_eigrls.M = M;   

% subplot(2,2,4);
eigrls_classifier = experiment_moon(x,y1,xt,yt,'eigrls',gamma_I_eigrls,gamma_A_eigrls, options_eigrls);
alpha_eigrls = eigrls_classifier.alpha;
XTrain_eigrls = eigrls_classifier.xtrain;
% hold on;  
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')



%% Plot classifier using retreived alpahas
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);

X=[XX1(:) XX2(:)];

%--------------------------------------------------------------------------
K = calckernel(options_laprls.Kernel,options_laprls.KernelParam, XTrain_laprls, X);
Ka = K*alpha_laprls;
mKa = reshape(Ka,length(x1),length(x2));

figure;
sgtitle(sprintf('LapRLS - Laplacian and f from the same kernel\n gamma_A = %.4f, gamma_I = %.4f\nTest error = %d', laprlsc_classifier.gammas(1), laprlsc_classifier.gammas(2), laprlsc_classifier.test_error));
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
%--------------------------------------------------------------------------
% dim = 2;
% nxPoints = length(X);
% nyPoints = length(XTrain_laprls);
% 
% 
% tPhi_x = zeros(nxPoints, options.M, dim);
% tPhi_y = zeros(nyPoints, options.M, dim);
% vLambda = zeros(1, options.M);
% for m = 0:options.M-1    
%     [tPhi_x(:,m+1,1), vLambda(m+1)] = SqExpEig(a_k, b_k, m, X(:,1));
%     [tPhi_y(:,m+1,1), ~] = SqExpEig(a_k, b_k, m, XTrain_laprls(:,1));
%     [tPhi_x(:,m+1,2), ~] = SqExpEig(a_k, b_k, m, X(:,2));
%     [tPhi_y(:,m+1,2), ~] = SqExpEig(a_k, b_k, m, XTrain_laprls(:,2));
% 
% end
% 
% 
% k_mercer = zeros(nxPoints, nyPoints);
% lhs =  zeros(nxPoints, nyPoints);
% 
% for i = 1:length(X)
%     for j = 1:length(XTrain_laprls)
%         vLambda_Phix_Phiy = sum(vLambda.*tPhi_x(i,:,:).*tPhi_y(j,:,:), 2);
%         k_mercer(i,j) = prod(vLambda_Phix_Phiy, 3);
%         k_d = kernel(X(i,:), XTrain_laprls(j,:), kernel_sigma);
%         lhs(i,j) = prod(k_d, 2);
%     end
% end
% 
% figure;
% subplot(2,1,1);
% imagesc(lhs); colorbar;
% title('kernel');
% 
% subplot(2,1,2);
% imagesc(k_mercer); colorbar;
% title('mercer');
% 
% % assert(all(all(k_mercer - lhs < 10e-8)));
% % assert(all(all(k_mercer - K < 10e-8)));
% 
% %--------------------------------------------------------------------------
% Ka = k_mercer*alpha_laprls;
% mKa = reshape(Ka,length(x1),length(x2));
% 
% figure;
% title('LapRLS - mercer instead of kernel');
% subplot(2,1,1)
%     surf(XX1,XX2,mKa, 'edgecolor', 'none')
%     xlabel('$x_1$', 'Interpreter', 'latex')
%     ylabel('$x_2$', 'Interpreter', 'latex')
%     zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
%     colorbar;
% subplot(2,1,2)
%     contourf(XX1,XX2,mKa,[0 0]);shading flat;
%     hold on;
%     plot2D(x,y1,15);

%--------------------------------------------------------------------------


Phi_c = zeros(length(x1),length(x2));
Phi_c_test = zeros(length(xt),1);
for i = 0:options_eigrls.M-1 
    m = OneDim2TwoDimIndex(i, sParams.dim);
    vPhi_m_x = phi(sParams,m,X);
    vPhi_m_x = reshape(vPhi_m_x, length(x1),length(x2));
    Phi_c = Phi_c + alpha_eigrls(i+1) * vPhi_m_x ;
    
    vPhi_m_xt = phi(sParams,m,xt);
    Phi_c_test = Phi_c_test + alpha_eigrls(i+1) * vPhi_m_xt;
end


figure;
sgtitle(sprintf('EigRLS\n gamma_A = %.4f, gamma_I = %.4f\nTest error = %d', eigrls_classifier.gammas(1), eigrls_classifier.gammas(2), eigrls_classifier.test_error));

subplot(2,1,1)
    surf(XX1,XX2,Phi_c, 'edgecolor', 'none')
    hold on;
    scatter3(xt(:,1), xt(:,2), Phi_c_test, 'filled');
    xlabel('$x_1$', 'Interpreter', 'latex')
    ylabel('$x_2$', 'Interpreter', 'latex')
    zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
    colorbar;
subplot(2,1,2)
    contourf(XX1,XX2,Phi_c,[0 0]);shading flat;
    hold on;
    plot2D(x,y1,15);