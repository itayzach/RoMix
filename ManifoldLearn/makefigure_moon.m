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

%% Set kernel sigma

% kernel_sigma = 0.35; 

% a_k = 1;
% kernel_sigma = 1/sqrt(2); 
% scale_factor = 10;
% a_k = 0.00095;
a_k = 0.0001;
kernel_sigma = 1/sqrt(2)/10; 
scale_factor = 1;


b_k = 1/(2*kernel_sigma^2);
M = 50;        

%% Load dataset
load 2moons.mat;

l=1; % number of labeled examples

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
options=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma, ...0.35, ...
                   'GraphWeightParam', 1, ...
                   'GraphWeights', 'binary', ...
                   'GraphNormalize', true);
% subplot(2,2,1); 
experiment_moon(x,y1,xt,yt,'rlsc',0.03125,0, options);
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
options = ml_options('gamma_A', 0.1, ...
                     'NN',6, ...
                     'Kernel','rbf',...
                     'KernelParam',kernel_sigma, ...0.35...
                     'GraphWeightParam',1, ...
                     'GraphWeights', 'binary', ...
                     'GraphNormalize',true);    
             
% subplot(2,2,2);
experiment_moon(x,y1,xt,yt,'laprlsc',0.03125,1, options);
% hold on;  
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')
%% LapRLS - Laplacian and f from the same kernel
options=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma, ...
                   'GraphWeightParam', kernel_sigma, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);

options.a_k = a_k;
options.b_k = b_k;
options.M = M;                
               
% subplot(2,2,3);
[alpha_laprls, XTrain_laprls] = experiment_moon(x,y1,xt,yt,'laprlsc',0.03125,1, options);
% hold on;  
% plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions
% options=ml_options('gamma_A',0.1, ...
%                    'NN',6, ...
%                    'Kernel', 'rbf', ...
%                    'KernelParam', kernel_sigma, ...
%                    'GraphWeightParam', kernel_sigma, ...
%                    'GraphWeights', 'my_heat', ...
%                    'GraphNormalize', false);
%                
% options.a_k = a_k;
% options.b_k = b_k;
% options.M = M;   
% 
% % subplot(2,2,4);
% [alpha_eigrls, ~] = experiment_moon(x,y1,xt,yt,'eigrls',0.03125, 1, options);
% % hold on;  
% % plot2D(x,y1,15);
% fprintf('---------------------------------------------------\n')



%% Plot classifier using retreived alpahas
step = (xMax - xMin)/100;
x1 = xMin:step:xMax;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);

X=[XX1(:) XX2(:)];
K=calckernel(options.Kernel,options.KernelParam, XTrain_laprls, X);
Ka = K*alpha_laprls;
mKa = reshape(Ka,length(x1),length(x2));

figure;
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

% mPhi_xy_plane = zeros(length(X), options.M);
% mPhi_training = zeros(length(XTrain_laprls), options.M);
% for m = 0:options.M-1 
%     [vPhi_m_x1, lambda] = SqExpEig(options.a_k, options.b_k, m, X(:,1));
%     [vPhi_m_x2, ~] = SqExpEig(options.a_k, options.b_k, m, X(:,2));
% 
%     mPhi_xy_plane(:, m+1) = vPhi_m_x1.*vPhi_m_x2;   
%     vLambda(m+1) = lambda;
% 
%     [vPhi_m_xt1, ~] = SqExpEig(options.a_k, options.b_k, m, XTrain_laprls(:,1));
%     [vPhi_m_xt2, ~] = SqExpEig(options.a_k, options.b_k, m, XTrain_laprls(:,2));
% 
%     mPhi_training(:, m+1) = vPhi_m_xt1.*vPhi_m_xt2;
% end
% 
% for m = 0:options.M-1 
%     k_mercer = k_mercer + vLambda(m+1)*repmat(mPhi_xy_plane(:,m+1), 1, length(XTrain_laprls)).*repmat(mPhi_training(:,m+1).',length(X),1); 
% end
dim = 2;
nxPoints = length(X);
nyPoints = length(XTrain_laprls);


tPhi_x = zeros(nxPoints, options.M, dim);
tPhi_y = zeros(nyPoints, options.M, dim);
vLambda = zeros(1, options.M);
for m = 0:options.M-1    
    [tPhi_x(:,m+1,1), vLambda(m+1)] = SqExpEig(a_k, b_k, m, X(:,1));
    [tPhi_y(:,m+1,1), ~] = SqExpEig(a_k, b_k, m, XTrain_laprls(:,1));
    [tPhi_x(:,m+1,2), ~] = SqExpEig(a_k, b_k, m, X(:,2));
    [tPhi_y(:,m+1,2), ~] = SqExpEig(a_k, b_k, m, XTrain_laprls(:,2));

end


k_mercer = zeros(nxPoints, nyPoints);
lhs =  zeros(nxPoints, nyPoints);

for i = 1:length(X)
    for j = 1:length(XTrain_laprls)
        vLambda_Phix_Phiy = sum(vLambda.*tPhi_x(i,:,:).*tPhi_y(j,:,:), 2);
        k_mercer(i,j) = prod(vLambda_Phix_Phiy, 3);
        k_d = kernel(X(i,:), XTrain_laprls(j,:), kernel_sigma);
        lhs(i,j) = prod(k_d, 2);
    end
end

    figure;
    subplot(2,1,1);
    imagesc(lhs); colorbar;
    title('kernel');
    
    subplot(2,1,2);
    imagesc(k_mercer); colorbar;
    title('mercer');

% assert(all(all(k_mercer - lhs < 10e-8)));
% assert(all(all(k_mercer - K < 10e-8)));

Ka = k_mercer*alpha_laprls;
mKa = reshape(Ka,length(x1),length(x2));

figure;
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


% Phi_c = zeros(length(x1),length(x2));
% for m = 0:options.M-1 
%     [vPhi_m_x1, lambda_m1] = SqExpEig(options.a_k, options.b_k, m, x1);
%     [vPhi_m_x2, lambda_m2] = SqExpEig(options.a_k, options.b_k, m, x2);
%     
%     Phi_c = Phi_c + alpha_eigrls(m+1) * vPhi_m_x1.' * vPhi_m_x2 ;
% end
% 
% 
% figure;
% surf(XX1,XX2,(Phi_c), 'edgecolor', 'none')
% % view(2)
% hold on;
% plot2D(x,y1,15);
% xlabel('$x_1$', 'Interpreter', 'latex')
% ylabel('$x_2$', 'Interpreter', 'latex')
% zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
% colorbar;
