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

%% Load dataset
load 2moons.mat;

l=1; % number of labeled examples

% x = 10*x;
% xt = 10*xt;

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
fig = figure();
fig_left_loc = -1500;
fig_bottom_loc = 100;
fig_width = 750;
fig_height = 750;
set(fig,'position',[fig_left_loc,fig_bottom_loc,fig_width,fig_height])

%% Set kernel sigma

% kernel_sigma = 0.35; 
kernel_sigma = 1/sqrt(2)/10; 



%% RLS
options=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma, ...0.35, ...
                   'GraphWeightParam', 1, ...
                   'GraphWeights', 'binary', ...
                   'GraphNormalize', true);
subplot(2,2,1); 
experiment_moon(x,y1,xt,yt,'rlsc',0.03125,0, options);
plot2D(x,y1,15);
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
             
subplot(2,2,2);
experiment_moon(x,y1,xt,yt,'laprlsc',0.03125,1, options);
hold on;  
plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')
%% LapRLS - Laplacian and f from the same kernel
options=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma, ...
                   'GraphWeightParam', kernel_sigma, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
subplot(2,2,3);
[alpha_laprls, XTrain_laprls] = experiment_moon(x,y1,xt,yt,'laprlsc',0.03125,1, options);
hold on;  
plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')

%% EigRLS - Laplacian and f from the same kernel using eigenfunctions


options=ml_options('gamma_A',0.1, ...
                   'NN',6, ...
                   'Kernel', 'rbf', ...
                   'KernelParam', kernel_sigma, ...
                   'GraphWeightParam', kernel_sigma, ...
                   'GraphWeights', 'my_heat', ...
                   'GraphNormalize', false);
               
options.a_k = 105;
options.b_k = 1/(2*kernel_sigma^2);

options.M = 15; 

subplot(2,2,4);
[alpha_eigrls, ~] = experiment_moon(x,y1,xt,yt,'eigrls',0.03125, 1, options);
hold on;  
plot2D(x,y1,15);
fprintf('---------------------------------------------------\n')

%% Plot classifier using retreived alpahas
x1 = -2:0.1:3;
x2 = x1;
[XX1,XX2] = meshgrid(x1,x2);

X=[XX1(:) XX2(:)];
K=calckernel(options.Kernel,options.KernelParam, XTrain_laprls, X);
Ka = K*alpha_laprls;
mKa = reshape(Ka,length(x1),length(x2));

figure;
surf(XX1,XX2,mKa)
view(2)
hold on;
plot2D(x,y1,15);
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
colorbar;




Phi_c = zeros(length(x1),length(x2));
for m = 0:options.M-1 
    [vPhi_m_x1, lambda_m1] = SqExpEig(options.a_k, options.b_k, m, x1);
    [vPhi_m_x2, lambda_m2] = SqExpEig(options.a_k, options.b_k, m, x2);
    
    Phi_c = Phi_c + alpha_eigrls(m+1) * vPhi_m_x1.' * vPhi_m_x2 ;
end


figure;
surf(XX1,XX2,(Phi_c))
view(2)
hold on;
plot2D(x,y1,15);
xlabel('$x_1$', 'Interpreter', 'latex')
ylabel('$x_2$', 'Interpreter', 'latex')
zlabel('$f(x_1,x_2)$', 'Interpreter', 'latex')
colorbar;