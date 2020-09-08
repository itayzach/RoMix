clc; clear; close all; rng('default')

%% Parameters
N = 2500;
nEigs = 30;
dataDim = 2;
nComponents = 1;
omega = 0.3;

GMMRegVal = 0;
dataDist = 'Gaussian';

sSimParams.PlotEigenFuncsM = min(nEigs, 20);
sSimParams.outputFolder = 'figs';
%% Generate data
sDataset = GenerateDataset(dataDist, dataDim, nComponents, N, 0);
v = sDataset.sData.x;
if dataDim == 1
    v_for_graph = [v zeros(size(v))];
else
    v_for_graph = v;
end

%% Estimate distribution and get kernel parameters
sDistParams = EstimateDistributionParameters(sDataset, nComponents, GMMRegVal);
sKernelParams = GetKernelParams(sDataset, sDistParams, omega);
[sKernelParams.vLambdaAnaytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, sDataset.dim, nComponents);


%% Adjacency (gaussian kernel)
W = sparse(CalcAdjacency(sKernelParams, v));
%% Graph
G1 = gsp_graph(W, v_for_graph);
G = G1;
G.lap_type = 'normalized';
G = gsp_graph_default_parameters(G);

% G = gsp_sensor(N);
% G = gsp_2dgrid(N);
% G = gsp_two_moons(0.05,N,0.1);
%% Plot the graph
param.show_edges = false;
% figure; gsp_plot_graph(G,param);

% figure; 
%     scatter(G.coords(:,1), G.coords(:,2), 'filled')
%     set(gca,'FontSize', 14);
% figure; 
%     hist3(G.coords,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
%     colormap(gca, 'hot')
%     colorbar()
%     view(2)
%     title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
%     set(gca,'FontSize', 14);
figure; 
subplot(131)
    imagesc(G.W); 
    colorbar; 
    title('Kerel (weights) matrix', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);
subplot(132); 
    imagesc(diag(G.d)); 
    colorbar; 
    title('Degree matrix', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);      
subplot(133); 
    imagesc(G.L); 
    colorbar; 
    title('Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
    set(gca,'FontSize', 14);    
set(gcf,'Position', [100 200 1800 400])

%% Eigenvectors
[V, vLambdaNumeric] = CalcNumericEigenvectors(nEigs, sKernelParams, v);
PlotEigenfuncvecScatter(sSimParams, dataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [])
% figure;
% for i=1:nEigs
%     subplot(2,5,i)
%     gsp_plot_signal(G,V(:,i),param); 
%     title(['v_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaNumeric(i), '%.4f')]);
% end
% set(gcf,'Position', [200 200 1600 600])


%% Eigenfunctions
[mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
PlotEigenfuncvecScatter(sSimParams, dataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [])
% figure;
% for i=1:nEigs
%     subplot(2,5,i)
%     gsp_plot_signal(G,mPhiAnalytic(:,i),param); 
%     title(['\phi_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaAnalytic(i), '%.4f')]);
% end
% set(gcf,'Position', [200 200 1600 600])
%% Signal on the graph
G = gsp_compute_fourier_basis(G);
figure;
for i=1:10
    subplot(2,5,i)
    gsp_plot_signal(G,G.U(:,i),param); 
    title(['u_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(G.e(i), '%.4f')]);
end
set(gcf,'Position', [200 200 1600 600])
paramf.log = 1;
Nf = 5;
g = gsp_design_warped_translates(G, Nf,paramf);  
s = sign(G.U(:,2));
sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
orig_signal = sf(:,2);

% orig_c = [0 0 100 0 0 0 0 0 0 0]';
% orig_signal = G.U(:,1:nEigs)*orig_c;
%% Sample the signal
randPerm = randperm(N)';
sample_ind = sort(randPerm(1:floor(0.5*N)));
sampled_signal = orig_signal(sample_ind);
sampled_signal_padded = zeros(size(orig_signal));
sampled_signal_padded(sample_ind) = orig_signal(sample_ind);
%% GSP interpolation
gsp_interp_signal = gsp_interpolate(G,sampled_signal,sample_ind);

%% My interpolation
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0;0.01;  
c = eigrls(sampled_signal_padded, mPhiAnalytic, vLambdaAnalytic, gamma_A_eigrls, gamma_I_eigrls, G1.L);
fprintf('c =\n\t');
fprintf('%f\n\t', c);
fprintf('\n');
my_interp_signal = mPhiAnalytic*c;

%% Plot interpolation
figure;
subplot(131); 
    gsp_plot_signal(G,orig_signal,param); 
    title('Original signal');
subplot(132); 
    gsp_plot_signal(G,gsp_interp_signal,param); 
    title('gsp interpolate');
subplot(133); 
    gsp_plot_signal(G,my_interp_signal,param); 
    title('Mine!');    
set(gcf,'Position', [400 400 1500 400])

fprintf('GSP interpolation error: %.3f\n', norm(orig_signal - gsp_interp_signal)/norm(orig_signal));
fprintf('Our interpolation error: %.3f\n', norm(orig_signal - my_interp_signal)/norm(orig_signal));
