clc; clear; close all; rng('default')

%% Parameters
N = 2500;
nEigs = 30;
nComponents = 1;
omega = 0.3;

sSimParams.PlotEigenFuncsM = min(nEigs, 20);
sSimParams.outputFolder = 'figs';
b_showFigures = true;
b_gaussianDataset = false;
b_normlizedLaplacian = true;
estDataDist = 'Gaussian';
%% Generate dataset
if b_gaussianDataset
    dataDim = 2;
    sDataset = GenerateDataset('Gaussian', dataDim, nComponents, N, 0);
    GMMRegVal = 0;
else
    G_tmp = gsp_sensor(N);
%     G_tmp = gsp_bunny(); N = G_tmp.N;
%     G_tmp = gsp_minnesota(); N = G_tmp.N;
%     G_tmp = gsp_david_sensor_network(N);
    sDataset.sData.x = G_tmp.coords;
    sDataset.estDataDist = estDataDist;
    sDataset.dim = size(sDataset.sData.x,2);
    GMMRegVal = 0.1;
end
v = sDataset.sData.x; % For convenience
%% Estimate distribution and get kernel parameters
sDistParams = EstimateDistributionParameters(sDataset, nComponents, GMMRegVal);
sKernelParams = GetKernelParams(sDataset, sDistParams, omega);
[sKernelParams.vLambdaAnaytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, sDataset.dim, nComponents);

%% Adjacency (gaussian kernel)
W = sparse(CalcAdjacency(sKernelParams, v));
%% Graph
G = gsp_graph(W, v);
if b_normlizedLaplacian
    G.lap_type = 'normalized';
    G = gsp_graph_default_parameters(G);
end

%% Plot the graph
param.show_edges = false;
if b_showFigures
    if sDataset.dim == 2
        figure; 
            gsp_plot_graph(G,param);
        figure; 
            hist3(G.coords,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
            colormap(gca, 'hot')
            colorbar()
            view(2)
            title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
            set(gca,'FontSize', 14);
    end
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
end
%% Eigenvectors
[V, vLambdaNumeric] = CalcNumericEigenvectors(nEigs, sKernelParams, v);
if b_showFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [])
    figure;    
    for i=1:10
        subplot(2,5,i)
        gsp_plot_signal(G,V(:,i),param); 
        view(0,90);
        title(['v_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaNumeric(i), '%.4f')]);
    end
    set(gcf,'Position', [200 200 1600 600])
end

set(gcf,'Position', [200 200 1600 600])%% Eigenfunctions
[mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
if b_showFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [])
    figure;   
    for i=1:10
        subplot(2,5,i)
        gsp_plot_signal(G,mPhiAnalytic(:,i),param); 
        view(0,90)
        title(['\phi_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaAnalytic(i), '%.4f')]);
    end
    set(gcf,'Position', [200 200 1600 600])
end

%% Signal on the graph
G = gsp_compute_fourier_basis(G);
figure;
for i=1:10
    subplot(2,5,i)
    gsp_plot_signal(G,G.U(:,i),param); 
    view(0,90)
    title(['u_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(G.e(i), '%.4f')]);
end
set(gcf,'Position', [200 200 1600 600])

% attemp #1:
paramf.log = 1;
Nf = 5;
g = gsp_design_warped_translates(G, Nf,paramf);  
s = sign(G.U(:,2));
sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
f = sf*[0 1 0 0 0]';
% orig_c = [0 0 100 0 0 0 0 0 0 0]';
% f = G.U(:,1:nEigs)*orig_c;

% attempt #2:
% K = floor(samplingRatio*N);
% f_hat = randn(K,1);
% f_hat_full = [f_hat; zeros(N-K,1)];
% f = G.U*f_hat_full;

% attempt #3:
% f_tilde = randn(samplingRatio*N,1);
% f = G.U*[ f_tilde; zeros((1-samplingRatio)*N,1)];

%% Sample the signal
samplingRatio = 0.1;

% attemp #1:
randPerm = randperm(N)';
sample_ind = sort(randPerm(1:floor(samplingRatio*N)));

% attempt #2,#3:
% sample_ind = (1:samplingRatio*N)';

sampled_signal = f(sample_ind);
sampled_signal_padded = zeros(size(f));
sampled_signal_padded(sample_ind) = f(sample_ind);

%% GSP interpolation
gsp_interp_signal = gsp_interpolate(G,sampled_signal,sample_ind);

%% My interpolation using eigenFUNCTIONS
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0;0.1;  
c = eigrls(sampled_signal_padded, mPhiAnalytic, vLambdaAnalytic, gamma_A_eigrls, gamma_I_eigrls, G.L);
% fprintf('c =\n\t');
% fprintf('%f\n\t', c);
% fprintf('\n');
my_interp_signal = mPhiAnalytic*c;

%% My interpolation using eigenVECTORS
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0;0.1;  
c_eigenvecs = eigrls(sampled_signal_padded, V, vLambdaNumeric, gamma_A_eigrls, gamma_I_eigrls, G.L);
% fprintf('c =\n\t');
% fprintf('%f\n\t', c_eigenvecs);
% fprintf('\n');
eigenvecs_interp_signal = V*c_eigenvecs;

%% Plot interpolation
figure;
subplot(221); 
    gsp_plot_signal(G,f,param); 
%     hold on; scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
    view(0,90)
    title('Original signal');
subplot(222); 
    gsp_plot_signal(G,gsp_interp_signal,param); 
%     hold on; scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
    view(0,90)
    title(sprintf('gsp interpolate\n error: %.3f\n', norm(f - gsp_interp_signal)/norm(f)));
subplot(223); 
    gsp_plot_signal(G,eigenvecs_interp_signal,param); 
    title(sprintf('Using V\n error: %.3f\n', norm(f - eigenvecs_interp_signal)/norm(f)));   
    view(0,90)
subplot(224); 
    gsp_plot_signal(G,my_interp_signal,param); 
    title(['Using \Phi. ' sprintf('\nerror: %.3f\n', norm(f - my_interp_signal)/norm(f))]);
    view(0,90)
set(gcf,'Position', [400 100 800 800])

fprintf('GSP interpolation error: %.3f\n', norm(f - gsp_interp_signal)/norm(f));
fprintf('V   interpolation error: %.3f\n', norm(f - eigenvecs_interp_signal)/norm(f));
fprintf('Our interpolation error: %.3f\n', norm(f - my_interp_signal)/norm(f));
