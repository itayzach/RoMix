clc; clear; close all; rng('default')

%% Parameters
N = 2500;
nEigs = 30;
dataDim = 2;
nComponents = 1;
omega = 0.3;

sSimParams.PlotEigenFuncsM = min(nEigs, 20);
sSimParams.outputFolder = 'figs';
b_showFigures = true;
b_gaussianDataset = false;
estDataDist = 'Gaussian';
%% Generate dataset
if b_gaussianDataset
    sDataset = GenerateDataset('Gaussian', dataDim, nComponents, N, 0);
    GMMRegVal = 0;
else
    G_tmp = gsp_sensor(N);
    sDataset.sData.x = G_tmp.coords;
    sDataset.estDataDist = estDataDist;
    sDataset.dim = dataDim;
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
G.lap_type = 'normalized';
G = gsp_graph_default_parameters(G);

%% Plot the graph
param.show_edges = false;
if b_showFigures
    figure; 
        gsp_plot_graph(G,param);
    figure; 
        scatter(G.coords(:,1), G.coords(:,2), 'filled')
        set(gca,'FontSize', 14);
    figure; 
        hist3(G.coords,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        title('$\hat{p}({\bf x})$', 'Interpreter', 'latex', 'FontSize', 14)
        set(gca,'FontSize', 14);
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
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [])
end

%% Eigenfunctions
[mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
if b_showFigures
    PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [])
end

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
f = sf*[0 1 0 0 0]';
% orig_c = [0 0 100 0 0 0 0 0 0 0]';
% f = G.U(:,1:nEigs)*orig_c;

%% Sample the signal
samplingRatio = 0.1;

randPerm = randperm(N)';
sample_ind = sort(randPerm(1:floor(samplingRatio*N)));
sampled_signal = f(sample_ind);
sampled_signal_padded = zeros(size(f));
sampled_signal_padded(sample_ind) = f(sample_ind);
%% GSP interpolation
gsp_interp_signal = gsp_interpolate(G,sampled_signal,sample_ind);

%% My interpolation
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0;0.1;  
c = eigrls(sampled_signal_padded, mPhiAnalytic, vLambdaAnalytic, gamma_A_eigrls, gamma_I_eigrls, G.L);
fprintf('c =\n\t');
fprintf('%f\n\t', c);
fprintf('\n');
my_interp_signal = mPhiAnalytic*c;

%% Plot interpolation
figure;
subplot(131); 
    gsp_plot_signal(G,f,param); 
    hold on; scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
    title('Original signal');
subplot(132); 
    gsp_plot_signal(G,gsp_interp_signal,param); 
    hold on; scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
    title('gsp interpolate');
subplot(133); 
    gsp_plot_signal(G,my_interp_signal,param); 
    title('Mine!');    
set(gcf,'Position', [400 400 1500 400])

fprintf('GSP interpolation error: %.3f\n', norm(f - gsp_interp_signal)/norm(f));
fprintf('Our interpolation error: %.3f\n', norm(f - my_interp_signal)/norm(f));
