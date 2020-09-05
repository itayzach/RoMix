clc; clear; close all; rng('default')

%% Graph
N = 1000;
G = gsp_sensor(N);

%% Vertices
v = G.coords;
figure; 
subplot(211); scatter(v(:,1), v(:,2), 'filled')
subplot(212); hist3(v,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');

%% Weights
dists = squareform(pdist(v, 'euclidean'));
omega = 0.1;
G.W = sparse(exp(-dists.^2/omega^2)); % Gaussian kernel
figure; imagesc(G.W); colorbar; title('Kerel (weights) matrix')
G = gsp_graph_default_parameters(G);

%% Plot the graph
param.show_edges = false;
figure; gsp_plot_graph(G,param);

%% Signal on the graph
G = gsp_compute_fourier_basis(G);
paramf.log = 1;
Nf = 5;
g = gsp_design_warped_translates(G, Nf,paramf);  
s = sign(G.U(:,2));
sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
signal = sf(:,2);

%% Sample the signal
sample_ind = randi(N, 0.1*N, 1);
sampled_signal = signal(sample_ind);

%% Interpolate
interpolated_signal = gsp_interpolate(G,sampled_signal,sample_ind);

figure;
subplot(121); 
    gsp_plot_signal(G,signal,param); 
    title('Original signal');
subplot(122); 
    gsp_plot_signal(G,interpolated_signal,param); 
    title('gsp interpolate (caxis [-1,1])'); caxis([-1 1])
set(gcf,'Position', [400 400 1100 400])

figure;
for i=1:10
    subplot(2,5,i)
    [V, mLambdaNumeric] = eigs(G.W,10);
    [vLambdaNumeric, idx] = sort(diag(mLambdaNumeric), 'descend');
    V = V(:,idx);
    gsp_plot_signal(G,V(:,i),param); 
    title(['v_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaNumeric(i), '%.4f')]);
end
set(gcf,'Position', [200 200 1600 600])