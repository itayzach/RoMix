clc; clear; close all; rng('default')

%% Parameters
nEigs = 30;
nComponents = 1;
omega = 0.3;
estDataDist = 'Gaussian';
sSimParams.PlotEigenFuncsM = min(nEigs, 20);
sSimParams.outputFolder = 'figs';

b_showFigures = true;
b_showEigenFigures = false;
b_showGraphMatricesFigures = false;

b_distributionDataset = false;
b_normlizedLaplacian = true;
verticesTransformation = 'mds'; % 'mds' / 'rp' / 'randOrth'
samplingRatio = 0.1;

%% Generate dataset
if b_distributionDataset
    N = 1000;
    dataDim = 1;
    graphName = 'Uniform';
    sDataset = GenerateDataset(graphName, dataDim, nComponents, N, 0);
    GMMRegVal = 0;
else
    N = 500;
%     N = round(sqrt(N))^2; G_tmp = gsp_2dgrid(sqrt(N),sqrt(N)); graphName = 'twodgrid';
%     G_tmp = gsp_sensor(N);  graphName = 'sensor';
%     G_tmp = gsp_bunny(); N = G_tmp.N; graphName = 'bunny';
    G_tmp = gsp_minnesota(); N = G_tmp.N; graphName = 'minnesota';
%     G_tmp = gsp_david_sensor_network(); N = G_tmp.N; graphName = 'david_seonsor';
%     G_tmp = gsp_swiss_roll(N); graphName = 'swiss_roll';
%     G_tmp = gsp_random_ring(N); graphName = 'random_ring';
    sDataset.sData.x = G_tmp.coords;
    sDataset.estDataDist = estDataDist;
    sDataset.dim = size(sDataset.sData.x,2);
    GMMRegVal = 0.1;
    dataDim = size(G_tmp.coords,2);
end
v = sDataset.sData.x; % For convenience
%% Estimate distribution and get kernel parameters
sDistParams = EstimateDistributionParameters(sDataset, nComponents, GMMRegVal);
sKernelParams = GetKernelParams(sDataset, sDistParams, omega);
[sKernelParams.vLambdaAnaytic, sKernelParams.vComponentIndex, sKernelParams.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParams, sDataset.dim, nComponents);

%% Adjacency (gaussian kernel)
[W, dist] = CalcAdjacency(sKernelParams, v);
W = sparse(W);
%% Graph
G = gsp_graph(W, v);
if b_normlizedLaplacian
    G.lap_type = 'normalized';
    G = gsp_graph_default_parameters(G);
end

%% Plot the graph
param.show_edges = false;
if b_showGraphMatricesFigures
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
        if b_normlizedLaplacian
            title('Normalized Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
        else
            title('Laplacian matrix', 'Interpreter', 'latex', 'FontSize', 14)
        end
        set(gca,'FontSize', 14);    
    set(gcf,'Position', [100 200 1800 400])
end
%% Eigenvectors
[V, vLambdaNumeric] = CalcNumericEigenvectors(N, sKernelParams, v);
if b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, V, vLambdaNumeric, 'Numeric', [])
    figure;    
    for i=1:min(10,nEigs)
        subplot(2,5,i)
        gsp_plot_signal(G,V(:,i),param); 
        view(0,90);
        title(['v_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaNumeric(i), '%.4f')]);
    end
    sgtitle('Adjacency numeric eigenvectors');
    set(gcf,'Position', [200 200 1600 600])

end
%% Eigenfunctions
[mPhiAnalytic, vLambdaAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParams, v, true);
if b_showEigenFigures
%     PlotEigenfuncvecScatter(sSimParams, estDataDist, v, [], 0, sSimParams.PlotEigenFuncsM-1, mPhiAnalytic, vLambdaAnalytic, 'Analytic', [])
    figure;   
    for i=1:min(10,nEigs)
        subplot(2,5,i)
        gsp_plot_signal(G,mPhiAnalytic(:,i),param); 
        view(0,90)
        title(['\phi_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(vLambdaAnalytic(i), '%.4f')]);
    end
    sgtitle('Adjacency analytic eigenfunctions');
    set(gcf,'Position', [200 200 1600 600])
end

%% Eigenvalues
if b_showEigenFigures
    PlotSpectrum(sSimParams, sDataset, [], vLambdaAnalytic, vLambdaNumeric, []);
end
%% Signal on the graph
if samplingRatio < 0.2
    b_plotSamplingPointsMarkers = true;
else
    b_plotSamplingPointsMarkers = false;
end

G = gsp_compute_fourier_basis(G);
if b_showEigenFigures
    figure;
    for i=1:min(10,nEigs)
        subplot(2,5,i)
        gsp_plot_signal(G,G.U(:,i),param); 
        view(0,90)
        title(['u_{' num2str(i) '}, \lambda_{' num2str(i) '} = ' num2str(G.e(i), '%.4f')]);
    end
    sgtitle('Laplacian eigenvectors');
    set(gcf,'Position', [200 200 1600 600])
end
% this:
% paramf.log = 1;
% Nf = 5;
% g = gsp_design_warped_translates(G, Nf,paramf);  
% s = sign(G.U(:,2));
% sf = gsp_vec2mat(gsp_filter_analysis(G,g,s),Nf);
% f = sf*[0 1 0 0 0]';
% f_hat = G.U.'*f;

% or this:
% orig_c = [0 0 100 0 0 0 0 0 0 0]';
% f = G.U(:,1:nEigs)*orig_c;

% or this:
% K = floor(samplingRatio*N);
% f_hat_K = 100*exp(-0.5*(1:K))';
% f_hat = [f_hat_K; zeros(N-K,1)];
% f = G.U*(f_hat);

% or this:
% f_hat = 5*exp(-0.5*(1:N))';
% f = G.U*(f_hat);

% or this:
k0 = round(0.01*N);
f_hat = zeros(N,1);
f_hat(1:k0) = 5*sort(abs(randn(k0,1)), 'descend');
f = G.U*f_hat;
% f = V.'*f_hat;

figure; 
subplot(1,2,1)
gsp_plot_signal_spectral(G,f_hat,param);
xlabel('$\lambda_k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
subplot(1,2,2)
stem(f_hat(1:k0+3));
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
set(gca,'FontSize', 14);
sgtitle('Graph signal spectrum', 'Interpreter', 'latex', 'FontSize', 14)
set(gcf,'Position', [100 200 1200 400])


%% Sample the signal
% this:
randPerm = randperm(N)';
sample_ind = sort(randPerm(1:round(samplingRatio*N)));
non_sample_ind = setdiff(1:N,sample_ind);
% or this:
% sample_ind = (1:samplingRatio*N)';

f_sampled = f(sample_ind);
f_sampled_padded = zeros(size(f));
f_sampled_padded(sample_ind) = f(sample_ind);

%% GSP interpolation
gsp_interp_signal = gsp_interpolate(G,f_sampled,sample_ind);

%% My interpolation using eigenFUNCTIONS
gamma_A_eigrls = 0;0.01;
gamma_I_eigrls = 0;0.1; 
c = eigrls(f_sampled_padded, mPhiAnalytic, vLambdaAnalytic, gamma_A_eigrls, gamma_I_eigrls, G.L);
f_int_no_transform = mPhiAnalytic*c;

%% My interpolation using eigenVECTORS
% gamma_A_eigrls = 0;0.01;
% gamma_I_eigrls = 0;0.1;  
% c_eigenvecs = eigrls(sampled_signal_padded, V, vLambdaNumeric, gamma_A_eigrls, gamma_I_eigrls, G.L);
% % fprintf('c =\n\t');
% % fprintf('%f\n\t', c_eigenvecs);
% % fprintf('\n');
% eigenvecs_interp_signal = V*c_eigenvecs;

%% Plot interpolation
fig = figure;
subplot(131);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline ...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(132); 
    if dataDim > 1
        gsp_plot_signal(G,gsp_interp_signal,param); 
    else
        plot(G.coords, gsp_interp_signal, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Pesenson''s interpolation' newline ...
        '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(gsp_interp_signal'*G.L*gsp_interp_signal, '%.3f') newline ...
        'error: ' num2str(norm(f - gsp_interp_signal)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(133); 
    if dataDim > 1
        gsp_plot_signal(G,f_int_no_transform,param);
    else
        plot(G.coords, f_int_no_transform, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    set(gca,'FontSize', 14);
    title(['Our interpolation (no transformation)' newline '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(f_int_no_transform'*G.L*f_int_no_transform, '%.3f') newline ...
        'error: ' num2str(norm(f - f_int_no_transform)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    view(0,90)
set(gcf,'Position', [400 100 1200 400])   
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_interp_no_transform_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');


%% Transform
% "Random projection"
if strcmp(verticesTransformation, 'rp')
    newDim = dataDim;
    R = randn(dataDim,newDim); % Random projection matrix
    for i = 1:dataDim
       R(i,:) = R(i,:)/norm(R(i,:)); % Normalization
    end
    v_tilde = v*R; % Projection
%     Rf = f;
elseif strcmp(verticesTransformation, 'randOrth')
    % % Linear combination of Gaussians
    % R = randn(N, N);
    % % for i = 1:N
    % %    R(:,i) = R(:,i)/norm(R(:,i)); % Normalization
    % % end
    % for i = 1:N
    %    R(i,:) = R(i,:)/norm(R(i,:)); % Normalization
    % end
    % v_tilde = R*v;

    % Random orthonormal matrix
    R = RandOrthMat(N);
    v_tilde = R*v;
%     Rf = R*f;
    f_tilde_sampled_padded = R*f_sampled_padded;

    % % Whitening matrix
    % [v_tilde, R, invR] = whiten(v,1e-5);
elseif strcmp(verticesTransformation, 'mds')
    [v_tilde, mdsLambda] = cmdscale(dist);
%     Rf = f;
    f_tilde_sampled_padded = f_sampled_padded;
end
if dataDim == 1
    histPlotRows = 1;
else
    histPlotRows = 3;
end
fig = figure; 
subplot(histPlotRows,2,1)
histogram(v(:,1),100);
set(gca,'FontSize', 14);
title('$v_1$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
subplot(histPlotRows,2,2)
histfit(v_tilde(:,1),100);
set(gca,'FontSize', 14);
title('$\tilde{v}_1$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
if dataDim == 2 || dataDim == 3
    subplot(histPlotRows,2,3)
    histogram(v(:,2),100);
    title('$v_2$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'FontSize', 14);
    subplot(histPlotRows,2,4)
    histfit(v_tilde(:,2),100);
    title('$\tilde{v}_2$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca,'FontSize', 14);
    
    if dataDim == 2
        subplot(histPlotRows,2,5)
        hist3(v,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        xlim([ min(v(:,1)) max(v(:,1))])
        ylim([ min(v(:,2)) max(v(:,2))])
        title('$v$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);

        subplot(histPlotRows,2,6)
        hist3(v_tilde,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
        colormap(gca, 'hot')
        colorbar()
        view(2)
        xlim([ min(v_tilde(:,1)) max(v_tilde(:,1))])
        ylim([ min(v_tilde(:,2)) max(v_tilde(:,2))])
        title('$\tilde{v}$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
    else
        subplot(histPlotRows,2,5)
        histogram(v(:,3),100);
        title('$v_3$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
        subplot(histPlotRows,2,6)
        histfit(v_tilde(:,3),100);
        title('$\tilde{v}_3$ histogram', 'Interpreter', 'latex', 'FontSize', 14);
        set(gca,'FontSize', 14);
    end    
end
set(gcf,'Position', [100 200 1000 histPlotRows*250])
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_hists_transformed'), 'epsc');

% Distances
ii = 5; jj = 60;
[pdist2(v(ii,:),v(jj,:)) pdist2(v_tilde(ii,:),v_tilde(jj,:))]

sDatasetTilde.sData.x = v_tilde;
sDatasetTilde.estDataDist = 'Gaussian';
sDatasetTilde.dim = size(v_tilde,2);
GMMRegVal = 0;
nTildeComponents = nComponents;
% Estimate distribution and get kernel parameters
sDistParamsTilde = EstimateDistributionParameters(sDatasetTilde, nTildeComponents, GMMRegVal);
sKernelParamsTilde = GetKernelParams(sDatasetTilde, sDistParamsTilde, omega);
[sKernelParamsTilde.vLambdaAnaytic, sKernelParamsTilde.vComponentIndex, sKernelParamsTilde.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParamsTilde, sDatasetTilde.dim, nTildeComponents);
% Adjacency (gaussian kernel)
[W_tilde, dist_tilde] = CalcAdjacency(sKernelParamsTilde, v_tilde);
W_tilde = sparse(W_tilde);


% Graph
G_tilde = gsp_graph(W_tilde, v_tilde);
if b_normlizedLaplacian
    G_tilde.lap_type = 'normalized';
    G_tilde = gsp_graph_default_parameters(G_tilde);
end

[mPhiTildeAnalytic, vLambdaTildeAnalytic] = CalcAnalyticEigenfunctions(nEigs, sKernelParamsTilde, v_tilde, true);
% [~, vLambdaTildeNumeric] = CalcNumericEigenvectors(nEigs, sKernelParamsTilde, v_tilde);

% Spectrum
% PlotSpectrum(sSimParams, sDatasetTilde, [], vLambdaTildeAnalytic, vLambdaTildeNumeric, []);

%% Find coeffs for graph-signal - EigRLS
tilde_gamma_A_eigrls = 0;0;
tilde_gamma_I_eigrls = 0;0.1;
c_tilde = eigrls(f_tilde_sampled_padded, mPhiTildeAnalytic, vLambdaTildeAnalytic, tilde_gamma_A_eigrls, tilde_gamma_I_eigrls, G.L);
% c_tilde = lsqminnorm(mPhiTildeAnalytic,R*f);
% c_tilde = (mPhiTildeAnalytic.'*mPhiTildeAnalytic)\(mPhiTildeAnalytic.'*(R*f));
% c_tilde = pinv(mPhiTildeAnalytic)*R*f;
% f_tilde = mPhiTildeAnalytic*c_tilde;
f_tilde_int_eigrls = mPhiTildeAnalytic*c_tilde;
if strcmp(verticesTransformation, 'randOrth')
    f_int_eigrls = R'*f_tilde_int_eigrls;
else
    f_int_eigrls = f_tilde_int_eigrls;
end

%% Find coeffs for graph-signal - fminsearch
fun = @(c)f_sampled'*G.L(sample_ind,sample_ind)*f_sampled - c'*mPhiTildeAnalytic.'*G_tilde.L*mPhiTildeAnalytic*c;
c0 = randn(nEigs,1);
options = optimset('MaxIter',100);
c_tilde = fminsearch(fun,c0,options);
f_tilde_int_fminsearch = mPhiTildeAnalytic*c_tilde;
if strcmp(verticesTransformation, 'randOrth')
    f_int_fminsearch = R'*f_tilde_int_fminsearch;
else
    f_int_fminsearch = f_tilde_int_fminsearch;
end


%% Reconstruction
% if strcmp(verticesTransformation, 'rp') || strcmp(verticesTransformation, 'mds')
%     f_rec = f;
%     f_int = f;
%     v_rec = v;
% elseif strcmp(verticesTransformation, 'randOrth')
%     %This is the same: [mPhiTildeAnalytic.'*f_tilde mPhiTildeAnalytic.'*R*f]
%     % keyboard;
% 
% 
%     % c_rec = (mPhiTildeAnalytic.'*mPhiTildeAnalytic)\(mPhiTildeAnalytic.'*(R'*f_tilde));
%     % f_reconstructed = mPhiTildeAnalytic*c_rec;
%     % f_reconstructed = pinv(mPhiTildeAnalytic.'*R)*mPhiTildeAnalytic.'*f_tilde;
%     % invR = R';
%     % f_reconstructed = (invR.'*invR)\(invR.'*f_tilde);
%     % f_reconstructed = lsqminnorm(R',f_tilde);
%     f_rec = R'*Rf;
%     f_int = R'*f_tilde_int_eigrls;
%     v_rec = R'*v_tilde;
% end
% 
% 
% sDatasetRec.sData.x = v_rec;
% sDatasetRec.estDataDist = 'Gaussian';
% sDatasetRec.dim = size(v_rec,2);
% GMMRegVal = 0;
% nComponents = 1;
% 
% % Estimate distribution and get kernel parameters
% sDistParamsRec = EstimateDistributionParameters(sDatasetRec, nComponents, GMMRegVal);
% sKernelParamsRec = GetKernelParams(sDatasetRec, sDistParamsRec, omega);
% [sKernelParamsRec.vLambdaAnaytic, sKernelParamsRec.vComponentIndex, sKernelParamsRec.vEigIndex] ...
%     = CalcAnalyticEigenvalues(nEigs, sKernelParamsRec, sDatasetRec.dim, nComponents);
% 
% % Adjacency (gaussian kernel)
% W_rec = sparse(CalcAdjacency(sKernelParamsRec, v_rec));
% 
% 
% % Graph
% G_rec = gsp_graph(W_rec, v_rec);
% if b_normlizedLaplacian
%     G_rec.lap_type = 'normalized';
%     G_rec = gsp_graph_default_parameters(G_rec);
% end

%% Plot EigRLS
fig = figure;
subplot(1,3,1);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,3,2); 
    if dataDim > 1
        gsp_plot_signal(G_tilde,f_tilde_int_eigrls,param); 
    else
        plot(G_tilde.coords, f_tilde_int_eigrls, '.');
        xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G_tilde.coords(sample_ind), f_tilde_sampled_padded(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end    
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Interpolated graph-signal on $\tilde{G}$' newline ...
           'using EigRLS' newline...
           '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_int_eigrls'*G_tilde.L*f_tilde_int_eigrls, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,3,3); 
    if dataDim > 1
        gsp_plot_signal(G,f_tilde_int_eigrls,param); 
    else
        plot(G.coords, f_tilde_int_eigrls, '.');
        xlabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Our interpolation (with transformation)' newline ...
           '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(f_tilde_int_eigrls'*G.L*f_tilde_int_eigrls, '%.3f') newline ...
           'error: ' num2str(norm(f - f_tilde_int_eigrls)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
set(gcf,'Position', [400 100 1200 400])   
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_eigrls_interp_with_transform_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');


fig = figure;
subplot(1,2,1);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,2,2); 
    if dataDim > 1
        gsp_plot_signal(G_tilde,f_tilde_int_eigrls,param); 
    else
        plot(G_tilde.coords, f_tilde_int_eigrls, '.');
        xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G_tilde.coords(sample_ind), f_tilde_sampled_padded(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end    
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Interpolated graph-signal on $\tilde{G}$' newline ...
           'using EigRLS' newline...
           '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_int_eigrls'*G_tilde.L*f_tilde_int_eigrls, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
set(gcf,'Position', [400 100 800 400])   
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_eigrls_interp_with_transform_no_return_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');


%% Plot fminsearch
fig = figure;
subplot(1,3,1);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,3,2); 
    if dataDim > 1
        gsp_plot_signal(G_tilde,f_tilde_int_fminsearch,param); 
    else
        plot(G_tilde.coords, f_tilde_int_fminsearch, '.');
        xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G_tilde.coords(sample_ind), f_tilde_sampled_padded(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end    
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Interpolated graph-signal on $\tilde{G}$' newline ...
           'using fminsearch' newline...
           '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_int_fminsearch'*G_tilde.L*f_tilde_int_fminsearch, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,3,3); 
    if dataDim > 1
        gsp_plot_signal(G,f_tilde_int_fminsearch,param); 
    else
        plot(G.coords, f_tilde_int_fminsearch, '.');
        xlabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Our interpolation (with transformation)' newline ...
           '$f_{{\bf int}}^T L f_{{\bf int}}$ = ' num2str(f_tilde_int_fminsearch'*G.L*f_tilde_int_fminsearch, '%.3f') newline ...
           'error: ' num2str(norm(f - f_tilde_int_fminsearch)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
set(gcf,'Position', [400 100 1200 400])   
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_fminsearch_interp_with_transform_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');


fig = figure;
subplot(1,2,1);
    if dataDim > 1
        gsp_plot_signal(G,f,param); 
    else
        plot(G.coords, f, '.');
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G.coords(sample_ind), f(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Original graph-signal' newline ...
           '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline...
           'Sampling ratio = ' num2str(samplingRatio, '%.2f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(1,2,2); 
    if dataDim > 1
        gsp_plot_signal(G_tilde,f_tilde_int_fminsearch,param); 
    else
        plot(G_tilde.coords, f_tilde_int_fminsearch, '.');
        xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
    end
    if b_plotSamplingPointsMarkers
        hold on; 
        if dataDim == 1
            scatter(G_tilde.coords(sample_ind), f_tilde_sampled_padded(sample_ind), 'ko')
        elseif dataDim == 2
            scatter(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), 'ko')
        elseif dataDim == 3
            scatter3(G_tilde.coords(sample_ind,1), G_tilde.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
        end    
    end
    view(0,90)
    set(gca,'FontSize', 14);
    title(['Interpolated graph-signal on $\tilde{G}$' newline ...
           'using fminsearch' newline...
           '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_int_fminsearch'*G_tilde.L*f_tilde_int_fminsearch, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
set(gcf,'Position', [400 100 800 400])   
saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_', verticesTransformation, '_fminsearch_interp_with_transform_no_return_', num2str(round(samplingRatio*100), '%d'), 'prec_sampling'), 'epsc');




% fig = figure;
% subplot(2,3,1); 
%     if dataDim > 1
%         gsp_plot_signal(G,f,param);
%     else
%         plot(G.coords, f, '.');
%         xlabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Original graph-signal' newline ...
%            '$f(v)$' newline ...
%            '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
% %     hold on
% %     if dataDim == 3
% %         text(v(plot_v_indexes,1), v(plot_v_indexes,2), v(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
% %     else
% %         text(v(plot_v_indexes,1), v(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
% %     end     
%     set(gca,'FontSize', 14);
% subplot(2,3,2);    
%     if dataDim > 1
%         gsp_plot_signal(G_tilde,f_tilde_int_eigrls,param); 
%     else
%         plot(G_tilde.coords, f_tilde_int_eigrls, '.');
%         xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Interpolated graph-signal on $\tilde{G}$' newline ...
%            '$\tilde{f}_{\bf int}(\tilde{v}) = \tilde{{\bf \Phi}} c \approx {\bf R} f(v)$' newline ...
%            '$\tilde{f}_{\bf int}^T \tilde{L} \tilde{f}_{\bf int}$ = ' num2str(f_tilde_int_eigrls'*G_tilde.L*f_tilde_int_eigrls, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
% %     hold on
% %     if dataDim == 3
% %         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_tilde(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
% %     else
% %         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
% %     end
%     set(gca,'FontSize', 14);
% subplot(2,3,5); 
%     if dataDim > 1
%         gsp_plot_signal(G,f_tilde_int_eigrls,param); 
%     else
%         plot(G.coords, f_tilde_int_eigrls, '.');
%         xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Interpolated graph-signal on G (EigRLS)' newline ...
%            '$\tilde{f}(\tilde{v}) = \tilde{{\bf \Phi}} c \approx f(v)$' newline ...
%            '$\tilde{f}^T L \tilde{f}$ = ' num2str(f_tilde_int_eigrls'*G.L*f_tilde_int_eigrls, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
%     set(gca,'FontSize', 14);
% subplot(2,3,6);     
%     if dataDim > 1
%         gsp_plot_signal(G_tilde,f_tilde_int_fmin,param); 
%     else
%         plot(G_tilde.coords, f_tilde_int_fmin, '.');
%         xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['fminsearch on $\tilde{G}$' newline ...
%            '$\tilde{f}(\tilde{v}) = \tilde{{\bf \Phi}} c$      $[f^T L f \approx c^T \tilde{{\bf \Phi}}^T \tilde{L} \tilde{{\bf \Phi}} c]$' newline ...
%            '$\tilde{f}^T \tilde{L} \tilde{f}$ = ' num2str(f_tilde_int_fmin'*G_tilde.L*f_tilde_int_fmin, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
%     set(gca,'FontSize', 14);    
% set(gcf,'Position', [400 100 1400 800])  
% saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_transformation'), 'epsc');



% fig = figure;
% subplot(1,3,1); 
%     if dataDim > 1
%         gsp_plot_signal(G,f,param);
%     else
%         plot(G.coords, f, '.');
%         xlabel('$v$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Original graph-signal' newline '$f(v)$'], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)   
%     set(gca,'FontSize', 14);
% subplot(1,3,2); 
%     if dataDim > 1
%         gsp_plot_signal(G_tilde,Rf,param); 
%     else
%         plot(G_tilde.coords, Rf, '.');
%         xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Transformed graph-signal' newline '$\tilde{f}(\tilde{v}) = {\bf R} f(v)$'], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
%     set(gca,'FontSize', 14);
% subplot(1,3,3); 
%     if dataDim > 1
%         gsp_plot_signal(G_rec,f_rec,param); 
%     else
%         plot(G_rec.coords, f_rec, '.');
%         xlabel('$\tilde{v}$', 'Interpreter', 'latex', 'FontSize', 14);
%     end
%     title(['Reconstructed graph-signal' newline '$\hat{f}(v) = {\bf R}^T \tilde{f}(\tilde{v})$'], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
%     set(gca,'FontSize', 14);
% set(gcf,'Position', [400 100 1400 400])  
% saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_reconstruction'), 'epsc');




% plot_v_indexes = (1:3:20)';
% v_numers_cells = cellstr(num2str(plot_v_indexes))';
% figure;
% subplot(2,3,1); 
%     gsp_plot_signal(G,f,param); 
%     title(['Original graph-signal' newline '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     hold on;
%     if dataDim == 3
%         text(v(plot_v_indexes,1), v(plot_v_indexes,2), v(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
%     else
%         text(v(plot_v_indexes,1), v(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
%     end
%     view(0,90)
% subplot(2,3,2); 
%     gsp_plot_signal(G_tilde,f_tilde_sampled_padded,param); 
%     title(['Transformed graph-signal' newline '$\tilde{f}^T \tilde{L} \tilde{f}$ = ' num2str(f_tilde_sampled_padded'*G_tilde.L*f_tilde_sampled_padded, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     hold on;
%     if dataDim == 3
%         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_tilde(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
%     else
%         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
%     end
%     view(0,90)
% subplot(2,3,3); 
%     gsp_plot_signal(G_rec,f_rec,param); 
%     title(['Reconstructed graph-signal' newline '$\hat{f}^T L \hat{f}$ = ' num2str(f_rec'*G.L*f_rec, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
% subplot(2,3,5); 
%     gsp_plot_signal(G_tilde,f_tilde_int,param); 
%     title(['Transformed graph-signal 2' newline '$\tilde{f}^T \tilde{L} \tilde{f}$ = ' num2str(f_tilde_int'*G_tilde.L*f_tilde_int, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     hold on;
%     if dataDim == 3
%         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_tilde(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
%     else
%         text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
%     end
%     view(0,90)
% subplot(2,3,6); 
%     gsp_plot_signal(G_rec,f_int,param); 
%     title(['Reconstructed graph-signal 2' newline '$\hat{f}^T L \hat{f}$ = ' num2str(f_int'*G.L*f_int, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
%     view(0,90)
% set(gcf,'Position', [400 100 1200 800])   
