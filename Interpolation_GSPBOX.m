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

b_gaussianDataset = false;
b_normlizedLaplacian = true;
b_useRandomProjections = true;

%% Generate dataset
if b_gaussianDataset
    N = 1000;
    dataDim = 2;
    sDataset = GenerateDataset('Gaussian', dataDim, nComponents, N, 0);
    GMMRegVal = 0;
else
    N = 500;
%     N = round(sqrt(N))^2; G_tmp = gsp_2dgrid(sqrt(N),sqrt(N));
%     G_tmp = gsp_sensor(N);
    G_tmp = gsp_bunny(); N = G_tmp.N;
%     G_tmp = gsp_minnesota(); N = G_tmp.N;
%     G_tmp = gsp_david_sensor_network(); N =G_tmp.N;
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
W = sparse(CalcAdjacency(sKernelParams, v));
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
if b_showEigenFigures
    [V, vLambdaNumeric] = CalcNumericEigenvectors(nEigs, sKernelParams, v);
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
samplingRatio = 0.1;

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
f = G.U*(f_hat);

figure; 
subplot(1,2,1)
gsp_plot_signal_spectral(G,f_hat,param);
xlabel('$\lambda_k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
subplot(1,2,2)
stem(f_hat(1:k0+3));
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('$\hat{f}$', 'Interpreter', 'latex', 'FontSize', 16);
sgtitle('Graph signal spectrum')
set(gcf,'Position', [100 200 1200 400])


%% Sample the signal
% this:
randPerm = randperm(N)';
sample_ind = sort(randPerm(1:round(samplingRatio*N)));

% or this:
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
    gsp_plot_signal(G,f,param); 
    hold on; 
    if dataDim == 2
        scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
    elseif dataDim == 3
        scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
    end
    view(0,90)
    title(['Original graph-signal' newline '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f') newline], 'Interpreter', 'latex', 'FontSize', 14);
subplot(132); 
    gsp_plot_signal(G,gsp_interp_signal,param); 
    hold on; 
    if dataDim == 2
        scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
    elseif dataDim == 3
        scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
    end
    view(0,90)
    title(['Pesenson''s interpolation' newline '$f^T L f$ = ' num2str(gsp_interp_signal'*G.L*gsp_interp_signal, '%.3f') newline ...
        'error: ' num2str(norm(f - gsp_interp_signal)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
subplot(133); 
    gsp_plot_signal(G,my_interp_signal,param);
    hold on; 
    if dataDim == 2
        scatter(G.coords(sample_ind,1), G.coords(sample_ind,2), 'ko')
    elseif dataDim == 3
        scatter3(G.coords(sample_ind,1), G.coords(sample_ind,2), G.coords(sample_ind,3), 'ko')
    end
    title(['Our interpolation (no transformation)' newline '$f^T L f$ = ' num2str(my_interp_signal'*G.L*my_interp_signal, '%.3f') newline ...
        'error: ' num2str(norm(f - my_interp_signal)/norm(f), '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    view(0,90)
set(gcf,'Position', [400 100 1200 400])   

saveas(fig,strcat(sSimParams.outputFolder, filesep, graphName, '_interp_no_transform'), 'epsc');

% subplot(223); 
%     gsp_plot_signal(G,eigenvecs_interp_signal,param); 
%     title(sprintf('Using V\n error: %.3f\n', norm(f - eigenvecs_interp_signal)/norm(f)));   
%     view(0,90)
% set(gcf,'Position', [400 100 800 800])


% fprintf('GSP interpolation error: %.3f\n', norm(f - gsp_interp_signal)/norm(f));
% fprintf('V   interpolation error: %.3f\n', norm(f - eigenvecs_interp_signal)/norm(f));
% fprintf('Our interpolation error: %.3f\n', norm(f - my_interp_signal)/norm(f));

%% Transform
% "Random projection"
if b_useRandomProjections
    newDim = dataDim;
    R = randn(dataDim,newDim); % Random projection matrix
    for i = 1:dataDim
       R(i,:) = R(i,:)/norm(R(i,:)); % Normalization
    end
    v_tilde = v*R; % Projection

else
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

    % % Whitening matrix
    % [v_tilde, R, invR] = whiten(v,1e-5);
end

figure; 
subplot(dataDim+1,2,1)
histogram(v(:,1),100);
title('v1 histogram');
subplot(dataDim+1,2,2)
histfit(v_tilde(:,1),100);
title('v tilde1 histogram');
subplot(dataDim+1,2,3)
histogram(v(:,2),100);
title('v2 histogram');
subplot(dataDim+1,2,4)
histfit(v_tilde(:,2),100);
title('v tilde2 histogram');
if dataDim == 2
    subplot(dataDim+1,2,5)
    hist3(v,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
    colormap(gca, 'hot')
    colorbar()
    view(2)
    xlim([ min(v(:,1)) max(v(:,1))])
    ylim([ min(v(:,2)) max(v(:,2))])
    title('v histogram');
end
if size(v_tilde,2) == 2
    subplot(dataDim+1,2,6)
    hist3(v_tilde,'CdataMode','auto', 'Nbins', [50 50], 'edgecolor', 'flat');
    colormap(gca, 'hot')
    colorbar()
    view(2)
    xlim([ min(v_tilde(:,1)) max(v_tilde(:,1))])
    ylim([ min(v_tilde(:,2)) max(v_tilde(:,2))])
    title('v tilde3 histogram');
end
set(gcf,'Position', [100 200 1200 800])

% Distances
ii = 5; jj = 60;
[pdist2(v(ii,:),v(jj,:)) pdist2(v_tilde(ii,:),v_tilde(jj,:))]

sDatasetTilde.sData.x = v_tilde;
sDatasetTilde.estDataDist = 'Gaussian';
sDatasetTilde.dim = size(v_tilde,2);
GMMRegVal = 0;
nComponents = 1;
% Estimate distribution and get kernel parameters
sDistParamsTilde = EstimateDistributionParameters(sDatasetTilde, nComponents, GMMRegVal);
sKernelParamsTilde = GetKernelParams(sDatasetTilde, sDistParamsTilde, omega);
[sKernelParamsTilde.vLambdaAnaytic, sKernelParamsTilde.vComponentIndex, sKernelParamsTilde.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParamsTilde, sDatasetTilde.dim, nComponents);
% Adjacency (gaussian kernel)
W_tilde = sparse(CalcAdjacency(sKernelParamsTilde, v_tilde));


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

%% Find coeffs for graph-signal
if b_useRandomProjections
    f_tilde = f; % Assuming the distances are preserved
else
    f_tilde = R*f;
end

tilde_gamma_A_eigrls = 0;0;
tilde_gamma_I_eigrls = 0;0.1;
c_tilde = eigrls(f_tilde, mPhiTildeAnalytic, vLambdaTildeAnalytic, tilde_gamma_A_eigrls, tilde_gamma_I_eigrls, G.L);
% c_tilde = lsqminnorm(mPhiTildeAnalytic,R*f);
% c_tilde = (mPhiTildeAnalytic.'*mPhiTildeAnalytic)\(mPhiTildeAnalytic.'*(R*f));
% c_tilde = pinv(mPhiTildeAnalytic)*R*f;
% f_tilde = mPhiTildeAnalytic*c_tilde;

f_tilde2 = mPhiTildeAnalytic*c_tilde;


%% Reconstruction
if b_useRandomProjections
    f_rec = f;
    f_rec2 = f;
    v_rec = v;
else
    %This is the same: [mPhiTildeAnalytic.'*f_tilde mPhiTildeAnalytic.'*R*f]
    % keyboard;


    % c_rec = (mPhiTildeAnalytic.'*mPhiTildeAnalytic)\(mPhiTildeAnalytic.'*(R'*f_tilde));
    % f_reconstructed = mPhiTildeAnalytic*c_rec;
    % f_reconstructed = pinv(mPhiTildeAnalytic.'*R)*mPhiTildeAnalytic.'*f_tilde;
    % invR = R';
    % f_reconstructed = (invR.'*invR)\(invR.'*f_tilde);
    % f_reconstructed = lsqminnorm(R',f_tilde);
    f_rec = R'*f_tilde;
    f_rec2 = R'*f_tilde2;
    v_rec = R'*v_tilde;
end


sDatasetRec.sData.x = v_rec;
sDatasetRec.estDataDist = 'Gaussian';
sDatasetRec.dim = size(v_rec,2);
GMMRegVal = 0;
nComponents = 1;

% Estimate distribution and get kernel parameters
sDistParamsRec = EstimateDistributionParameters(sDatasetRec, nComponents, GMMRegVal);
sKernelParamsRec = GetKernelParams(sDatasetRec, sDistParamsRec, omega);
[sKernelParamsRec.vLambdaAnaytic, sKernelParamsRec.vComponentIndex, sKernelParamsRec.vEigIndex] ...
    = CalcAnalyticEigenvalues(nEigs, sKernelParamsRec, sDatasetRec.dim, nComponents);

% Adjacency (gaussian kernel)
W_rec = sparse(CalcAdjacency(sKernelParamsRec, v_rec));


% Graph
G_rec = gsp_graph(W_rec, v_rec);
if b_normlizedLaplacian
    G_rec.lap_type = 'normalized';
    G_rec = gsp_graph_default_parameters(G_rec);
end

%% Plot
plot_v_indexes = (1:3:20)';
v_numers_cells = cellstr(num2str(plot_v_indexes))';
figure;
subplot(2,3,1); 
    gsp_plot_signal(G,f,param); 
    title(['Original graph-signal' newline '$f^T L f$ = ' num2str(f'*G.L*f, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    hold on;
    if dataDim == 3
        text(v(plot_v_indexes,1), v(plot_v_indexes,2), v(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
    else
        text(v(plot_v_indexes,1), v(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
    end
    view(0,90)
subplot(2,3,2); 
    gsp_plot_signal(G_tilde,f_tilde,param); 
    title(['Transformed graph-signal' newline '$\tilde{f}^T \tilde{L} \tilde{f}$ = ' num2str(f_tilde'*G_tilde.L*f_tilde, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    hold on;
    if dataDim == 3
        text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_tilde(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
    else
        text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
    end
    view(0,90)
subplot(2,3,3); 
    gsp_plot_signal(G_rec,f_rec,param); 
    title(['Reconstructed graph-signal' newline '$\hat{f}^T L \hat{f}$ = ' num2str(f_rec'*G.L*f_rec, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    view(0,90)
subplot(2,3,5); 
    gsp_plot_signal(G_tilde,f_tilde2,param); 
    title(['Transformed graph-signal 2' newline '$\tilde{f}^T \tilde{L} \tilde{f}$ = ' num2str(f_tilde2'*G_tilde.L*f_tilde2, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    hold on;
    if dataDim == 3
        text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_tilde(plot_v_indexes,3), v_numers_cells,'FontWeight','bold','Color', 'black')
    else
        text(v_tilde(plot_v_indexes,1), v_tilde(plot_v_indexes,2), v_numers_cells,'FontWeight','bold','Color', 'black')
    end
    view(0,90)
subplot(2,3,6); 
    gsp_plot_signal(G_rec,f_rec2,param); 
    title(['Reconstructed graph-signal 2' newline '$\hat{f}^T L \hat{f}$ = ' num2str(f_rec2'*G.L*f_rec2, '%.3f')], 'Interpreter', 'latex', 'FontSize', 14);
    view(0,90)
set(gcf,'Position', [400 100 1200 800])   
